# JEMRIS vs. MRzeroCore

JEMRIS and MRzeroCore simulate the same underlying MR physics using different
numerical representations. JEMRIS integrates the Bloch equations for a spin
ensemble. MRzeroCore evolves Phase Distribution Graph (PDG) states and exposes
the calculation through PyTorch.

Framework-specific details remain on the [JEMRIS Inputs](Inputs) and
[MRzeroCore Inputs](MRzero-Inputs) pages. This page is the dedicated comparison.

## Detailed comparison

| Aspect | JEMRIS | MRzeroCore | Practical implication |
|---|---|---|---|
| Fundamental method | Numerically integrates the Bloch ODE for each live spin using CVODE | Evolves analytical PDG transverse and longitudinal states | JEMRIS emphasizes general forward physics; MRzeroCore emphasizes efficient state evolution |
| Simulated state | Three-component magnetization for each spin or compartment | Collections of phase-distribution states associated with voxels and coherence pathways | PDG compresses phase behavior that otherwise requires many explicitly sampled spins |
| Time evolution | Adaptive numerical integration through sequence waveforms | Event-level RF, relaxation, off-resonance, diffusion, and gradient operators | JEMRIS can resolve detailed continuous waveforms; MRzeroCore can skip unnecessary fine time stepping |
| Spatial discretization | Cartesian sample grid or explicit spin ensemble; every nonzero-`M0` entry becomes a live spin | Voxel grid reduced to signal-producing voxels, or explicit `CustomVoxelPhantom` entries | Runtime scales with different state counts in the two methods |
| Primary interface | Compiled simulator configured by XML and HDF5 | Python package using PyTorch tensors and objects | JEMRIS experiments are declarative file bundles; MRzeroCore experiments are readily scripted and modified |
| Experiment wiring | `simu.xml` connects the sample, sequence, TX array, RX array, and model | Python passes a sequence and phantom to `mr0.util.simulate()` | JEMRIS makes the virtual scanner components explicit; MRzeroCore minimizes orchestration code |
| Sequence representation | Hierarchical JEMRIS sequence XML with analytic or HDF5 waveform elements | Native `mr0.Sequence`, PyPulseq sequence object, or Pulseq `.seq` file | Both can participate in Pulseq workflows, but their internal sequence models differ |
| Pulseq support | This repository converts `.seq` files into XML and referenced waveform files | High-level simulation accepts PyPulseq objects and `.seq` paths | Pulseq is a useful exchange format for running related experiments in both tools |
| Static B0 input | Per-spin/grid `DB` value plus optional model-specific field contributions | Scalar `phantom.B0[nx,ny,nz]` off-resonance tensor in Hz | Both normally need a rotating-frame offset rather than the absolute field in tesla |
| Three-dimensional B0 map | Represented through spatially varying sample/model values evaluated for individual spins | Represented as a scalar value at every voxel | A 3D scalar map is distinct from a three-component vector field |
| Arbitrary vector field | Direct Bloch integration is a natural basis for spatially varying effective-field directions, provided the field model supplies them | Standard PDG free precession uses scalar off-resonance about a common longitudinal axis | Strong transverse static components require more than MRzeroCore's standard scalar B0 map |
| Ideal linear gradients | Sequence gradients generate position-dependent spin evolution | Gradient moments shift PDG states and advance the k-space trajectory | Both reproduce ordinary linear spatial encoding |
| Gradient moment | Can be derived from the integrated ideal gradient waveform | Stored directly in native sequence events as `rep.gradm` | MRzeroCore's native representation is naturally event/moment based |
| Nonlinear gradients | Spatial gradient fields can be evaluated locally for each spin; configured models can include nonlinear effects | Native PDG gradient operator assumes a linear spatial phase ramp | Nonlinear encoding cannot generally be reduced to one global PDG k-space shift |
| Concomitant fields | Supported through the corresponding JEMRIS field/model configuration | Not part of the standard `rep.gradm` operator | Detailed hardware-field validation favors a direct local-field calculation |
| RF representation | RF waveform plus TX coil-array sensitivity is applied during Bloch integration | RF events mix PDG states; transmit sensitivity is stored in `phantom.B1` | The PDG operator is efficient, while numerical integration naturally represents detailed time-domain rotations |
| Shaped and off-resonant RF | Numerical integration follows the configured waveform and local effective field | Supported workflows depend on the MRzeroCore pulse/event representation and current PDG RF model | Validate specialized pulses against the physics supported by the chosen MRzeroCore version |
| Receive coils | RX coil-array XML defines sensitivities and signal channels | `phantom.coil_sens[n_rx,nx,ny,nz]` defines complex receive maps | Both form coil-weighted transverse signals, but use different input containers |
| Equilibrium magnetization | `M0` for each live spin/grid point | `PD` provides proton-density/equilibrium signal weighting | Scaling and voxel-volume normalization must be aligned before comparing amplitudes |
| T1 input | Stored as `R1 = 1/T1` in the standard sample | Stored directly as `T1` in seconds | Convert rates and time constants carefully when recreating a phantom |
| T2 input | Stored as `R2 = 1/T2` | Stored directly as `T2` in seconds | Time/rate units are a common source of mismatched results |
| T2-star/dephasing input | `R2S = 1/T2*` and spatial `DB` contribute to dephasing | `T2dash` and `B0` represent reversible dephasing/off-resonance behavior | Parameters are related physically but should not be copied field-for-field without conversion |
| Diffusion | Available through configured JEMRIS physical models | Included in the PDG simulation framework through voxel diffusion properties | Confirm that both runs use equivalent diffusion models and units |
| Motion | Supported through JEMRIS trajectories and model configuration | Phantom and voxel motion can be attached to MRzeroCore phantoms | Motion parameterization and numerical treatment differ |
| Multiple tissue pools | Multi-pool sample layout and Bloch-McConnell model | Availability depends on the selected MRzeroCore phantom and simulation workflow | Exchange-model equivalence must be checked explicitly |
| Differentiability | Primarily a forward simulator; parameter gradients require an external method | PyTorch automatic differentiation can propagate through supported simulation and reconstruction operations | MRzeroCore is suited to gradient-based sequence and reconstruction optimization |
| CPU parallelism | MPI distributes spin ensembles across processes | PyTorch and compiled components provide vectorized CPU execution | Scaling behavior depends on spin count versus PDG state count |
| GPU execution | GPU-JEMRIS exists as a research extension rather than the classic default workflow | PyTorch tensors and supported simulation paths can run on CUDA | MRzeroCore integrates more directly with common Python GPU workflows |
| Output format | Native signal HDF5, optional ISMRMRD, and sequence diagnostics | PyTorch tensors for signal and k-space | JEMRIS outputs are portable files; MRzeroCore outputs enter Python reconstruction immediately |
| Reconstruction workflow | Usually downstream through ISMRMRD/BART or custom readers | FFT, NUFFT, loss functions, and learned reconstruction can remain in one tensor graph | MRzeroCore is convenient for end-to-end experiments |
| Configuration complexity | More files and explicit component definitions | Smaller Python-first setup | JEMRIS exposes more scanner structure at the cost of setup overhead |
| Reproducibility style | Fixed XML/HDF5 inputs provide a declarative experiment record | Python environment, code, tensors, and random state define the experiment | Both are reproducible when their complete input conventions and versions are recorded |
| Typical use | Detailed forward simulation, hardware/field effects, and validation | Rapid simulation, parameter sweeps, reconstruction, ML, and optimization | The intended scientific question should determine the simulator |

## Mathematical distinction

### Direct Bloch evolution

JEMRIS evolves magnetization according to

```math
\frac{d\mathbf M(\mathbf r,t)}{dt}
=
\gamma\,\mathbf M(\mathbf r,t)\times\mathbf B_{\mathrm{eff}}(\mathbf r,t)
-
\begin{bmatrix}
M_x/T_2\\
M_y/T_2\\
(M_z-M_0)/T_1
\end{bmatrix}
```

The local field can include RF, linear or nonlinear gradient fields,
off-resonance, and other configured terms.

### PDG state evolution

MRzeroCore summarizes event-level evolution as

```math
\mathbf s_{n+1}
=
\mathcal D
\!\left(\Delta t_n,T_1,T_2,T_2',\Delta f,\Delta\mathbf k_n\right)
\mathcal R(\alpha_n,\phi_n)\mathbf s_n
+\mathbf b_n
```

For scalar off-resonance, a transverse state evolves as

```math
F^+(t+\Delta t)
=
F^+(t)
\exp\!\left(-\frac{\Delta t}{T_2}\right)
\exp\!\left(-i2\pi\Delta f\,\Delta t\right)
```

For an ideal linear gradient,

```math
\Delta\mathbf k_n
=
\frac{\gamma}{2\pi}
\int_{t_{n-1}}^{t_n}\mathbf G(\tau)\,d\tau,
\qquad
\mathbf k_n=\mathbf k_{n-1}+\Delta\mathbf k_n
```

The analytical state representation is the source of much of MRzeroCore's
speed, but it also defines which spatial-field effects fit its standard
operators.

## Recommended combined workflow

1. Define consistent tissue, field, sequence, coil, and unit conventions.
2. Prototype or optimize the sequence in MRzeroCore.
3. Export the sequence through Pulseq.
4. Convert Pulseq to JEMRIS XML with `tools/seq_to_jemris_xml.py`.
5. Configure the detailed JEMRIS sample, coil, and field models.
6. Compare signal phase, magnitude, k-space position, and reconstruction using
   identical ADC times and normalization.
7. Attribute remaining differences to model assumptions rather than only to
   numerical tolerances.

## References

- [JEMRIS documentation](https://www.jemris.org/)
- [MRzeroCore documentation](https://mrsources.github.io/MRzero-Core/)
- [MRzeroCore source](https://github.com/MRsources/MRzero-Core)
- [Phase distribution graphs paper](https://doi.org/10.1002/mrm.30055)

