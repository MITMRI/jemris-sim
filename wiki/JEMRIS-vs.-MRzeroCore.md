# JEMRIS vs. MRzeroCore

JEMRIS and MRzeroCore simulate the same underlying MR physics using different
numerical representations. JEMRIS integrates the Bloch equations for a spin
ensemble. MRzeroCore evolves Phase Distribution Graph (PDG) states and exposes
the calculation through PyTorch.

Framework-specific details remain on the [JEMRIS Inputs](Inputs) and
[MRzeroCore Inputs](MRzero-Inputs) pages. This page is the dedicated comparison.

## Detailed comparison

The **availability** column is deliberately explicit:

- **Both** means both standard simulation paths contain the feature, although
  their models may differ.
- **JEMRIS only** or **MRzeroCore only** means the other simulator has no
  equivalent feature in its standard simulation path. It may still be possible
  to add one by extending the source.
- **Neither standard** means that both require an extension.
- **Different models** means that both expose a related feature, but the two
  calculations are not physically interchangeable without validation.

### Numerical representation and execution

| Capability | JEMRIS | MRzeroCore | Availability | Consequence |
|---|---|---|---|---|
| Fundamental state | One three-component magnetization vector per live spin and, for exchange, per pool | Longitudinal and transverse phase-distribution states attached to voxels and coherence pathways | **Different models** | JEMRIS samples magnetization directly in real space; MRzeroCore compresses common spatial phase patterns into PDG states |
| Bloch evolution | CVODE adaptively integrates the Bloch ODE through time | Analytical RF, relaxation, phase, diffusion, and motion operators evolve retained PDG states | **Different models** | Agreement is expected only where both physical models describe the same experiment |
| Continuously sampled local field | The ODE evaluates sequence and field values at the solver's internal time points | Standard evolution uses event-level parameters and integrated gradient moments | **JEMRIS only** | JEMRIS can directly resolve simultaneous waveform effects; MRzeroCore needs an appropriate event-level operator |
| Adaptive local-error control | CVODE changes its step size to satisfy integration tolerances | No ODE step-size tolerance; accuracy is controlled by event representation, voxel sampling, and graph thresholds | **JEMRIS only** | Tightening a PDG pruning threshold is not equivalent to tightening an ODE tolerance |
| Coherence-pathway graph | Not used; all spin vectors are evolved directly | `compute_graph()` discovers, merges, and prunes longitudinal and transverse pathways | **MRzeroCore only** | MRzeroCore avoids evolving pathways predicted to have negligible latent or emitted signal |
| Automatic differentiation | The classic C++/CVODE simulator does not expose derivatives of the signal with respect to inputs | Supported PyTorch operations can propagate gradients from signal or reconstruction losses back to sequence and phantom tensors | **MRzeroCore only** | Native gradient-based sequence or reconstruction optimization is a distinguishing MRzeroCore capability |
| Spatial discretization | Cartesian sample grid or explicit spin ensemble; each nonzero-`M0` entry is a live spin | Signal-producing voxel grid or `CustomVoxelPhantom` entries, each carrying amplitudes for retained PDG states | **Both, implemented differently** | JEMRIS cost scales mainly with spins, pools, and ODE work; MRzeroCore cost scales with voxels, retained states, and events |
| CPU parallel execution | MPI distributes spin ensembles across processes | Vectorized PyTorch and compiled components execute over voxels and states | **Different models; native MPI is JEMRIS only** | Large jobs are partitioned differently and do not share a common scaling curve |
| Standard GPU path | Classic JEMRIS is CPU/MPI; GPU-JEMRIS is a separate research extension | Supported tensors and simulation paths run on CUDA through PyTorch | **MRzeroCore only in the standard workflows compared here** | MRzeroCore fits more directly into a Python CUDA pipeline |

### Sequence, RF, and magnetic fields

| Capability | JEMRIS | MRzeroCore | Availability | Consequence |
|---|---|---|---|---|
| Native sequence representation | Hierarchical XML sequence tree with analytic pulses and optional sampled HDF5 waveforms | Repetitions containing an RF event, gradient moments, durations, ADC flags, and phases | **Different models** | JEMRIS preserves waveform structure; MRzeroCore is naturally event- and moment-oriented |
| Pulseq `.seq` input | This repository supplies a converter that produces JEMRIS XML and waveform files; classic JEMRIS does not directly execute `.seq` | `Sequence.import_file()` and the high-level interface directly accept Pulseq input | **MRzeroCore only for direct import; both through this repository's workflow** | Pulseq can drive both, but JEMRIS has an explicit conversion step |
| Ideal linear gradients | Gradient waveforms generate a position-dependent longitudinal field at every spin | Integrated moments translate transverse PDG states in k-space | **Both** | Both model ordinary linear encoding when moments, timing, gyromagnetic-ratio conventions, and units match |
| Arbitrary gradient waveform shape between RF/ADC events | Sampled or analytic waveform is evaluated during ODE integration | The standard free-precession operator primarily uses duration and integrated gradient moment; additional event subdivision can retain more trajectory detail | **JEMRIS only for fully waveform-resolved Bloch evolution** | Equal total moment can still produce different diffusion or simultaneous-field behavior |
| Nonlinear gradient fields | `NLG_field` expressions are evaluated at each spin | Standard `rep.gradm` produces one linear spatial phase ramp | **JEMRIS only** | A nonlinear encoding field cannot generally be represented by one global PDG k-space shift |
| Concomitant fields | Built-in concomitant-field contribution can be enabled in the simulation configuration | No standard concomitant-field term in the PDG gradient operator | **JEMRIS only** | Maxwell-field validation requires JEMRIS or a new MRzeroCore operator |
| Scalar three-dimensional B0/off-resonance map | Per-spin or per-grid `DB` supplies a spatially varying longitudinal offset | `phantom.B0[nx,ny,nz]` supplies one longitudinal offset in hertz per voxel | **Both** | “Three-dimensional” describes the map domain, not three components of the B vector |
| Arbitrary three-component static magnetic field | The standard model accepts RF in the transverse plane and gradients/off-resonance along its longitudinal axis; adding static Bx and By requires a field-model/source extension | The standard PDG model assumes free precession about one common longitudinal axis; transverse static components require a new or hybrid propagator | **Neither standard** | JEMRIS's direct vector ODE is the more natural extension point, but this feature is not already exposed as a standard vector-B0 input in either simulator |
| Time-varying off-resonance/tissue trajectories | JEMRIS has dynamic-variable trajectories, including T2-prime/off-resonance-related evolution | Standard `phantom.B0` is fixed for a simulation; motion can change positions but does not provide an equivalent general time-varying B0 tensor | **JEMRIS only as a built-in model** | Dynamic field experiments require a custom MRzeroCore extension or explicit sequence/time segmentation |
| RF waveform resolved during excitation | RF amplitude and phase waveforms, gradients, and local offsets are evaluated together by the ODE solver | The standard PDG path mixes states using effective flip angle and phase; imported waveform details are generally reduced to pulse parameters, with specialized handling for selected off-resonant/selective cases | **JEMRIS only for general waveform-resolved integration** | Long, strongly off-resonant, adiabatic, or simultaneous-gradient pulses should not be assumed equivalent without validation |
| Spatial B1+ and parallel transmit | TX coil-array sensitivities and channel waveforms enter the local RF field | Complex `phantom.B1` maps and RF shim weights modify voxelwise flip angle and phase | **Both** | Match map normalization, channel phases, and voltage-to-field conventions |
| Receive coil sensitivities | RX coil-array definitions produce separate receive channels | Complex `phantom.coil_sens[n_rx,nx,ny,nz]` maps weight voxel signals | **Both** | Both sum coil-weighted transverse magnetization, but their normalization conventions must be aligned |

### Sample physics and dynamics

| Capability | JEMRIS | MRzeroCore | Availability | Consequence |
|---|---|---|---|---|
| Proton density/equilibrium magnetization | `M0` is stored for every live spin or pool | `PD` weights equilibrium signal for each voxel | **Both** | Voxel volume and amplitude normalization must be matched before comparing absolute signal |
| T1 and T2 relaxation | Standard samples store the rates `R1 = 1/T1` and `R2 = 1/T2` | Phantoms store `T1` and `T2` time constants in seconds | **Both** | Convert rates and time constants instead of copying numbers directly |
| Reversible off-resonance and intravoxel dephasing | `DB`, spin sampling, and configured T2-star-related terms contribute to phase dispersion | `B0`, `T2dash`, voxel shape, and dephasing functions contribute analytical phase and attenuation factors | **Different models** | Matching a single nominal T2-star value does not guarantee the same echo shape |
| Multiple non-exchanging tissues | Different spins or voxels can carry different tissue parameters | Every voxel can carry different `PD`, `T1`, `T2`, `T2dash`, `D`, `B0`, and B1 values | **Both** | Heterogeneous anatomy is supported in both |
| Chemical exchange between pools | The Bloch–McConnell CVODE model evolves multiple coupled pools and exchange rates | The standard PDG simulation has no equivalent coupled-pool exchange operator | **JEMRIS only** | A multi-pool phantom without exchange is not a substitute for Bloch–McConnell dynamics |
| Free Gaussian diffusion attenuation | Can be modeled through diffusion trajectories/models | Analytical attenuation is computed from gradient history and voxel diffusion coefficient `D` | **Both, implemented differently** | Equivalent scalar-D experiments can be compared after checking units and gradient timing |
| Restricted, anisotropic, or microstructure trajectory diffusion | Dedicated trajectory code includes fibers, crossings, kissing fibers, compartments, and box geometries | Standard PDG diffusion uses the supplied voxel diffusion property and does not perform stochastic restricted random walks through microstructure | **JEMRIS only** | Microstructure/restriction studies need JEMRIS or a new MRzeroCore diffusion model |
| Prescribed object/voxel motion | Motion and spatial trajectories update spin positions | A voxel-position function can supply moving voxel coordinates and corresponding phase | **Both** | Match whether field maps move with the object or remain fixed in scanner coordinates |
| Dedicated flow and respiration models | Built-in trajectory classes include flow and respiratory motion | No equivalent dedicated physiology modules in the standard API; prescribed motion can be implemented through a custom position function | **JEMRIS only as built-in models** | MRzeroCore can reproduce selected trajectories, but it does not provide the same ready-made model set |
| Time-varying T1, T2, M0, and T2-prime | Dynamic-variable trajectories can change these properties during a run | Standard main-pass tissue tensors are fixed during a run | **JEMRIS only** | Dynamic contrast or relaxation-changing experiments require custom MRzeroCore work |

### Inputs, outputs, and intended workflow

| Capability | JEMRIS | MRzeroCore | Availability | Consequence |
|---|---|---|---|---|
| Experiment assembly | `simu.xml` explicitly links sample HDF5, sequence XML, TX array, RX array, model, and options | Python passes sequence and phantom objects to `mr0.util.simulate()` | **Different interfaces** | JEMRIS produces a declarative scanner bundle; MRzeroCore favors concise programmatic experiments |
| Immediate signal result | Native signal HDF5 and optional downstream formats | Complex PyTorch signal tensor returned in memory | **Both** | JEMRIS favors portable files; MRzeroCore favors immediate numerical processing |
| Explicit k-space trajectory result | Can be derived or exported through sequence diagnostics and downstream tooling | High-level simulation returns k-space alongside signal | **MRzeroCore only as a direct standard return value** | MRzeroCore makes simulation-to-reconstruction handoff especially short |
| ISMRMRD interoperability | JEMRIS can write/export ISMRMRD for downstream tools | MRzeroCore includes signal-to-ISMRMRD conversion utilities | **Both** | ISMRMRD is a useful neutral output when comparing reconstruction pipelines |
| End-to-end differentiable reconstruction | Reconstruction is normally a separate downstream program and is not connected to CVODE derivatives | Simulation, FFT/NUFFT, learned reconstruction, and loss evaluation can remain in one PyTorch computation graph where operations support autograd | **MRzeroCore only** | This is a major reason to select MRzeroCore for sequence/reconstruction co-optimization |
| Declarative self-contained experiment record | XML/HDF5 inputs naturally capture the configured virtual scanner | Reproducibility requires preserving Python code, tensor data, package versions, device behavior, and random state | **JEMRIS only as the native configuration style** | Either can be reproducible, but the artifacts that must be archived differ |
| Typical strength | Detailed forward simulation with waveform, hardware-field, exchange, and trajectory models | Rapid differentiable simulation, parameter sweeps, reconstruction research, and optimization | **Complementary** | Use the simulator whose standard physics contains the effect under study; cross-validate where their supported models overlap |

### Features present in only one standard framework

The most important asymmetric capabilities are therefore:

| JEMRIS standard feature with no standard MRzeroCore equivalent | MRzeroCore standard feature with no classic JEMRIS equivalent |
|---|---|
| Nonlinear gradient fields | PyTorch automatic differentiation through supported simulation inputs |
| Concomitant-field calculation | End-to-end differentiable simulation and reconstruction |
| General waveform-resolved RF/gradient/B0 interaction through an adaptive ODE | PDG coherence-pathway construction, merging, and pruning |
| Bloch–McConnell multi-pool chemical exchange | Standard CUDA tensor workflow |
| Restricted and microstructure-based diffusion trajectories | Signal and k-space returned directly as tensors to Python |
| Built-in flow, respiration, and time-varying tissue-property trajectories | Direct Pulseq import without an XML conversion stage |

This table describes the standard paths documented here, not a claim that an
absent feature is impossible to implement. In particular, direct Bloch
integration makes some new local-field terms comparatively natural to add to
JEMRIS, while PyTorch and the PDG operator structure make new differentiable
event operators comparatively natural to add to MRzeroCore.

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
