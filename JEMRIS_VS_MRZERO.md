# JEMRIS vs. MRzeroCore

JEMRIS and MRzeroCore are complementary MRI simulation frameworks. JEMRIS
prioritizes general Bloch-equation modeling of spins and scanner fields;
MRzeroCore prioritizes fast, differentiable simulation for sequence design,
reconstruction, and optimization.

## Comparison

| Aspect | JEMRIS | MRzeroCore |
|---|---|---|
| Core simulation method | Numerically integrates the Bloch equations for individual spins, using CVODE | Uses Phase Distribution Graphs (PDG), a spatially encoded extension of phase-graph methods |
| Main design goal | General-purpose, high-fidelity forward simulation | Fast, differentiable sequence simulation and optimization |
| Typical workflow | Sample HDF5 + coil descriptions + sequence XML + simulation XML → signal HDF5 | Python/PyTorch phantom + native or Pulseq sequence → tensors for signal, k-space, and reconstruction |
| Primary interface | C++ executable configured through XML; MATLAB tools are available in upstream JEMRIS | Python package with PyTorch tensors |
| Sequence input | JEMRIS XML; this repository also converts Pulseq `.seq` files to XML | Native MRzeroCore sequences, PyPulseq objects, or Pulseq `.seq` files |
| Physical model | Direct time evolution under RF, gradients, relaxation, off-resonance, and optional physical effects | Analytical PDG state evolution representing spatial and frequency distributions |
| Differentiability | Not inherently differentiable | Differentiable through PyTorch automatic differentiation |
| Computational profile | More expensive per spin; designed to parallelize large spin ensembles with MPI | Often substantially faster for supported sequence physics; can use PyTorch CPU or GPU execution |
| Static-field input | Per-spin off-resonance (`DB`) and field-related models; this fork documents planned/general B0 map workflows | A scalar three-dimensional off-resonance map `B0[x,y,z]` expressed in Hz |
| Full vector fields | Better foundation for evolving spins under detailed effective fields, though the model and input maps must still be configured | Standard PDG evolution assumes a common longitudinal axis; arbitrary three-component fields require a new or hybrid propagator |
| Gradient modeling | Supports sequence gradients, nonlinear gradient fields, and concomitant-field modeling | Primarily sequence encoding gradients represented through PDG state shifts |
| RF modeling | Numerically follows RF waveforms and spatial transmit-coil fields | Efficient RF event/state transformations; recent PDG work extends shaped-pulse and off-resonant response modeling |
| Transmit/receive coils | Explicit TX and RX coil-array descriptions, including multiple channels | B1 and coil-sensitivity maps stored with the phantom |
| Motion and changing parameters | Supports arbitrary spatiotemporal parameter variation and motion-oriented models | Supports phantom and voxel motion within the MRzero simulation model |
| Diffusion | Available through JEMRIS physical models | Supported within the PDG framework |
| Reconstruction integration | Writes native signal HDF5 and can export ISMRMRD for downstream reconstruction | Naturally integrates with PyTorch, NumPy, NUFFT, and differentiable reconstruction pipelines |
| Machine-learning integration | Usually used to generate or validate data before a separate ML stage | Designed for end-to-end optimization with PyTorch models and loss functions |
| Scaling strategy | CPU and MPI parallelism across spin ensembles; GPU-JEMRIS research extensions also exist | Vectorized PyTorch operations on CPU/GPU plus a compiled Rust simulation component |
| Setup and configuration | More components and configuration files, but exposes detailed experiment structure | Simpler Python-first setup for rapid experiments |
| Best use cases | Hardware-aware validation, complex RF/coil behavior, nonlinear fields, detailed 3D experiments | Rapid sequence prototyping, parameter sweeps, reconstruction research, and gradient-based optimization |

## Input model and file formats

The most visible practical difference is how an experiment is assembled.
JEMRIS uses a set of files connected by a simulation-driver XML document.
MRzeroCore normally constructs the same concepts as Python objects containing
PyTorch tensors.

| Experimental input | JEMRIS representation | MRzeroCore representation |
|---|---|---|
| Simulation driver | `simu.xml` connects the sample, sequence, TX array, RX array, and physical model | A Python script passes a sequence and phantom to `mr0.util.simulate()` |
| Sequence | JEMRIS sequence-tree XML; waveform nodes may reference HDF5 arrays | `mr0.Sequence`, a PyPulseq `Sequence`, or a Pulseq `.seq` file |
| Pulseq conversion | This repository converts `.seq` into JEMRIS XML plus referenced waveform files | Pulseq is accepted directly by the high-level simulation interface |
| Tissue phantom | Sample HDF5 containing spin locations and properties such as $M_0$, $R_1$, $R_2$, $R_2^*$, and `DB` | `VoxelGridPhantom` tensors such as `PD`, `T1`, `T2`, `T2dash`, `D`, and `B0` |
| Static-field error | Per-spin `DB` data in the sample; additional field models can be configured by the JEMRIS model | `phantom.B0`, a scalar tensor in Hz with spatial shape `(nx, ny, nz)` |
| Transmit field | TX coil-array XML and associated spatial sensitivity data | Complex `phantom.B1` tensor, conventionally shaped `(n_tx, nx, ny, nz)` |
| Receive field | RX coil-array XML and sensitivity maps | Complex `phantom.coil_sens` tensor shaped `(n_rx, nx, ny, nz)` |
| Geometry | Spin coordinates and resolution stored with the sample | Phantom arrays plus a $3\times4$ affine transformation in millimeters |
| Output | Native signal HDF5 and optional ISMRMRD; sequence diagnostics are separate | PyTorch tensors, typically `(signal, kspace)`, immediately available to Python |

### JEMRIS input bundle

A typical JEMRIS experiment has this structure:

```text
experiment/
├── simu.xml            # connects every input and selects Bloch/CVODE
├── sequence.xml        # RF, gradients, timing, and ADC events
├── sample.h5           # spins, tissue properties, and DB offsets
├── tx_coils.xml        # transmit-coil definitions
├── rx_coils.xml        # receive-coil definitions
└── waveforms/          # optional RF/gradient HDF5 arrays
```

The simulation driver explicitly connects those files:

```xml
<simulate name="JEMRIS">
  <sample name="phantom" uri="sample.h5"/>
  <TXcoilarray uri="tx_coils.xml"/>
  <RXcoilarray uri="rx_coils.xml"/>
  <sequence name="sequence" uri="sequence.xml"/>
  <model name="Bloch" type="CVODE"/>
</simulate>
```

This separation makes each part of the virtual scanner explicit. It is useful
for exchanging fixed scanner, coil, and phantom descriptions, although it
requires more file management. In this repository, a Pulseq sequence can be
converted with

```bash
python tools/seq_to_jemris_xml.py input.seq -o converted -n sequence
```

which produces `sequence.xml` and any RF or gradient waveform files referenced
by that XML.

### MRzeroCore input objects

The corresponding MRzeroCore setup is normally expressed in one Python data
flow:

```python
import MRzeroCore as mr0

sequence = mr0.Sequence.import_file("sequence.seq")
phantom = mr0.VoxelGridPhantom.load("phantom.npz")

# Scalar off-resonance map: shape (nx, ny, nz), unit Hz
phantom.B0 = measured_delta_f_hz

signal, kspace = mr0.util.simulate(sequence, phantom)
```

For a voxel-grid phantom, the principal arrays are conceptually

$$
\begin{aligned}
\mathrm{PD},T_1,T_2,T_2',D,B_0
    &\in \mathbb R^{N_x\times N_y\times N_z},\\
B_1
    &\in \mathbb C^{N_{\mathrm{TX}}\times N_x\times N_y\times N_z},\\
C
    &\in \mathbb C^{N_{\mathrm{RX}}\times N_x\times N_y\times N_z}.
\end{aligned}
$$

This representation is convenient for programmatically modifying inputs,
placing them on a GPU, or making selected tensors trainable. It is less focused
on a self-contained, declarative description of an entire scanner experiment.

### Input-flow summary

```text
JEMRIS
sample.h5 ─┐
TX/RX XML ─┼─> simu.xml ─> Bloch/CVODE executable ─> signal HDF5
sequence  ─┘

MRzeroCore
phantom tensors ─┐
sequence object ─┴─> mr0.util.simulate(...) ─> signal and k-space tensors
```

## Mathematical view

### JEMRIS: direct Bloch-equation integration

JEMRIS evolves the magnetization vector of each simulated spin by numerically
integrating the Bloch equation:

$$
\frac{d\mathbf M(\mathbf r,t)}{dt}
=
\gamma\,\mathbf M(\mathbf r,t)\times\mathbf B_{\mathrm{eff}}(\mathbf r,t)
-
\begin{bmatrix}
M_x/T_2 \\
M_y/T_2 \\
(M_z-M_0)/T_1
\end{bmatrix}.
$$

The effective field can contain static-field error, applied gradients, and RF
fields:

$$
\mathbf B_{\mathrm{eff}}(\mathbf r,t)
=
\mathbf B_0(\mathbf r)
+ \sum_{j\in\{x,y,z\}} g_j(t)\,\mathbf G_j(\mathbf r)
+ \mathbf B_1(\mathbf r,t).
$$

This formulation naturally permits rotations about spatially varying axes.
Its cost is that the ODE must be evaluated over many time steps and spins. A
receive channel approximately measures

$$
S_c(t)
=
\int_V C_c(\mathbf r)
\left[M_x(\mathbf r,t)+iM_y(\mathbf r,t)\right]d\mathbf r,
$$

where $C_c(\mathbf r)$ is the complex sensitivity of receive coil $c$.

### MRzeroCore: phase-distribution state evolution

MRzeroCore represents transverse and longitudinal magnetization using a set of
phase-distribution states rather than integrating every isochromat through
every small time step. At the event level, its evolution can be summarized as

$$
\mathbf s_{n+1}
=
\mathcal D\!\left(\Delta t_n,T_1,T_2,T_2',\Delta f,\Delta\mathbf k_n\right)
\mathcal R\!\left(\alpha_n,\phi_n\right)\mathbf s_n
+ \mathbf b_n.
$$

Here, $\mathcal R$ mixes transverse and longitudinal states for an RF event,
while $\mathcal D$ applies relaxation, off-resonance phase, dephasing, and the
gradient-induced phase-state shift. For example, scalar off-resonance evolves
a transverse component as

$$
F^+(t+\Delta t)
=
F^+(t)\,
e^{-\Delta t/T_2}
e^{-i2\pi\Delta f\,\Delta t}.
$$

Gradient encoding accumulates a spatial phase

$$
\Delta\mathbf k
=
\frac{\gamma}{2\pi}\int_t^{t+\Delta t}\mathbf G(\tau)\,d\tau,
\qquad
M_{xy}(\mathbf r)\mapsto
M_{xy}(\mathbf r)e^{-i2\pi\Delta\mathbf k\cdot\mathbf r}.
$$

The analytical state representation avoids densely sampling all intra-voxel
spin phases. It is therefore fast and deterministic for physics represented by
the PDG operators, but arbitrary vector-field rotations do not reduce to a
single scalar phase factor.

### Why differentiability matters

If sequence parameters are collected in $\boldsymbol\theta$, MRzeroCore can
differentiate a loss through simulation and reconstruction:

$$
\boldsymbol\theta^*
=
\operatorname*{arg\,min}_{\boldsymbol\theta}
\mathcal L\!\left(
\mathcal R_{\mathrm{recon}}
\left[
\mathcal S_{\mathrm{PDG}}(\boldsymbol\theta)
\right],
I_{\mathrm{target}}
\right),
$$

with PyTorch evaluating $\partial\mathcal L/\partial\boldsymbol\theta$ by
automatic differentiation. JEMRIS is principally a forward simulator, so the
same optimization generally requires finite differences, surrogate models, or
an external optimization strategy.

## Choosing between them

Use **MRzeroCore** when iteration speed or gradients with respect to sequence,
phantom, or reconstruction parameters are important. It is especially useful
for optimization and machine-learning workflows.

Use **JEMRIS** when the experiment depends on detailed time-domain Bloch
evolution, explicit coil configurations, nonlinear or concomitant fields, or
other hardware and physical effects that do not fit the standard PDG model.

A useful combined workflow is:

1. Prototype or optimize a sequence in MRzeroCore.
2. Export the sequence as Pulseq.
3. Convert Pulseq to JEMRIS XML with this repository's
   `tools/seq_to_jemris_xml.py` utility.
4. Validate selected cases in JEMRIS using the detailed sample, field, and coil
   configuration.

## Important B0 distinction

Both frameworks can represent a scalar off-resonance value over three spatial
dimensions. In MRzeroCore this is a tensor such as `B0[x,y,z]`, normally in Hz.
That is different from evolving magnetization under an arbitrary vector field
`(Bx, By, Bz)` whose direction changes over space. Direct Bloch integration is
the more natural model for the latter.

For a measured field magnitude or longitudinal component in tesla, the usual
MRzeroCore conversion is

$$
\Delta f(x,y,z)
=
\frac{\gamma}{2\pi}
\left[B_{\parallel}(x,y,z)-B_{\mathrm{reference}}\right].
$$

For hydrogen, $\gamma/(2\pi)$ is approximately
$42.57747892\times10^6\ \mathrm{Hz/T}$.

## References

- [This repository's JEMRIS simulation and Pulseq workflow](README.md)
- [JEMRIS documentation](https://www.jemris.org/ug_intro.html)
- [JEMRIS command-line and XML configuration](https://www.jemris.org/ug_command_line.html)
- [MRzeroCore documentation](https://mrsources.github.io/MRzero-Core/)
- [MRzeroCore source](https://github.com/MRsources/MRzero-Core)
- [Phase distribution graphs for fast, differentiable, and spatially encoded Bloch simulations](https://doi.org/10.1002/mrm.30055)
