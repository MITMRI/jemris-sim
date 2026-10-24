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

```text
delta_f(x,y,z) = gamma_hz_per_t * (B_parallel(x,y,z) - B_reference)
```

where `gamma_hz_per_t` is approximately `42.57747892e6` Hz/T for hydrogen.

## References

- [This repository's JEMRIS simulation and Pulseq workflow](README.md)
- [JEMRIS documentation](https://www.jemris.org/ug_intro.html)
- [JEMRIS command-line and XML configuration](https://www.jemris.org/ug_command_line.html)
- [MRzeroCore documentation](https://mrsources.github.io/MRzero-Core/)
- [MRzeroCore source](https://github.com/MRsources/MRzero-Core)
- [Phase distribution graphs for fast, differentiable, and spatially encoded Bloch simulations](https://doi.org/10.1002/mrm.30055)

