# JEMRIS vs. MRzeroCore

JEMRIS and MRzeroCore solve related MRI simulation problems with different
numerical representations. JEMRIS numerically integrates the Bloch equations
for spin ensembles. MRzeroCore uses Phase Distribution Graphs (PDG) to provide
fast, differentiable sequence simulations.

The maintained, detailed comparison includes:

- simulation methods and mathematical formulas;
- sequence, phantom, field, and coil input formats;
- JEMRIS XML/HDF5 and MRzeroCore Python/tensor examples;
- static-field and full-vector-field distinctions;
- performance, differentiation, reconstruction, and ML considerations; and
- a combined MRzeroCore-to-JEMRIS validation workflow.

**[Read the complete comparison in the repository](https://github.com/MITMRI/jemris-sim/blob/master/JEMRIS_VS_MRZERO.md).**

## At a glance

| Aspect | JEMRIS | MRzeroCore |
|---|---|---|
| Method | Direct Bloch ODE integration | Phase Distribution Graphs |
| Input style | XML and HDF5 experiment bundle | Python objects and PyTorch tensors |
| Sequence input | JEMRIS XML; Pulseq conversion is provided | Native object, PyPulseq object, or `.seq` file |
| Static-field input | Per-spin `DB` and configured field models | Scalar `B0[nx,ny,nz]` tensor in Hz |
| Differentiable | Not inherently | Yes, through PyTorch |
| Best fit | Detailed forward simulation and validation | Rapid simulation and gradient-based optimization |

