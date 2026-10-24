# jemris-sim Wiki

`jemris-sim` is the MIT MRI group's JEMRIS-based MRI physics simulation
environment. It connects field, tissue, coil, and pulse-sequence inputs to a
Bloch/CVODE simulation and produces receive-coil signal as a function of time.

## Documentation sections

### 1. [JEMRIS Documentation](JEMRIS-Documentation)

Set up complete JEMRIS experiments using sample HDF5, sequence XML, transmit
and receive coils, and a simulation driver. This section contains dedicated
pages for [Inputs](Inputs), [Simulation](Simulation), and [Outputs](Outputs).

### 2. [MRzero Documentation](MRzero-Documentation)

Build differentiable MRI simulations in Python with MRzeroCore, PyTorch phantom
tensors, native sequences, and Pulseq `.seq` files.

### 3. [JEMRIS vs. MRzeroCore](JEMRIS-vs.-MRzeroCore)

Compare their numerical methods, equations, data formats, field models,
performance characteristics, and recommended use cases.

## Simulation pipeline

```text
sample.h5 ─┐
TX/RX XML ─┼─> simu.xml ─> Bloch/CVODE ─> signal HDF5 / ISMRMRD
sequence  ─┘
```

The simulation driver `simu.xml` connects the sample, transmit and receive
coils, sequence, and physical model. See the checked-in reference files under
`share/examples/approved/` and their outputs under `examples/outputs/`.

## Related projects

- [Upstream JEMRIS](https://github.com/JEMRIS/jemris)
- [JEMRIS documentation](https://www.jemris.org/)
- [MRzeroCore](https://github.com/MRsources/MRzero-Core)
