# jemris-sim Wiki

`jemris-sim` is the MIT MRI group's JEMRIS-based MRI physics simulation
environment. It connects field, tissue, coil, and pulse-sequence inputs to a
Bloch/CVODE simulation and produces receive-coil signal as a function of time.

## Start here

- [Inputs](Inputs) — samples, fields, coils, and sequence files.
- [Simulation](Simulation) — connect inputs in `simu.xml` and run JEMRIS.
- [Outputs](Outputs) — native signal HDF5, ISMRMRD, and sequence diagnostics.
- [Repository README](https://github.com/MITMRI/jemris-sim#readme) — installation,
  the reference EPI example, output layout, and Pulseq conversion.
- [JEMRIS vs. MRzeroCore](JEMRIS-vs.-MRzeroCore) — methods, equations, input
  formats, and guidance on choosing a simulator.
- [Pulseq to JEMRIS](https://github.com/MITMRI/jemris-sim#pulseq--jemris-xml-sequence-conversion)
  — convert a Pulseq `.seq` file into JEMRIS XML and waveform files.

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
