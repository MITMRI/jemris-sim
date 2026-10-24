# JEMRIS Inputs

JEMRIS assembles an experiment from a sample, magnetic-field and coil models,
a sequence, and a simulation-driver XML file.

## 1. Static magnetic field B0

### Representation

The standard sample stores one `DB` value per grid point or spin. `DB` supplies
the local static off-resonance used by the Bloch model, including chemical-shift
or susceptibility contributions. Additional field terms may be supplied by the
selected physical model.

### Units and reference field

JEMRIS normally operates in a rotating frame. Convert an experimental field to
an offset relative to the RF reference; do not store a large absolute field in
tesla directly as `DB`.

```math
\Delta f(\mathbf r)
=
\frac{\gamma}{2\pi}
\left[B_{\parallel}(\mathbf r)-B_{\mathrm{reference}}\right]
```

Keep `DB`, sequence time, and the selected model's frequency units consistent.
The checked-in reference sample has `DB = 0`.

## 2. Gradients

### Representation

Gradient pulses and waveforms are elements of the JEMRIS sequence XML tree.
Analytic pulses can be defined directly in XML; sampled gradient waveforms can
be referenced from HDF5 files.

### Spatial fields

The sequence waveform multiplies the configured spatial gradient field. Models
may add nonlinear-gradient and concomitant-field effects.

```math
B_z(\mathbf r,t)
=B_0(\mathbf r)
+g_x(t)G_x(\mathbf r)
+g_y(t)G_y(\mathbf r)
+g_z(t)G_z(\mathbf r)
```

### Pulseq input

Convert Pulseq into JEMRIS XML and referenced waveform files with:

```bash
python tools/seq_to_jemris_xml.py input.seq -o converted -n sequence
```

## 3. Sample, voxels, and spins

### Geometry and spin positions

A standard gridded sample uses:

```text
/sample/data         property maps on the spatial grid
/sample/resolution   voxel spacing in x, y, z
/sample/offset       spatial origin offset in x, y, z
```

JEMRIS generates a coordinate for every grid entry from `resolution` and
`offset`. Entries with `M0 > 0` become live spins; entries with `M0 = 0` do not
contribute signal. An explicit, nongridded spin ensemble can also be used.

### Tissue and magnetization properties

For a standard single-pool sample, the property axis of `/sample/data` is:

| Index | Field | Meaning |
|---:|---|---|
| 0 | `M0` | Equilibrium magnetization and signal-density scale |
| 1 | `R1` | Longitudinal relaxation rate, `1/T1` |
| 2 | `R2` | Transverse relaxation rate, `1/T2` |
| 3 | `R2S` | Effective transverse relaxation rate, `1/T2*` |
| 4 | `DB` | Static off-resonance term |

```math
R_1=\frac{1}{T_1},
\qquad R_2=\frac{1}{T_2},
\qquad R_2^*=\frac{1}{T_2^*}
```

Internally, each live spin carries:

```text
x, y, z, M0, R1, R2, R2S, DB, spin ID
```

Magnetization is normally initialized along the longitudinal direction from
`M0`. Dynamic trajectories and alternate physical models can modify tissue or
magnetization behavior during a simulation.

```math
\mathbf M(\mathbf r,0)
=
\begin{bmatrix}0 & 0 & M_0(\mathbf r)\end{bmatrix}^{\mathsf T}
```

### Multiple pools and diffusion

A five-dimensional HDF5 layout can add a pool/compartment dimension. The
Bloch-McConnell model handles exchanging pools. Diffusion and motion require
the corresponding JEMRIS model/configuration rather than extra columns in the
basic five-property sample.

## 4. RF, transmit fields, and receive coils

### Transmit fields

TX coil-array XML defines transmit channels and their spatial sensitivities.
RF amplitude and phase waveforms come from sequence XML or referenced waveform
files.

```math
B_1^+(\mathbf r,t)=\sum_c C_c^+(\mathbf r)u_c(t)
```

### Receive fields

RX coil-array XML defines receive sensitivities. Each receive coil produces an
output signal channel.

```math
S_c(t)
=
\int_V C_c^-(\mathbf r)
\left[M_x(\mathbf r,t)+iM_y(\mathbf r,t)\right]d\mathbf r
```

The minimal uniform-coil definition is:

```xml
<CoilArray>
  <IdealCoil/>
</CoilArray>
```

## 5. Sequence and ADC schedule

JEMRIS sequence XML is a tree of RF pulses, gradient pulses, delays, loops, and
ADC events. The simulation driver connects every input:

```xml
<simulate name="JEMRIS">
  <sample name="phantom" uri="sample.h5"/>
  <TXcoilarray uri="tx_coils.xml"/>
  <RXcoilarray uri="rx_coils.xml"/>
  <sequence name="sequence" uri="sequence.xml"/>
  <model name="Bloch" type="CVODE"/>
</simulate>
```

## Input checklist

- Confirm `DB` frequency and sequence time units.
- Confirm `R1`, `R2`, and `R2S` are rates, not time constants.
- Check sample resolution, offset, and axis ordering.
- Check RF and gradient waveform units.
- Check transmit and receive sensitivity normalization.
- Run JEMRIS from a directory where relative XML/HDF5 paths resolve.

Next: [Simulation](Simulation) or [Outputs](Outputs).
