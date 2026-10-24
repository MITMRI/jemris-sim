# Inputs

This page describes JEMRIS and MRzeroCore inputs using the same organization:
static field, gradients, sample properties, RF/coil fields, and sequence events.

## Unified overview

| Input category | JEMRIS | MRzeroCore |
|---|---|---|
| Static field | Per-spin `DB` in `sample.h5`; optional model-specific field terms | Scalar `phantom.B0[nx, ny, nz]` tensor in Hz |
| Gradients | Waveforms in sequence XML multiplied by configured spatial gradient fields | Gradient moments in `mr0.Sequence`, or physical waveforms imported from Pulseq |
| Sample geometry | Cartesian HDF5 grid or explicit spin ensemble; coordinates come from `resolution` and `offset` | Voxel tensors plus a `3 × 4` affine, or explicit `CustomVoxelPhantom` positions |
| Tissue properties | `M0`, `R1`, `R2`, `R2*`, and `DB` per grid point/spin | `PD`, `T1`, `T2`, `T2dash`, `D`, and `B0` per voxel |
| Transmit field | TX coil-array XML and sensitivity data | Complex `B1[n_tx, nx, ny, nz]` tensor |
| Receive field | RX coil-array XML; one output channel per receive coil | Complex `coil_sens[n_rx, nx, ny, nz]` tensor |
| Sequence | JEMRIS sequence-tree XML; waveform HDF5 may be referenced | Native `mr0.Sequence`, PyPulseq object, or Pulseq `.seq` file |
| Simulation wiring | `simu.xml` connects sample, coils, sequence, and model | Python passes sequence and phantom to `mr0.util.simulate()` |

## 1. Static magnetic field B0

Both simulators normally work in a rotating frame. The useful input is the
local off-resonance relative to the RF carrier, rather than a large absolute
field such as 3 T. For a measured longitudinal field map,

```math
\Delta f(\mathbf r)
=
\frac{\gamma}{2\pi}
\left[B_{\parallel}(\mathbf r)-B_{\mathrm{reference}}\right],
\qquad
\frac{\gamma}{2\pi}
\approx 42.57747892\times10^6\ \mathrm{Hz/T}
```

for hydrogen.

### JEMRIS representation

- The standard sample contains one `DB` value per grid point or spin.
- `DB` supplies the local static off-resonance contribution used by the Bloch
  model, including effects such as chemical shift or susceptibility.
- The reference sample has `DB = 0`, representing a flat off-resonance map.
- Additional field contributions can be supplied by the selected physical
  model.

Keep the numerical convention consistent with the JEMRIS model and sequence
time units. Do not put an absolute field in tesla directly into `DB` without
converting it to the model's off-resonance convention.

### MRzeroCore representation

- `phantom.B0` is a real tensor with shape `(nx, ny, nz)`.
- Values are scalar off-resonance frequencies in Hz.
- A uniform absolute field is implicit in the rotating-frame carrier.
- A 3D scalar map `B0(x,y,z)` is supported; this is not an arbitrary vector
  field `(Bx, By, Bz)` with spatially varying direction.

```python
gamma_hz_per_t = 42.57747892e6
phantom.B0 = gamma_hz_per_t * (measured_b_parallel_t - reference_field_t)
```

## 2. Gradient representation

The encoding field can be written conceptually as

```math
B_z(\mathbf r,t)
=
B_0(\mathbf r)
+ g_x(t)G_x(\mathbf r)
+ g_y(t)G_y(\mathbf r)
+ g_z(t)G_z(\mathbf r)
```

where `g_j(t)` is a sequence waveform and `G_j(r)` is its spatial field.

### JEMRIS representation

- Gradient pulses and waveforms are sequence XML elements.
- Analytic pulses can be defined directly in the XML tree.
- External waveforms can be referenced from HDF5 files.
- Nonlinear-gradient and concomitant-field effects can be included when the
  corresponding models are configured.
- Convert Pulseq to XML and referenced waveforms with:

```bash
python tools/seq_to_jemris_xml.py input.seq -o converted -n sequence
```

### MRzeroCore representation

- Native repetitions contain gradient moments in `rep.gradm`.
- With normalized gradients, moment steps describe encoding and are scaled to
  the phantom field of view during simulation.
- With `sequence.normalized_grads = False`, physical SI gradient moments are
  used instead.
- Pulseq input supplies physical waveforms that MRzeroCore converts into
  simulation events and k-space evolution.

Gradient moments determine the accumulated k-space displacement:

```math
\Delta\mathbf k
=
\frac{\gamma}{2\pi}
\int_t^{t+\Delta t}\mathbf G(\tau)\,d\tau
```

## 3. Sample, voxels, and spins

The sample is distinct from the sequence and field waveforms. It defines which
spins exist, where they are located, equilibrium magnetization, relaxation, and
dephasing properties.

### Shared physical properties

| Physical property | Meaning | JEMRIS | MRzeroCore |
|---|---|---|---|
| Position | Spatial location of a spin or voxel | HDF5 grid `resolution` and `offset`, or explicit ensemble | Voxel affine, or `CustomVoxelPhantom.pos` |
| Equilibrium magnetization | Available magnetization before excitation | `M0` | `PD` proton-density weighting |
| Longitudinal relaxation | Recovery toward equilibrium | `R1 = 1/T1` | `T1` in seconds |
| Transverse relaxation | Irreversible transverse decay | `R2 = 1/T2` | `T2` in seconds |
| Reversible dephasing | Additional intravoxel/static dephasing | `R2* = 1/T2*` plus `DB` | `T2dash` and `B0` |
| Diffusion | Motion-driven attenuation | Model-dependent | `D` map |
| Initial magnetization | State at simulation start | Normally initialized from `M0`; trajectories/models may alter it | Normally equilibrium longitudinal magnetization from `PD`; custom initial states are possible in simulation data |
| Multiple compartments | Several tissue pools | Multi-pool sample dimension and Bloch-McConnell model | Support depends on the selected phantom/model workflow |

### JEMRIS sample format

A gridded sample is stored in HDF5:

```text
/sample/data         property maps on the spatial grid
/sample/resolution   voxel spacing in x, y, z
/sample/offset       spatial origin offset in x, y, z
```

For the standard single-pool sample, the property axis of `/sample/data` is:

| Property index | Field | Description |
|---:|---|---|
| 0 | `M0` | Equilibrium magnetization |
| 1 | `R1` | Longitudinal relaxation rate, `1/T1` |
| 2 | `R2` | Transverse relaxation rate, `1/T2` |
| 3 | `R2S` | Effective transverse relaxation rate, `1/T2*` |
| 4 | `DB` | Static off-resonance term |

JEMRIS generates a coordinate for every grid entry from `resolution` and
`offset`. Entries with `M0 > 0` become live spins; zero-`M0` entries do not
contribute signal. Internally, each live spin carries:

```text
x, y, z, M0, R1, R2, R2S, DB, spin ID
```

A five-dimensional HDF5 layout can add a pool/compartment dimension for
multi-pool simulations.

### MRzeroCore sample format

`VoxelGridPhantom` stores one tensor per property:

```math
\mathrm{PD},T_1,T_2,T_2',D,B_0
\in \mathbb R^{N_x\times N_y\times N_z}
```

and coil fields as

```math
B_1
\in\mathbb C^{N_{\mathrm{TX}}\times N_x\times N_y\times N_z},
\qquad
C
\in\mathbb C^{N_{\mathrm{RX}}\times N_x\times N_y\times N_z}
```

The `3 × 4` affine maps voxel indices to physical coordinates. Calling
`phantom.build()` selects signal-producing voxels (by default, non-negligible
`PD`) and creates compact simulation data for the PDG engine.

For a small explicit experiment, `CustomVoxelPhantom` accepts positions and
property lists directly:

```python
phantom = mr0.CustomVoxelPhantom(
    pos=[[-0.01, 0.0, 0.0], [0.01, 0.0, 0.0]],
    PD=[1.0, 0.7],
    T1=[0.9, 0.3],
    T2=[0.08, 0.07],
    T2dash=[0.03, 0.05],
    D=[0.0, 0.0],
    B0=[0.0, -120.0],
    voxel_size=0.005,
)
```

## 4. RF, transmit fields, and receive coils

### JEMRIS representation

- TX coil-array XML defines transmit channels and spatial sensitivities.
- RX coil-array XML defines receive sensitivities and output channels.
- RF amplitude/phase waveforms come from sequence XML or referenced HDF5.
- The minimal `<IdealCoil/>` has uniform sensitivity.

```xml
<CoilArray>
  <IdealCoil/>
</CoilArray>
```

### MRzeroCore representation

- `phantom.B1` is the complex transmit sensitivity/relative field tensor.
- `phantom.coil_sens` is the complex receive-sensitivity tensor.
- RF angle, phase, usage, and timing are sequence-event properties.
- Leading coil dimensions represent multiple channels.

## 5. Sequence and ADC schedule

### JEMRIS representation

JEMRIS sequence XML is a tree of RF pulses, gradients, delays, loops, and ADC
events. A separate `simu.xml` connects the sequence, sample, and coil files:

```xml
<simulate name="JEMRIS">
  <sample name="phantom" uri="sample.h5"/>
  <TXcoilarray uri="tx_coils.xml"/>
  <RXcoilarray uri="rx_coils.xml"/>
  <sequence name="sequence" uri="sequence.xml"/>
  <model name="Bloch" type="CVODE"/>
</simulate>
```

### MRzeroCore representation

MRzeroCore accepts a native sequence, PyPulseq object, or `.seq` path:

```python
sequence = mr0.Sequence.import_file("sequence.seq")
signal, kspace = mr0.util.simulate(sequence, phantom)
```

Signal samples correspond to ADC-marked events; the k-space tensor contains
the corresponding trajectory.

## Unit and consistency checklist

Before comparing the two simulators, verify:

- the same rotating-frame reference and off-resonance sign convention;
- the field/frequency unit expected by JEMRIS `DB` versus MRzeroCore Hz;
- spatial coordinates and field-of-view units;
- `T1`/`T2` time units versus JEMRIS `R1`/`R2` rate units;
- gradient waveform and moment units;
- RF phase and signal sign conventions;
- voxel volume or spin-density normalization; and
- TX/RX sensitivity normalization.

Next: [JEMRIS simulation](Simulation), [MRzero Documentation](MRzero-Documentation),
or [JEMRIS vs. MRzeroCore](JEMRIS-vs.-MRzeroCore).
