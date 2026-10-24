# MRzeroCore Inputs

MRzeroCore assembles an experiment from a phantom object and a sequence object
or Pulseq file. Inputs are normally PyTorch tensors held by Python objects.

## Input sections

1. [Static magnetic field B0](#1-static-magnetic-field-b0)
2. [Gradients](#2-gradients)
3. [Sample, voxels, and spins](#3-sample-voxels-and-spins)
   - [Geometry and voxel positions](#geometry-and-voxel-positions)
   - [Tissue and magnetization properties](#tissue-and-magnetization-properties)
   - [Explicit voxel example](#explicit-voxel-example)
4. [RF, transmit fields, and receive coils](#4-rf-transmit-fields-and-receive-coils)
5. [Sequence and ADC schedule](#5-sequence-and-adc-schedule)
6. [Input checklist](#input-checklist)

## 1. Static magnetic field B0

### Representation

`phantom.B0` is a real tensor with shape `(nx, ny, nz)`. Each value is the
scalar off-resonance frequency for one voxel. A three-dimensional scalar map
`B0(x,y,z)` is supported; it is not a three-component vector field.

### Units and reference field

`phantom.B0` is expressed in Hz relative to the rotating-frame reference. For
an experimental longitudinal field map:

```math
\Delta f(\mathbf r)
=
\frac{\gamma}{2\pi}
\left[B_{\parallel}(\mathbf r)-B_{\mathrm{reference}}\right],
\qquad
\frac{\gamma}{2\pi}
\approx42.57747892\times10^6\ \mathrm{Hz/T}
```

```python
gamma_hz_per_t = 42.57747892e6
phantom.B0 = gamma_hz_per_t * (measured_b_parallel_t - reference_field_t)
```

## 2. Gradients

### Representation

Native `mr0.Sequence` repetitions store gradient moments in `rep.gradm`.
Gradient moments accumulate k-space displacement:

```math
\Delta\mathbf k
=
\frac{\gamma}{2\pi}
\int_t^{t+\Delta t}\mathbf G(\tau)\,d\tau
```

For physical gradients, the k-space increment has units of cycles per meter
when the gyromagnetic ratio divided by two pi is in hertz per tesla and the
gradient is in tesla per meter. The full trajectory is the cumulative sum of
the event increments:

```math
\mathbf k_n
=
\mathbf k_{n-1}+\Delta\mathbf k_n
```

Here, `n` indexes a sequence event or time interval. The increment during that
event and the position after it are, respectively,

```math
\Delta\mathbf k_n
\quad\text{and}\quad
\mathbf k_n.
```

### Normalized and physical gradients

When `sequence.normalized_grads = True`, values in `rep.gradm` are dimensionless
encoding steps rather than fixed physical k-space increments. During
simulation, MRzeroCore reads the physical phantom extent

```math
\mathbf L=
\begin{bmatrix}L_x & L_y & L_z\end{bmatrix}^{\mathsf T}
```

and applies the axis-wise scale

```math
\mathbf q=
\begin{bmatrix}
1/L_x & 1/L_y & 1/L_z
\end{bmatrix}^{\mathsf T}.
```

Let the normalized gradient-moment vector for event `n` be denoted by

```math
\Delta\widetilde{\mathbf k}_n.
```

Its physical k-space increment is

```math
\Delta\mathbf k_n
=
\begin{bmatrix}
1/L_x & 0 & 0\\
0 & 1/L_y & 0\\
0 & 0 & 1/L_z
\end{bmatrix}
\Delta\widetilde{\mathbf k}_n.
```

Consequently, a normalized step of one along an axis is one reciprocal-FOV
step:


```math
\Delta k_x=\frac{1}{\mathrm{FOV}_x},
\qquad
\Delta k_y=\frac{1}{\mathrm{FOV}_y},
\qquad
\Delta k_z=\frac{1}{\mathrm{FOV}_z}
```

For example, if

```math
\mathrm{FOV}_x=0.20\ \mathrm{m},
```

then a normalized x-gradient step of one becomes

```math
\Delta k_x
=
\frac{1}{0.20\ \mathrm m}
=
5\ \mathrm{cycles/m}.
```

For `N_x` Cartesian samples separated by this increment, the nominal spatial
resolution is

```math
\Delta x
=
\frac{1}{N_x\Delta k_x}
=
\frac{\mathrm{FOV}_x}{N_x}.
```

The purpose of normalization is to let the same dimensionless sequence encode
the same matrix layout for phantoms with different physical field-of-view
sizes. Changing the phantom size changes the physical gradient moment and
k-space spacing automatically. An axis whose phantom extent is zero is assigned
a zero scale to avoid division by zero.

When `sequence.normalized_grads = False`, MRzeroCore does not apply reciprocal-
FOV scaling: `rep.gradm` already contains physical SI-based k-space increments.
Sequences imported from Pulseq use this physical convention.

```math
M_{xy}(\mathbf r)
\longmapsto
M_{xy}(\mathbf r)
\exp\!\left(-i2\pi\,\Delta\mathbf k\cdot\mathbf r\right)
```

### Linear-gradient assumption

The standard PDG gradient operator assumes an ideal spatially linear encoding
field:

```math
B_{\mathrm{grad}}(\mathbf r,t)
=
G_x(t)x+G_y(t)y+G_z(t)z
```

This linearity is what turns a gradient event into a global k-space shift in
the phase-distribution state. A nonlinear gradient field would instead produce

```math
\Delta\phi(\mathbf r)
=
\gamma\int
\left[
g_x(t)G_x(\mathbf r)
+g_y(t)G_y(\mathbf r)
+g_z(t)G_z(\mathbf r)
\right]dt
```

which generally cannot be represented by one global k-space shift. The native
`rep.gradm` operator therefore represents ideal linear encoding gradients.

## 3. Sample, voxels, and spins

### Geometry and voxel positions

`VoxelGridPhantom` stores property tensors on an `(nx, ny, nz)` grid. Its
`3 × 4` affine maps voxel indices to physical coordinates. For explicit small
experiments, `CustomVoxelPhantom` accepts a list of positions and voxel sizes.

### Tissue and magnetization properties

| Field | Meaning | Typical shape/unit |
|---|---|---|
| `PD` | Proton density and equilibrium signal weighting | `(nx, ny, nz)` |
| `T1` | Longitudinal relaxation time | `(nx, ny, nz)`, seconds |
| `T2` | Transverse relaxation time | `(nx, ny, nz)`, seconds |
| `T2dash` | Reversible intravoxel dephasing time | `(nx, ny, nz)`, seconds |
| `D` | Diffusion property | `(nx, ny, nz)` |
| `B0` | Scalar off-resonance | `(nx, ny, nz)`, Hz |

The principal scalar maps are therefore:

```math
\mathrm{PD},T_1,T_2,T_2',D,B_0
\in\mathbb R^{N_x\times N_y\times N_z}
```

```math
E_1(\Delta t)=\exp\!\left(-\frac{\Delta t}{T_1}\right),
\qquad
E_2(\Delta t)=\exp\!\left(-\frac{\Delta t}{T_2}\right)
```

Calling `phantom.build()` selects signal-producing voxels—by default voxels
with non-negligible `PD`—and creates compact simulation data. Magnetization is
normally initialized at longitudinal equilibrium from `PD`; lower-level
simulation data can carry custom initial magnetization.

```math
\mathbf M(\mathbf r,0)
=
\begin{bmatrix}0 & 0 & \mathrm{PD}(\mathbf r)\end{bmatrix}^{\mathsf T}
```

### Explicit voxel example

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

### Transmit fields

`phantom.B1` is a complex tensor with shape
`(n_tx, nx, ny, nz)`. It contains relative transmit-field sensitivity. RF
angle, phase, usage, and timing are stored in sequence events.

```math
B_1^+(\mathbf r,t)=\sum_c B_{1,c}(\mathbf r)u_c(t)
```

### Receive fields

`phantom.coil_sens` is a complex tensor with shape
`(n_rx, nx, ny, nz)`. The leading dimension identifies receive channels.

```math
B_1\in\mathbb C^{N_{\mathrm{TX}}\times N_x\times N_y\times N_z},
\qquad
C\in\mathbb C^{N_{\mathrm{RX}}\times N_x\times N_y\times N_z}
```

```math
S_c(t)
=
\sum_v C_{c,v}\left[M_{x,v}(t)+iM_{y,v}(t)\right]
```

## 5. Sequence and ADC schedule

MRzeroCore accepts:

- a native `mr0.Sequence`;
- a PyPulseq sequence object; or
- a Pulseq `.seq` path.

```python
sequence = mr0.Sequence.import_file("sequence.seq")
signal, kspace = mr0.util.simulate(sequence, phantom)
```

Signal samples correspond to ADC-marked events. The returned k-space tensor
contains their trajectory.

## Input checklist

- Confirm `B0` is in Hz relative to the intended carrier.
- Confirm `T1`, `T2`, and `T2dash` are in seconds.
- Check affine orientation, field of view, and voxel dimensions.
- Check normalized versus physical gradient moments.
- Check RF phase and signal sign conventions.
- Check transmit/receive tensor channel order and normalization.
- Check voxel-volume or proton-density normalization.

Next: [MRzero Documentation](MRzero-Documentation).
