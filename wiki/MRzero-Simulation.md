# MRzeroCore Simulation

MRzeroCore does not normally integrate the Bloch ODE with tiny numerical time
steps. It builds a Phase Distribution Graph (PDG), applies analytical event
operators to the retained states, evaluates voxel phase/relaxation factors, and
sums the receive-coil-weighted transverse magnetization at ADC events.

## Simulation sections

1. [Build simulation data](#1-build-simulation-data)
2. [Construct the phase graph](#2-construct-the-phase-graph)
3. [Initialize magnetization](#3-initialize-magnetization)
4. [Apply RF and free-evolution operators](#4-apply-rf-and-free-evolution-operators)
5. [Evaluate voxel magnetization](#5-evaluate-voxel-magnetization)
6. [Form the net signal](#6-form-the-net-signal)
7. [Accuracy, performance, and differentiation](#7-accuracy-performance-and-differentiation)

## 1. Build simulation data

The high-level interface converts the phantom into compact `SimData` and the
sequence into MRzeroCore events:

```python
signal, kspace = mr0.util.simulate(sequence, phantom)
```

The simulation data contains signal-producing voxel positions, proton density,
relaxation times, diffusion, B0/B1 maps, receive sensitivities, voxel shape,
and dephasing information. Voxels below the proton-density threshold are
discarded before execution.

## 2. Construct the phase graph

`compute_graph()` performs a prepass through the sequence. It determines which
longitudinal and transverse phase-distribution states can be created by RF and
gradient events and records their ancestor relationships.

Conceptually, every retained state is a node, while an RF transformation is an
edge such as

```math
Z\rightarrow Z,
\quad
Z\rightarrow F^+,
\quad
F^+\rightarrow Z,
\quad
F^-\rightarrow F^+.
```

Signal estimates from the prepass allow states with negligible latent or
emitted signal to be pruned. This prevents the number of coherence pathways
from growing without bound.

## 3. Initialize magnetization

The default graph begins with a fully relaxed longitudinal state:

```math
Z_0(\mathbf r,0)=1,
\qquad
F^+(\mathbf r,0)=F^-(\mathbf r,0)=0.
```

Proton density is applied during receive-signal formation. A custom complex
initial-magnetization tensor can replace the relaxed state for steady-state or
specialized experiments.

## 4. Apply RF and free-evolution operators

### RF mixing

For each repetition, transmit sensitivity modifies the nominal flip angle and
phase. The RF operator analytically mixes longitudinal and transverse ancestor
states. In abstract form,

```math
\mathbf s^+
=
\mathcal R(\alpha B_1,\phi+\arg B_1)\mathbf s^-.
```

The implementation evaluates trigonometric mixing coefficients such as

```math
\cos\alpha,
\qquad
\cos^2\!\left(\frac{\alpha}{2}\right),
\qquad
\sin\alpha\,e^{i\phi},
```

and combines all simulated ancestors of each graph state.

### Gradient and time trajectory

Within a repetition, MRzeroCore cumulatively sums gradient moments and event
durations:

```math
\mathbf k_n
=
\mathbf k_{n-1}+\Delta\mathbf k_n,
\qquad
t_n=t_{n-1}+\Delta t_n.
```

Normalized gradients are first divided by the phantom extent along each axis;
physical Pulseq gradients already use the physical convention.

### Relaxation and off-resonance

Instead of numerically stepping through free precession, MRzeroCore evaluates
closed-form factors at the required event/ADC times:

```math
E_1(\Delta t)=\exp\!\left(-\frac{\Delta t}{T_1}\right),
\qquad
E_2(\Delta t)=\exp\!\left(-\frac{\Delta t}{T_2}\right),
```

```math
P(\Delta t,\Delta f)
=
\exp\!\left(-i2\pi\Delta f\,\Delta t\right).
```

The longitudinal equilibrium state recovers according to

```math
Z_0(t+\Delta t)
=
E_1Z_0(t)+\left(1-E_1\right).
```

MRzeroCore also applies T2-prime dephasing, voxel-shape dephasing, diffusion
attenuation, and motion phase using their corresponding analytical factors.

## 5. Evaluate voxel magnetization

For a transverse PDG state `q` at ADC time `t_a`, the voxelwise transverse
magnetization has the general form

```math
m_{q,v}(t_a)
=
A_{q,v}
E_{2,v}(t_a)
E_{2,v}'(t_a)
D_{q,v}(t_a)
V_q(t_a)
\exp\!\left[
-i2\pi
\left(
\mathbf k_q(t_a)\cdot\mathbf r_v
+t_q(t_a)\Delta f_v
+\phi_{\mathrm{motion},v}(t_a)
\right)
\right].
```

Here, `A` is the RF-mixed state amplitude, `D` is diffusion attenuation, and
`V` is the intravoxel shape/dephasing factor. The time-like PDG coordinate
tracks off-resonance dephasing, while the spatial coordinates track gradient
encoding.

## 6. Form the net signal

For every ADC sample and receive coil, MRzeroCore sums all emitted transverse
states and all signal-producing voxels:

```math
S_c(t_a)
=
e^{i\phi_{\mathrm{ADC}}(t_a)}
\sum_{q\in\mathcal E(t_a)}
\sum_{v=1}^{N_{\mathrm{voxel}}}
m_{q,v}(t_a)\,
\mathrm{PD}_v\,
C_{c,v}.
```

In the implementation, proton density is multiplied into the coil-sensitivity
matrix, and each state's voxel magnetization is evaluated with a tensor matrix
product:

```math
\underbrace{\mathbf m_q}_{N_{\mathrm{ADC}}\times N_{\mathrm{voxel}}}
\underbrace{\mathbf C}_{N_{\mathrm{voxel}}\times N_{\mathrm{coil}}}
=
\underbrace{\mathbf S_q}_{N_{\mathrm{ADC}}\times N_{\mathrm{coil}}}.
```

State contributions are added to obtain the final complex signal tensor. The
ADC phase is applied uniformly to the receive channels, and the corresponding
k-space coordinates are returned with the signal.

## 7. Accuracy, performance, and differentiation

The principal accuracy controls are the minimum latent-signal and
emitted-signal thresholds used to retain and measure PDG states. Smaller
thresholds preserve more pathways but require more memory and computation.

Runtime depends on the number of repetitions, retained graph states, voxels,
ADC samples, and coils. PyTorch vectorizes voxel and coil calculations on CPU
or supported GPU paths. Because event operators and signal formation are tensor
operations, automatic differentiation can propagate a loss through supported
sequence, phantom, and reconstruction parameters.

Previous: [MRzeroCore Inputs](MRzero-Inputs). See also:
[JEMRIS vs. MRzeroCore](JEMRIS-vs.-MRzeroCore).

