# MRzeroCore Simulation

MRzeroCore does not normally integrate the Bloch ODE with tiny numerical time
steps. It builds a Phase Distribution Graph (PDG), applies analytical event
operators to the retained states, evaluates voxel phase/relaxation factors, and
sums the receive-coil-weighted transverse magnetization at ADC events.

## Simulation sections

1. [Prepare the simulation](#1-prepare-the-simulation)
2. [Initialize the simulation state](#2-initialize-the-simulation-state)
3. [Evaluate the sequence and fields](#3-evaluate-the-sequence-and-fields)
4. [Evolve the magnetization](#4-evolve-the-magnetization)
5. [Evaluate magnetization at ADC samples](#5-evaluate-magnetization-at-adc-samples)
6. [Form the net receive signal](#6-form-the-net-receive-signal)
7. [Control accuracy and performance](#7-control-accuracy-and-performance)
8. [Run the simulation](#8-run-the-simulation)

## 1. Prepare the simulation

The high-level interface converts the phantom into compact `SimData` and the
sequence into MRzeroCore events:

```python
signal, kspace = mr0.util.simulate(sequence, phantom)
```

The simulation data contains signal-producing voxel positions, proton density,
relaxation times, diffusion, B0/B1 maps, receive sensitivities, voxel shape,
and dephasing information. Voxels below the proton-density threshold are
discarded before execution.

### Construct the phase graph

`compute_graph()` performs a lightweight pass through the complete sequence
before the voxelwise simulation begins. The purpose of this pass is easier to
understand after changing the object being followed.

#### From spin vectors to phase distributions

An explicit-spin simulator stores a magnetization vector for every spin and
lets different spins acquire different phases. MRzeroCore instead groups
magnetization that has the same spatial phase pattern. Schematically, the
transverse and longitudinal magnetization may be expanded as

```math
M_+(\mathbf r)
=M_x(\mathbf r)+iM_y(\mathbf r)
=\sum_q F_q^+ e^{-i2\pi\mathbf k_q\cdot\mathbf r},
```

```math
M_z(\mathbf r)
=\sum_q Z_q e^{-i2\pi\mathbf k_q\cdot\mathbf r}.
```

The sign of the exponential is a convention; the important idea is that one
coefficient represents an entire phase pattern. A phase-distribution state is
therefore **not a spin or a voxel**. It is a coherent packet of magnetization
with

1. a type: longitudinal Z, positive transverse F⁺, or its
   conjugate/oppositely rotating form F⁻;
2. a phase coordinate **k**, describing how rapidly its phase winds
   across space, plus a time-like coordinate for accumulated off-resonance
   phase; and
3. a complex amplitude, which the main pass later evaluates for every voxel.

For example, an F⁺ state at zero spatial frequency is in phase throughout the
object and can produce a large net signal. A state with a large spatial
frequency winds through many cycles across the object, so the
voxel contributions mostly cancel even though the local transverse
magnetization may be large.

#### What gradients and RF pulses do to a state

A gradient adds a position-dependent phase. In the phase-distribution picture,
this translates a transverse state through phase space:

```math
F^+(\mathbf k)
\xrightarrow{\,\Delta\mathbf k\,}
F^+(\mathbf k+\Delta\mathbf k),
\qquad
\Delta\mathbf k
=\frac{\gamma}{2\pi}\int \mathbf G(t)\,dt.
```

The conjugate transverse component follows the opposite sign convention. Thus,
a gradient need not be applied separately to thousands of spins: MRzeroCore
updates the phase coordinate carried by the state.

An RF pulse does something different. It rotates magnetization locally and
therefore mixes the three state types:

```math
\begin{bmatrix}
F^+\\ F^-\\ Z
\end{bmatrix}_{\!\mathrm{after}}
=
\mathcal R(\alpha,\phi)
\begin{bmatrix}
F^+\\ F^-\\ Z
\end{bmatrix}_{\!\mathrm{before}}.
```

Consequently, one incoming state can have several descendants. Some of a
longitudinal state may remain longitudinal, while some is tipped into
transverse states; some transverse magnetization may remain transverse, return
to Z, or enter the conjugate transverse pathway. The actual edge weights
are elements of the RF rotation matrix and depend on flip angle and RF phase.

#### Why this becomes a graph

The prepass records these possible transfers as a directed graph:

- A **node** is one retained phase-distribution state at a particular point in
  sequence history, including its type and phase coordinate.
- An **edge** says that an RF transformation transfers some complex amplitude
  from a parent state into a child state.
- Several parent edges may end at the same child. Their complex amplitudes are
  added, so constructive and destructive interference are preserved.
- A complete route from the initial equilibrium state Z₀ to a measured
  transverse state is a **coherence pathway**.

Typical edge types include

```math
Z\rightarrow Z,
\qquad
Z\rightarrow F^+,
\qquad
F^+\rightarrow Z,
\qquad
F^-\rightarrow F^+.
```

These arrows do not mean that every arrow is followed as a separate full
simulation. They are ancestry instructions. During the main pass, the value of
a child node is assembled from all of its parents:

```math
m_{\mathrm{child}}(\mathbf r)
=
\sum_{p\in\operatorname{parents}(\mathrm{child})}
c_{p\rightarrow\mathrm{child}}\,
m_p(\mathbf r),
```

after which relaxation, off-resonance, diffusion, motion, and gradient phase
are applied.

#### A spin-echo intuition

MRzeroCore combines the three spatial-frequency coordinates and the time-like
off-resonance coordinate in its four-component `kt_vec`. Denote that general
phase coordinate by **κ**. A useful pathway in a simple spin echo then looks
schematically like

```math
Z_0
\xrightarrow{90^\circ}
F^+(\boldsymbol\kappa=\mathbf 0)
\xrightarrow{\text{dephase}}
F^+(\Delta\boldsymbol\kappa)
\xrightarrow{180^\circ}
F^+(-\Delta\boldsymbol\kappa)
\xrightarrow{\text{continued evolution}}
F^+(\boldsymbol\kappa\approx\mathbf 0).
```

The middle RF operation represents the conjugate pathway: internally, an edge
can tell the main pass to conjugate the parent amplitude and reverse its phase
coordinate. This is why the implementation need not store a completely
separate array of F⁻ states even though F⁻ is useful in the physical notation.

The intermediate transverse state can contribute almost no net signal because
it is strongly dephased. It must nevertheless be retained because a later RF
pulse and gradient evolution can refocus it. The graph is a compact record of
that echo history: *where did this state come from, and can it eventually reach
a signal-producing state?*

#### Why a prepass and pruning are needed

Repeated RF pulses continually branch states. If every possible branch were
retained, the number of coherence pathways would grow combinatorially.
`compute_graph()` therefore runs an inexpensive representative-tissue
simulation to discover the reachable states and estimate two importance
measures:

- **latent signal** estimates whether a state carries magnetization worth
  propagating, including magnetization that is currently longitudinal or
  dephased but may later be refocused;
- **emitted signal** estimates whether a transverse state can contribute
  appreciably at an ADC event.

States below the configured importance thresholds are omitted, and
`max_state_count` bounds the number retained. The resulting graph is a sparse
execution plan. The main pass then places actual voxelwise magnetization on
those nodes, uses the real tissue and field maps, combines every retained
ancestor coherently, and evaluates the emitted states at the ADC times.

This separation is the central computational idea: the prepass answers
**which coherence pathways are worth evaluating?**; the main pass answers
**what complex signal do those pathways produce in this phantom?**

## 2. Initialize the simulation state

The default graph begins with a fully relaxed longitudinal state:

```math
Z_0(\mathbf r,0)=1,
\qquad
F^+(\mathbf r,0)=F^-(\mathbf r,0)=0.
```

Proton density is applied during receive-signal formation. A custom complex
initial-magnetization tensor can replace the relaxed state for steady-state or
specialized experiments.

## 3. Evaluate the sequence and fields

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

## 4. Evolve the magnetization

For every retained graph node, MRzeroCore gathers its simulated ancestors,
applies the RF mixing coefficients, and then propagates the resulting state
through the event trajectory. The event-level update can be summarized as

```math
\mathbf s_{n+1}
=
\mathcal D_n
\mathcal P_n
\mathcal R_n
\mathbf s_n
+\mathbf b_n,
```

where the RF operator mixes states, the phase operator applies gradient,
off-resonance, and motion phase, the decay operator applies relaxation and
diffusion, and the recovery term restores longitudinal equilibrium. These are
analytical tensor operations rather than adaptive ODE steps.

## 5. Evaluate magnetization at ADC samples

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

## 6. Form the net receive signal

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

## 7. Control accuracy and performance

The principal accuracy controls are the minimum latent-signal and
emitted-signal thresholds used to retain and measure PDG states. Smaller
thresholds preserve more pathways but require more memory and computation.

Runtime depends on the number of repetitions, retained graph states, voxels,
ADC samples, and coils. PyTorch vectorizes voxel and coil calculations on CPU
or supported GPU paths. Because event operators and signal formation are tensor
operations, automatic differentiation can propagate a loss through supported
sequence, phantom, and reconstruction parameters.

## 8. Run the simulation

The high-level interface prepares the graph, executes it, and returns the
signal and trajectory:

```python
signal, kspace = mr0.util.simulate(sequence, phantom)
```

Previous: [MRzeroCore Inputs](MRzero-Inputs). See also:
[JEMRIS vs. MRzeroCore](JEMRIS-vs.-MRzeroCore).
