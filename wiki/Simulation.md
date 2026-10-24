# JEMRIS Simulation

JEMRIS numerically integrates the Bloch equations for every live spin and then
sums the receive-coil-weighted magnetization at each ADC time.

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

The simulation driver selects the sample, coils, sequence, and physical model:

```xml
<simulate name="JEMRIS">
  <sample name="phantom" uri="sample.h5"/>
  <RXcoilarray uri="uniform.xml"/>
  <TXcoilarray uri="uniform.xml"/>
  <parameter PositionRandomness="0"/>
  <sequence name="epi" uri="epi.xml"/>
  <model name="Bloch" type="CVODE"/>
</simulate>
```

JEMRIS expands the sequence tree into atomic modules and time points of
interest. These include waveform boundaries, ADC samples, and other times at
which the solver must stop or expose magnetization.

## 2. Initialize the simulation state

For every sample entry with nonzero equilibrium magnetization, JEMRIS loads its
position and tissue properties and initializes longitudinal equilibrium:

```math
\mathbf M_j(0)
=
\begin{bmatrix}
0 & 0 & M_{0,j}
\end{bmatrix}^{\mathsf T}.
```

The sequence tree is then traversed independently for spin `j`. With multiple
compartments, each spin carries one magnetization state per compartment.

## 3. Evaluate the sequence and fields

Whenever the ODE solver requests a derivative at time `t`, JEMRIS evaluates the
current sequence waveforms at the spin's current position. In the rotating
frame, the implementation constructs

```math
B_x(t)=B_1(t)\cos\phi_{\mathrm{RF}}(t),
\qquad
B_y(t)=B_1(t)\sin\phi_{\mathrm{RF}}(t),
```

and

```math
B_z(\mathbf r,t)
=
xG_x(t)+yG_y(t)+zG_z(t)
+\Delta B(\mathbf r,t)
+B_{\mathrm{concomitant}}(\mathbf r,t)
+B_{\mathrm{nonlinear}}(\mathbf r,t).
```

Dynamic-variable modules may first update position, flow activation, motion,
respiration, diffusion, off-resonance, relaxation rates, or equilibrium
magnetization. Thus the right-hand side is evaluated using the state and field
values at the solver's requested time, not merely at ADC times.

## 4. Evolve the magnetization

For a single pool, the physical equation is

```math
\frac{d\mathbf M}{dt}
=
\gamma\,\mathbf M\times\mathbf B_{\mathrm{eff}}
-
\begin{bmatrix}
M_x/T_2\\
M_y/T_2\\
(M_z-M_0)/T_1
\end{bmatrix}.
```

Equivalently, using relaxation rates,

```math
\begin{aligned}
\dot M_x &= \gamma(M_yB_z-M_zB_y)-R_2M_x,\\
\dot M_y &= \gamma(M_zB_x-M_xB_z)-R_2M_y,\\
\dot M_z &= \gamma(M_xB_y-M_yB_x)+R_1(M_0-M_z).
\end{aligned}
```

The current implementation stores transverse magnetization in cylindrical
form—magnitude, phase, and longitudinal component—but evaluates the Bloch
derivative through the equivalent Cartesian components.

### CVODE integration

JEMRIS uses SUNDIALS CVODE with the Adams method, an adaptive multistep method
suited to nonstiff ODEs. For each spin and atomic sequence interval, JEMRIS:

1. reinitializes CVODE from the current magnetization;
2. provides the Bloch derivative callback;
3. sets the next sequence boundary or time point as a stop time;
4. lets CVODE select internal time steps to satisfy relative and componentwise
   absolute error tolerances; and
5. retrieves magnetization at the requested stop time.

The numerical approximation advances from

```math
\mathbf M(t_n)
\longrightarrow
\mathbf M(t_{n+1})
```

using as many internal derivative evaluations as required by the tolerance.
The solver is therefore not restricted to one fixed step per RF, gradient, or
ADC event. If an integration attempt fails, JEMRIS can repeat the current atom
with tighter accuracy.

For multi-pool exchange, the Bloch-McConnell CVODE model integrates the coupled
magnetization system instead.

## 5. Evaluate magnetization at ADC samples

At an ADC time, the solver has produced a magnetization state for spin `j`:

```math
M_{xy,j}(t)=M_{x,j}(t)+iM_{y,j}(t).
```

For receive coil `c`, JEMRIS evaluates the coil sensitivity magnitude and phase
at the current spin position. Including receiver phase, the spin contributes

```math
\Delta S_{c,j}(t_a)
=
C_c^-(\mathbf r_j(t_a))
M_{xy,j}(t_a)
\exp\!\left(-i\phi_{\mathrm{ADC}}(t_a)\right).
```

The implementation accumulates real and imaginary transverse components into
the selected ADC row for every coil and compartment. It can also retain a
longitudinal component for diagnostic output.

## 6. Form the net receive signal

The net channel signal is the discrete sum over all live spins and, when
present, tissue compartments:

```math
S_c(t_a)
=
\sum_{j=1}^{N_{\mathrm{spin}}}
\sum_{p=1}^{N_{\mathrm{pool}}}
C_c^-(\mathbf r_j(t_a))
M_{xy,j,p}(t_a)
\exp\!\left(-i\phi_{\mathrm{ADC}}(t_a)\right).
```

This is the discrete counterpart of

```math
S_c(t)
=
\int_V C_c^-(\mathbf r)M_{xy}(\mathbf r,t)\,d\mathbf r.
```

Optional output normalization divides accumulated samples by the total spin
count. Optional Gaussian noise is added after accumulation with dwell-time
dependent scaling. The resulting arrays are written under
`/signal/channels/NN`, with ADC timestamps in `/signal/times`; ISMRMRD export
places the complex transverse samples into acquisition channels.

## 7. Control accuracy and performance

Runtime generally increases with spin count, sequence duration and waveform
complexity, solver tolerances, ADC/evolution time points, motion/physical
models, compartments, and coils. MPI partitions the spin ensemble across
processes; partial coil signals are combined to obtain the global result.

Accuracy is primarily controlled by CVODE's relative and componentwise
absolute tolerances. Tighter tolerances force smaller or more numerous internal
steps. The number of output ADC samples does not by itself set the ODE step
size: CVODE chooses additional internal evaluations as required.

## 8. Run the simulation

Run the checked-in example from the repository root:

```bash
./run_jemris.sh
```

Previous: [JEMRIS Inputs](Inputs). Next: [Outputs](Outputs).
