# Simulation

JEMRIS reads a simulation-driver XML file, loads its referenced inputs, and
numerically evolves every spin with the Bloch equations.

## Connect the inputs

A minimal `simu.xml` looks like this:

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

All relative paths are resolved from the working directory used to run JEMRIS.

## Physical evolution

For each simulated spin, JEMRIS integrates

$$
\frac{d\mathbf M}{dt}
=
\gamma\mathbf M\times\mathbf B_{\mathrm{eff}}
-
\begin{bmatrix}
M_x/T_2\\
M_y/T_2\\
(M_z-M_0)/T_1
\end{bmatrix}.
$$

The effective field combines the static field, gradient fields, and transmit RF
field:

$$
\mathbf B_{\mathrm{eff}}(\mathbf r,t)
=
\mathbf B_0(\mathbf r)
+\sum_j g_j(t)\mathbf G_j(\mathbf r)
+\mathbf B_1(\mathbf r,t).
$$

At each ADC time, the receive coils weight and sum transverse magnetization to
form the measured channels.

## Run the reference example

From the repository root:

```bash
./run_jemris.sh
```

The script runs the checked-in reference configuration. For another experiment,
invoke the installed JEMRIS executable with that experiment's simulation XML,
following the same directory and relative-path conventions.

## What controls runtime

Runtime generally increases with:

- the number of simulated spins;
- sequence duration and waveform complexity;
- the number of ADC and evolution time points;
- solver accuracy requirements;
- motion or additional physical models; and
- the number and complexity of coil fields.

MPI execution distributes the spin ensemble across processes for larger jobs.

Previous: [prepare inputs](Inputs). Next: [inspect outputs](Outputs).

