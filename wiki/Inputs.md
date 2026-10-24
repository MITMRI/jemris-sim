# Inputs

A JEMRIS experiment combines four main inputs: a sample, a pulse sequence,
transmit coils, and receive coils. A simulation-driver XML file identifies the
files used in one run.

```text
sample.h5 ─────── tissue, spin positions, relaxation, off-resonance
sequence.xml ─── RF, gradients, timing, and ADC schedule
tx_coils.xml ─── transmit sensitivities
rx_coils.xml ─── receive sensitivities and output channels
```

## Sample and tissue

The sample HDF5 describes the simulated spins and their physical properties.
The reference example uses `sample.h5` with spin locations and tissue values
such as $M_0$, $R_1$, $R_2$, $R_2^*$, and `DB`.

| Quantity | Meaning | Main signal effect |
|---|---|---|
| $M_0$ | Equilibrium magnetization / proton-density scale | Signal amplitude |
| $R_1=1/T_1$ | Longitudinal relaxation rate | Recovery between excitations |
| $R_2=1/T_2$ | Transverse relaxation rate | Transverse signal decay |
| $R_2^*$ | Effective transverse decay rate | Additional dephasing |
| `DB` | Per-spin static off-resonance | Phase accrual and signal cancellation |
| Position | Spin coordinate | Gradient encoding and coil weighting |

For a longitudinal field deviation $\Delta B_\parallel$ in tesla, the
corresponding proton off-resonance is

$$
\Delta f(\mathbf r)
=
\frac{\gamma}{2\pi}\Delta B_\parallel(\mathbf r),
\qquad
\frac{\gamma}{2\pi}\approx42.57747892\times10^6\ \mathrm{Hz/T}.
$$

## Sequence

JEMRIS sequences are XML trees containing RF pulses, gradient events, delays,
loops, and ADC events. The repository includes hand-written reference sequences
under `share/examples/approved/`.

Pulseq `.seq` files can be converted into JEMRIS XML:

```bash
python tools/seq_to_jemris_xml.py input.seq -o converted -n sequence
```

The converter writes `sequence.xml` and any referenced RF or gradient waveform
HDF5 files. Run JEMRIS from the directory containing these files so relative
paths resolve correctly.

## Transmit and receive coils

TX coil definitions determine the spatial RF excitation field $B_1^+$. RX
coils determine the spatial receive weighting $B_1^-$ and the number of output
signal channels. The minimal ideal-coil definition is

```xml
<CoilArray>
  <IdealCoil/>
</CoilArray>
```

## Reference inputs

The checked-in EPI example starts with:

- `share/examples/approved/simu.xml`
- `share/examples/approved/sample.h5`
- `share/examples/approved/epi.xml`
- `share/examples/approved/uniform.xml`

Next: [run the simulation](Simulation).

