# jemris-sim

MRI physics simulation based on [JEMRIS](https://github.com/JEMRIS/jemris), set up for the MIT MRI group.

Pipeline in one line:

**field / tissue / sequence inputs → Bloch (CVODE) simulation → coil signal vs time**

---

## How inputs become the signal

At every time step, each spin feels an effective field roughly:

$$B_x, B_y \leftarrow \text{RF waveform} \times B_1^+(\mathbf{r})$$

$$B_z \leftarrow B_0(\mathbf{r}) + G_x(\mathbf{r})\,g_x(t) + G_y(\mathbf{r})\,g_y(t) + G_z(\mathbf{r})\,g_z(t) + \cdots$$

Magnetization $\mathbf{M}(\mathbf{r},t)$ evolves with the Bloch equations. At each ADC time, coils form the receive signal from transverse magnetization weighted by $B_1^-(\mathbf{r})$:

$$S(t) = \sum_{\text{spins}} B_1^-(\mathbf{r})\, M_{xy}(\mathbf{r},t)\, e^{i\phi}$$

So **every input below exists to shape $B$, $\mathbf{M}$, or the sampling times that produce $S(t)$.**

> Today’s checked-in example uses a **uniform** $B_0$ / linear gradients / ideal coil. The table also lists the **general** inputs we want (B0 map, G maps, …) and how they would correlate to the same signal outputs.

---

## Input ↔ simulation ↔ output (correlated)

| You provide | What it controls in the sim | What you see in the output |
|-------------|-----------------------------|----------------------------|
| **Tissue / sample** (`sample.h5`) | Spin positions; $M_0$, $R_1$, $R_2$, $R_2^*$; per-spin off-resonance **DB** (stand-in for a $B_0$ error map) | Overall signal scale, decay / $T_2$ effects, dephasing from ΔB |
| **$B_0(x,y,z)$ map** *(desired / general)* | Spatial main-field / off-resonance → enters $B_z$ | Phase accrual, ghosting-like dephasing, peak $\|S\|$ drops if spins disagree |
| **$G_x,G_y,G_z(x,y,z)$ maps** *(desired / general)* | Spatial encoding fields (architecture-agnostic); with waveforms $g_i(t)$ set $B_z$ | k-space trajectory / encoding; distorted maps → distorted $S(t)$ / images |
| **Gradient + RF waveforms** (sequence XML) | Time courses: RF amp/phase, $g_x(t),g_y(t),g_z(t)$, ADC schedule | *When* samples appear in `/signal/times`; encoding structure of $S(t)$ |
| **TX coil / $B_1^+$ map** | Flip-angle / excitation pattern in space | How much $M_{xy}$ is created → amplitude of $S(t)$ |
| **RX coil / $B_1^-$ map** | Spatial weight on each spin’s contribution | Channel signals in `/signal/channels/NN`; multi-coil if multiple maps |
| **Sim driver** (`simu.xml`) | Wires sample + coils + sequence + Bloch/CVODE | One run → one (or more) signal files |

**Outputs (always correlated back to those inputs):**

| Output | Contents | Tied to |
|--------|----------|---------|
| `/signal/times` | ADC timestamps (ms) | Sequence ADC schedule |
| `/signal/channels/00` (…`NN`) | $(S_x, S_y, S_z)$ vs time | RX maps × evolved $M$; length = number of ADCs |
| `*_ismrmrd.h5` | Same acquisition, ISMRMRD layout | Sequence headers + signal (for recon) |
| Sequence dump (`*_seq` / `seqdiag/*`) | Waveforms TX/GX/GY/GZ vs time | Sequence XML only (no tissue) |

$S(t)$ is **magnetization × sensitivity** (sim units), not SI probe volts—unless you apply a calibration factor later.

---

## Diagram

```
 GENERAL INPUTS                         TODAY'S EXAMPLE FILES
 ────────────────                       ─────────────────────
 B0(x,y,z)          ─┐                  (uniform / DB=0 in sample.h5)
 Gx,Gy,Gz(x,y,z)    ─┼─►  B(r,t)       
 gx,gy,gz(t) + RF   ─┤                  epi.xml  (waveforms + ADCs)
 B1+(r), B1-(r)     ─┤                  uniform.xml  (IdealCoil)
 tissue M0,R1,R2,…  ─┘                  sample.h5

              │
              ▼
        Bloch + CVODE  (jemris)
              │
              ▼
     S(t) per RX channel
              │
              ├─► *.h5          /signal/times , /signal/channels/*
              └─► *_ismrmrd.h5  recon-ready copy
```

---

## Reference example (what is checked in)

All specs under `share/examples/approved/`. Outputs under `examples/outputs/`.

### Driver — `simu.xml`

```xml
<?xml version="1.0" encoding="utf-8"?>
<simulate name="JEMRIS">
   <sample name="sphere" uri="sample.h5"/>
   <RXcoilarray uri="uniform.xml"/>
   <TXcoilarray uri="uniform.xml"/>
   <parameter PositionRandomness="0"/>
   <sequence name="epi" uri="epi.xml"/>
   <model name="Bloch" type="CVODE"/>
</simulate>
```

| This line | Feeds into signal how |
|-----------|------------------------|
| `sample.h5` | Which spins exist and their relaxation / ΔB → shape of $S(t)$ |
| `uniform.xml` (RX) | Weight = 1 for every spin → single channel `/signal/channels/00` |
| `uniform.xml` (TX) | Ideal excitation |
| `epi.xml` | 4096 ADC times and EPI encoding → length and structure of $S(t)$ |
| Bloch / CVODE | Physics engine that turns the above into $M$ then $S$ |

### Coil — `uniform.xml`

```xml
<?xml version="1.0" encoding="utf-8"?>
<CoilArray>
  <IdealCoil/>
</CoilArray>
```

Ideal coil ⇒ $B_1^\pm = 1$ everywhere ⇒ one channel, no spatial RX weighting.

### Sequence — `epi.xml`

```xml
<?xml version="1.0" encoding="utf-8"?>
<Parameters FOVx="256" FOVy="256" GradMaxAmpl="2" GradSlewRate="10" Name="P" Nx="64" Ny="64" TE="50" TR="100">
   <ConcatSequence Name="C1">
      <AtomicSequence Name="A1">
         <HARDRFPULSE Axis="RF" Duration="0.1" FlipAngle="90" Name="P1"/>
      </AtomicSequence>
      <AtomicSequence Name="A2">
         <TrapGradPulse Area="-0.5*abs(A)" Axis="GX" Name="P2" Observe="A=P4.Area"/>
         <TrapGradPulse Area="KMY" Axis="GY" Name="P3" Observe="KMY=P.KMAXy"/>
      </AtomicSequence>
      <DelayAtomicSequence Delay="TE" DelayType="C2C" Name="D1" Observe="TE=P.TE" StartSeq="A1" StopSeq="C2"/>
      <ConcatSequence LoopFlag="2" Name="C2" Observe="NY=P.Ny" Repetitions="NY">
         <AtomicSequence Name="A3">
            <TrapGradPulse ADCFlag="2" ADCs="NX" Axis="GX" FlatTopArea="2*KMX*(-1)^C" FlatTopTime="1" Name="P4" Observe="KMX=P.KMAXx, C=C2.Counter, NX=P.Nx"/>
         </AtomicSequence>
         <AtomicSequence Name="A4">
            <TrapGradPulse Area="ite(1+C,R,0,-DKY)" Axis="GY" Name="P5" Observe="DKY=P.DKy, C=C2.Counter, R=C2.Repetitions"/>
         </AtomicSequence>
      </ConcatSequence>
      <DelayAtomicSequence Delay="TR" DelayType="B2E" Name="D2" Observe="TR=P.TR" StartSeq="A1"/>
   </ConcatSequence>
</Parameters>
```

Plain terms: FOV 256×256, matrix 64×64, hard 90° RF, prephasers, TE = 50 ms, EPI readout (64 ADCs × 64 lines = **4096** samples), TR = 100 ms.

Those 4096 ADCs are exactly why the output arrays have length 4096.

### Sample — `sample.h5`

| Field | Value in this example |
|--------|------------------------|
| Resolution | 1 mm × 1 mm × 1 mm |
| Offset | (0, 0, 0) |
| Live spins | 25 (“sphere”) |
| Tissue props | $M_0=1$, $R_1=0.001$, $R_2=0.01$, $R_2^*=0.01$, $\mathrm{DB}=0$ |

`DB = 0` here means **flat / uniform off-resonance** (no extra $B_0$ map). To study $B_0$ inhomogeneity with today’s engine, put a spatial ΔB into the sample’s `DB` field (or later: a dedicated $B_0$ volume).

---

## Example outputs (from this same input set)

| File | What it is | Correlates to |
|------|------------|---------------|
| `examples/outputs/approved_sim.h5` | `/signal/times` shape `(4096,)`; `/signal/channels/00` shape `(4096, 3)` | `epi.xml` ADCs × `sample.h5` × ideal RX |
| `examples/outputs/approved_sim_ismrmrd.h5` | ISMRMRD copy of the acquisition | same run, recon tools |
| `examples/outputs/gre_seq.h5` | Sequence waveforms only | from `gre.xml` (no tissue signal) |

From the EPI run:

- Times ≈ **4.83 → 95.27 ms** (set by TE / readout timing in `epi.xml`)
- **4096** samples (= $N_x \times N_y$)
- Runtime ≈ **0.21 s** for this tiny phantom

```python
import h5py
with h5py.File("examples/outputs/approved_sim.h5") as h:
    t = h["signal/times"][:]           # from sequence ADCs
    s = h["signal/channels/00"][:]     # (n_adc, 3) = (Sx, Sy, Sz)
    # peak |Mxy|-like magnitude from this input set:
    import numpy as np
    peak = np.max(np.hypot(s[:, 0], s[:, 1]))
```

Change an input → expect a matching change in $S(t)$ (examples):

- Stronger / denser sample $M_0$ → larger $\|S\|$
- Nonzero $B_0$ / `DB` map → faster dephasing, usually lower peak $\|S\|$
- Different TE / ADC grid in the sequence → different `/signal/times` and contrast
- Non-uniform RX map → channel-weighted $S$; multiple coils → `channels/00`, `01`, …

---

## Run the reference example

```bash
cd share/examples/approved
../../../run_jemris.sh simu.xml -o ../../../examples/outputs -f approved_sim
```

Sequence diagram only (no sample signal):

```bash
./run_jemris.sh share/examples/gre.xml -o examples/outputs -f gre_seq
```

---

## Build (cluster / Linux)

Dependencies: CMake, MPI (`mpicxx`), HDF5 (C+C++), Xerces-C, Sundials CVODE (~6.x), GiNaC+CLN, Boost, ISMRMRD.

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpicxx
cmake --build . --parallel
```

Includes a small HDF5 / `NDData` fix when `hsize_t` and `size_t` are the same type.

`run_jemris.sh` sets library paths if local deps live under `deps/` (gitignored).

---

## Layout

```
src/                         # simulator core
cmake/
share/examples/approved/     # reference EPI inputs (simu / sample / coils / epi)
examples/outputs/            # signals produced from those inputs
run_jemris.sh
```

Upstream: [JEMRIS/jemris](https://github.com/JEMRIS/jemris).
