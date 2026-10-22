# jemris-sim

MRI physics simulation based on [JEMRIS](https://github.com/JEMRIS/jemris), set up for the MIT MRI group.

This repo contains the simulator source (`src/`), build files, and the **approved EPI example** we use as the reference pipeline:

**inputs → Bloch / CVODE simulation → signal outputs**

---

## Inputs → simulation → outputs

```
sample.h5          tissue / spin phantom
uniform.xml        coil (ideal, sensitivity = 1)
epi.xml            pulse sequence (64×64 EPI)
        │
        ▼
   simu.xml  ──►  jemris  ──►  signal HDF5 (+ optional ISMRMRD)
        │
   Bloch equations
   (CVODE solver)
```

| Stage | What it is | In this example |
|--------|------------|-----------------|
| **Input: sample** | Grid of spins with \(M_0\), \(T_1/T_2\) rates, off-resonance, positions | `share/examples/approved/sample.h5` — “sphere”, 25 live spins, 1 mm resolution |
| **Input: coils** | Transmit / receive sensitivity | `uniform.xml` — single ideal coil |
| **Input: sequence** | RF + gradients + ADC timing | `epi.xml` — 64×64 EPI, TE=50 ms, TR=100 ms |
| **Input: sim driver** | Wires the pieces together + chooses the model | `simu.xml` — Bloch + CVODE |
| **Simulation** | Each spin evolves under RF / \(G\cdot r\) / ΔB; coils sample \(M_{xy}\) at ADC times | ~0.2 s for this tiny phantom |
| **Output: signal** | Complex receive signal vs time | `examples/outputs/approved_sim.h5` |
| **Output: ISMRMRD** | Same acquisition in ISMRMRD form (for recon tools) | `examples/outputs/approved_sim_ismrmrd.h5` |

The receive “signal” is magnetization × coil sensitivity (physics units), **not** SI probe volts.

---

## Example spec we used

All files live under `share/examples/approved/`.

### `simu.xml` (top-level driver)

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

### `uniform.xml` (coil)

```xml
<?xml version="1.0" encoding="utf-8"?>
<CoilArray>
  <IdealCoil/>
</CoilArray>
```

### `epi.xml` (pulse sequence)

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

**In plain terms:** FOV 256×256, matrix 64×64, hard 90° RF, prephasers, wait to TE, EPI readout loop (Nx ADCs per line × Ny lines), wait out TR.

### `sample.h5` (phantom)

| Field | Value |
|--------|--------|
| Resolution | 1 mm × 1 mm × 1 mm |
| Offset | (0, 0, 0) |
| Live spins | 25 (“sphere”) |
| Typical tissue props | \(M_0=1\), \(R_1=0.001\), \(R_2=0.01\), \(R_2^*=0.01\), \(\Delta B=0\) |

---

## Example outputs

Committed under `examples/outputs/` from a successful run:

| File | Contents |
|------|-----------|
| `approved_sim.h5` | `/signal/times` (4096,) in ms; `/signal/channels/00` (4096, 3) = \((S_x, S_y, S_z)\) |
| `approved_sim_ismrmrd.h5` | ISMRMRD dataset for reconstruction tooling |
| `gre_seq.h5` | Sequence-diagram dump from the GRE XML (waveform QA) |

Observed from the EPI sim:

- ADC times ≈ **4.83 → 95.27 ms**
- **4096** samples
- Wall time ≈ **0.21 s** on the reference run

Inspect with Python:

```python
import h5py
with h5py.File("examples/outputs/approved_sim.h5") as h:
    print(h["signal/times"].shape)
    print(h["signal/channels/00"].shape)  # (n_adc, 3)
```

---

## Run the example

After building `jemris` (see below):

```bash
cd share/examples/approved
# ensure epi.xml / sample.h5 / uniform.xml sit next to simu.xml
../../../run_jemris.sh simu.xml -o ../../../examples/outputs -f approved_sim
```

Or dump a sequence diagram only:

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

This tree includes a small HDF5/`NDData` compatibility fix for platforms where `hsize_t` and `size_t` are the same type.

`run_jemris.sh` sets library paths when deps were installed under `deps/` locally (that directory is gitignored).

---

## Layout

```
src/                         # JEMRIS simulator core (Bloch, sequence, coils, IO)
cmake/                       # Find* modules
share/examples/              # Sequences, coils, samples
share/examples/approved/     # Reference EPI sim (simu.xml + friends)
examples/outputs/            # Checked-in signal outputs from that example
run_jemris.sh                # Convenience launcher
```

Upstream project: [JEMRIS/jemris](https://github.com/JEMRIS/jemris).
