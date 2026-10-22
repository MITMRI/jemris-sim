# jemris-sim

MRI physics simulation based on [JEMRIS](https://github.com/JEMRIS/jemris), set up for the MIT MRI group.

**Pipeline:** field / tissue / sequence inputs → Bloch (CVODE) simulation → coil signal vs time

**Quick links:** [convert Pulseq → XML](#pulseq--jemris-xml-sequence-conversion) · [parse outputs](#how-to-parse-it-python) · [run the EPI example](#run-the-reference-example) · [inputs ↔ outputs](#input--simulation--output-correlated)

---

## Pulseq → JEMRIS XML (sequence conversion)

JEMRIS sequences are **XML**. Scanner / PyPulseq files are usually **Pulseq `.seq`**. Use the repo converter to go either way into a sim-ready sequence:

```text
your.seq  →  python tools/seq_to_jemris_xml.py  →  name.xml + rf_*.h5 / grad_*.h5
                                                    └─ use as <sequence uri="…"/> in simu.xml
```

```bash
# deps once (Python >= 3.10)
pip install -r tools/requirements-seq.txt

# convert any Pulseq file
python tools/seq_to_jemris_xml.py path/to/seq.seq -o out_dir -n myseq
# writes: out_dir/myseq.xml  (+ waveform HDF5 files next to it)
```

Checked-in demo (our EPI, exported then converted):

| File | Role |
|------|------|
| `examples/sequences/epi_example.seq` | Pulseq input |
| `examples/sequences/epi_from_pulseq/epi_from_pulseq.xml` | JEMRIS XML output |
| `examples/sequences/epi_from_pulseq/rf_1.h5` | RF waveform referenced by that XML |

```bash
python tools/seq_to_jemris_xml.py examples/sequences/epi_example.seq \
  -o examples/sequences/epi_from_pulseq -n epi_from_pulseq
```

Hook it into a simulation driver:

```xml
<sequence name="epi_from_pulseq" uri="epi_from_pulseq.xml"/>
```

Run `jemris` **from the directory that contains the XML and `.h5` files** so relative waveform paths resolve.

**Notes**

- Pulseq is unrolled → converted XML is flat (no recovered loops); can be slower than hand-written looped XML like `share/examples/approved/epi.xml`.
- Other direction (XML → Pulseq): `jemris sequence.xml -x -o out_dir -f outname`
- Converter adapted from [py2jemris](https://github.com/imr-framework/py2jemris); code lives in `tools/seq_to_jemris_xml.py` and `tools/seq2xml/`.

---

## How inputs become the signal

At each time and position **r**, the spin feels an effective field:

```text
Bx, By  ←  RF waveform(t)  ×  B1+(r)          # transmit
Bz      ←  B0(r)
         + Gx(r)·gx(t) + Gy(r)·gy(t) + Gz(r)·gz(t)   # maps × waveforms
         + …                                          # optional extras
```

Magnetization **M(r,t)** evolves with the Bloch equations. At each ADC time, each receive coil forms:

```text
S(t)  =  Σ_spins   B1-(r) · Mxy(r,t) · exp(i φ)
```

So **every input exists to shape B, M, or the ADC times that produce S(t).**

Today’s checked-in example uses a **uniform** B0, **linear** gradients, and an **ideal** coil. The tables below also list the **general** inputs we care about (B0 map, G maps, …) and how each one shows up in the **same** signal outputs.

---

## Input ↔ simulation ↔ output (correlated)

| You provide | What it controls | What changes in the output |
|-------------|------------------|----------------------------|
| **Tissue / sample** (`sample.h5`) | Spin positions; M0, R1, R2, R2\*; per-spin off-resonance **DB** (B0-error stand-in) | Scale, decay, dephasing of `S(t)` |
| **B0(x,y,z) map** *(general / desired)* | Spatial main field / off-resonance → enters **Bz** | Phase accrual; peak \|S\| drops if spins disagree |
| **Gx, Gy, Gz (x,y,z) maps** *(general / desired)* | Spatial encoding fields (architecture-agnostic), times waveforms gx,gy,gz(t) | Encoding / trajectory baked into `S(t)` |
| **RF + gradient waveforms** (sequence XML) | RF amp/phase, gx(t), gy(t), gz(t), **ADC schedule** | `/signal/times` length & values; structure of `S(t)` |
| **TX / B1+ map** | Excitation pattern | How much Mxy is created → amplitude of `S` |
| **RX / B1− map** | Spatial weight per spin / coil | `/signal/channels/00`, `01`, … |
| **Sim driver** (`simu.xml`) | Wires sample + coils + sequence + Bloch/CVODE | One run → one signal (+ optional ISMRMRD) |

**Outputs (always tied back to those inputs):**

| Output path | Contents | Driven by |
|-------------|----------|-----------|
| `/signal/times` | ADC timestamps (ms) | Sequence ADC schedule |
| `/signal/channels/NN` | `(Sx, Sy, Sz)` vs time | RX maps × evolved M; length = #ADCs |
| `*_ismrmrd.h5` | Same acquisition, ISMRMRD layout | Sequence + signal (for recon) |
| Sequence dump (`seqdiag/*`) | TX / GX / GY / GZ vs time | Sequence XML only (no tissue) |

`S(t)` is magnetization × sensitivity (sim units), **not** SI probe volts unless you apply a calibration later.

---

## Diagram

```text
 GENERAL INPUTS                         TODAY'S EXAMPLE FILES
 ────────────────                       ─────────────────────
 B0(x,y,z)          ─┐                  (flat; DB=0 in sample.h5)
 Gx,Gy,Gz(x,y,z)    ─┼─►  B(r,t)
 gx,gy,gz(t) + RF   ─┤                  epi.xml   (waveforms + ADCs)
 B1+(r), B1-(r)     ─┤                  uniform.xml (IdealCoil)
 tissue M0,R1,R2,…  ─┘                  sample.h5

              │
              ▼
        Bloch + CVODE  (jemris)
              │
              ▼
     S(t) per RX channel
              │
              ├─► *.h5           /signal/times , /signal/channels/*
              └─► *_ismrmrd.h5   recon-ready copy
```

---

## Reference example (checked in)

Specs: `share/examples/approved/`  
Outputs from **this same input set**: `examples/outputs/`

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

| Line in `simu.xml` | How it shows up in the signal |
|--------------------|-------------------------------|
| `sample.h5` | Which spins exist + relaxation / DB → shape of `S(t)` |
| `uniform.xml` (RX) | Weight = 1 everywhere → single channel `/signal/channels/00` |
| `uniform.xml` (TX) | Ideal excitation |
| `epi.xml` | 4096 ADC times + EPI encoding → length/structure of `S(t)` |
| Bloch / CVODE | Engine that turns the above into M, then S |

### Coil — `uniform.xml`

```xml
<?xml version="1.0" encoding="utf-8"?>
<CoilArray>
  <IdealCoil/>
</CoilArray>
```

Ideal coil ⇒ B1± = 1 ⇒ one channel, no spatial RX weighting.

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

Plain terms: FOV 256×256, matrix 64×64, hard 90° RF, prephasers, TE = 50 ms, EPI readout (**64 ADCs × 64 lines = 4096 samples**), TR = 100 ms.

Those 4096 ADCs are why the output arrays have length **4096**.

### Sample — `sample.h5`

| Field | Value in this example |
|--------|------------------------|
| Resolution | 1 mm × 1 mm × 1 mm |
| Offset | (0, 0, 0) |
| Live spins | 25 (“sphere”) |
| Tissue props | M0 = 1, R1 = 0.001, R2 = 0.01, R2\* = 0.01, **DB = 0** |

`DB = 0` ⇒ flat off-resonance (no extra B0 map). To study B0 inhomogeneity with today’s engine, write a spatial ΔB into the sample’s `DB` field (later: a dedicated B0 volume).

---

## Example outputs (from this same input set)

| File | What it is | Correlates to |
|------|------------|---------------|
| `examples/outputs/approved_sim.h5` | JEMRIS native signal HDF5 | `epi.xml` ADCs × `sample.h5` × ideal RX |
| `examples/outputs/approved_sim_ismrmrd.h5` | Same acquisition in ISMRMRD layout | recon tools |
| `examples/outputs/gre_seq.h5` | Sequence waveforms only | `gre.xml` (no tissue signal) |

From the EPI run:

- Times ≈ **4.83 → 95.27 ms** (median Δt ≈ 0.0156 ms) ← TE / readout in `epi.xml`
- **4096** samples ← Nx × Ny
- Peak \|Mxy\|-like magnitude ≈ **0.65**
- Runtime ≈ **0.21 s** for this tiny phantom

### Output layout (`approved_sim.h5`)

```text
approved_sim.h5
├── signal/
│   ├── times              float64 (N,)      ADC times in milliseconds
│   └── channels/
│       └── 00             float64 (N, 3)    channel 0: columns [Sx, Sy, Sz]
│       └── 01             …                 (only if multi-coil RX)
└── kspace                 float64 (…,)      optional k-space helper data
```

| Dataset | Shape (this example) | Meaning |
|---------|----------------------|---------|
| `/signal/times` | `(4096,)` | ADC clock (ms), set by the sequence |
| `/signal/channels/00` | `(4096, 3)` | Coil 0 signal: **Sx, Sy, Sz** (mag × sensitivity units) |
| `/signal/channels/NN` | `(N, 3)` | Extra RX channels if you used a multi-coil map |

Complex transverse signal is `Sx + 1j·Sy`; `Sz` is the longitudinal component contributed at ADC times.

### How to parse it (Python)

```bash
pip install h5py numpy matplotlib   # if needed
```

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

path = "examples/outputs/approved_sim.h5"

with h5py.File(path, "r") as h:
    # list everything
    h.visit(lambda name: print(name))

    t = h["signal/times"][:]                 # (N,) ms
    s = h["signal/channels/00"][:]           # (N, 3)

sx, sy, sz = s[:, 0], s[:, 1], s[:, 2]
mxy = sx + 1j * sy
mag = np.abs(mxy)
phase = np.angle(mxy)

print(f"N = {len(t)}, t = [{t[0]:.3f}, {t[-1]:.3f}] ms")
print(f"peak |Mxy| = {mag.max():.4f}")

# optional: all coil channels
with h5py.File(path, "r") as h:
    chans = sorted(h["signal/channels"].keys())
    stack = np.stack([h[f"signal/channels/{c}"][:, :2] for c in chans], axis=0)
    # stack shape: (n_coils, N, 2) as (Sx, Sy)

plt.figure(figsize=(8, 3))
plt.plot(t, mag)
plt.xlabel("time (ms)")
plt.ylabel("|Sxy|")
plt.title("Receive signal magnitude vs time")
plt.tight_layout()
plt.show()
```

### ISMRMRD file (`approved_sim_ismrmrd.h5`)

Same acquisition, packaged for reconstruction toolchains:

```text
dataset/
├── xml     experiment / encoding header (string)
└── data    structured acquisitions (headers + traj + samples)
```

Prefer the native `approved_sim.h5` for quick analysis of `S(t)`. Use the ISMRMRD file with [ismrmrd](https://github.com/ismrmrd/ismrmrd) / BART / etc. when you want image recon.

### What changes when you change inputs

- Larger sample M0 → larger \|S\|
- Nonzero B0 / DB map → more dephasing, usually lower peak \|S\|
- Different TE / ADC grid → different `/signal/times` and contrast
- Non-uniform / multi-coil RX → weighted or multiple `/signal/channels/*`

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

```text
src/                         # simulator core
cmake/
share/examples/approved/     # reference EPI inputs (simu / sample / coils / epi)
examples/outputs/            # signals produced from those inputs
examples/sequences/          # sample .seq + converted XML demo
tools/seq_to_jemris_xml.py   # Pulseq → JEMRIS XML CLI
tools/seq2xml/               # converter library
run_jemris.sh
```

Upstream: [JEMRIS/jemris](https://github.com/JEMRIS/jemris).
