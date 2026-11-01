# sim6480: the 6.4800 tabletop MRI, simulated from field maps

Give it **field maps** (B0, the three gradient channels, the RF coil) and a **pulse sequence**. It returns the
**voltage on the receive coil**, the ADC counts and FPGA bytes the lab would capture, and the 1-D image the lab's
processing makes from them. The hardware enters only through the maps. Swap a map for a measured one, or for one
generated from a different magnet or coil, and nothing else changes. That makes it easy to compare strategies: a
higher peak echo voltage means the spins are better aligned, as in Irene's thesis.

```bash
pip install -r requirements.txt
python -m sim6480 run --plot shot.png                       # one shot of the default design
python -m sim6480 run design/examples/aubert_b0.yaml        # the same shot in an Aubert-pair magnet
python -m sim6480 run --set gradients.axis=y --set phantom.liquid=pure_water
python -m sim6480 numbers                                   # scaling: B0, M0, EMF per mm^3, noise, gain
python -m sim6480 export-maps maps/                         # the maps a shot uses, as files you can edit
python -m sim6480 jemris out/                               # the same shot as JEMRIS input, for the cluster
python -m pytest tests -q
```

![one simulated shot](docs/shot_z.png)

*Default design, readout along z, three-slice phantom. Top: coil voltage, with RF pulses in red and the readout
gradient in green. Middle: ch2 counts as captured. Bottom: the image the lab's processing makes, against the
phantom.*

## What it does (Kevin's questions)

**What is it imaging?** The real 3D-printed phantoms. Their STL files, from
[6480-pulse-sequence](https://github.com/MITMRI/6480-pulse-sequence), are voxelized at 0.5 mm where they sit in
the scanner and filled with a liquid from the design: doped water by default, or pure water or mineral oil. The
default is the three-slice phantom: three 2 mm slots, 4 mm apart, 10 mm deep, 853 mm³ of water.

**How does it simulate pulse sequences?** The sequence is the firmware's own spin echo:
- **Timing:** read from the Teensy firmware by `pulseseq.timeline`, which counts its GPIO writes. The pulses are
  90° at 27–60 µs and 180° at 780–848 µs; the gradient lobes run dephase, readout, rewind; the echo is at about
  1.58 ms.
- **Format:** written as a standard Pulseq `.seq` file. MRzero reads that file directly, and JEMRIS reads it after
  `tools/seq_to_jemris_xml.py` converts it.
- **Your own sequence:** design key `sequence.file`.

**What is the kernel?** [MRzero](https://github.com/MRsources/MRzero-Core) (MRzeroCore 1.3), a *phase
distribution graph* (PDG) simulator:
- **Per pulse:** each RF pulse splits each voxel's magnetization into a few states. Between pulses every state
  relaxes (T1, T2, T2′) and gains phase from the voxel's off-resonance and from the gradients.
- **Speed:** about 2 s per shot on a laptop CPU, so parameter sweeps are cheap.
- **Arbitrary gradient maps:** each voxel's gradient phase is set from the gradient field it actually sits in.
  Its encoding coordinate becomes γ̄·B_grad,∥(r) / 10⁶ Hz/m, so a nonlinear gradient map distorts the image the
  way the real coil would.
- **Cross-check:** for a high-fidelity Bloch-equation check, the same experiment exports to
  [JEMRIS](https://www.jemris.org) and runs on the cluster. See `../JEMRIS_VS_MRZERO.md` for why both.

**Inputs, concretely:**

| input | form | default |
|---|---|---|
| B0 map | vector field (T) on a grid, `.npz`/`.h5` | uniform 107.3 mT (4.567 MHz) + 20 Hz/mm along z; or permanent-magnet rings (`rings:`, e.g. an **Aubert pair**) |
| gx, gy, gz maps | vector field (T) per unit of drive | linear, 3000 Hz/mm per drive; or coils by their dimensions (`loops:`, e.g. a **Maxwell / anti-Helmholtz pair**) |
| RF coil map | vector field (T per amp) | 20-turn solenoid on the 19 mm coil form, axis x |
| sequence | the firmware's timings, or a `.seq` | firmware spin echo, 128 × 1024 samples at 50 MS/s |
| phantom + liquid | STL + T1, T2, proton density | three-slice phantom, doped water (T1 0.2 s, T2 0.15 s) |
| receive chain | gains, noise, ADC | AFE5808A: LNA 24 dB + PGA 30 dB, 14-bit, 2 Vpp, 1 nV/√Hz |

**Outputs, concretely:**

| output | value for the default design |
|---|---|
| coil voltage vs time | complex baseband, V; peak echo **8.4 µV** |
| noise at the coil | 3.2 µV rms in the 10 MHz bandwidth |
| ADC counts | about ±35 counts of echo on ch2 (8.4 µV × 54 dB ÷ 0.12 mV per count); RF pulses saturate |
| FPGA bytes | exactly what the laptop receives (6480-signal-processing's packing) |
| 1-D image | along the readout axis, 0.52 mm pixels: three slices for z, one block for y |

**Scale check.** Water at 107 mT and 22 °C has M0 = 3.5×10⁻⁴ A/m. With the coil's 1.05 mT/A, reciprocity gives
10.6 nV per mm³ of water tipped 90°. For the phantom's 853 mm³ that is about 9 µV. The full simulation gives
8.4 µV, a little less because of T2′ and the coil's falloff (`python -m sim6480 numbers`).

## Field maps: the only way hardware enters

The simulator never sees "a magnet" or "a gradient coil". It sees what they produce:

```yaml
fields:
  b0:   {file: maps/b0.npz}        # or rings: [...] / offset_hz + linear + quadratic
  gz:   {file: maps/gz.npz}        # or loops: [...] + amps_per_drive / hz_per_mm_per_drive
  coil: {file: maps/coil.npz}      # or loops: [...] / solenoid: {...}
```

**File format.** A map file has `x_mm, y_mm, z_mm` (1-D grid axes) and either `Bx, By, Bz` (vector field) or `B`
(a scalar field along `b0_direction`), in tesla. Gradient maps are per unit of drive; the coil map is per amp.
`python -m sim6480 export-maps maps/` writes the current maps in this format, so you can start from them. Relative
paths are read from `design/`.

**What each map becomes:**

| map | becomes | per voxel |
|---|---|---|
| B0 | off-resonance | γ̄·\|B0(r)\| − f_rf |
| gradient maps | encoding | the field component along the local B0 |
| coil (transmit) | flip-angle map | \|B1⊥(r)\| relative to the centre |
| coil (receive) | sensitivity | \|B1⊥(r)\| in T/A; volts by reciprocity, V = ω·(B1/I)·M0·V_voxel |

**Built-in generators.** Each is checked against a closed form in `tests/`:
- **Permanent-magnet rings, `generated_rings`:** axial or radial magnetization (dipole sum over the ring volume).
  The Aubert pair, a radial-out ring plus a radial-in ring, is in `design/examples/aubert_b0.yaml`.
  `match_larmor: true` scales the map so its centre sits at `magnet.larmor_hz`.
- **Current loops, `generated_loops`:** coil pairs and solenoids, via Biot–Savart with elliptic integrals.
  `design/examples/maxwell_gz.yaml` is an anti-Helmholtz z-gradient pair.

> The ring and coil dimensions in both examples are **placeholders**. Put the real ones in.

## Changing the design

`design/default.yaml` holds every number, each tagged with where it came from:
- **[firmware]** the Teensy firmware
- **[fpga]** the FPGA driver
- **[afe]** the AFE5808A settings the firmware programs
- **[schematic]** the main board, newMRI_26_08A
- **[stl]** the 3D-printed parts
- **[guess]** not measured yet

Three ways to change it:
1. Edit the file.
2. Layer a file over it: `python -m sim6480 run my_magnet.yaml`. It only needs the keys it changes.
3. Override one key: `--set key.path=value`.

Unknown keys are an error, so typos don't pass silently.

## How it fits together

```
design/default.yaml (+ your layers)
  │
  ├─ fieldmaps.py   B0, gx, gy, gz, coil maps (files or generated) → per-voxel off-resonance, encoding, B1, sensitivity
  ├─ phantom.py     STL phantom → voxels (6480-pulse-sequence)
  ├─ sequence.py    firmware timings → Pulseq .seq (6480-pulse-sequence)
  ├─ kernel_mrzero.py   MRzero PDG → coil voltage (V)
  ├─ receive.py     RF at 4.567 MHz, noise, AFE gain, RF saturation, 14-bit ADC → FPGA words → bytes
  ├─ process.py     the lab's own processing: unpack, downconvert, echo row, image (6480-signal-processing)
  └─ jemris_export.py   the same experiment as JEMRIS files (sample.h5, sequence.xml, simu.xml, coils)
```

## Known limits

- **One shot, 1-D images.** The firmware's spin echo with a single readout gradient. Any Pulseq sequence can be
  given instead with `sequence.file`; phase encoding for 2-D works in principle but hasn't been tried here.
- **B0 direction.** MRzero assumes one common B0 direction. The maps' local direction sets each voxel's
  frequency and gradient component, but a field that rotates a lot across the sample is beyond the PDG model.
  JEMRIS handles that better.
- **JEMRIS export.**
  - Gradients are made linear (each map's centre slope); `export.json` reports the worst position error.
  - Coils are uniform.
  - Reading JEMRIS's signal back (`load_signal`) hasn't been run against a real cluster output yet.
- **Not measured yet** (the [guess] values): gradient strength per drive, coil turns and resistance, the
  preamp gain, noise density, ADC full scale, T2′. The magnet and gradient-coil dimensions for the ring and
  loop generators are also missing.
- **Orientation.** With `gradients.sign: -1` (the firmware's default) the image matches the phantom's
  orientation. Flipping the sign mirrors it, as it would on the bench.
