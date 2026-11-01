"""The same experiment as JEMRIS input files, for a high-fidelity Bloch check on the cluster (JEMRIS isn't installed
locally; see ../run_jemris.sh).

    out/
    ├── simu.xml        sample + sequence + coils + Bloch/CVODE
    ├── sample.h5       the phantom's water voxels: M0, R1, R2, R2*, DB (off-resonance from the B0 map)
    ├── sequence.xml    the 6.4800 spin echo, converted from Pulseq by ../tools/seq_to_jemris_xml.py (+ .h5 waveforms)
    ├── sequence.seq    the Pulseq file it came from
    ├── coils.xml       IdealCoil, for both transmit and receive
    └── export.json     what was approximated, and the factor that turns JEMRIS's signal into volts

What carries over exactly: phantom, relaxation, T2', the B0 map (as DB), the timing, flips and RF phases.
What JEMRIS gets approximately: the gradients (linear, at each map's centre slope: JEMRIS places spins at their
true positions; export.json says how far each map is from linear over the phantom) and the coils (uniform; the
coil map's centre value sets the volt scale).

    cd out && ../../../run_jemris.sh simu.xml       # writes signals.h5
    shot = jemris_export.load_signal("out/signals.h5", plan)   # back through receive.capture and process.run
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from sim6480 import fieldmaps, physics, sequence
from sim6480.physics import GAMMA_BAR

AXES = ("x", "y", "z")
REPO_ROOT = Path(__file__).resolve().parents[2]
CONVERTER = REPO_ROOT / "tools" / "seq_to_jemris_xml.py"


def _centre_slopes(maps, direction, points_mm, d_mm=0.5):
    """Each gradient map's slope at the centre (Hz/m per drive), and its worst deviation from that line over the
    phantom (mm of position error)."""
    out = {}
    di = AXES.index(direction)
    for i, a in enumerate(AXES):
        m = maps[f"g{a}"]
        p = np.zeros((2, 3))
        p[0, i], p[1, i] = -d_mm, d_mm
        b = m.sample(p)[di]
        slope = GAMMA_BAR * (b[1] - b[0]) / (2 * d_mm * 1e-3)           # Hz/m per drive
        f = GAMMA_BAR * m.sample(points_mm)[di]
        f0 = GAMMA_BAR * m.sample(np.zeros((1, 3)))[di][0]
        err = (f - f0) / slope * 1e3 - points_mm[:, i] if abs(slope) > 0 else np.zeros(len(f))
        out[a] = (float(slope), float(np.abs(err).max()))
    return out


def write_sample(path, phantom, fields):
    """sample.h5 on the phantom's grid. /sample/data is (nz, ny, nx, 5): M0, R1, R2, R2* (1/ms), DB (rad/s)."""
    import h5py

    pd = phantom.pd
    m = phantom.mask
    t1, t2, t2d = phantom.liquid["t1_s"] * 1e3, phantom.liquid["t2_s"] * 1e3, phantom.t2dash_s * 1e3
    db = np.zeros(pd.shape)
    db[m] = 2 * np.pi * fields.off_resonance_hz
    data = np.zeros(pd.shape + (5,))
    data[..., 0] = pd
    data[m, 1] = 1 / t1
    data[m, 2] = 1 / t2
    data[m, 3] = 1 / t2 + 1 / t2d
    data[..., 4] = db
    v = phantom.voxel_mm
    centre = [float((ax[0] + ax[-1]) / 2) for ax in phantom.axes_mm]
    with h5py.File(path, "w") as f:
        f["/sample/data"] = np.ascontiguousarray(data.transpose(2, 1, 0, 3))
        f["/sample/resolution"] = np.full((3, 1), v)
        f["/sample/offset"] = np.array(centre).reshape(3, 1)


SIMU = """<?xml version="1.0" encoding="utf-8"?>
<simulate name="JEMRIS">
   <sample name="sim6480" uri="sample.h5"/>
   <TXcoilarray uri="coils.xml"/>
   <RXcoilarray uri="coils.xml"/>
   <parameter PositionRandomness="0"/>
   <sequence name="{name}" uri="{name}.xml"/>
   <model name="Bloch" type="CVODE"/>
</simulate>
"""

COILS = """<?xml version="1.0" encoding="utf-8"?>
<CoilArray>
  <IdealCoil/>
</CoilArray>
"""


def export(design, out_dir, shot=None):
    """Write the JEMRIS bundle for a design to out_dir; returns the export.json contents."""
    from sim6480 import phantom as phantom_mod
    from sim6480.config import PACKAGE_ROOT
    from sim6480.simulate import map_extent

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ph = shot.phantom if shot else phantom_mod.build(design)
    step = design["simulation"]["field_map_step_mm"]
    maps = shot.maps if shot else fieldmaps.build_maps(design, map_extent(ph, step), step, PACKAGE_ROOT / "design")
    direction = design["fields"]["b0_direction"]
    slopes = _centre_slopes(maps, direction, ph.points_mm)
    axis = design["gradients"]["axis"]

    plan = sequence.build(design, hz_per_m_per_drive=slopes[axis][0])
    fields = fieldmaps.voxel_fields(maps, ph.points_mm, plan.f_rf_hz)
    write_sample(out / "sample.h5", ph, fields)
    name = "sequence"
    plan.seq.write(str(out / f"{name}.seq"))
    subprocess.run([sys.executable, str(CONVERTER), str(out / f"{name}.seq"), "-o", str(out), "-n", name],
                   check=True, capture_output=True, text=True)
    (out / "simu.xml").write_text(SIMU.format(name=name))
    (out / "coils.xml").write_text(COILS)

    b0 = physics.b0_tesla(design["magnet"]["larmor_hz"])
    m0 = physics.magnetization(b0, design["magnet"]["temperature_c"])
    info = {
        "f_rf_hz": plan.f_rf_hz,
        "readout_axis": axis,
        "gradient_hz_per_m_per_drive": {a: s for a, (s, _) in slopes.items()},
        "gradient_nonlinearity_mm": {a: e for a, (_, e) in slopes.items()},
        "coil_t_per_a_used": fields.coil_centre_t_per_a,
        "coil_variation_over_phantom": float(fields.coil_t_per_a.max() / fields.coil_t_per_a.min() - 1),
        "b1_variation_over_phantom": float(fields.b1_relative.max() / fields.b1_relative.min() - 1),
        # JEMRIS's signal is magnetization summed over spins with M0 = PD; this makes it volts.
        "volts_per_unit": 2 * np.pi * plan.f_rf_hz * m0 * (ph.voxel_mm * 1e-3) ** 3 * fields.coil_centre_t_per_a,
        "adc_samples": int(len(plan.adc_times_s)),
    }
    (out / "export.json").write_text(json.dumps(info, indent=2))
    return info


def load_signal(path, plan, volts_per_unit):
    """JEMRIS's signals.h5 as complex baseband coil voltage at plan.adc_times_s (feed to receive.capture).

    Not yet checked against a cluster run: check the sample count and the sign of the imaginary part there.
    """
    import h5py

    with h5py.File(path, "r") as f:
        ch = np.asarray(f["/signal/channels/00"])
    m = ch.reshape(-1, 3) if ch.shape[-1] == 3 else ch.reshape(3, -1).T
    sig = m[:, 0] + 1j * m[:, 1]
    if len(sig) != len(plan.adc_times_s):
        raise RuntimeError(f"JEMRIS wrote {len(sig)} samples; the sequence has {len(plan.adc_times_s)}")
    return sig * volts_per_unit
