"""Sweep the spacing between the two magnet holders (the adjustable "hub" separation) and report what it does to the
field over the phantom and, with --shots, to the echo.

    python tools/sweep_gap.py                        # field only: a few seconds per spacing
    python tools/sweep_gap.py --shots                # also one simulated shot per spacing (the echo peak)
    python tools/sweep_gap.py --deltas -4 -2 0 2 4   # spacing changes from the fitted design, mm

Each spacing moves the two holders symmetrically along z by delta/2. The RF is retuned to the magnet's new centre
field for each shot (magnet.larmor_hz), as the lab does with the RF knob.
"""

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sim6480 import fieldmaps, phantom as phantom_mod                 # noqa: E402
from sim6480.config import DEFAULT_DESIGN, apply_override, load_design  # noqa: E402
from sim6480.physics import GAMMA_BAR                                  # noqa: E402
from sim6480.simulate import map_extent                                # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--deltas", type=float, nargs="+", default=[-6, -4, -2, 0, 2, 4, 6],
                help="change in holder spacing from aubert_6480.yaml, mm")
ap.add_argument("--step", type=float, default=1.0, help="field-map grid step, mm")
ap.add_argument("--shots", action="store_true", help="also simulate one shot per spacing")
ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="design overrides for the shots")
args = ap.parse_args()

base = load_design([DEFAULT_DESIGN, ROOT / "design" / "examples" / "aubert_6480.yaml"], args.set)
rings = base["fields"]["b0"]["rings"]
ph = phantom_mod.build(base)
pts = ph.points_mm
axes = fieldmaps.grid_axes(map_extent(ph, args.step), args.step)
spacing0 = max(r["centre_mm"]["z"] for r in rings) - min(r["centre_mm"]["z"] for r in rings)

print(f"{'spacing':>8s} {'delta':>6s} {'B0 centre':>10s} {'f0':>9s} {'spread over water':>18s} {'rms':>7s}"
      + (f" {'peak echo':>10s}" if args.shots else ""))
print(f"{'mm':>8s} {'mm':>6s} {'mT':>10s} {'MHz':>9s} {'kHz (ppm)':>18s} {'Hz':>7s}" + (f" {'µV':>10s}" if args.shots else ""))
with tempfile.TemporaryDirectory() as tmp:
    for delta in args.deltas:
        moved = []
        for r in rings:
            z = r["centre_mm"]["z"]
            moved.append({**r, "centre_mm": {**r["centre_mm"], "z": z + np.sign(z) * delta / 2}})
        m = fieldmaps.generated_rings(axes, moved)
        f0 = GAMMA_BAR * np.linalg.norm(m.sample(np.zeros((1, 3)))[:, 0])
        f = GAMMA_BAR * np.linalg.norm(m.sample(pts), axis=0) - f0
        line = (f"{spacing0 + delta:8.1f} {delta:+6.1f} {f0 / GAMMA_BAR * 1e3:10.2f} {f0 / 1e6:9.4f} "
                f"{(f.max() - f.min()) / 1e3:8.2f} ({(f.max() - f.min()) / f0 * 1e6:5.0f})  {f.std():7.0f}")
        if args.shots:
            from sim6480.simulate import simulate
            path = Path(tmp) / f"b0_{delta:+.1f}.npz"
            fieldmaps.save_map(path, m)
            d = load_design(overrides=args.set)
            for o in (f"fields.b0.file={path}", "fields.b0.match_larmor=false", f"magnet.larmor_hz={f0}",
                      "receive.noise_nv_per_rthz=0"):
                apply_override(d, o)
            line += f" {simulate(d).peak_coil_uv:10.2f}"
        print(line, flush=True)
