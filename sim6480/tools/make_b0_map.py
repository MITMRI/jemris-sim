"""Precompute design/maps/b0_aubert.npz from design/examples/aubert_6480.yaml (the 6.4800 magnet).

    python tools/make_b0_map.py [--extent 15] [--step 1.0]
"""

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sim6480 import fieldmaps                                  # noqa: E402
from sim6480.config import DEFAULT_DESIGN, load_design         # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--extent", type=float, default=15.0, help="half-width of the cubic grid, mm")
ap.add_argument("--step", type=float, default=1.0, help="grid step, mm")
args = ap.parse_args()

d = load_design([DEFAULT_DESIGN, ROOT / "design" / "examples" / "aubert_6480.yaml"])
m = fieldmaps.generated_rings(fieldmaps.grid_axes(args.extent, args.step), d["fields"]["b0"]["rings"])
out = ROOT / "design" / "maps" / "b0_aubert.npz"
np.savez_compressed(out, x_mm=m.x_mm, y_mm=m.y_mm, z_mm=m.z_mm, Bx=m.B[0].astype(np.float32),
                    By=m.B[1].astype(np.float32), Bz=m.B[2].astype(np.float32))
c = np.linalg.norm(m.sample(np.zeros((1, 3)))[:, 0])
print(f"wrote {out}: centre {c * 1e3:.2f} mT ({c * 42.577e6 / 1e6:.4f} MHz)")
