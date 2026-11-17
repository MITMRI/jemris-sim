"""Precompute a magnet's B0 map from its design file: by default design/maps/b0_spokes.npz from
design/examples/spokes_6480.yaml (the 6.4800 magnet, Jacob White's design).

    python tools/make_b0_map.py [--design design/examples/aubert_6480.yaml --out design/maps/b0_aubert.npz]
                                [--extent 15] [--step 1.0]
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
ap.add_argument("--design", type=Path, default=ROOT / "design" / "examples" / "spokes_6480.yaml")
ap.add_argument("--out", type=Path, default=ROOT / "design" / "maps" / "b0_spokes.npz")
ap.add_argument("--extent", type=float, default=15.0, help="half-width of the cubic grid, mm")
ap.add_argument("--step", type=float, default=1.0, help="grid step, mm")
args = ap.parse_args()

d = load_design([DEFAULT_DESIGN, args.design])
m = fieldmaps.generated_rings(fieldmaps.grid_axes(args.extent, args.step), d["fields"]["b0"]["rings"])
np.savez_compressed(args.out, x_mm=m.x_mm, y_mm=m.y_mm, z_mm=m.z_mm, Bx=m.B[0].astype(np.float32),
                    By=m.B[1].astype(np.float32), Bz=m.B[2].astype(np.float32))
c = np.linalg.norm(m.sample(np.zeros((1, 3)))[:, 0])
print(f"wrote {args.out}: centre {c * 1e3:.2f} mT ({c * 42.577e6 / 1e6:.4f} MHz)")
