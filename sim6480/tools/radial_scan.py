"""Simulate the 2026-11-13 radial phantom scan (design/examples/lab10_radial.yaml, 64 readout angles) and compare it
with the capture, shot by shot, through the same processing.

    python tools/radial_scan.py [--capture ../../datasets/6480_phantom_2026-11-13] [--angles 64] [--out radial.png]

Both go through 6480-signal-processing the same way: ch1, offset binary, demodulated at magnet.larmor_hz, 64 samples
at 200 kS/s around the echo, then |FFT|. Compared: the echo and FID in counts, the noise in counts, and the
projections (frequency x angle). With the capture's E_64_64.mat and index_inside.mat it also reconstructs both,
the way the lab's generate_image.m does (regularized least squares, lambda 800).
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sim6480 import repos                                                 # noqa: E402
from sim6480.config import DEFAULT_DESIGN, apply_override, load_design    # noqa: E402
from sim6480.simulate import simulate                                     # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--capture", type=Path, default=ROOT.parents[1] / "datasets" / "6480_phantom_2026-11-13")
ap.add_argument("--angles", type=int, default=64, help="readout angles over 180°")
ap.add_argument("--out", type=Path, default=Path("radial.png"))
ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
args = ap.parse_args()

base = load_design([DEFAULT_DESIGN, ROOT / "design" / "examples" / "lab10_radial.yaml"], args.set)
sp = repos.sigproc(base)
r, s = base["receive"], base["sequence"]
f_demod, rate, n = base["magnet"]["larmor_hz"], s["out_rate_hz"], s["echo_row_samples"]
angles = np.arange(args.angles) * 180.0 / args.angles


def measure(x, t_echo_s):
    """What both sides are scored on, from one shot's counts (signed, coil channel)."""
    x = np.asarray(x, dtype=float)
    bb, fs = sp.downconvert.downconvert(x, r["fs_hz"], f_demod, rate)
    t = np.arange(len(bb)) / fs
    row = sp.segment.segment_shot(x, r["fs_hz"], f_demod, rate, t_echo_s, n=n)
    late = bb[(t > 5.0e-3) & (t < 6.0e-3)]
    return dict(proj=np.abs(sp.recon.projection(row, 1.0)[1]), echo=float(np.abs(row.samples).max()),
                fid=float(np.abs(bb[(t > 0.16e-3) & (t < 0.25e-3)]).mean()),
                noise=float(np.sqrt(np.mean(np.abs(late - late.mean()) ** 2))) if len(late) else np.nan,
                t_echo=row.t_centre_s)


sim = []
for i, a in enumerate(angles):
    d = load_design([DEFAULT_DESIGN, ROOT / "design" / "examples" / "lab10_radial.yaml"], args.set)
    apply_override(d, f"gradients.angle_deg={a}")
    shot = simulate(d, rng=np.random.default_rng(i))
    sim.append(measure(shot.process.counts, shot.plan.spin_echo_s))
    print(f"\rsimulated {i + 1}/{len(angles)}", end="", flush=True)
print()

real = None
scan_file = args.capture / "mri-data" / "scan.json"
if scan_file.is_file():
    scan = json.loads(scan_file.read_text())
    if len(scan["angles_deg"]) == len(angles):
        real = []
        for i in range(len(angles)):
            raw = np.load(args.capture / "mri-data" / f"projection_{i:02d}.npz")[r["channel"]].astype(np.int32)
            x = (raw & 0x3FFF) - 8192 if r["number_format"] == "offset_binary" else sp.unpack.to_signed14(raw)
            real.append(measure(x, 2.945e-3))
    else:
        print(f"the capture has {len(scan['angles_deg'])} angles; run with --angles {len(scan['angles_deg'])} to compare")
else:
    print(f"no capture in {args.capture}: simulation only")


def corr(a, b):
    a, b = a - a.mean(), b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def summary(rows):
    return {k: np.nanmedian([x[k] for x in rows]) for k in ("echo", "fid", "noise", "t_echo")}


P_sim = np.array([x["proj"] for x in sim])
rows = [("sim", summary(sim))] + ([("capture", summary(real))] if real else [])
print(f"{'':8s} {'echo (counts)':>14s} {'FID 0.16-0.25 ms':>17s} {'noise rms':>10s} {'echo at (ms)':>13s}")
for name, m in rows:
    print(f"{name:8s} {m['echo']:14.0f} {m['fid']:17.0f} {m['noise']:10.1f} {m['t_echo'] * 1e3:13.3f}")
recon = {}
if real:
    P_real = np.array([x["proj"] for x in real])
    print(f"projections (frequency x angle): correlation {corr(P_sim, P_real):.3f}")
    try:
        from scipy.io import loadmat
        E = loadmat(args.capture / "E_64_64.mat")["E"]
        idx = loadmat(args.capture / "index_inside.mat")["index_inside"].ravel().astype(int) - 1
        A = np.vstack([E, np.sqrt(800.0) * np.eye(E.shape[1])])
        for name, P in (("sim", P_sim), ("capture", P_real)):
            p = P.T * (P_real.sum() / P.sum())                 # (position, angle), the capture's overall scale
            m = np.linalg.lstsq(A, np.r_[p.ravel(order="F"), np.zeros(E.shape[1])], rcond=None)[0]
            img = np.zeros(64 * 64)
            img[idx] = m
            recon[name] = img.reshape(64, 64, order="F")
        print(f"images (the lab's E matrix): correlation {corr(recon['sim'], recon['capture']):.3f}")
    except (FileNotFoundError, KeyError, ValueError) as e:
        print(f"no 2-D reconstruction: {e}")

import matplotlib                                                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                            # noqa: E402
panels = [("sim projections", P_sim.T)] + ([("capture projections", P_real.T)] if real else []) + \
    [(f"{k} image", np.abs(v)) for k, v in recon.items()]
fig, ax = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 4))
for a, (title, im) in zip(np.atleast_1d(ax), panels):
    a.imshow(im, origin="lower" if "projections" in title else "upper", aspect="auto" if "projections" in title else "equal")
    a.set_title(title, fontsize=9)
    a.set_xticks([])
    a.set_yticks([])
fig.tight_layout()
fig.savefig(args.out, dpi=110)
print(f"wrote {args.out}")
