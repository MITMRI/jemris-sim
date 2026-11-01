"""Command line.

    python -m sim6480 run [design.yaml ...] [--set key.path=value ...] [--plot shot.png] [--save shot.npz]
    python -m sim6480 numbers [design.yaml ...] [--set ...]          derived scaling numbers (B0, M0, EMF, noise)
    python -m sim6480 jemris OUT_DIR [design.yaml ...] [--set ...]   the same shot as JEMRIS input files
    python -m sim6480 export-maps OUT_DIR [design.yaml ...]          the four field maps, as files fields.*.file reads

Design files layer on top of design/default.yaml, then --set overrides apply, e.g.
    python -m sim6480 run design/examples/aubert_b0.yaml --set gradients.axis=y --set phantom.liquid=pure_water
"""

import argparse
import json
from pathlib import Path

import numpy as np

from sim6480.config import DEFAULT_DESIGN, PACKAGE_ROOT, load_design


def _design(args):
    return load_design([DEFAULT_DESIGN] + [Path(f) for f in args.designs], args.set)


def plot(shot, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p, ax = shot.plan, shot.design["gradients"]["axis"]
    fig, axs = plt.subplots(3, 1, figsize=(9, 9))
    t = p.adc_times_s * 1e6
    axs[0].plot(t, np.abs(shot.coil_v) * 1e6, lw=0.8)
    for t0, t1 in p.rf_windows_s:
        axs[0].axvspan(t0 * 1e6, t1 * 1e6, color="tab:red", alpha=0.2)
    axs[0].axvspan(*(np.array(p.readout_s) * 1e6), color="tab:green", alpha=0.1, label="readout gradient")
    axs[0].set(xlabel="time from trigger (µs)", ylabel="|coil voltage| (µV)",
               title=f"coil voltage (peak echo {shot.peak_coil_uv:.2f} µV)")
    axs[0].legend(loc="upper right")
    c = shot.capture
    axs[1].plot(c.t_s * 1e6, c.counts[shot.design["receive"]["channel"]], lw=0.3)
    axs[1].set(xlabel="time from trigger (µs)", ylabel="ADC counts",
               title=f"{shot.design['receive']['channel']} as captured (RF pulses saturate the receiver)")
    pos, img = shot.process.positions_mm, np.abs(shot.process.image)
    tx, tw = shot.phantom.profile(ax)
    axs[2].plot(pos, img / img.max(), label="simulated image")
    axs[2].plot(tx, tw / tw.max(), "k--", lw=1, label="phantom (truth)")
    axs[2].set(xlabel=f"{ax} (mm)", ylabel="normalized", title=f"1-D image along {ax}")
    axs[2].legend()
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m sim6480", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "numbers", "jemris", "export-maps"):
        sp = sub.add_parser(name)
        if name == "jemris":
            sp.add_argument("out_dir", type=Path)
        if name == "export-maps":
            sp.add_argument("out", type=Path)
        sp.add_argument("designs", nargs="*", help="design YAML files, layered over design/default.yaml")
        sp.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
        if name == "run":
            sp.add_argument("--plot", type=Path)
            sp.add_argument("--save", type=Path)
    args = ap.parse_args(argv)
    design = _design(args)

    if args.cmd == "numbers":
        from sim6480.physics import summary
        for k, v in summary(design).items():
            print(f"{k:32s} {v:.4g}")
    elif args.cmd == "run":
        from sim6480.simulate import simulate
        shot = simulate(design)
        o = shot.fields.off_resonance_hz
        print(f"simulated in {shot.seconds:.1f} s: {int(shot.phantom.mask.sum())} water voxels "
              f"({shot.phantom.water_mm3:.0f} mm^3), readout {shot.hz_per_mm:.0f} Hz/mm along {design['gradients']['axis']}")
        print(f"off-resonance over the water: {o.min():.0f} .. {o.max():.0f} Hz")
        print(f"peak echo at the coil: {shot.peak_coil_uv:.2f} µV  (noise {shot.capture.noise_v_rms * 1e6:.2f} µV rms)")
        if args.plot:
            plot(shot, args.plot)
            print(f"wrote {args.plot}")
        if args.save:
            np.savez(args.save, adc_times_s=shot.plan.adc_times_s, coil_v=shot.coil_v,
                     positions_mm=shot.process.positions_mm, image=shot.process.image, raw=np.frombuffer(shot.capture.raw, np.uint8))
            print(f"wrote {args.save}")
    elif args.cmd == "jemris":
        from sim6480 import jemris_export
        info = jemris_export.export(design, args.out_dir)
        print(json.dumps(info, indent=2))
        print(f"wrote {args.out_dir}; run it with: cd {args.out_dir} && <jemris-sim>/run_jemris.sh simu.xml")
    elif args.cmd == "export-maps":
        from sim6480 import fieldmaps, phantom
        from sim6480.simulate import map_extent
        step = design["simulation"]["field_map_step_mm"]
        maps = fieldmaps.build_maps(design, map_extent(phantom.build(design), step), step, PACKAGE_ROOT / "design")
        args.out.mkdir(parents=True, exist_ok=True)
        for k, m in maps.items():
            fieldmaps.save_map(args.out / f"{k}.npz", m)
        print(f"wrote {', '.join(maps)} (.npz) to {args.out}: Bx, By, Bz in T on x_mm, y_mm, z_mm "
              f"(gradients per unit drive, coil per amp)")


if __name__ == "__main__":
    main()
