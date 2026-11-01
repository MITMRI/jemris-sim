"""Field maps + pulse sequences in, coil voltage out.

    python demo/demo.py                  # ~20 s; writes demo/out/ and docs/demo.png

1. Write the field maps (B0, gx, gy, gz, RF coil) to plain .npz files: arrays of (Bx, By, Bz) in tesla on a grid.
   From here on the simulator sees only those files, never a magnet or coil.
2. Write pulse sequences as standard Pulseq .seq files, built here with pypulseq.
3. Read out the coil voltage for each (maps, sequence) pair, and compare strategies by the voltage:
   a. RF pulse length sweep: the peak FID voltage is highest when the pulse tips the spins a full 90°.
   b. Swapping only the B0 map: the fitted magnet, the holders 2 mm further apart, and an ideal uniform field.
   c. Swapping only the gz map: a linear gradient against an anti-Helmholtz coil pair, seen in the 1-D image
      (with an ideal uniform B0 for reference: the magnet, not the gradient, limits this image).
"""

import math
import sys
from pathlib import Path

import numpy as np
import pypulseq as pp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sim6480 import fieldmaps, phantom as phantom_mod                         # noqa: E402
from sim6480.config import DEFAULT_DESIGN, load_design                        # noqa: E402
from sim6480.fieldmaps import ENCODING_HZ_PER_M                               # noqa: E402
from sim6480.physics import GAMMA_BAR, b0_tesla                               # noqa: E402
from sim6480.simulate import map_extent, readout_hz_per_mm, simulate_sequence  # noqa: E402

OUT = ROOT / "demo" / "out"
EXAMPLES = ROOT / "design" / "examples"
T90_US = 33.0                 # the firmware's 90° pulse length (27-60 µs); the RF amplitude is set so this is 90°
DWELL = 1e-6


# ---- 1. field maps, as files -------------------------------------------------------------------------------------

def write_maps():
    """Every map the demo uses, as .npz files. Returns {name: path}."""
    (OUT / "maps").mkdir(parents=True, exist_ok=True)
    d = load_design()
    ph = phantom_mod.build(d)
    step = d["simulation"]["field_map_step_mm"]
    ext = map_extent(ph, step)
    paths = {}
    for k, m in fieldmaps.build_maps(d, ext, step, ROOT / "design").items():
        paths[k] = OUT / "maps" / f"{k}.npz"
        fieldmaps.save_map(paths[k], m)

    coarse = fieldmaps.grid_axes(map_extent(ph, 1.0), 1.0)      # the generated magnets: 1 mm grid is plenty
    rings = load_design([DEFAULT_DESIGN, EXAMPLES / "aubert_6480.yaml"])["fields"]["b0"]["rings"]
    wider = [{**r, "centre_mm": {"z": r["centre_mm"]["z"] + np.sign(r["centre_mm"]["z"])}} for r in rings]
    extra = {
        "b0_holders_2mm_apart": fieldmaps.generated_rings(coarse, wider),
        "b0_uniform": fieldmaps.generated_b0(coarse, d["magnet"]["larmor_hz"], 0, {}, {}, "z"),
        "gz_maxwell_pair": fieldmaps.build_maps(load_design([DEFAULT_DESIGN, EXAMPLES / "maxwell_gz.yaml"]),
                                                ext, step, ROOT / "design")["gz"],
    }
    for k, m in extra.items():
        paths[k] = OUT / "maps" / f"{k}.npz"
        fieldmaps.save_map(paths[k], m)
    return paths


def design_from_maps(paths, b0="b0", gz="gz"):
    """A design whose hardware is only these map files. The RF is tuned to the B0 map's centre, as on the bench."""
    m = fieldmaps.load_map(paths[b0])
    f0 = GAMMA_BAR * np.linalg.norm(m.sample(np.zeros((1, 3)))[:, 0])
    return load_design(overrides=[
        f"fields.b0.file={paths[b0]}", "fields.b0.match_larmor=false", f"magnet.larmor_hz={f0}",
        f"fields.gx.file={paths['gx']}", f"fields.gy.file={paths['gy']}", f"fields.gz.file={paths[gz]}",
        f"fields.coil.file={paths['coil']}", "receive.noise_nv_per_rthz=0"])


# ---- 2. pulse sequences, as Pulseq ----------------------------------------------------------------------------------

def _system(larmor_hz):
    return pp.Opts(B0=b0_tesla(larmor_hz), max_grad=1e12, grad_unit="Hz/m", max_slew=1e18, slew_unit="Hz/m/s",
                   rf_raster_time=1e-7, grad_raster_time=1e-6, adc_raster_time=1e-7, block_duration_raster=1e-6,
                   rf_ringdown_time=0, rf_dead_time=0, adc_dead_time=0)


def _rf(seq, sysopts, dur_us, phase_deg=0.0):
    seq.add_block(pp.make_block_pulse(math.radians(90 * dur_us / T90_US), duration=dur_us * 1e-6,
                                      phase_offset=math.radians(phase_deg), system=sysopts))


def _adc(seq, sysopts, dur_s, grad=None):
    adc = pp.make_adc(num_samples=int(round(dur_s / DWELL)), dwell=DWELL, system=sysopts)
    seq.add_block(*([grad] if grad else []), adc, pp.make_delay(dur_s))


def fid(larmor_hz, pulse_us, read_s=200e-6):
    """One RF pulse at fixed amplitude, then read the FID."""
    s = _system(larmor_hz)
    seq = pp.Sequence(system=s)
    seq.add_block(pp.make_delay(20e-6))
    _rf(seq, s, pulse_us)
    _adc(seq, s, read_s)
    return seq


def spin_echo(larmor_hz, te_s=2e-3, drive=0.0, axis="z", read_s=600e-6):
    """90° - TE/2 - 180° - echo at TE. With drive > 0: a dephase lobe after the 180° and a readout lobe centred on
    the echo (drive in the gradient maps' units; 1 = drive 1 of that map)."""
    s = _system(larmor_hz)
    seq = pp.Sequence(system=s)
    seq.add_block(pp.make_delay(20e-6))
    _rf(seq, s, T90_US)
    t90_end = 20e-6 + T90_US * 1e-6
    t180 = 20e-6 + T90_US / 2 * 1e-6 + te_s / 2 - T90_US * 1e-6            # 180° centred at TE/2
    _adc(seq, s, round((t180 - t90_end) / DWELL) * DWELL)
    _rf(seq, s, 2 * T90_US, 90.0)
    t_after = t180 + 2 * T90_US * 1e-6
    echo = 20e-6 + T90_US / 2 * 1e-6 + te_s
    if drive:
        amp = drive * ENCODING_HZ_PER_M
        lobe = lambda a, d: pp.make_trapezoid(axis, amplitude=a, flat_time=d - 2e-6, rise_time=1e-6, system=s)
        gap = round((echo - read_s / 2 - read_s / 2 - t_after) / DWELL) * DWELL
        _adc(seq, s, gap)
        _adc(seq, s, read_s / 2, lobe(-amp, read_s / 2))      # dephase (the 180° already flipped the phase)
        _adc(seq, s, read_s, lobe(amp, read_s))               # readout: gradient echo at the spin echo
    else:
        _adc(seq, s, round((echo + read_s / 2 - t_after) / DWELL) * DWELL)
    return seq


# ---- 3. coil voltage for each strategy -----------------------------------------------------------------------------

def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = write_maps()
    print("field maps (the only hardware the simulator sees):")
    for k, p in paths.items():
        print(f"  {p.relative_to(ROOT)}")
    (OUT / "seqs").mkdir(parents=True, exist_ok=True)

    d = design_from_maps(paths)
    ph = phantom_mod.build(d)
    f0 = d["magnet"]["larmor_hz"]

    # a. pulse length sweep
    pulses = np.arange(4, 128, 4.0)                           # up to ~350°
    peaks = []
    for p in pulses:
        seq = fid(f0, p)
        seq.write(str(OUT / "seqs" / f"fid_{p:.0f}us.seq"))
        peaks.append(simulate_sequence(d, seq, phantom=ph).peak_coil_uv)
    best = pulses[int(np.argmax(peaks))]
    print(f"\na. pulse length sweep: peak FID {max(peaks):.2f} µV at {best:.0f} µs (90° is {T90_US:.0f} µs)")

    # b. swap only the B0 map
    print("\nb. same spin echo (TE 2 ms), only the B0 map changes:")
    traces = {}
    for name in ("b0", "b0_holders_2mm_apart", "b0_uniform"):
        db = design_from_maps(paths, b0=name)
        seq = spin_echo(db["magnet"]["larmor_hz"])
        seq.write(str(OUT / "seqs" / f"spin_echo_{name}.seq"))
        sh = simulate_sequence(db, seq, phantom=ph)
        t, v = sh.t_s, np.abs(sh.coil_v) * 1e6
        fid_v = v[t < 0.5e-3]
        low = fid_v < fid_v.max() / math.e
        t2s = f"{(t[np.argmax(low)] - t[0]) * 1e3:.2f} ms" if low.any() else "> 0.5 ms"
        traces[name] = (t, v)
        print(f"  {name:22s} FID decays to 1/e in {t2s}; echo peak {v[t > 1.5e-3].max():.2f} µV")

    # c. swap only the gz map
    print("\nc. spin echo with a readout gradient, swapping the gz map (and, for reference, the B0 map):")
    images = {}
    for name, b0, gz in (("linear gz map", "b0", "gz"), ("Maxwell-pair gz map", "b0", "gz_maxwell_pair"),
                         ("linear gz, ideal uniform B0", "b0_uniform", "gz")):
        dg = design_from_maps(paths, b0=b0, gz=gz)
        seq = spin_echo(dg["magnet"]["larmor_hz"], drive=1.0)
        seq.write(str(OUT / "seqs" / f"spin_echo_read_{gz}_{b0}.seq"))
        sh = simulate_sequence(dg, seq, phantom=ph)
        hz_mm = readout_hz_per_mm(fieldmaps.build_maps(dg, map_extent(ph, 0.5), 0.5, ROOT / "design"), "z", 1.0, "z")
        n = 600
        ro = sh.coil_v[-n:]
        spec = np.fft.fftshift(np.fft.fft(ro, 8 * n))
        x = np.fft.fftshift(np.fft.fftfreq(8 * n, DWELL)) / abs(hz_mm)
        if hz_mm > 0:
            x = -x                                     # MRzero's signal is conjugate to the lab convention
        images[name] = (x, np.abs(spec) / np.abs(spec).max())
        print(f"  {name:28s} {hz_mm:6.0f} Hz/mm at the centre; echo peak {np.abs(sh.coil_v[-n:]).max() * 1e6:.2f} µV")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(pulses, peaks, "o-")
    ax[0].axvline(T90_US, color="k", ls=":", label="90° (33 µs)")
    ax[0].set(xlabel="RF pulse length (µs), fixed amplitude", ylabel="peak FID voltage at the coil (µV)",
              title="a. Pulse length sweep: max voltage = best tip")
    ax[0].legend()
    labels = {"b0": "fitted 6.4800 magnet", "b0_holders_2mm_apart": "holders 2 mm further apart",
              "b0_uniform": "ideal uniform B0"}
    for k, (t, v) in traces.items():
        ax[1].plot(t * 1e3, v, lw=0.9, label=labels[k])
    ax[1].set(xlabel="time (ms)", ylabel="|coil voltage| (µV)", title="b. Same spin echo, only the B0 map swapped")
    ax[1].legend()
    tx, tw = ph.profile("z")
    for k, (x, im) in images.items():
        ax[2].plot(x, im, label=k, **({"ls": "--", "lw": 1} if "uniform" in k else {}))
    ax[2].fill_between(tx, tw / tw.max(), color="0.85", label="phantom (truth)")
    ax[2].set(xlim=(-8, 8), xlabel="z (mm)", ylabel="normalized", title="c. Swapping the gz map: 1-D image")
    ax[2].legend()
    fig.tight_layout()
    out = ROOT / "docs" / "demo.png"
    fig.savefig(out, dpi=120)
    print(f"\nwrote {out.relative_to(ROOT)}, and the maps and .seq files under {OUT.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
