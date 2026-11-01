"""From coil voltage to what the laptop receives: the 6.4800's receive chain.

    coil voltage (complex baseband at the simulator's ADC times)
      -> real voltage at 50 MS/s: Re{v(t) e^(i 2 pi f_rf t)}, zero during the RF pulses
      -> + noise (white, `noise_nv_per_rthz`, low-passed at `lpf_mhz`)
      -> x gain (preamp + LNA + PGA)
      -> RF pulses: the receiver saturates (clipped bursts, `rf_leak_counts`)
      -> 14-bit ADC (rounding and clipping), x digital gain
      -> FPGA words, coil on its channel and noise on the other (6480-signal-processing's packing)
      -> bytes
"""

from dataclasses import dataclass

import numpy as np
from scipy import signal as sps

from sim6480 import physics, repos


@dataclass
class Capture:
    t_s: np.ndarray            # sample times, 50 MS/s from the trigger
    v_coil: np.ndarray         # real coil voltage (V), before noise
    counts: dict               # {"ch1": int16 counts, "ch2": ...} after the ADC
    words: np.ndarray          # FPGA words
    raw: bytes                 # what the laptop receives
    noise_v_rms: float         # noise at the coil, in the AFE's bandwidth


def _quantize(v, vpp, bits):
    lsb = vpp / 2**bits
    return np.clip(np.round(v / lsb), -(2 ** (bits - 1)), 2 ** (bits - 1) - 1).astype(np.int16)


def capture(plan, baseband_v, design, rng=None):
    rng = rng if rng is not None else np.random.default_rng(design["simulation"]["seed"])
    r = design["receive"]
    fs = r["fs_hz"]
    n = int(design["sequence"]["capture_blocks"] * 1024)
    t = np.arange(n) / fs

    # Baseband coil voltage on the 50 MS/s clock: interpolate between the simulator's samples, zero in RF pulses.
    bb = np.interp(t, plan.adc_times_s, baseband_v.real, left=0, right=0) + \
        1j * np.interp(t, plan.adc_times_s, baseband_v.imag, left=0, right=0)
    in_rf = np.zeros(n, dtype=bool)
    for t0, t1 in plan.rf_windows_s:
        in_rf |= (t >= t0) & (t < t1)
    bb[in_rf] = 0
    # MRzero's signal is conjugate to the lab-frame convention the processing demodulates with.
    v_coil = np.real(np.conj(bb) * np.exp(2j * np.pi * plan.f_rf_hz * t))

    def noise():
        density = r["noise_nv_per_rthz"] * 1e-9
        white = rng.normal(0, density * np.sqrt(fs / 2), n)
        sos = sps.butter(4, r["lpf_mhz"] * 1e6, fs=fs, output="sos")
        return sps.sosfilt(sos, white)

    gain = physics.db_to_gain(r["preamp_db"] + r["lna_db"] + r["pga_db"])
    lsb = r["adc_vpp"] / 2 ** r["adc_bits"]
    digital = physics.db_to_gain(r["digital_db"])
    leak = np.sign(np.cos(2 * np.pi * plan.f_rf_hz * t)) * r["rf_leak_counts"] * lsb

    def to_counts(v_in):
        v_adc = v_in * gain
        v_adc = np.where(in_rf, leak, v_adc)
        c = _quantize(v_adc, r["adc_vpp"], r["adc_bits"]).astype(float)
        return np.clip(np.round(c * digital), -(2 ** 13), 2**13 - 1).astype(np.int16)

    n_coil = noise()
    other = "ch1" if r["channel"] == "ch2" else "ch2"
    counts = {r["channel"]: to_counts(v_coil + n_coil), other: _quantize(noise() * gain, r["adc_vpp"], r["adc_bits"])}

    sp = repos.sigproc(design)
    words = sp.transfer.pack_words(counts["ch1"], counts["ch2"])     # select mode: ch1 low, ch2 high
    return Capture(t, v_coil, counts, words, sp.transfer.words_to_bytes(words), float(np.std(n_coil)))
