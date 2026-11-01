"""The 6.4800 spin echo as a Pulseq sequence, timed from the Teensy firmware.

The event times come from pulseseq.timeline (6480-pulse-sequence), which counts GPIO writes in the firmware, so
this follows the firmware when it changes. The design sets what the firmware can't tell us: the flip angles the
RF knob gives, the gradient drive, the readout axis (the gradient strength comes from the field maps).

One Pulseq file serves both kernels: MRzero imports it directly, and tools/seq_to_jemris_xml.py turns it into
JEMRIS XML. The ADC samples every gap after the 90° between RF pulses (the FID and the echo) at `simulation.adc_dwell_us`; the
receive model resamples that to the 50 MS/s capture.
"""

import math
from dataclasses import dataclass, field

import numpy as np
import pypulseq as pp

from sim6480 import repos
from sim6480.fieldmaps import ENCODING_HZ_PER_M
from sim6480.physics import b0_tesla


@dataclass
class SequencePlan:
    seq: object                      # pypulseq.Sequence
    f_rf_hz: float                   # transmit and demodulation frequency
    adc_times_s: np.ndarray          # time of every simulated ADC sample, from the trigger
    rf_windows_s: list               # (start, end) of each RF pulse
    spin_echo_s: float
    gradient_echo_s: float
    capture_s: float
    rf_phases_deg: list = field(default_factory=list)
    readout_s: tuple = (0.0, 0.0)    # (start, end) of the readout gradient


def _us(x):
    return int(round(x))


def build(design, firmware=None, hz_per_m_per_drive=None):
    """Build the Pulseq sequence for one shot of the design (readout along gradients.axis).

    hz_per_m_per_drive: give the gradient in real units instead (for JEMRIS, which places spins at their true
    positions and so needs the true, linear, gradient); default: the field-map encoding units MRzero uses.
    """
    timeline, _ = repos.pulseseq(design)
    fw_seq = timeline.build_sequence(firmware, design["sequence"]["write_ns"])
    step = {s.name: s for s in fw_seq.timed_steps if s.name != "Clamp"}

    f_rf = design["magnet"]["larmor_hz"] + design["rf"]["frequency_offset_hz"]
    dwell_us = design["simulation"]["adc_dwell_us"]
    capture_us = design["sequence"]["capture_blocks"] * 1024 / design["receive"]["fs_hz"] * 1e6
    g = design["gradients"]
    # Gradient amplitudes are drive x ENCODING_HZ_PER_M: the field maps (fieldmaps.py) turn that into each voxel's
    # actual frequency shift, so the sequence carries only the drive, never a gradient strength.
    amp = g["drive"] * (ENCODING_HZ_PER_M if hz_per_m_per_drive is None else hz_per_m_per_drive)
    sign = g["sign"]

    sysopts = pp.Opts(B0=b0_tesla(design["magnet"]["larmor_hz"]), max_grad=1e12, grad_unit="Hz/m", max_slew=1e18,
                      slew_unit="Hz/m/s", rf_raster_time=1e-7, grad_raster_time=1e-6, adc_raster_time=1e-7,
                      block_duration_raster=1e-6, rf_ringdown_time=0, rf_dead_time=0, adc_dead_time=0)
    seq = pp.Sequence(system=sysopts)

    # Block edges, in whole microseconds (the firmware's times rounded; well under one 10 µs processing sample).
    e = {
        "s90": _us(step["90° RF pulse"].start_us), "e90": _us(step["90° RF pulse"].end_us),
        "s180": _us(step["180° RF pulse"].start_us), "e180": _us(step["180° RF pulse"].end_us),
        "sdeph": _us(step["Dephase gradient"].start_us), "sro": _us(step["Readout gradient"].start_us),
        "srw": _us(step["Rewind gradient"].start_us), "erw": _us(step["Rewind gradient"].end_us),
        "cap": int(math.floor(capture_us)),
    }

    adc_times, rf_windows, phases = [], [], []

    def adc_block(t0, t1, grad=None):
        n = int((t1 - t0) // dwell_us)
        if n <= 0:
            return
        adc = pp.make_adc(num_samples=n, dwell=dwell_us * 1e-6, system=sysopts)
        if grad is None:
            seq.add_block(pp.make_delay((t1 - t0) * 1e-6), adc)
        else:
            seq.add_block(grad, adc, pp.make_delay((t1 - t0) * 1e-6))
        adc_times.extend((t0 + (np.arange(n) + 0.5) * dwell_us) * 1e-6)

    def rf_block(t0, t1, flip_deg, default_phase_deg):
        phase = (360.0 * f_rf * t0 * 1e-6) % 360.0 if design["rf"]["phase_from_timing"] else default_phase_deg
        rf = pp.make_block_pulse(math.radians(flip_deg), duration=(t1 - t0) * 1e-6,
                                 phase_offset=math.radians(phase), system=sysopts)
        seq.add_block(rf)
        rf_windows.append((t0 * 1e-6, t1 * 1e-6))
        phases.append(phase)

    def lobe(t0, t1, amplitude):
        return pp.make_trapezoid(g["axis"], amplitude=amplitude, flat_time=(t1 - t0 - 2) * 1e-6, rise_time=1e-6,
                                 system=sysopts)

    seq.add_block(pp.make_delay(e["s90"] * 1e-6))      # before the 90°: no spins tipped, and MRzero drops ADC here
    rf_block(e["s90"], e["e90"], design["rf"]["flip90_deg"], 0.0)
    adc_block(e["e90"], e["s180"])
    rf_block(e["s180"], e["e180"], design["rf"]["flip180_deg"], 90.0)
    adc_block(e["e180"], e["sdeph"])
    adc_block(e["sdeph"], e["sro"], lobe(e["sdeph"], e["sro"], sign * amp))      # dephase
    adc_block(e["sro"], e["srw"], lobe(e["sro"], e["srw"], -sign * amp))         # readout
    adc_block(e["srw"], e["erw"], lobe(e["srw"], e["erw"], sign * amp))          # rewind
    adc_block(e["erw"], e["cap"])

    ok, report = seq.check_timing()
    if not ok:
        raise RuntimeError(f"the Pulseq sequence has timing errors: {report[:3]}")
    seq.set_definition("Name", "sim6480_spin_echo")
    seq.set_definition("ReadoutAxis", g["axis"])

    return SequencePlan(seq, f_rf, np.asarray(adc_times), rf_windows, fw_seq.spin_echo_us * 1e-6,
                        fw_seq.gradient_echo_us * 1e-6, capture_us * 1e-6, phases,
                        (e["sro"] * 1e-6, e["srw"] * 1e-6))


def plan_from_seq(seq, f_rf_hz):
    """A SequencePlan for any Pulseq sequence (a pypulseq.Sequence, or the path of a .seq file).

    Gradient amplitudes are read as drive x ENCODING_HZ_PER_M: a lobe of 1e6 Hz/m in the file is drive 1 of that
    channel's field map. ADC samples before the first RF pulse are left out, as MRzero leaves them out.
    """
    if not isinstance(seq, pp.Sequence):
        path, seq = str(seq), pp.Sequence()
        seq.read(path)
    t, adc_times, rf_windows, phases = 0.0, [], [], []
    for i in seq.block_events:
        b = seq.get_block(i)
        if b.rf is not None:
            rf_windows.append((t + b.rf.delay, t + b.rf.delay + b.rf.shape_dur))
            phases.append(math.degrees(b.rf.phase_offset))
        if b.adc is not None and rf_windows:
            adc_times.extend(t + b.adc.delay + (np.arange(b.adc.num_samples) + 0.5) * b.adc.dwell)
        t += seq.block_durations[i]
    if not rf_windows:
        raise ValueError("the sequence has no RF pulse")
    return SequencePlan(seq, f_rf_hz, np.asarray(adc_times), rf_windows, float("nan"), float("nan"), t, phases,
                        (rf_windows[0][1], t))
