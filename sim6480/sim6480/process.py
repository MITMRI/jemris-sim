"""What the lab's processing makes of the capture: 6480-signal-processing, run exactly as on real data.

    bytes -> words -> channels (unpack) -> downconvert -> echo row (segment_shot) -> 1-D image (projection)
"""

from dataclasses import dataclass

import numpy as np

from sim6480 import repos


@dataclass
class Result:
    counts: np.ndarray        # the coil's channel, unpacked
    baseband: np.ndarray      # downconverted, at out_rate_hz
    rate_hz: float
    row: object               # sigproc EchoRow
    positions_mm: np.ndarray
    image: np.ndarray         # complex 1-D image


def readout_label(design):
    """The row's axis label: x, y or z, or the angle in the y-z plane for a turned readout."""
    g = design["gradients"]
    return g["axis"] if g["angle_deg"] is None else f"yz{g['angle_deg']:g}"


def run(raw, plan, design, hz_per_mm):
    """Process a capture the way the lab does. hz_per_mm: the calibration used to turn frequency into position."""
    sp = repos.sigproc(design)
    r, s = design["receive"], design["sequence"]
    ch = sp.unpack.unpack(sp.unpack.bytes_to_words(raw), sel_0=1, sel_1=2)[r["channel"]]
    if r["number_format"] == "offset_binary":     # unpack read two's complement: undo that, then 8192 is zero
        ch = ((ch.astype(np.int32) & 0x3FFF) - 8192).astype(np.int16)
    bb, rate = sp.downconvert.downconvert(ch, r["fs_hz"], plan.f_rf_hz, s["out_rate_hz"])
    row = sp.segment.segment_shot(ch, r["fs_hz"], plan.f_rf_hz, s["out_rate_hz"], plan.spin_echo_s,
                                  n=s["echo_row_samples"], labels={"axis": readout_label(design)})
    pos, img = sp.recon.projection(row, hz_per_mm)
    return Result(np.asarray(ch), bb, rate, row, np.asarray(pos), np.asarray(img))
