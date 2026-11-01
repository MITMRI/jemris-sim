"""One shot of the 6.4800, end to end:

    design -> field maps -> sequence (firmware timings) -> phantom (STL) -> MRzero -> coil voltage
           -> receive chain (AFE, ADC, FPGA words) -> lab processing (echo row, 1-D image)

    from sim6480.config import load_design
    from sim6480.simulate import simulate
    shot = simulate(load_design(overrides=["gradients.axis=y"]))
    shot.peak_coil_uv, shot.process.positions_mm, abs(shot.process.image)
"""

import time
from dataclasses import dataclass

import numpy as np

from sim6480 import fieldmaps, kernel_mrzero, phantom as phantom_mod, process, receive, sequence
from sim6480.config import PACKAGE_ROOT
from sim6480.physics import GAMMA_BAR


@dataclass
class Shot:
    design: dict
    plan: object               # sequence.SequencePlan
    phantom: object            # phantom.Phantom
    maps: dict                 # fieldmaps.FieldMap by name
    fields: object             # fieldmaps.VoxelFields at the water voxels
    coil_v: np.ndarray         # complex baseband coil voltage at plan.adc_times_s
    capture: object            # receive.Capture
    process: object            # process.Result
    hz_per_mm: float           # readout calibration the processing used
    seconds: float

    @property
    def peak_coil_uv(self):
        """Peak coil voltage during the readout gradient, in microvolts (the echo)."""
        t0, t1 = self.plan.readout_s
        m = (self.plan.adc_times_s >= t0) & (self.plan.adc_times_s <= t1)
        return float(np.abs(self.coil_v[m]).max() * 1e6)


def readout_hz_per_mm(maps, axis, drive, b0_direction):
    """The readout gradient's strength at the centre (Hz/mm at this drive), from its map: what a lab calibration
    with a known phantom would measure."""
    d = 0.5
    pts = np.zeros((2, 3))
    i = "xyz".index(axis)
    pts[0, i], pts[1, i] = -d, d
    g = maps[f"g{axis}"].sample(pts)[("xyz").index(b0_direction)]
    return float(GAMMA_BAR * (g[1] - g[0]) / (2 * d) * drive)


def map_extent(phantom, step_mm):
    """Half-width (mm) of the cubic field-map grid: the phantom's grid plus two steps."""
    return max(np.abs(np.concatenate(phantom.axes_mm)).max(), 1.0) + 2 * step_mm


def simulate(design, rng=None):
    start = time.time()
    ph = phantom_mod.build(design)
    maps = fieldmaps.build_maps(design, map_extent(ph, design["simulation"]["field_map_step_mm"]), design["simulation"]["field_map_step_mm"], PACKAGE_ROOT / "design")
    plan = sequence.build(design)
    fields = fieldmaps.voxel_fields(maps, ph.points_mm, plan.f_rf_hz)
    data = phantom_mod.to_mrzero(ph, fields, design["simulation"]["voxel_shape"])
    coil_v, _ = kernel_mrzero.run(plan, ph, data, design)
    cap = receive.capture(plan, coil_v, design, rng)
    g = design["gradients"]
    hz_mm = readout_hz_per_mm(maps, g["axis"], g["drive"], design["fields"]["b0_direction"])
    proc = process.run(cap.raw, plan, design, abs(hz_mm))
    return Shot(design, plan, ph, maps, fields, coil_v, cap, proc, hz_mm, time.time() - start)
