"""Run MRzero: the sequence + phantom + field maps -> coil voltage (complex baseband, volts) at the ADC samples.

MRzero simulates in a frame rotating at the RF frequency (`larmor_hz` on import), using phase distribution graphs:
each RF pulse splits the magnetization into a few states, and between pulses every state's phase grows with the
voxel's off-resonance (B0 map) and the gradient field it sits in (gradient maps x drive). Its signal is
sum over voxels of PD x coil sensitivity x transverse magnetization; with the coil sensitivity in T/A this becomes
volts by reciprocity:

    V(t) = omega_rf x M0_water x V_voxel x signal(t)
"""

import tempfile
from pathlib import Path

import numpy as np

from sim6480 import physics


def run(plan, phantom, data, design, seq_path=None):
    """Simulate; returns (complex baseband coil voltage at plan.adc_times_s, volts per signal unit)."""
    import MRzeroCore as mr0

    sim = design["simulation"]
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(seq_path) if seq_path else Path(tmp) / "sequence.seq"
        if not seq_path:
            plan.seq.write(str(path))
        seq = mr0.Sequence.import_file(str(path), larmor_hz=plan.f_rf_hz)
    graph = mr0.compute_graph(seq, data, sim["max_states"], sim["min_state_mag"])
    signal = mr0.execute_graph(graph, seq, data, min_emitted_signal=sim["min_state_mag"],
                               min_latent_signal=sim["min_state_mag"], print_progress=False)
    signal = signal.detach().cpu().numpy()[:, 0]
    if len(signal) != len(plan.adc_times_s):
        raise RuntimeError(f"MRzero returned {len(signal)} samples; the sequence has {len(plan.adc_times_s)}")

    b0 = physics.b0_tesla(design["magnet"]["larmor_hz"])
    m0 = physics.magnetization(b0, design["magnet"]["temperature_c"])
    volts_per_unit = 2 * np.pi * plan.f_rf_hz * m0 * (phantom.voxel_mm * 1e-3) ** 3
    return signal * volts_per_unit, volts_per_unit
