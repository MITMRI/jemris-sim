# Outputs

JEMRIS produces time-domain receive signals. This repository includes native
JEMRIS HDF5 outputs, optional ISMRMRD output, and sequence-only diagnostic data.

## Native signal HDF5

The primary datasets are:

| HDF5 path | Shape | Meaning |
|---|---:|---|
| `/signal/times` | `(n_adc,)` | ADC timestamps in milliseconds |
| `/signal/channels/00` | `(n_adc, 3)` | Signal components `S_x`, `S_y`, and `S_z` for receive channel 0 |
| `/signal/channels/NN` | `(n_adc, 3)` | Additional receive channels when configured |

The transverse complex signal is commonly formed as

```math
S_{xy}(t)=S_x(t)+iS_y(t).
```

A minimal Python reader is:

```python
import h5py

with h5py.File("simulation.h5", "r") as h5:
    time_ms = h5["/signal/times"][:]
    components = h5["/signal/channels/00"][:]

signal = components[:, 0] + 1j * components[:, 1]
```

These values are simulated magnetization weighted by receive sensitivity. They
are not calibrated SI probe voltages unless an external calibration is applied.

## ISMRMRD

An optional `*_ismrmrd.h5` file stores the acquisition in a reconstruction-ready
ISMRMRD layout. Use this form when passing simulated data into reconstruction
software that understands ISMRMRD.

## Sequence diagnostics

Sequence diagnostic output records RF and gradient waveforms versus time. It is
determined by the sequence alone and does not contain the tissue-dependent
receive signal.

## Reference outputs

Checked-in examples are under `examples/outputs/`:

| File | Contents |
|---|---|
| `approved_sim.h5` | Native JEMRIS receive signal for the EPI example |
| `approved_sim_ismrmrd.h5` | The same acquisition in ISMRMRD layout |
| `gre_seq.h5` | GRE sequence waveform diagnostics |

The reference EPI sequence has `64 × 64 = 4096` ADC samples, so its time and
signal arrays contain 4096 rows.

Previous: [run the simulation](Simulation). Return to [Home](Home).
