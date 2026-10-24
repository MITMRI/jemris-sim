# MRzero Documentation

MRzeroCore is a Python/PyTorch MRI simulation framework based on Phase
Distribution Graphs (PDG). Its differentiable simulation is designed for rapid
sequence experiments, reconstruction research, and gradient-based optimization.

## Installation

Create and activate a virtual environment, then install the package:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install MRzeroCore
```

## Inputs

MRzeroCore normally receives two Python objects:

- a sequence: `mr0.Sequence`, a PyPulseq sequence, or a Pulseq `.seq` path;
- a phantom: commonly `mr0.VoxelGridPhantom` or `mr0.CustomVoxelPhantom`.

A voxel-grid phantom contains tensors such as

```math
\mathrm{PD},T_1,T_2,T_2',D,B_0
\in\mathbb R^{N_x\times N_y\times N_z},
```

plus complex transmit and receive sensitivity arrays. `phantom.B0` is a scalar
off-resonance map in Hz, not an absolute main-field value in tesla.

For hydrogen, convert a measured longitudinal field deviation using

```math
\Delta f(\mathbf r)
=
\frac{\gamma}{2\pi}
\left[B_\parallel(\mathbf r)-B_{\mathrm{reference}}\right],
\qquad
\frac{\gamma}{2\pi}\approx42.57747892\times10^6\ \mathrm{Hz/T}.
```

## Simulation

The high-level workflow is concise:

```python
import MRzeroCore as mr0

sequence = mr0.Sequence.import_file("sequence.seq")
phantom = mr0.VoxelGridPhantom.load("phantom.npz")

signal, kspace = mr0.util.simulate(sequence, phantom)
```

PDG tracks transverse and longitudinal phase-distribution states. At the event
level, the evolution can be summarized as

```math
\mathbf s_{n+1}
=
\mathcal D(\Delta t,T_1,T_2,T_2',\Delta f,\Delta\mathbf k)
\mathcal R(\alpha,\phi)\mathbf s_n+\mathbf b_n.
```

This representation avoids integrating a dense set of isochromats over every
small time step. PyTorch automatic differentiation can propagate gradients
through supported simulation and reconstruction operations.

## Outputs

`mr0.util.simulate()` returns PyTorch tensors:

- `signal`: complex ADC samples, including receive-coil channels;
- `kspace`: the sampled trajectory, used for reconstruction.

Because outputs remain tensors, they can be passed directly to FFT, NUFFT,
loss, optimization, or machine-learning code.

## Important scope distinction

MRzeroCore supports a scalar off-resonance map distributed over three spatial
dimensions, `B_0(x,y,z)`. That is not the same as direct Bloch evolution under
an arbitrary vector field `(B_x,B_y,B_z)` whose direction varies spatially.
See the [comparison](JEMRIS-vs.-MRzeroCore) for the mathematical distinction.

## External references

- [Official MRzeroCore documentation](https://mrsources.github.io/MRzero-Core/)
- [MRzeroCore source](https://github.com/MRsources/MRzero-Core)
- [MRzeroCore FID example](https://github.com/MRsources/MRzero-Core/blob/main/documentation/playground_mr0/mr0_FID_seq.ipynb)
- [PDG method paper](https://doi.org/10.1002/mrm.30055)

