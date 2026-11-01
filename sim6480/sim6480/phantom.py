"""The phantom: the real 3D-printed phantoms, voxelized where they sit in the scanner.

The pockets come from the STL files in 6480-pulse-sequence (pulseseq.phantoms), filled with one of the design's
liquids. Placement, in scanner axes:

    phantom.axis       the cylinder's axis (default x, along the RF coil form)
    rotation_deg       turn about that axis. At 0, the slots' narrow side lies along the first of z, y, x that
                       isn't the cylinder axis (z by default): z shots then show three peaks and y shots one block
    centre_mm          where the middle of the water sits (moves the phantom through the field maps)

The voxel grid is regular and aligned with the scanner axes, which is what MRzero (VoxelGridPhantom) and JEMRIS
(sample.h5) both expect.
"""

from dataclasses import dataclass

import numpy as np

from sim6480 import repos

AXES = ("x", "y", "z")


@dataclass
class Phantom:
    name: str
    liquid: dict                 # pd, t1_s, t2_s
    t2dash_s: float
    voxel_mm: float
    axes_mm: tuple               # (x, y, z) 1D coordinate arrays of the grid, true positions in mm
    pd: np.ndarray               # (nx, ny, nz) proton density (0 outside the water)
    water_mm3: float

    @property
    def mask(self):
        return self.pd > 0

    @property
    def points_mm(self):
        """(n, 3) positions of the water voxels, in the order MRzero builds them (C order of the mask)."""
        X, Y, Z = np.meshgrid(*self.axes_mm, indexing="ij")
        m = self.mask
        return np.column_stack([X[m], Y[m], Z[m]])

    def profile(self, axis):
        """The true 1-D image along a scanner axis: (positions in mm, water in mm^2 per mm)."""
        i = AXES.index(axis)
        other = tuple(j for j in range(3) if j != i)
        water = self.pd.sum(axis=other) * self.voxel_mm**2
        return self.axes_mm[i], water


def _plane_axes(cyl_axis):
    """The two scanner axes across the cylinder: (where the phantom's x goes at rotation 0, where its y goes)."""
    return tuple(a for a in ("z", "y", "x") if a != cyl_axis)


def build(design):
    """Voxelize design['phantom'] on a grid aligned with the scanner axes."""
    _, phantoms = repos.pulseseq(design)
    p = design["phantom"]
    v = float(p["voxel_mm"])
    ph = phantoms.read_phantom(p["stl"])
    liquid = design["liquids"][p["liquid"]]

    cyl = p["axis"]
    pa, pb = _plane_axes(cyl)
    centre = np.array([p["centre_mm"][a] for a in AXES], dtype=float)

    # Grid: covers the cylinder across, and the water along the axis, with a voxel of margin.
    half_plane = ph.diameter_mm / 2 + v
    half_axis = ph.depth_mm / 2 + v
    half = {a: (half_axis if a == cyl else half_plane) for a in AXES}
    axes = []
    for i, a in enumerate(AXES):
        n = int(np.ceil(half[a] / v))
        axes.append(centre[i] + np.arange(-n, n + 1) * v)
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    P = {"x": X - centre[0], "y": Y - centre[1], "z": Z - centre[2]}

    # Across the cylinder: rotate scanner coordinates back into the phantom's own frame and test the pockets.
    th = np.radians(p["rotation_deg"])
    u, w = P[pa], P[pb]
    xp = np.cos(th) * u + np.sin(th) * w
    yp = -np.sin(th) * u + np.cos(th) * w
    pts = np.column_stack([xp.ravel(), yp.ravel()])
    wet = np.zeros(len(pts), dtype=bool)
    for pocket in ph.pockets:
        wet |= phantoms.inside(pts, pocket.outline)
    wet = wet.reshape(X.shape)

    # Along the cylinder: the pockets' depth, centred on centre_mm.
    along = np.abs(P[cyl]) <= ph.depth_mm / 2
    pd = np.where(wet & along, float(liquid["pd"]), 0.0)
    return Phantom(p["stl"], liquid, design["magnet"]["t2dash_s"], v, tuple(axes), pd,
                   float((pd > 0).sum() * v**3))


def to_mrzero(phantom, voxel_fields, voxel_shape="sinc"):
    """MRzero SimData for the phantom, with the field maps' per-voxel values put in.

    voxel_fields: fieldmaps.VoxelFields sampled at phantom.points_mm (same order).
    """
    import MRzeroCore as mr0
    import torch

    t = lambda a: torch.tensor(np.asarray(a), dtype=torch.float32)
    pd = phantom.pd.astype(np.float32)
    shape = pd.shape
    ones = np.ones(shape, dtype=np.float32)
    v = phantom.voxel_mm
    vox = mr0.VoxelGridPhantom(
        PD=t(pd), T1=t(ones * phantom.liquid["t1_s"]), T2=t(ones * phantom.liquid["t2_s"]),
        T2dash=t(ones * phantom.t2dash_s), D=t(np.zeros(shape)), B0=t(np.zeros(shape)),
        B1=t(np.ones((1,) + shape)), coil_sens=t(np.ones((1,) + shape)),
        affine=torch.tensor([[v, 0, 0, 0], [0, v, 0, 0], [0, 0, v, 0]], dtype=torch.float32))
    data = vox.build(voxel_shape=voxel_shape)
    n = int(data.PD.shape[0])
    if n != len(voxel_fields.off_resonance_hz):
        raise RuntimeError(f"MRzero built {n} voxels but the field maps were sampled at {len(voxel_fields.off_resonance_hz)}")
    data.B0 = t(voxel_fields.off_resonance_hz)
    data.voxel_pos = t(voxel_fields.encoding_m)
    data.B1 = torch.tensor(voxel_fields.b1_relative[None, :], dtype=torch.complex64)
    data.coil_sens = torch.tensor(voxel_fields.coil_t_per_a[None, :], dtype=torch.complex64)
    return data
