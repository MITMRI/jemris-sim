"""Field maps: the simulator's view of the hardware. Only fields, never the parts that make them.

A field map is a 3D grid of magnetic field values, in tesla:

    x_mm, y_mm, z_mm    1D grid axes (mm), increasing
    B                   (nx, ny, nz) field magnitude, or
    Bx, By, Bz          (nx, ny, nz) each: the field vector

stored in a .npz or .h5 file with those names (`save_map` writes one). Four maps describe the scanner:

    b0      the static field (T)
    gx, gy, gz   the field each gradient channel makes per unit of drive (T per unit); the sequence supplies the
            drive over time (in whatever units the map uses: amps, DAC steps, ...)
    coil    the RF coil's field per amp (T/A): sets each voxel's flip angle (transmit) and how well the coil hears it
            (receive, by reciprocity)

How they enter the simulation, per voxel at position r:

    off-resonance   f(r) - f_rf, with f(r) = gamma_bar |B0(r)|
    encoding        a gradient's effect on precession is its field along the local B0 direction; MRzero's linear
                    -gradient math is fed that as the voxel's encoding coordinate, so any (nonlinear) map works
    transmit        relative flip angle |B1_perp(r)| / |B1_perp(centre)|
    receive         |B1_perp(r)| in T/A: EMF = omega * B1_perp/I * M * V  (reciprocity)

where B1_perp is the coil field's component perpendicular to the local B0.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import ellipe, ellipk

from sim6480.physics import GAMMA_BAR, MU0

AXES = ("x", "y", "z")
ENCODING_HZ_PER_M = 1e6   # Pulseq gradient amplitude per unit of drive (Hz/m); fixes the encoding-coordinate scale


@dataclass
class FieldMap:
    x_mm: np.ndarray
    y_mm: np.ndarray
    z_mm: np.ndarray
    B: np.ndarray            # (3, nx, ny, nz) vector field in T (a scalar map is stored along `scalar_axis`)
    name: str = ""

    def sample(self, points_mm):
        """The field vector (3, n) at each point (n, 3) in mm. Points outside the map are an error."""
        pts = np.asarray(points_mm, dtype=float)
        lo = np.array([self.x_mm[0], self.y_mm[0], self.z_mm[0]])
        hi = np.array([self.x_mm[-1], self.y_mm[-1], self.z_mm[-1]])
        outside = np.any((pts < lo - 1e-9) | (pts > hi + 1e-9), axis=1)
        if outside.any():
            raise ValueError(f"{self.name} map covers {lo}..{hi} mm, but {outside.sum()} sample points fall outside it")
        grid = (self.x_mm, self.y_mm, self.z_mm)
        return np.stack([RegularGridInterpolator(grid, self.B[c])(pts) for c in range(3)])


def load_map(path, name="", scalar_axis="z"):
    """Read a field map file (.npz or .h5). A scalar `B` is treated as a field along `scalar_axis`."""
    path = Path(path)
    if path.suffix == ".npz":
        f = dict(np.load(path))
    elif path.suffix in (".h5", ".hdf5"):
        import h5py
        with h5py.File(path, "r") as h:
            f = {k: np.array(h[k]) for k in h.keys()}
    else:
        raise ValueError(f"{path}: field maps are .npz or .h5")
    axes = [np.asarray(f[f"{a}_mm"], dtype=float) for a in AXES]
    shape = tuple(len(a) for a in axes)
    if "B" in f:
        B = np.zeros((3,) + shape)
        B[AXES.index(scalar_axis)] = f["B"]
    elif all(k in f for k in ("Bx", "By", "Bz")):
        B = np.stack([f["Bx"], f["By"], f["Bz"]]).astype(float)
    else:
        raise ValueError(f"{path}: needs B, or Bx, By and Bz")
    if B.shape[1:] != shape:
        raise ValueError(f"{path}: field shape {B.shape[1:]} doesn't match the grid {shape}")
    return FieldMap(*axes, B, name or path.stem)


def save_map(path, fmap):
    """Write a FieldMap as .npz (x_mm, y_mm, z_mm, Bx, By, Bz)."""
    np.savez(path, x_mm=fmap.x_mm, y_mm=fmap.y_mm, z_mm=fmap.z_mm, Bx=fmap.B[0], By=fmap.B[1], Bz=fmap.B[2])


def grid_axes(extent_mm, step_mm):
    """Symmetric 1D axes covering +-extent_mm."""
    n = int(np.ceil(extent_mm / step_mm))
    a = np.arange(-n, n + 1) * step_mm
    return a, a.copy(), a.copy()


# ---- Default maps, generated from the design when no file is given ------------------------------------------------

def generated_b0(axes, larmor_hz, offset_hz, linear_hz_per_mm, quadratic_hz_per_mm2, direction="z"):
    """B0 = (larmor + offset + linear·r + quadratic·r²) / gamma_bar, along `direction`."""
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    f = np.full(X.shape, larmor_hz + offset_hz, dtype=float)
    for a, R in zip(AXES, (X, Y, Z)):
        f += linear_hz_per_mm.get(a, 0) * R + quadratic_hz_per_mm2.get(a, 0) * R**2
    B = np.zeros((3,) + X.shape)
    B[AXES.index(direction)] = f / GAMMA_BAR
    return FieldMap(*axes, B, "b0")


def generated_gradient(axes, channel, hz_per_mm_per_drive, along="z"):
    """A linear gradient: field along `along` (the B0 direction) = (Hz/mm per drive) · r_channel / gamma_bar."""
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    R = {"x": X, "y": Y, "z": Z}[channel]
    B = np.zeros((3,) + X.shape)
    B[AXES.index(along)] = hz_per_mm_per_drive * R / GAMMA_BAR
    return FieldMap(*axes, B, f"g{channel}")


def _loop_field(rho, zz, a, current=1.0):
    """Field (B_rho, B_z) of a circular loop of radius a at the origin (axis z), at cylindrical (rho, z). SI units."""
    rho = np.maximum(rho, 1e-12)
    k2 = 4 * a * rho / ((a + rho) ** 2 + zz**2)
    K, E = ellipk(k2), ellipe(k2)
    c = MU0 * current / (2 * np.pi * np.sqrt((a + rho) ** 2 + zz**2))
    denom = (a - rho) ** 2 + zz**2
    bz = c * (K + (a**2 - rho**2 - zz**2) / denom * E)
    brho = c * zz / rho * (-K + (a**2 + rho**2 + zz**2) / denom * E)
    return brho, bz


def generated_loops(axes, loops, name, scale=1.0):
    """Field (T per amp, x scale) of any set of circular current loops (Biot-Savart): coil pairs, solenoids.

    Each loop: {radius_mm, axis, centre_mm: {x, y, z}, turns (default 1), sign (+1/-1, current direction)}.
    An anti-Helmholtz (Maxwell) pair along z, radius R and separation d:
        [{radius_mm: R, axis: z, centre_mm: {z: +d/2}, sign: +1}, {radius_mm: R, axis: z, centre_mm: {z: -d/2}, sign: -1}]
    """
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    B = np.zeros((3,) + X.shape)
    for loop in loops:
        axis = loop["axis"]
        c = loop.get("centre_mm", {}) or {}
        P = {"x": X - c.get("x", 0), "y": Y - c.get("y", 0), "z": Z - c.get("z", 0)}
        others = [a for a in AXES if a != axis]
        u, v, w = P[axis] * 1e-3, P[others[0]] * 1e-3, P[others[1]] * 1e-3
        rho = np.hypot(v, w)
        current = loop.get("turns", 1) * loop.get("sign", 1) * scale
        b_rho, b_ax = _loop_field(rho, u, loop["radius_mm"] * 1e-3, current)
        B[AXES.index(axis)] += b_ax
        with np.errstate(invalid="ignore", divide="ignore"):
            B[AXES.index(others[0])] += np.where(rho > 0, b_rho * v / rho, 0.0)
            B[AXES.index(others[1])] += np.where(rho > 0, b_rho * w / rho, 0.0)
    return FieldMap(*axes, B, name)


def generated_rings(axes, rings, name="b0", n_r=6, n_ax=6, n_phi=96):
    """Field (T) of permanent-magnet rings, as a sum of point dipoles over each ring's volume: an Aubert pair, a
    Halbach-free ring magnet, any stack of rings. Valid away from the magnet material (in the bore).

    Each ring: {inner_radius_mm, outer_radius_mm, thickness_mm, axis, centre_mm: {x, y, z}, br_t (remanence),
    magnetization: axial | -axial | radial_out | radial_in}.
    An Aubert pair: two coaxial rings, one magnetized radially outward and one inward, which gives a field along
    the axis at the centre; their spacing sets the homogeneity.
    """
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]) * 1e-3
    B = np.zeros_like(pts)
    for ring in rings:
        ai = AXES.index(ring["axis"])
        u, v = [i for i in range(3) if i != ai]
        r1, r2 = ring["inner_radius_mm"] * 1e-3, ring["outer_radius_mm"] * 1e-3
        t = ring["thickness_mm"] * 1e-3
        c = ring.get("centre_mm", {}) or {}
        centre = np.array([c.get(a, 0) for a in AXES]) * 1e-3
        dr, dz, dphi = (r2 - r1) / n_r, t / n_ax, 2 * np.pi / n_phi
        R, H, PHI = np.meshgrid(r1 + (np.arange(n_r) + 0.5) * dr, -t / 2 + (np.arange(n_ax) + 0.5) * dz,
                                (np.arange(n_phi) + 0.5) * dphi, indexing="ij")
        R, H, PHI = R.ravel(), H.ravel(), PHI.ravel()
        src = np.tile(centre, (len(R), 1))
        src[:, ai] += H
        src[:, u] += R * np.cos(PHI)
        src[:, v] += R * np.sin(PHI)
        moment = ring["br_t"] / MU0 * R * dr * dz * dphi          # A m^2 per element (M = Br / mu0)
        m = np.zeros_like(src)
        kind = ring["magnetization"]
        if kind in ("axial", "-axial"):
            m[:, ai] = moment * (1 if kind == "axial" else -1)
        elif kind in ("radial_out", "radial_in"):
            sgn = 1 if kind == "radial_out" else -1
            m[:, u], m[:, v] = sgn * moment * np.cos(PHI), sgn * moment * np.sin(PHI)
        else:
            raise ValueError(f"unknown magnetization {kind!r}")
        for k in range(0, len(pts), 4096):                        # chunked: points x dipoles
            d = pts[k:k + 4096, None, :] - src[None, :, :]
            r = np.linalg.norm(d, axis=-1)
            mr = np.einsum("pnc,nc->pn", d, m)
            B[k:k + 4096] += MU0 / (4 * np.pi) * np.einsum(
                "pnc->pc", 3 * d * (mr / r**5)[..., None] - m[None] / (r**3)[..., None])
    return FieldMap(*axes, B.T.reshape((3,) + X.shape), name)


def solenoid_loops(turns, radius_mm, length_mm, axis="x", centre_mm=None):
    """An N-turn solenoid as loops, evenly spaced over its length."""
    c = dict(centre_mm or {})
    zs = np.linspace(-length_mm / 2, length_mm / 2, int(turns)) if turns > 1 else [0.0]
    return [{"radius_mm": radius_mm, "axis": axis, "centre_mm": {**c, axis: c.get(axis, 0) + z}} for z in zs]


def generated_solenoid(axes, turns, radius_mm, length_mm, axis="x"):
    """Field per amp (T/A) of an N-turn solenoid centred at the origin, its axis along `axis` (Biot-Savart)."""
    return generated_loops(axes, solenoid_loops(turns, radius_mm, length_mm, axis), "coil")


# ---- From maps to what the kernel needs ---------------------------------------------------------------------------

@dataclass
class VoxelFields:
    off_resonance_hz: np.ndarray     # (n,) f(r) - f_rf
    encoding_m: np.ndarray           # (n, 3) encoding coordinate per gradient channel x, y, z
    b1_relative: np.ndarray          # (n,) transmit flip-angle scale
    coil_t_per_a: np.ndarray         # (n,) receive sensitivity |B1_perp| per amp
    b0_tesla: np.ndarray             # (n,) |B0|
    coil_centre_t_per_a: float


def voxel_fields(maps, points_mm, f_rf_hz):
    """Sample the four maps at the voxel centres (n, 3, in mm) and turn them into per-voxel kernel inputs."""
    b0 = maps["b0"].sample(points_mm)
    b0_mag = np.linalg.norm(b0, axis=0)
    b0_hat = b0 / np.maximum(b0_mag, 1e-30)

    encoding = np.zeros((len(points_mm), 3))
    for i, a in enumerate(AXES):
        g = maps[f"g{a}"].sample(points_mm)
        # A gradient changes the precession frequency through its field along the local B0 (the secular part).
        g_par = np.einsum("cn,cn->n", g, b0_hat)
        encoding[:, i] = GAMMA_BAR * g_par / ENCODING_HZ_PER_M

    coil = maps["coil"].sample(points_mm)
    coil_par = np.einsum("cn,cn->n", coil, b0_hat)
    coil_perp = np.linalg.norm(coil - coil_par * b0_hat, axis=0)
    centre = maps["coil"].sample(np.zeros((1, 3)))
    b0_c_vec = maps["b0"].sample(np.zeros((1, 3)))[:, 0]
    b0_c = b0_c_vec / np.linalg.norm(b0_c_vec)
    centre_perp = float(np.linalg.norm(centre[:, 0] - np.dot(centre[:, 0], b0_c) * b0_c))

    return VoxelFields(GAMMA_BAR * b0_mag - f_rf_hz, encoding, coil_perp / centre_perp, coil_perp, b0_mag,
                       centre_perp)


def build_maps(design, extent_mm, step_mm, files_root=Path(".")):
    """The four maps for a design: from files where given, generated from the design otherwise."""
    fdesign = design["fields"]
    direction = fdesign["b0_direction"]
    axes = grid_axes(extent_mm, step_mm)
    maps = {}

    b0 = fdesign["b0"]
    if b0.get("file"):
        maps["b0"] = load_map(files_root / b0["file"], "b0", scalar_axis=direction)
    elif b0.get("rings"):           # the magnet by its dimensions, e.g. an Aubert ring pair
        maps["b0"] = generated_rings(axes, b0["rings"])
    else:
        maps["b0"] = generated_b0(axes, design["magnet"]["larmor_hz"], b0["offset_hz"], b0["linear_hz_per_mm"],
                                  b0["quadratic_hz_per_mm2"], direction)
    if b0.get("match_larmor"):      # retune: scale the map so its centre sits at magnet.larmor_hz (Br varies)
        f0 = GAMMA_BAR * np.linalg.norm(maps["b0"].sample(np.zeros((1, 3)))[:, 0])
        maps["b0"] = FieldMap(maps["b0"].x_mm, maps["b0"].y_mm, maps["b0"].z_mm,
                              maps["b0"].B * design["magnet"]["larmor_hz"] / f0, "b0")
    for a in AXES:
        g = fdesign[f"g{a}"]
        if g.get("file"):
            maps[f"g{a}"] = load_map(files_root / g["file"], f"g{a}", scalar_axis=direction)
        elif g.get("loops"):        # a real coil (e.g. an anti-Helmholtz pair), driven at amps_per_drive
            maps[f"g{a}"] = generated_loops(axes, g["loops"], f"g{a}", g["amps_per_drive"])
        else:
            maps[f"g{a}"] = generated_gradient(axes, a, g["hz_per_mm_per_drive"], direction)
    c = fdesign["coil"]
    if c.get("file"):
        maps["coil"] = load_map(files_root / c["file"], "coil")
    elif c.get("loops"):
        maps["coil"] = generated_loops(axes, c["loops"], "coil")
    else:
        s = c["solenoid"]
        maps["coil"] = generated_solenoid(axes, s["turns"], s["radius_mm"], s["length_mm"], s["axis"])
    return maps
