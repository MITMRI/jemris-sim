"""Checks of sim6480 against things known independently of it.

    cd sim6480 && python -m pytest tests -q        (needs the sibling repos 6480-pulse-sequence, 6480-signal-processing)
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sim6480 import fieldmaps, physics                                   # noqa: E402
from sim6480.config import PACKAGE_ROOT, load_design                     # noqa: E402

EXAMPLES = PACKAGE_ROOT / "design" / "examples"


# ---- design ----------------------------------------------------------------------------------------------------

def test_overrides_and_numbers():
    d = load_design(overrides=["gradients.axis=y", "receive.pga_db=24"])
    assert d["gradients"]["axis"] == "y" and d["receive"]["pga_db"] == 24
    assert isinstance(d["magnet"]["larmor_hz"], float)          # YAML reads 4.567e6 as a string; we coerce


def test_unknown_override_is_an_error():
    with pytest.raises(KeyError):
        load_design(overrides=["gradients.axes=y"])


def test_examples_load():
    for f in EXAMPLES.glob("*.yaml"):
        load_design([PACKAGE_ROOT / "design" / "default.yaml", f])


# ---- field maps against closed forms ---------------------------------------------------------------------------

def test_solenoid_centre_matches_formula():
    axes = fieldmaps.grid_axes(2, 1.0)
    m = fieldmaps.generated_solenoid(axes, 200, 9.5, 14.5, "x")
    b = m.sample(np.zeros((1, 3)))[0, 0]
    assert b == pytest.approx(physics.solenoid_field_per_amp(200, 9.5e-3, 14.5e-3), rel=0.01)


def test_maxwell_pair_gradient():
    """Anti-Helmholtz (Maxwell) pair: matches the on-axis field of two opposed loops, mu0 N I R^2 / 2 (R^2 + z^2)^1.5."""
    R, d, N = 15e-3, 26e-3, 20
    loops = [{"radius_mm": 15, "axis": "z", "centre_mm": {"z": 13}, "turns": N, "sign": 1},
             {"radius_mm": 15, "axis": "z", "centre_mm": {"z": -13}, "turns": N, "sign": -1}]
    m = fieldmaps.generated_loops(fieldmaps.grid_axes(2, 0.5), loops, "gz")
    z = np.array([-0.5, 0.5])
    bz = m.sample(np.column_stack([0 * z, 0 * z, z]))[2]
    on_axis = lambda zz: sum(s * physics.MU0 * N * R**2 / (2 * (R**2 + (zz * 1e-3 - c) ** 2) ** 1.5)
                             for c, s in ((d / 2, 1), (-d / 2, -1)))
    assert bz == pytest.approx([on_axis(v) for v in z], rel=1e-3)


def test_axial_ring_matches_closed_form():
    ring = dict(inner_radius_mm=20, outer_radius_mm=50, thickness_mm=20, axis="z", centre_mm={}, br_t=1.3,
                magnetization="axial")
    m = fieldmaps.generated_rings(fieldmaps.grid_axes(4, 2.0), [ring])
    z = np.array([0.0, 4.0])
    f = lambda z, r: (z + 10) / np.hypot(z + 10, r) - (z - 10) / np.hypot(z - 10, r)
    exact = 1.3 / 2 * (f(z, 50) - f(z, 20))
    assert m.sample(np.column_stack([0 * z, 0 * z, z]))[2] == pytest.approx(exact, rel=0.005)


def test_aubert_pair_field_is_axial():
    ring = dict(inner_radius_mm=20, outer_radius_mm=50, thickness_mm=20, axis="z", br_t=1.3)
    pair = [dict(ring, magnetization="radial_out", centre_mm={"z": 30}),
            dict(ring, magnetization="radial_in", centre_mm={"z": -30})]
    b = fieldmaps.generated_rings(fieldmaps.grid_axes(2, 2.0), pair).sample(np.zeros((1, 3)))[:, 0]
    assert abs(b[2]) > 0.1 and np.hypot(b[0], b[1]) < 1e-6 * abs(b[2])


def test_6480_magnet_map():
    """The shipped Aubert-pair map: ~107 mT along +z at the centre, and matches its design file."""
    m = fieldmaps.load_map(PACKAGE_ROOT / "design" / "maps" / "b0_aubert.npz")
    b = m.sample(np.zeros((1, 3)))[:, 0]
    assert b[2] == pytest.approx(physics.b0_tesla(4.567e6), rel=0.01)
    assert np.hypot(b[0], b[1]) < 1e-6 * b[2]
    d = load_design([PACKAGE_ROOT / "design" / "default.yaml", EXAMPLES / "aubert_6480.yaml"])
    p = np.array([[3.0, -2.0, 4.0]])
    gen = fieldmaps.generated_rings((np.array([2., 3, 4]), np.array([-3., -2, -1]), np.array([3., 4, 5])),
                                    d["fields"]["b0"]["rings"])
    np.testing.assert_allclose(gen.sample(p), m.sample(p), rtol=2e-3)


def test_map_file_round_trip(tmp_path):
    axes = fieldmaps.grid_axes(3, 1.0)
    m = fieldmaps.generated_gradient(axes, "y", 2500, "z")
    fieldmaps.save_map(tmp_path / "gy.npz", m)
    back = fieldmaps.load_map(tmp_path / "gy.npz")
    np.testing.assert_allclose(back.B, m.B)


def test_encoding_follows_the_map():
    """A gradient map 1.5x stronger encodes positions 1.5x further out (what makes the maps architecture-free)."""
    axes = fieldmaps.grid_axes(3, 1.0)
    maps = {"b0": fieldmaps.generated_b0(axes, 4.567e6, 0, {}, {}, "z"),
            "gx": fieldmaps.generated_gradient(axes, "x", 1500, "z"),
            "gy": fieldmaps.generated_gradient(axes, "y", 1000, "z"),
            "gz": fieldmaps.generated_gradient(axes, "z", 1000, "z"),
            "coil": fieldmaps.generated_solenoid(axes, 20, 9.5, 14.5, "x")}
    pts = np.array([[2.0, 2.0, 2.0]])
    enc = fieldmaps.voxel_fields(maps, pts, 4.567e6).encoding_m[0]
    assert enc[0] == pytest.approx(1.5 * enc[1])


# ---- scaling -----------------------------------------------------------------------------------------------------

def test_magnetization_of_water():
    """Curie law: water at 107 mT and 22 °C has M0 ~ 3.5e-4 A/m."""
    assert physics.magnetization(physics.b0_tesla(4.567e6), 22) == pytest.approx(3.5e-4, rel=0.02)


# ---- end to end --------------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def shot_z():
    from sim6480.simulate import simulate
    return simulate(load_design(overrides=["gradients.axis=z", "receive.noise_nv_per_rthz=0"]))


def test_echo_voltage_scale(shot_z):
    """The echo is about the EMF of the water fully tipped: (nV per mm^3) x water volume, a few µV."""
    expected = physics.summary(shot_z.design)["EMF per mm^3 of water (nV)"] * shot_z.phantom.water_mm3 * 1e-3
    assert 0.5 * expected < shot_z.peak_coil_uv < 1.1 * expected


def test_three_slices_where_they_are(shot_z):
    """Three separate regions above half height, spanning the phantom's water to within an image pixel."""
    pos, img = shot_z.process.positions_mm, np.abs(shot_z.process.image)
    above = img > 0.5 * img.max()
    assert np.count_nonzero(np.diff(above.astype(int)) == 1) == 3
    tx, tw = shot_z.phantom.profile("z")
    th = tx[tw > tw.max() / 2]
    assert pos[above].min() == pytest.approx(th.min(), abs=0.6)
    assert pos[above].max() == pytest.approx(th.max(), abs=0.6)


def test_image_moves_with_the_phantom():
    from sim6480.simulate import simulate
    s = simulate(load_design(overrides=["phantom.centre_mm.z=2", "receive.noise_nv_per_rthz=0"]))
    pos, img = s.process.positions_mm, np.abs(s.process.image)
    h = pos[img > img.max() / 2]
    assert (h.min() + h.max()) / 2 == pytest.approx(2.0, abs=0.6)       # one image pixel


def test_jemris_bundle(tmp_path, shot_z):
    import h5py
    from sim6480 import jemris_export
    info = jemris_export.export(shot_z.design, tmp_path, shot=shot_z)
    for f in ("simu.xml", "sample.h5", "sequence.xml", "coils.xml", "export.json"):
        assert (tmp_path / f).exists()
    nx, ny, nz = shot_z.phantom.pd.shape
    with h5py.File(tmp_path / "sample.h5") as h:
        assert h["/sample/data"].shape == (nz, ny, nx, 5)
    assert info["adc_samples"] == len(shot_z.plan.adc_times_s)


def test_any_pulseq_sequence(shot_z):
    """The firmware sequence, run through the generic Pulseq path, gives the same coil voltage."""
    from sim6480.simulate import simulate_sequence
    s = simulate_sequence(shot_z.design, shot_z.plan.seq, maps=shot_z.maps, phantom=shot_z.phantom)
    np.testing.assert_allclose(s.t_s, shot_z.plan.adc_times_s, atol=1e-9)
    np.testing.assert_allclose(s.coil_v, shot_z.coil_v, rtol=1e-4, atol=1e-4 * np.abs(shot_z.coil_v).max())
