"""The design: one nested dictionary, loaded from YAML files and `key.path=value` overrides.

    design = load_design(["design/default.yaml", "design/examples/pure_water.yaml"], ["gradients.axis=y"])
    design["gradients"]["strength_hz_per_mm"]

Later files and overrides replace only the keys they name, so an example file can be three lines long.
"""

import copy
import re
from pathlib import Path

import yaml

# PyYAML only reads exponents with a sign (4.567e+6) as numbers; 4.567e6 comes back as text. Accept both.
_NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


def _numbers(node):
    """Turn number-like strings anywhere in a loaded design into floats."""
    if isinstance(node, dict):
        return {k: _numbers(v) for k, v in node.items()}
    if isinstance(node, list):
        return [_numbers(v) for v in node]
    if isinstance(node, str) and _NUMBER.match(node.strip()):
        return float(node)
    return node

PACKAGE_ROOT = Path(__file__).resolve().parents[1]          # sim6480/ (holds design/, tests/)
DEFAULT_DESIGN = PACKAGE_ROOT / "design" / "default.yaml"

REQUIRED = {
    "magnet": ["larmor_hz", "temperature_c", "t2dash_s"],
    "fields": ["b0_direction", "b0", "gx", "gy", "gz", "coil"],
    "rf": ["frequency_offset_hz", "flip90_deg", "flip180_deg", "phase_from_timing"],
    "gradients": ["axis", "drive", "sign"],
    "receive": ["channel", "noise_nv_per_rthz", "preamp_db", "lna_db", "pga_db", "lpf_mhz", "digital_db",
                "adc_vpp", "adc_bits", "fs_hz", "rf_leak_counts"],
    "sequence": ["file", "capture_blocks", "write_ns", "out_rate_hz", "echo_row_samples"],
    "phantom": ["stl", "liquid", "rotation_deg", "axis", "centre_mm", "voxel_mm"],
    "liquids": [],
    "simulation": ["kernel", "adc_dwell_us", "voxel_shape", "max_states", "min_state_mag", "field_map_step_mm",
                   "seed"],
    "repos": ["pulse_sequence", "signal_processing"],
}


def deep_merge(base, extra):
    """A copy of base with extra's keys replacing it, recursing into nested dictionaries."""
    out = copy.deepcopy(base)
    for key, value in (extra or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def apply_override(design, assignment):
    """Apply one `a.b.c=value` override (value parsed as YAML: numbers, true/false, lists all work)."""
    if "=" not in assignment:
        raise ValueError(f"override {assignment!r} must look like key.path=value")
    path, raw = assignment.split("=", 1)
    keys = path.strip().split(".")
    node = design
    for key in keys[:-1]:
        if not isinstance(node.get(key), dict):
            raise KeyError(f"{path}: {key!r} isn't a section of the design")
        node = node[key]
    if keys[-1] not in node:
        raise KeyError(f"{path}: no such setting (check the spelling against design/default.yaml)")
    node[keys[-1]] = _numbers(yaml.safe_load(raw))
    return design


def validate(design):
    """Check every required setting is there and the choices are valid; returns the design."""
    for section, keys in REQUIRED.items():
        if section not in design:
            raise KeyError(f"the design has no `{section}` section")
        missing = [k for k in keys if k not in design[section]]
        if missing:
            raise KeyError(f"`{section}` is missing {missing}")
    if design["gradients"]["axis"] not in ("x", "y", "z"):
        raise ValueError("gradients.axis must be x, y or z")
    if design["phantom"]["axis"] not in ("x", "y", "z") or design["phantom"]["axis"] == design["gradients"]["axis"]:
        raise ValueError("phantom.axis must be x, y or z, and differ from the readout axis")
    if design["phantom"]["liquid"] not in design["liquids"]:
        raise KeyError(f"phantom.liquid {design['phantom']['liquid']!r} isn't one of {sorted(design['liquids'])}")
    if design["fields"]["b0_direction"] not in ("x", "y", "z"):
        raise ValueError("fields.b0_direction must be x, y or z")
    if design["simulation"]["voxel_shape"] not in ("sinc", "box", "point"):
        raise ValueError("simulation.voxel_shape must be sinc, box or point")
    if design["receive"]["adc_bits"] != 14:
        raise ValueError("receive.adc_bits: the FPGA packing carries 14-bit samples")
    if design["receive"]["channel"] not in ("ch1", "ch2"):
        raise ValueError("receive.channel must be ch1 or ch2")
    if design["simulation"]["kernel"] not in ("mrzero",):
        raise ValueError("simulation.kernel: only mrzero runs here (JEMRIS: use `sim6480 jemris` to write its inputs)")
    return design


def load_design(files=None, overrides=None):
    """The default design, with each file in `files` layered on top, then each `key=value` override."""
    files = list(files or [])
    design = _numbers(yaml.safe_load(DEFAULT_DESIGN.read_text()))
    for f in files:
        if Path(f).resolve() == DEFAULT_DESIGN:
            continue
        design = deep_merge(design, _numbers(yaml.safe_load(Path(f).read_text())))
    for assignment in overrides or []:
        apply_override(design, assignment)
    return validate(design)


def repo_path(design, name):
    """Absolute path of a sibling repo named in design['repos'] (relative paths are from the design folder)."""
    p = Path(design["repos"][name])
    return p if p.is_absolute() else (PACKAGE_ROOT / "design" / p).resolve()
