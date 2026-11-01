"""The sibling MITMRI repos this builds on, imported from the paths in design['repos'].

    pulse_sequence     6480-pulse-sequence: the firmware's timeline (pulseseq.timeline) and the phantom STLs
                       (pulseseq.phantoms)
    signal_processing  6480-signal-processing: FPGA packing, unpacking, downconversion, echo rows, 1-D images

Using their code (not copies of it) keeps the simulation in step with them: change the firmware there and the
simulated sequence changes too.
"""

import importlib
import sys

from sim6480.config import repo_path


def _add(design, name, marker):
    path = repo_path(design, name)
    if not (path / marker).exists():
        raise FileNotFoundError(f"repos.{name}: no {marker} under {path}; clone it there or fix the path in the design")
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
    return path


def pulseseq(design):
    """The pulseseq package from 6480-pulse-sequence (timeline, phantoms)."""
    _add(design, "pulse_sequence", "pulseseq")
    return importlib.import_module("pulseseq.timeline"), importlib.import_module("pulseseq.phantoms")


def sigproc(design):
    """The modules used from 6480-signal-processing, as a simple namespace."""
    _add(design, "signal_processing", "sigproc")

    class SP:
        transfer = importlib.import_module("sigproc.stage3_transfer")
        unpack = importlib.import_module("sigproc.stage4_unpack")
        downconvert = importlib.import_module("sigproc.downconvert")
        segment = importlib.import_module("sigproc.segment")
        recon = importlib.import_module("sigproc.recon")
        files = importlib.import_module("sigproc.files")

    return SP
