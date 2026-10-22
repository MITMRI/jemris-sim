# Pulseq (.seq) → JEMRIS (.xml + .h5 waveforms)
#
# Adapted from py2jemris seq2xml.py (Gehua Tong / imr-framework),
# updated for modern PyPulseq (1.4+/1.5) Sequence APIs.
# https://github.com/imr-framework/py2jemris
#
# Caveat: Pulseq is unrolled blocks; JEMRIS XML will be flat (no recovered
# high-level loops). Sims may be slower than hand-authored looped XML.

from __future__ import annotations

import os
from math import pi
from pathlib import Path
from typing import Dict, Tuple, Union
import xml.etree.ElementTree as ET

import h5py
import numpy as np
from pypulseq.Sequence.sequence import Sequence

# Pulseq → JEMRIS unit conversions
RF_CONST = 2 * pi / 1000          # Hz → rad/ms
G_CONST = 2 * pi / 1e6            # Hz/m → (rad/ms)/mm
SLEW_CONST = G_CONST / 1000       # Hz/(m·s) → (rad/ms)/(mm·ms)
GA_CONST = 2 * pi / 1000          # 1/m → 2π/mm
SEC2MS = 1000
RAD2DEG = 180 / pi
FREQ_CONST = 2 * pi / 1000        # Hz → rad/ms


def convert_seq_file(
    seq_path: Union[str, Path],
    out_dir: Union[str, Path],
    seq_name: str | None = None,
) -> Path:
    """Load a .seq file and write JEMRIS XML (+ RF/grad HDF5) into out_dir."""
    seq_path = Path(seq_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if seq_name is None:
        seq_name = seq_path.stem

    seq = Sequence()
    seq.read(str(seq_path))
    _, xml_path = seq2xml(seq, seq_name=seq_name, out_folder=str(out_dir))
    return Path(xml_path)


def seq2xml(seq: Sequence, seq_name: str, out_folder: str) -> Tuple[ET.ElementTree, str]:
    """
    Convert a loaded PyPulseq Sequence into JEMRIS XML.

    Returns (ElementTree, path_to_xml).
    """
    out_folder = str(out_folder)
    os.makedirs(out_folder, exist_ok=True)

    root = ET.Element("Parameters")
    root.set("GradMaxAmpl", str(seq.system.max_grad * G_CONST))
    root.set("GradSlewRate", str(seq.system.max_slew * SLEW_CONST))
    root.set("Name", seq_name)

    concat = ET.SubElement(root, "ConcatSequence")
    concat.set("Name", "C0")

    rf_shapes = _save_rf_library(seq, out_folder)
    grad_shapes = _save_grad_library(seq, out_folder)

    rf_i = grad_i = delay_i = adc_i = 1
    n_blocks = len(seq.block_events)

    for block_ind in range(1, n_blocks + 1):
        blk = seq.get_block(block_ind)
        has_adc = blk.adc is not None
        adc_added = False

        events = []
        if blk.rf is not None:
            events.append("rf")
        for ax in ("gx", "gy", "gz"):
            if getattr(blk, ax, None) is not None:
                events.append(ax)
        if not events and has_adc:
            events.append("adc_only")
        if not events and not has_adc:
            # pure delay
            delay_atom = ET.SubElement(concat, "DELAYATOMICSEQUENCE")
            delay_atom.set("Name", f"D{delay_i}")
            delay_i += 1
            delay_atom.set("Delay", str(float(blk.block_duration) * SEC2MS))
            delay_atom.set("DelayType", "B2E")
            continue

        atom = ET.SubElement(concat, "ATOMICSEQUENCE")
        atom.set("Name", f"C{block_ind}")

        if "rf" in events:
            rf = blk.rf
            rf_atom = ET.SubElement(atom, "EXTERNALRFPULSE")
            rf_atom.set("Name", f"R{rf_i}")
            rf_i += 1
            rf_atom.set("InitialDelay", str(float(rf.delay) * SEC2MS))
            rf_atom.set("InitialPhase", str(float(rf.phase_offset) * RAD2DEG))
            freq = float(getattr(rf, "freq_offset", 0.0) or 0.0)
            rf_atom.set("Frequency", str(freq * FREQ_CONST))
            rf_id = int(seq.block_events[block_ind][1])
            rf_atom.set("Filename", rf_shapes[rf_id])
            rf_atom.set("Scale", "1")
            rf_atom.set("Interpolate", "0")

        gnames_map = {"gx": 2, "gy": 3, "gz": 4}
        for key in ("gx", "gy", "gz"):
            g = getattr(blk, key, None)
            if g is None:
                continue

            if g.type == "trap":
                if float(g.amplitude) != 0.0:
                    g_atom = ET.SubElement(atom, "TRAPGRADPULSE")
                    g_atom.set("Name", f"G{grad_i}")
                    grad_i += 1
                    if float(g.flat_time) > 0:
                        g_atom.set("FlatTopArea", str(float(g.flat_area) * GA_CONST))
                        g_atom.set("FlatTopTime", str(float(g.flat_time) * SEC2MS))
                    else:
                        g_atom.set("MaxAmpl", str(abs(float(g.amplitude) * G_CONST)))
                        area = 0.5 * (float(g.rise_time) + float(g.fall_time)) * float(g.amplitude)
                        g_atom.set("Area", str(area * GA_CONST))
                    rise = max(float(g.rise_time), 1e-12)
                    fall = max(float(g.fall_time), 1e-12)
                    g_atom.set("SlewRate", str((float(g.amplitude) * G_CONST) / (fall * SEC2MS)))
                    g_atom.set("Asymmetric", str(fall / rise))
                    g_atom.set("Axis", key.upper())
                    g_atom.set("InitialDelay", str(float(g.delay) * SEC2MS))
            elif g.type == "grad":
                g_id = int(seq.block_events[block_ind][gnames_map[key]])
                g_atom = ET.SubElement(atom, "EXTERNALGRADPULSE")
                g_atom.set("Name", f"G{grad_i}")
                grad_i += 1
                g_atom.set("Axis", key.upper())
                g_atom.set("Filename", grad_shapes[g_id])
                g_atom.set("Scale", "1")
                g_atom.set("Interpolate", "0")
                g_atom.set("InitialDelay", str(float(g.delay) * SEC2MS))
            else:
                raise ValueError(f'Unsupported gradient type "{g.type}"')

            if has_adc and not adc_added:
                adc = blk.adc
                dwell = float(adc.dwell) * SEC2MS
                adc_delay = float(adc.delay) * SEC2MS
                nro = int(adc.num_samples)
                gzero = ET.SubElement(atom, "TRAPGRADPULSE")
                gzero.set("Name", f"S{adc_i}")
                adc_i += 1
                gzero.set("ADCs", str(nro))
                gzero.set("FlatTopTime", str(dwell * nro))
                gzero.set("FlatTopArea", "0")
                gzero.set("InitialDelay", str(adc_delay))
                adc_added = True

        if has_adc and not adc_added:
            # ADC-only block (no gradient to hang ADCs on)
            adc = blk.adc
            dwell = float(adc.dwell) * SEC2MS
            adc_delay = float(adc.delay) * SEC2MS
            nro = int(adc.num_samples)
            gzero = ET.SubElement(atom, "TRAPGRADPULSE")
            gzero.set("Name", f"S{adc_i}")
            adc_i += 1
            gzero.set("ADCs", str(nro))
            gzero.set("FlatTopTime", str(dwell * nro))
            gzero.set("FlatTopArea", "0")
            gzero.set("InitialDelay", str(adc_delay))

    tree = ET.ElementTree(root)
    xml_path = os.path.join(out_folder, f"{seq_name}.xml")
    _indent(root)
    tree.write(xml_path, encoding="utf-8", xml_declaration=True)
    return tree, xml_path


def _save_rf_library(seq: Sequence, out_folder: str) -> Dict[int, str]:
    paths: Dict[int, str] = {}
    for rf_id in list(seq.rf_library.data.keys()):
        rf_id = int(rf_id)
        partial = f"rf_{rf_id}.h5"
        full = os.path.join(out_folder, partial)
        rf = seq.rf_from_lib_data(seq.rf_library.data[rf_id])
        times = np.asarray(rf.t, dtype=float)
        mag = np.abs(np.asarray(rf.signal, dtype=complex))
        phase = np.angle(np.asarray(rf.signal, dtype=complex))
        times = times - times[0]
        n = len(mag)
        with h5py.File(full, "w") as f:
            ds = f.create_dataset("extpulse", (3, n), dtype="f")
            ds[0, :] = times * SEC2MS
            ds[1, :] = mag * RF_CONST
            ds[2, :] = phase
        paths[rf_id] = partial
    return paths


def _save_grad_library(seq: Sequence, out_folder: str) -> Dict[int, str]:
    paths: Dict[int, str] = {}
    processed = set()
    for nb in range(1, len(seq.block_events) + 1):
        gx_ind, gy_ind, gz_ind = [int(x) for x in seq.block_events[nb][2:5]]
        blk = seq.get_block(nb)
        for axis_ind, g_ind in enumerate((gx_ind, gy_ind, gz_ind)):
            if g_ind == 0 or g_ind in processed:
                continue
            # Arbitrary grads have 3 library entries in classic format; traps differ.
            gdata = seq.grad_library.data.get(g_ind)
            if gdata is None or len(gdata) != 3:
                continue
            axis = ("gx", "gy", "gz")[axis_ind]
            g = getattr(blk, axis, None)
            if g is None or getattr(g, "type", None) != "grad":
                continue
            partial = f"grad_{g_ind}.h5"
            full = os.path.join(out_folder, partial)
            t_points = np.asarray(g.t, dtype=float)
            g_shape = np.asarray(g.waveform, dtype=float)
            n = len(t_points)
            with h5py.File(full, "w") as f:
                ds = f.create_dataset("extpulse", (2, n), dtype="f")
                ds[0, :] = t_points * SEC2MS
                ds[1, :] = g_shape * G_CONST
            paths[g_ind] = partial
            processed.add(g_ind)
    return paths


def _indent(elem: ET.Element, level: int = 0) -> None:
    pad = "\n" + "  " * level
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = pad + "  "
        for child in elem:
            _indent(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = pad
    if level and (not elem.tail or not elem.tail.strip()):
        elem.tail = pad
