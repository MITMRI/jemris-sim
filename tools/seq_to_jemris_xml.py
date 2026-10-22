#!/usr/bin/env python3
"""
Convert a Pulseq / PyPulseq sequence (.seq) into JEMRIS XML + waveform HDF5.

Usage:
  python tools/seq_to_jemris_xml.py path/to/seq.seq -o out_dir
  ./tools/seq_to_jemris_xml.py examples/sequences/epi_example.seq -o /tmp/jemris_seq

Then point a simu.xml sequence uri at out_dir/<name>.xml (waveforms sit beside it).

Requires: Python >= 3.10, pypulseq, h5py, numpy
  pip install -r tools/requirements-seq.txt
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running as a script from repo root or tools/
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from seq2xml import convert_seq_file  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(
        description="Convert Pulseq .seq → JEMRIS XML (+ RF/grad .h5 waveforms)"
    )
    p.add_argument("seq", type=Path, help="Input Pulseq .seq file")
    p.add_argument(
        "-o",
        "--out-dir",
        type=Path,
        default=Path("jemris_seq_out"),
        help="Output directory for .xml and .h5 files (default: ./jemris_seq_out)",
    )
    p.add_argument(
        "-n",
        "--name",
        default=None,
        help="Base name for output XML (default: input stem)",
    )
    args = p.parse_args()

    if not args.seq.is_file():
        print(f"error: sequence file not found: {args.seq}", file=sys.stderr)
        return 1

    xml_path = convert_seq_file(args.seq, args.out_dir, seq_name=args.name)
    print(f"Wrote JEMRIS sequence: {xml_path}")
    print(f"Waveform HDF5 files (if any) are in: {args.out_dir.resolve()}")
    print()
    print("Wire it into a simu.xml like:")
    print(f'  <sequence name="{xml_path.stem}" uri="{xml_path.name}"/>')
    print("(run jemris from that output directory so relative .h5 paths resolve)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
