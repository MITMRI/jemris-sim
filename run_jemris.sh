#!/usr/bin/env bash
# Launcher for jemris (sets library paths for conda-built deps)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export PREFIX="$ROOT/deps/install"
export CONDA="$ROOT/deps/conda"
export LD_LIBRARY_PATH="$CONDA/lib:$PREFIX/lib:$PREFIX/lib64:${LD_LIBRARY_PATH:-}"
export LD_PRELOAD="$CONDA/lib/libstdc++.so.6${LD_PRELOAD:+:$LD_PRELOAD}"
export PATH="$ROOT/build/src:$PREFIX/bin:$PATH"

# Prefer cluster MPI if available
module load gcc/12.2.0 openmpi/4.1.4 2>/dev/null || true

exec jemris "$@"
