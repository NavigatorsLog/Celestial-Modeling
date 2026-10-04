#!/usr/bin/env bash
# setup.sh  [CL draft, 4 October 2026]: fetch the pinned toolchain and build formation_run.
# Pins: REBOUND git tag 4.4.10; fragmentation module github.com/annacrnn/rebound_fragmentation commit d67ab2b.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p build && cd build
[ -d rebound ] || git clone -q --depth 1 --branch 4.4.10 https://github.com/hannorein/rebound.git rebound
if [ ! -d rebound_fragmentation ]; then
  git clone -q https://github.com/annacrnn/rebound_fragmentation.git rebound_fragmentation
  git -C rebound_fragmentation checkout -q d67ab2b
fi
# the module is used unmodified except for its include path
sed 's|#include "rebound.h"|#include "rebound.h" /* path set by -I */|' rebound_fragmentation/fragmentation.c > frag.c
make -s -C rebound/src librebound.so OPENGL=0 SERVER=0 >/dev/null
cp rebound/src/librebound.so .
gcc -O3 -I rebound/src -I . -Wl,-rpath,'$ORIGIN' "$ROOT/src/formation_run.c" -L. -lrebound -lm -o formation_run
sha256sum frag.c "$ROOT/src/formation_run.c" formation_run | tee build_hashes.txt
echo "built: $ROOT/build/formation_run"
