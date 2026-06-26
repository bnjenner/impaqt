#!/bin/bash
set -euxo pipefail

# conda-build sets CC/CXX to the conda toolchain, PREFIX to the install root,
# CPU_COUNT to the available cores, and PKG_CONFIG_PATH to $PREFIX/lib/pkgconfig
# (so the host bamtools is found). Use the conda-provided bamtools and skip the
# unit tests (and thus GoogleTest) so the build needs no network.
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DUSE_SYSTEM_BAMTOOLS=ON \
    -DIMPAQT_BUILD_TESTS=OFF

cmake --build build -j"${CPU_COUNT}"
cmake --install build
