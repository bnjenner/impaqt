#!/usr/bin/env bash
#
# Build (and optionally install) impaqt.
#
# Third-party dependencies (bamtools, GoogleTest) are fetched automatically by
# CMake on the first configure, so an internet connection is required the first
# time you build.
#
# Usage:
#   ./install.sh                      # configure + build into ./build
#   ./install.sh --install            # ...then install (uses sudo for system prefix)
#   ./install.sh --install --prefix ~/.local   # install to a custom prefix (no sudo)
#   ./install.sh --build-type Debug --jobs 4
#
set -euo pipefail

# Repo root = directory containing this script
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${ROOT}/build"
BUILD_TYPE="Release"
JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
DO_INSTALL=0
PREFIX=""

usage() {
    sed -n '3,14p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --install)            DO_INSTALL=1; shift ;;
        --prefix)             PREFIX="${2:?--prefix needs a directory}"; shift 2 ;;
        --build-type)         BUILD_TYPE="${2:?--build-type needs a value}"; shift 2 ;;
        --jobs|-j)            JOBS="${2:?--jobs needs a number}"; shift 2 ;;
        -h|--help)            usage 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage 1 ;;
    esac
done

# Require cmake
if ! command -v cmake >/dev/null 2>&1; then
    echo "ERROR: cmake not found. Install it first:" >&2
    echo "  Linux: sudo apt install cmake zlib1g-dev" >&2
    echo "  macOS: brew install cmake zlib" >&2
    exit 1
fi

CONFIGURE_ARGS=( -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE" )
if [ -n "$PREFIX" ]; then
    CONFIGURE_ARGS+=( -DCMAKE_INSTALL_PREFIX="$PREFIX" )
fi

echo "==> Configuring (build type: ${BUILD_TYPE})"
cmake "${CONFIGURE_ARGS[@]}"

echo "==> Building (-j ${JOBS})"
cmake --build "$BUILD_DIR" -j "$JOBS"

if [ "$DO_INSTALL" -eq 1 ]; then
    if [ -n "$PREFIX" ]; then
        echo "==> Installing to ${PREFIX}"
        cmake --install "$BUILD_DIR"
    else
        echo "==> Installing to system prefix (using sudo)"
        sudo cmake --install "$BUILD_DIR"
    fi
fi

echo "==> Done. Binary: ${BUILD_DIR}/impaqt"
