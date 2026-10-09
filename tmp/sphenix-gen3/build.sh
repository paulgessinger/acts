#!/usr/bin/env bash
# Build the standalone Gen3 builder against an existing ACTS build tree.
#   ACTS_BUILD=/path/to/acts/build ./build.sh [sphenix_gen3.cpp]
# Needs ROOT (root-config on PATH) and Eigen (EIGEN_INCLUDE, default
# /usr/include/eigen3). ACTS must be built with ACTS_BUILD_PLUGIN_ROOT=ON.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../.." && pwd)
: "${ACTS_BUILD:?set ACTS_BUILD to the ACTS build directory}"
: "${EIGEN_INCLUDE:=/usr/include/eigen3}"
: "${CXX:=g++}"
in=${1:-$here/sphenix_gen3.cpp}
libdir=$(dirname "$(find "$ACTS_BUILD" -name 'libActsCore.so' | head -1)")
$CXX -std=c++23 -O1 -g \
  -I"$src/Plugins/Root/include" -I"$src/Core/include" -I"$ACTS_BUILD/Core/include" \
  -isystem "$EIGEN_INCLUDE" $(root-config --cflags | sed 's/-std=[^ ]*//') \
  "$in" -o "$here/$(basename "${in%.cpp}")" \
  -L"$libdir" -Wl,-rpath,"$libdir" -lActsPluginRoot -lActsCore \
  $(root-config --libs) -lGeom
