# This file is part of the Acts project.
#
# Copyright (C) 2021 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
include(FetchContent)

# Tell the user what's happening.
message(STATUS "Building Covfie as part of the Acts project")

set(COVFIE_VERSION "v${_acts_covfie_version}")

# Declare where to get covfie from.
FetchContent_Declare(covfie ${ACTS_COVFIE_SOURCE})

# Options used for covfie.
set(COVFIE_BUILD_EXAMPLES OFF CACHE BOOL "Build covfie examples")
set(COVFIE_BUILD_TESTS OFF CACHE BOOL "Build covfie tests")
set(COVFIE_BUILD_BENCHMARKS OFF CACHE BOOL "Build covfie benchmarks")

set(COVFIE_PLATFORM_CPU ON CACHE BOOL "Enable covfie CPU platform")
set(COVFIE_PLATFORM_CUDA
    ${ACTS_BUILD_PLUGIN_CUDA}
    CACHE BOOL
    "Enable covfie CUDA platform"
)

set(COVFIE_REQUIRE_CXX20 OFF CACHE BOOL "Enable covfie C++20 requirement")
set(COVFIE_QUIET ON CACHE BOOL "Quiet covfie feature warnings")

# Now set up its build.
FetchContent_MakeAvailable(covfie)
