cmake_minimum_required(VERSION 3.20)
project(ArrowZeroCopyPoc CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

# Arrow's CMakeLists.txt does string(TOLOWER ${CMAKE_BUILD_TYPE} ...) without
# guarding against an unset variable, which errors in CMake 4.x when no build
# type is specified. Force a default here before Arrow is included.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ── Arrow ──────────────────────────────────────────────────────────────────
# Disable every optional feature to keep the build fast.
set(ARROW_BUILD_SHARED      ON  CACHE BOOL   "" FORCE)
set(ARROW_BUILD_STATIC      OFF CACHE BOOL   "" FORCE)
set(ARROW_COMPUTE           OFF CACHE BOOL   "" FORCE)
set(ARROW_CSV               OFF CACHE BOOL   "" FORCE)
set(ARROW_DATASET           OFF CACHE BOOL   "" FORCE)
set(ARROW_FILESYSTEM        OFF CACHE BOOL   "" FORCE)
set(ARROW_FLIGHT            OFF CACHE BOOL   "" FORCE)
set(ARROW_IPC               OFF CACHE BOOL   "" FORCE)
set(ARROW_JSON              OFF CACHE BOOL   "" FORCE)
set(ARROW_PARQUET           OFF CACHE BOOL   "" FORCE)
set(ARROW_TESTING           OFF CACHE BOOL   "" FORCE)
set(ARROW_BUILD_TESTS       OFF CACHE BOOL   "" FORCE)
set(ARROW_BUILD_BENCHMARKS  OFF CACHE BOOL   "" FORCE)
set(ARROW_BUILD_EXAMPLES    OFF CACHE BOOL   "" FORCE)
set(ARROW_BUILD_UTILITIES   OFF CACHE BOOL   "" FORCE)
set(ARROW_WITH_BROTLI       OFF CACHE BOOL   "" FORCE)
set(ARROW_WITH_BZ2          OFF CACHE BOOL   "" FORCE)
set(ARROW_WITH_LZ4          OFF CACHE BOOL   "" FORCE)
set(ARROW_WITH_SNAPPY       OFF CACHE BOOL   "" FORCE)
set(ARROW_WITH_ZLIB         OFF CACHE BOOL   "" FORCE)
set(ARROW_WITH_ZSTD         OFF CACHE BOOL   "" FORCE)
# BUNDLED: Arrow downloads and builds its own deps — no system Arrow required.
set(ARROW_DEPENDENCY_SOURCE "BUNDLED" CACHE STRING "" FORCE)
# Arrow 23 leaves ARROW_SIMD_LEVEL empty on ARM, breaking io/CMakeLists.txt.
set(ARROW_SIMD_LEVEL        "NONE"    CACHE STRING "" FORCE)

FetchContent_Declare(
  Arrow
  GIT_REPOSITORY https://github.com/apache/arrow.git
  GIT_TAG        apache-arrow-23.0.1
  GIT_SHALLOW    TRUE
  GIT_PROGRESS   TRUE
  SOURCE_SUBDIR  cpp          # Arrow's CMakeLists.txt lives under cpp/
)
FetchContent_MakeAvailable(Arrow)

# ── pybind11 ────────────────────────────────────────────────────────────────
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v3.0.2
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(pybind11)

# ── pybind11 extension module ───────────────────────────────────────────────
pybind11_add_module(arrow_bridge src/arrow_bridge.cpp)
target_link_libraries(arrow_bridge PRIVATE arrow_shared)
# Arrow sets includes via directory-scoped include_directories() calls
# (cpp/CMakeLists.txt:559-563), not target_include_directories(), so they
# don't cross the FetchContent boundary.  We have to add the two dirs ourselves:
#   arrow_SOURCE_DIR/cpp/src  — public headers (SOURCE_SUBDIR=cpp, headers in src/)
#   arrow_BINARY_DIR/src      — generated headers (e.g. arrow/util/config.h)
target_include_directories(arrow_bridge SYSTEM PRIVATE
  ${arrow_SOURCE_DIR}/cpp/src
  ${arrow_BINARY_DIR}/src
)
