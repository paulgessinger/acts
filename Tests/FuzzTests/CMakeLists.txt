# Fuzz testing infrastructure using LibFuzzer
#
# Fuzz tests require Clang and LibFuzzer support.
# They are automatically instrumented with ASAN and UBSAN for optimal coverage.

# Check if we're using a compiler that supports LibFuzzer
if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(
        FATAL_ERROR
        "Fuzz tests require Clang compiler with LibFuzzer support. "
        "Current compiler: ${CMAKE_CXX_COMPILER_ID}"
    )
endif()

# Check if LibFuzzer is actually available by trying a test compilation
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-fsanitize=fuzzer")
check_cxx_source_compiles(
    "extern \"C\" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) { return 0; } int main() { return 0; }"
    LIBFUZZER_AVAILABLE
)
set(CMAKE_REQUIRED_FLAGS "")

if(NOT LIBFUZZER_AVAILABLE)
    message(
        FATAL_ERROR
        "LibFuzzer is not available in the current Clang installation.\n"
        "Apple Clang (from Xcode) does not include LibFuzzer.\n"
        "Please install LLVM/Clang via Homebrew or another package manager:\n"
        "  brew install llvm\n"
        "Then configure with:\n"
        "  cmake -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++ ...\n"
        "Or disable fuzz tests with:\n"
        "  cmake -DACTS_BUILD_FUZZTESTS=OFF ..."
    )
endif()

add_custom_target(fuzz_tests)

# Macro to add a fuzz test executable
#
# Usage: add_fuzztest(TestName source.cpp [EXTRA_LIBRARIES lib1 lib2 ...])
#
# This creates an executable with LibFuzzer instrumentation.
# Fuzz tests are NOT registered with CTest as they run continuously
# rather than as pass/fail tests.
macro(add_fuzztest _name _source)
    # Parse optional EXTRA_LIBRARIES argument
    set(_extra_libs "")
    set(_parsing_libs FALSE)
    foreach(_arg ${ARGN})
        if(_arg STREQUAL "EXTRA_LIBRARIES")
            set(_parsing_libs TRUE)
        elseif(_parsing_libs)
            list(APPEND _extra_libs ${_arg})
        endif()
    endforeach()

    # Automatically prefix the target name
    set(_target "ActsFuzzTest${_name}")
    add_executable(${_target} ${_source})

    # Add LibFuzzer instrumentation
    target_compile_options(
        ${_target}
        PRIVATE -fsanitize=fuzzer,address,undefined
    )
    target_link_options(${_target} PRIVATE -fsanitize=fuzzer,address,undefined)

    # Add common dependencies
    target_include_directories(${_target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(
        ${_target}
        PRIVATE Acts::TestsCommonHelpers ${_extra_libs}
    )

    # Add to the fuzz_tests target (but NOT to CTest)
    add_dependencies(fuzz_tests ${_target})

    # Print information about how to run the fuzz test
    message(STATUS "Fuzz test '${_name}' will be built as '${_target}'")
endmacro()

# Add fuzz tests
add_fuzztest(MaterialInteractions InteractionsFuzz.cpp EXTRA_LIBRARIES Acts::Core)
