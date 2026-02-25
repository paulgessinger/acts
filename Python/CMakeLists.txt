# Configure setup and install python files
set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

# List to collect built plugins
set(_plugins_built)

if(ACTS_BUILD_PYTHON_STUBS)
    execute_process(
        COMMAND
            ${Python_EXECUTABLE} -c "import pybind11_stubgen"
        RESULT_VARIABLE
            _acts_has_pybind11_stubgen
        OUTPUT_QUIET
        ERROR_QUIET
    )
    find_program(_acts_uvx_executable uvx)

    if(
        NOT _acts_has_pybind11_stubgen EQUAL 0
        AND NOT _acts_uvx_executable
    )
        if(ACTS_BUILD_PYTHON_WHEEL)
            message(
                FATAL_ERROR
                "ACTS_BUILD_PYTHON_STUBS=ON requires pybind11-stubgen. Install it in the wheel build environment."
            )
        else()
            message(
                WARNING
                "Disabling ACTS_BUILD_PYTHON_STUBS because neither pybind11-stubgen nor uvx was found."
            )
            set(ACTS_BUILD_PYTHON_STUBS OFF)
        endif()
    endif()
endif()

# Add post-build stub generation commands to a Python module target.
#
# Arguments:
#   TARGET  : CMake target name to attach POST_BUILD commands to.
#   MODULES : One or more Python module import paths for pybind11-stubgen.
#   DEPENDS : Optional target dependencies required before stub generation.
#   WORKING_DIRECTORY : Optional working directory for stubgen invocation.
#                       Defaults to `${CMAKE_BINARY_DIR}/python`.
#
# Behavior:
#   - No-op when ACTS_BUILD_PYTHON_STUBS is OFF.
#   - Validates required arguments and rejects unknown ones.
#   - Invokes pybind11-stubgen directly (via Python module or uvx fallback).
function(acts_python_add_stubgen_post_build)
    set(_stubgen_options)
    set(_stubgen_one_value_args TARGET WORKING_DIRECTORY)
    set(_stubgen_multi_value_args MODULES DEPENDS)
    cmake_parse_arguments(
        PARSE_ARGV 0
        _stubgen
        "${_stubgen_options}"
        "${_stubgen_one_value_args}"
        "${_stubgen_multi_value_args}"
    )

    if(NOT ACTS_BUILD_PYTHON_STUBS)
        return()
    endif()

    if(_stubgen_UNPARSED_ARGUMENTS)
        message(
            FATAL_ERROR
            "acts_python_add_stubgen_post_build got unexpected arguments: ${_stubgen_UNPARSED_ARGUMENTS}"
        )
    endif()

    if(NOT _stubgen_TARGET)
        message(
            FATAL_ERROR
            "acts_python_add_stubgen_post_build requires TARGET."
        )
    endif()

    if(NOT _stubgen_MODULES)
        message(
            FATAL_ERROR
            "acts_python_add_stubgen_post_build requires at least one MODULES entry."
        )
    endif()

    if(NOT _stubgen_WORKING_DIRECTORY)
        set(_stubgen_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/python")
    endif()

    if(_stubgen_DEPENDS)
        add_dependencies(${_stubgen_TARGET} ${_stubgen_DEPENDS})
    endif()

    set(_stubgen_commands)

    if(_acts_has_pybind11_stubgen EQUAL 0)
        set(_stubgen_cli ${Python_EXECUTABLE} -m pybind11_stubgen)
    elseif(_acts_uvx_executable)
        set(
            _stubgen_cli
            ${_acts_uvx_executable}
            -p
            ${Python_EXECUTABLE}
            -w
            numpy
            pybind11-stubgen
        )
    else()
        message(
            FATAL_ERROR
            "ACTS_BUILD_PYTHON_STUBS=ON requires pybind11-stubgen or uvx."
        )
    endif()

    foreach(_stubgen_module IN LISTS _stubgen_MODULES)
        list(
            APPEND
            _stubgen_commands
            COMMAND
                ${CMAKE_COMMAND}
                -E
                env
                "PYTHONPATH=${_stubgen_WORKING_DIRECTORY}:${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}"
                ${_stubgen_cli}
                ${_stubgen_module}
                -o
                .
                --ignore-invalid-expressions
                .*
                --ignore-invalid-identifiers
                "GeoModelTree::FpvConstLink"
                --ignore-unresolved-names
                "_abc\\._abc_data"
        )
    endforeach()

    list(
        APPEND
        _stubgen_commands
        COMMAND
            ${Python_EXECUTABLE}
            ${PROJECT_SOURCE_DIR}/Python/scripts/stubs/patch_pybind11_stub_signatures.py
            --root
            ${CMAKE_BINARY_DIR}/python/acts
    )

    add_custom_command(
        TARGET
            ${_stubgen_TARGET}
        POST_BUILD
        ${_stubgen_commands}
        WORKING_DIRECTORY
            ${_stubgen_WORKING_DIRECTORY}
        COMMENT
            "Generating stubs for ${_stubgen_TARGET}"
        VERBATIM
    )
endfunction()

# Add the bindings
add_subdirectory(Core)
add_subdirectory(Fatras)
add_subdirectory_if(Examples ACTS_BUILD_EXAMPLES)
add_subdirectory(Plugins)
add_subdirectory(Utilities)

if(ACTS_BUILD_PYTHON_STUBS)
    # Include generated .pyi files in install trees and wheels.
    install(
        DIRECTORY ${_python_dir}/acts/
        DESTINATION ${_acts_python_install_dir}
        FILES_MATCHING
        PATTERN "*.pyi"
    )
endif()
