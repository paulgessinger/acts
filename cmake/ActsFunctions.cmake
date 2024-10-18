function(acts_add_library name type)
    add_library(${name} ${type} ${ARGN})
    if(type STREQUAL "SHARED" OR type STREQUAL "STATIC")
        _add_options(${name} ${type})
    endif()
endfunction()

function(acts_add_executable name)
    add_executable(${name} ${ARGN})
    _add_options(${name} EXECUTABLE)
endfunction()

function(_add_options name type)
    string(REPLACE " " ";" cxx_flags ${ACTS_CXX_FLAGS})
    target_compile_options(${name} PUBLIC ${cxx_flags})

    if(ACTS_RUN_UB_SANITIZER)
        if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            target_compile_options(
                ${name}
                PUBLIC
                    "-fsanitize=undefined"
                    "-fno-omit-frame-pointer"
                    "-fno-sanitize-recover=all"
            )

            target_link_options(
                ${name}
                PUBLIC "SHELL:-fuse-ld=lld -fsanitize=undefined"
            )
        endif()
    endif()
endfunction()
