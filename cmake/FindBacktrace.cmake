include(FindPackageHandleStandardArgs)

# macro(done)
  # find_package_handle_standard_args(Backtrace DEFAULT_MSG
      # REQUIRED_VARS backtrace_INCLUDE_DIR)
# endmacro()

if(APPLE)
if(Backtrace_FIND_REQUIRED)
  message(SEND_ERROR "FindBacktrace.cmake currently not supported on macOS")
endif()
endif()


file(WRITE
  "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
  "
  #include \"backtrace.h\"
  int main() {}
  " )


if(Backtrace_FIND_QUIETLY)
message(CHECK_START "Does backtrace work with the default include")
endif()

try_compile(_backtrace_default_header "${CMAKE_BINARY_DIR}"
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
    OUTPUT_VARIABLE __OUTPUT)


if(_backtrace_default_header)
  if(Backtrace_FIND_QUIETLY)
    message(CHECK_PASS "yes")
  endif()
  set(Backtrace_HEADER "backtrace.h")
else()
  if(Backtrace_FIND_QUIETLY)
    message(CHECK_FAIL "no")
  endif()
  find_file(Backtrace_HEADER "backtrace.h" REQUIRED)
  # get_filename_component(Backtrace_INCLUDE_DIR ${Backtrace_HEADER} DIRECTORY)
  if(Backtrace_FIND_QUIETLY)
    message(CHECK_START "Does backtrace work with ${Backtrace_HEADER}")
  endif()

  file(WRITE
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
    "
    #include \"${Backtrace_HEADER}\"
    int main() {}
    " )


  try_compile(_backtrace_explicit_header "${CMAKE_BINARY_DIR}"
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
      OUTPUT_VARIABLE __OUTPUT)

  if(_backtrace_explicit_header)
    if(Backtrace_FIND_QUIETLY)
      message(CHECK_PASS "yes")
    endif()
  else()
    if(Backtrace_FIND_QUIETLY)
      message(CHECK_FAIL "no")
    endif()
    if(Backtrace_FIND_REQUIRED)
      message(SEND_ERROR "FindBacktrace.cmake did not manage to compile with default or explicit header")
    endif()
  endif()
endif()

file(WRITE
  "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
  "#include <${Backtrace_HEADER}>\n"
  "int main() { backtrace_state* state = backtrace_create_state(NULL, 0, NULL, NULL); }\n" )

if(Backtrace_FIND_QUIETLY)
  message(CHECK_START "Does backtrace work without linker flag")
endif()
try_compile(_backtrace_nolink "${CMAKE_BINARY_DIR}"
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
    OUTPUT_VARIABLE __OUTPUT)

if(_backtrace_nolink)
  if(Backtrace_FIND_QUIETLY)
    message(CHECK_PASS "yes")
  endif()
else()
  if(Backtrace_FIND_QUIETLY)
    message(CHECK_FAIL "no")
  endif()

  find_library(Backtrace_LIBRARY NAMES backtrace)
  if("${Backtrace_LIBRARY}" STREQUAL "Backtrace_LIBRARY-NOTFOUND")
    set(Backtrace_LIBRARY backtrace)
  endif()

file(WRITE
  "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
  "#include <${Backtrace_HEADER}>\n"
  "int main() { backtrace_state* state = backtrace_create_state(NULL, 0, NULL, NULL); }\n" )

  if(Backtrace_FIND_QUIETLY)
    message(CHECK_START "Does backtrace work with explicit library")
  endif()
  try_compile(_backtrace_link "${CMAKE_BINARY_DIR}"
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
      LINK_LIBRARIES ${Backtrace_LIBRARY}
      OUTPUT_VARIABLE __OUTPUT)

  if(_backtrace_link)
    if(Backtrace_FIND_QUIETLY)
      message(CHECK_PASS "yes")
    endif()
  else()
    if(Backtrace_FIND_QUIETLY)
      message(CHECK_FAIL "no")
    endif()
    if(Backtrace_FIND_REQUIRED)
      message(SEND_ERROR "FindBacktrace.cmake did not manage to compile with or without explicit library")
    endif()
  endif()

endif()

find_package_handle_standard_args(Backtrace
  REQUIRED_VARS Backtrace_HEADER)
