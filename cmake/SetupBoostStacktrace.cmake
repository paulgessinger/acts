include_guard(GLOBAL)

include(CheckCXXSourceCompiles)

find_library(dl_LIBRARY dl REQUIRED)

find_package(Backtrace)
find_program(addr2line_EXECUTABLE addr2line)

add_library(BoostStacktrace INTERFACE)
add_library(boost::stacktrace ALIAS BoostStacktrace)

set(_setup_complete FALSE)

if(APPLE)
  # @TODO: implement
  return()
else()
  target_link_libraries(BoostStacktrace INTERFACE ${dl_LIBRARY})
    
  if(Backtrace_FOUND)
    message(CHECK_START "Does stacktrace work with backtrace")

    file(WRITE
      "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
      "#include <boost/stacktrace.hpp>\n"
      "#include <iostream>\n"
      "int main() { std::cout << boost::stacktrace::stacktrace(); }\n" )

    try_compile(_stacktrace_backtrace "${CMAKE_BINARY_DIR}"
        "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/backtrace.cpp"
        LINK_LIBRARIES ${dl_LIBRARY} ${Backtrace_LIBRARY}
        COMPILE_DEFINITIONS 
          -DBOOST_STACKTRACE_USE_BACKTRACE
          -DBOOST_STACKTRACE_BACKTRACE_INCLUDE_FILE=\"${Backtrace_HEADER}\"
        OUTPUT_VARIABLE __OUTPUT)

    if(_stacktrace_backtrace)
      message(CHECK_PASS "yes")

      target_link_libraries(BoostStacktrace INTERFACE ${Backtrace_LIBRARY})
      target_compile_definitions(BoostStacktrace INTERFACE
        -DBOOST_STACKTRACE_USE_BACKTRACE
        -DBOOST_STACKTRACE_BACKTRACE_INCLUDE_FILE=\"${Backtrace_HEADER}\"
      )

    set(_setup_complete TRUE)

    else()
      message(CHECK_FAIL "no")
      if(NOT addr2line_FOUND)
        message(WARNING "Backtrace not available, and addr2line not found, symbols will not be resolved")
      endif()
    endif()
  endif()


  if(NOT _setup_complete)
    message(CHECK_START "Is addr2line available")
    if(addr2line_FOUND)
      message(CHECK_PASS "yes")
      target_compile_definitions(BoostStacktrace INTERFACE
        -DBOOST_STACKTRACE_USE_ADDR2LINE
        -DBOOST_STACKTRACE_ADDR2LINE_LOCATION=${addr2line_EXECUTABLE}
      )
    else()
      message(CHECK_FAIL "no")
    endif()
  endif()
endif()

