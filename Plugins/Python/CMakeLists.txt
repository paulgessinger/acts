# if(SKBUILD)
#   # Scikit-Build does not add your site-packages to the search path automatically,
#   # so we need to add it _or_ the pybind11 specific directory here.
#   execute_process(
#     COMMAND
#       "${PYTHON_EXECUTABLE}" -c
#       "import pybind11; print(pybind11.get_cmake_dir())"
#     OUTPUT_VARIABLE _tmp_dir
#     OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
#   list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
# endif()

add_library(ActsPluginPython INTERFACE)
#   ActsPluginJson SHARED
#   src/JsonGeometryConverter.cpp
#   src/AlgebraJsonConverter.cpp
#   src/SurfaceBoundsJsonConverter.cpp
#   src/SurfaceJsonConverter.cpp
#   src/UtilitiesJsonCon_actsverter.cpp)
# target_include_directories(
#   ActsPluginJson
#   PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
#     $<INSTALL_INTERFACE:include>)
# target_link_libraries(
#   ActsPluginJson PUBLIC ActsCore nlohmann_json::nlohmann_json)

install(
  TARGETS ActsPluginPython
  EXPORT ActsPluginPythonTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
# install(
#   DIRECTORY include/Acts
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

set(_python_dir "${CMAKE_BINARY_DIR}/python/")
# file(MAKE_DIRECTORY ${_python_dir})
# file(MAKE_DIRECTORY ${_python_dir}/acts)

pybind11_add_module(_acts 
  src/ActsModule.cpp
  src/Base.cpp
  src/ExamplesDetector.cpp
  src/Material.cpp
  src/Geometry.cpp
  src/ExampleAlgorithms.cpp
  src/MagneticField.cpp
  src/Output.cpp
  src/Propagation.cpp
  src/Generators.cpp
)


if(ACTS_BUILD_PLUGIN_DD4HEP AND ACTS_BUILD_EXAMPLES_DD4HEP)
  target_link_libraries(_acts PUBLIC ActsExamplesDetectorDD4hep)
  target_sources(_acts PRIVATE src/DD4hepComponent.cpp)
endif()

target_include_directories(_acts PRIVATE include)

set_target_properties(_acts PROPERTIES INSTALL_RPATH "\$ORIGIN/${CMAKE_INSTALL_LIBDIR}")

# set_target_properties(_acts PROPERTIES OUTPUT_NAME _acts)
set_target_properties(_acts PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts)
# set_target_properties(_acts PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${_python_dir}/acts)

target_link_libraries(_acts PUBLIC 
  ActsCore 
  ActsExamplesFramework
  ActsExamplesDetectorsCommon
  ActsExamplesDetectorGeneric
  ActsExamplesMagneticField
  ActsExamplesIoRoot
  ActsExamplesIoCsv
  ActsExamplesIoObj
  ActsExamplesIoJson
  ActsExamplesIoPerformance
  ActsExamplesGenerators
  ActsExamplesTrackFinding
  ActsExamplesFatras
)


add_custom_target(ActsPythonPyCode)

# configure_file(python/_version.py ${_python_dir}/acts/)

set(py_files 
  __init__.py 
  examples.py
)

foreach(f ${py_files})
add_custom_command(
        TARGET ActsPythonPyCode POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_SOURCE_DIR}/python/acts/${f}
                ${_python_dir}/acts/${f})
endforeach()

# add_custom_command(
#         TARGET ActsPythonPyCode POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E create_symlink
#                 ${CMAKE_CURRENT_SOURCE_DIR}/python/__init__.py
#                 ${_python_dir}/acts/__init__.py)
# add_custom_command(
#         TARGET ActsPythonPyCode POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E create_symlink
#                 ${CMAKE_CURRENT_SOURCE_DIR}/python/setup.py
#                 ${_python_dir}/setup.py)

add_dependencies(_acts ActsPythonPyCode)

# configure_file(python/__init__.py ${_python_dir}/acts/)
# configure_file(python/setup.py ${_python_dir}/)

install(TARGETS _acts DESTINATION .)
# install(TARGETS _acts
#   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#   RUNTIME DESTssddNATION ${CMAKE_INSTALL_BINDIR})
  