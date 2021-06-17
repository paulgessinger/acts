add_library(ActsPluginPython INTERFACE)

install(
  TARGETS ActsPluginPython
  EXPORT ActsPluginPythonTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

set(_python_dir "${CMAKE_BINARY_DIR}/python")
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

pybind11_add_module(ActsPythonBindings 
  src/ActsModule.cpp
  src/Base.cpp
  src/ExamplesDetector.cpp
  src/Material.cpp
  src/Geometry.cpp
  src/ExampleAlgorithms.cpp
  src/MagneticField.cpp
  src/Output.cpp
  src/Input.cpp
  src/Propagation.cpp
  src/Generators.cpp
  src/TruthTracking.cpp
  src/TrackFitting.cpp
  src/TrackFinding.cpp
)


if(ACTS_BUILD_PLUGIN_DD4HEP AND ACTS_BUILD_EXAMPLES_DD4HEP)
  target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesDetectorDD4hep)
  target_sources(ActsPythonBindings PRIVATE src/DD4hepComponent.cpp)
endif()

if(ACTS_BUILD_EXAMPLES_PYTHIA8)
  target_link_libraries(ActsPythonBindings PUBLIC ActsExamplesGeneratorsPythia8)
  target_sources(ActsPythonBindings PRIVATE src/Pythia8Component.cpp)
endif()

target_include_directories(ActsPythonBindings PRIVATE src)

set_target_properties(ActsPythonBindings PROPERTIES INSTALL_RPATH "\$ORIGIN/${CMAKE_INSTALL_LIBDIR}")

set_target_properties(ActsPythonBindings PROPERTIES OUTPUT_NAME _acts)
set_target_properties(ActsPythonBindings PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts)
# set_target_properties(ActsPythonBindings PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${_python_dir}/acts)

target_link_libraries(ActsPythonBindings PUBLIC 
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
  ActsExamplesTrackFitting
  ActsExamplesFatras
  ActsExamplesPrinters
)


add_custom_target(ActsPythonGlueCode)
configure_file(setup.sh.in ${_python_dir}/setup.sh)

set(py_files 
  __init__.py 
  examples.py
)

foreach(f ${py_files})
add_custom_command(
        TARGET ActsPythonGlueCode POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_SOURCE_DIR}/python/acts/${f}
                ${_python_dir}/acts/${f})
endforeach()

add_dependencies(ActsPythonBindings ActsPythonGlueCode)

install(TARGETS ActsPythonBindings DESTINATION .)
  