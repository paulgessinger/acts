add_library(ActsPluginPython INTERFACE)
#   ActsPluginJson SHARED
#   src/JsonGeometryConverter.cpp
#   src/AlgebraJsonConverter.cpp
#   src/SurfaceBoundsJsonConverter.cpp
#   src/SurfaceJsonConverter.cpp
#   src/UtilitiesJsonConverter.cpp)
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
file(MAKE_DIRECTORY ${_python_dir})
file(MAKE_DIRECTORY ${_python_dir}/acts)

pybind11_add_module(ActsPython 
  src/ActsModule.cpp
  src/ExamplesDetector.cpp
  src/Material.cpp
  src/Geometry.cpp
  src/ExampleAlgorithms.cpp
  src/MagneticField.cpp
  src/Output.cpp
)

if(ACTS_BUILD_PLUGIN_DD4HEP AND ACTS_BUILD_EXAMPLES_DD4HEP)
  target_compile_definitions(ActsPython PRIVATE -DPYBIND_DD4HEP)
  target_link_libraries(ActsPython PUBLIC ActsExamplesDetectorDD4hep)
endif()

set_target_properties(ActsPython PROPERTIES OUTPUT_NAME _acts)
set_target_properties(ActsPython PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${_python_dir}/acts)
set_target_properties(ActsPython PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${_python_dir}/acts)

target_link_libraries(ActsPython PUBLIC 
  ActsCore 
  ActsExamplesFramework
  ActsExamplesDetectorsCommon
  ActsExamplesDetectorGeneric
  ActsExamplesMagneticField
  ActsExamplesIoRoot
  ActsExamplesIoObj
)


add_custom_target(ActsPythonPyCode)

configure_file(python/_version.py ${_python_dir}/acts/)

add_custom_command(
        TARGET ActsPythonPyCode POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/python/__init__.py
                ${_python_dir}/acts/__init__.py)
add_custom_command(
        TARGET ActsPythonPyCode POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/python/setup.py
                ${_python_dir}/setup.py)

add_dependencies(ActsPython ActsPythonPyCode)

# configure_file(python/__init__.py ${_python_dir}/acts/)
# configure_file(python/setup.py ${_python_dir}/)
