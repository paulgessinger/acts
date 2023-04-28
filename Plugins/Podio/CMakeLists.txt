add_library(
  ActsPluginPodio SHARED
  src/PodioUtil.cpp
)

target_include_directories(
  ActsPluginPodio
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsPluginPodio
  PUBLIC ActsCore ActsPluginIdentification)

# target_link_libraries(ActsPluginEDM4hep PUBLIC EDM4HEP::edm4hep)


# message(STATUS "IO HANDLERS: ${PODIO_IO_HANDLERS}")

PODIO_GENERATE_DATAMODEL(
  ActsPodioEdm 
  ${CMAKE_CURRENT_LIST_DIR}/edm.yml 
  headers
  sources
  IO_BACKEND_HANDLERS ${PODIO_IO_HANDLERS}
)


# message(STATUS "HEADERS: ${headers}")
# message(STATUS "SOURCES: ${sources}")

PODIO_ADD_DATAMODEL_CORE_LIB(ActsPodioEdm "${headers}" "${sources}")
# set_target_properties(ActsPodioEdm PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(ActsPluginPodio PUBLIC 
  ActsPodioEdm 
  ROOT::Core 
  podio::podio
  podio::podioRootIO
)

PODIO_ADD_ROOT_IO_DICT(ActsPodioEdmDict ActsPodioEdm "${headers}" src/selection.xml)
add_library(Acts::ActsPodioEdmDict ALIAS ActsPodioEdmDict)

install(
  TARGETS ActsPodioEdm
  EXPORT ActsPodioEdmTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ActsPodioEdm)

install(
  DIRECTORY ActsPodioEdm
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ActsPodioEdm)

set(install_package_config_dir "${CMAKE_INSTALL_LIBDIR}/cmake/Acts")
install(
  EXPORT ActsPodioEdmTargets
  DESTINATION ${install_package_config_dir})

install(
  TARGETS ActsPluginPodio
  EXPORT ActsPluginPodioTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ActsPodioEdmDict
  EXPORT ActsPodioEdmTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ActsPodioEdm"
  COMPONENT dev)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/ActsPodioEdmDictDict.rootmap"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT dev)

if (${ROOT_VERSION} GREATER 6)
  install(FILES
      "${CMAKE_CURRENT_BINARY_DIR}/libActsPodioEdmDict_rdict.pcm"
      DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT dev)
endif()

# if (nlohmann_json_FOUND)
  # target_compile_definitions(edm4hep PUBLIC PODIO_JSON_OUTPUT)
  # target_link_libraries(edm4hep PUBLIC nlohmann_json::nlohmann_json)
# endif()


# # list(APPEND EDM4HEP_INSTALL_LIBS edm4hep edm4hepDict)

# PODIO_ADD_SIO_IO_BLOCKS(ActsPodioEdmDict "${headers}" "${sources}")
# IF(TARGET edm4hepSioBlocks)
  # message(STATUS "Building and installing the SioBlocks since podio supports it")
  # list(APPEND EDM4HEP_INSTALL_LIBS edm4hepSioBlocks)
# ENDIF()

# install(TARGETS ${EDM4HEP_INSTALL_LIBS}
  # EXPORT EDM4HEPTargets
  # RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  # LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  # PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/edm4hep"
  # COMPONENT dev)

# install(FILES
  # "${PROJECT_BINARY_DIR}/edm4hep/edm4hepDictDict.rootmap"
  # DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT dev)

# install(FILES
  # ../edm4hep.yaml
  # DESTINATION "${CMAKE_INSTALL_DATADIR}/edm4hep" COMPONENT dev)

# if (${ROOT_VERSION} GREATER 6)
  # install(FILES
      # "${PROJECT_BINARY_DIR}/edm4hep/libedm4hepDict_rdict.pcm"
      # DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT dev)
# endif()

