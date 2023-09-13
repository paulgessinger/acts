add_library(
	ActsPluginFpeMonitoring SHARED
	src/FpeMonitor.cpp)
target_include_directories(
  ActsPluginFpeMonitoring
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(
  ActsPluginFpeMonitoring PUBLIC ActsCore)


# find_package(Backtrace)
include(SetupBoostStacktrace)

target_link_libraries(ActsPluginFpeMonitoring PRIVATE boost::stacktrace)

install(
  TARGETS ActsPluginFpeMonitoring
  EXPORT ActsPluginFpeMonitoringTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
