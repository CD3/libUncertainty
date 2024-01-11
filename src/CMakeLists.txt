if(TARGET libUncertainty)
  return()
endif(TARGET libUncertainty)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_library(libUncertainty INTERFACE)
add_library(libUncertainty::libUncertainty ALIAS libUncertainty)
target_sources(
  libUncertainty
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libUncertainty/uncertain.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libUncertainty/propagate.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libUncertainty/correlation.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libUncertainty/tags.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libUncertainty/utils.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/libUncertainty/statistics.hpp>
)
target_include_directories(
  libUncertainty
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(libUncertainty INTERFACE)
target_compile_features(libUncertainty INTERFACE cxx_std_20)

# Install
install(
  TARGETS libUncertainty
  EXPORT libUncertaintyTargets)
install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/src/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
  PATTERN *.hpp
  PATTERN *.h)
install(
  DIRECTORY ${PROJECT_BINARY_DIR}/include/
  DESTINATION  ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
  PATTERN *.hpp
  PATTERN *.h)
install(
  EXPORT libUncertaintyTargets
  FILE libUncertaintyTargets.cmake
  NAMESPACE libUncertainty::
  DESTINATION cmake)
file(
  WRITE ${CMAKE_CURRENT_BINARY_DIR}/libUncertaintyConfig.cmake
  "include(CMakeFindDependencyMacro)
include(\${CMAKE_CURRENT_LIST_DIR}/libUncertaintyTargets.cmake)
")
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/libUncertaintyConfigVersion.cmake
  VERSION ${libUncertainty_VERSION}
  COMPATIBILITY SameMajorVersion)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/libUncertaintyConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/libUncertaintyConfigVersion.cmake
  DESTINATION cmake)

