cmake_minimum_required(VERSION 3.11)

project(libUncertainty VERSION 0.2)

if(TARGET libUncertainty)
  return()
endif(TARGET libUncertainty)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
set(INSTALL_PATH ${INSTALL_PATH}
                 "${CMAKE_CURRENT_SOURCE_DIR}/externals/install/")


set(LIB_TARGET libUncertainty)
add_library(${LIB_TARGET} INTERFACE)
add_library(${PROJECT_NAME}::${LIB_TARGET} ALIAS ${LIB_TARGET})
target_sources(
  ${LIB_TARGET}
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libUncertainty/uncertain.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libUncertainty/propagate.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libUncertainty/correlation.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libUncertainty/tags.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libUncertainty/utils.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libUncertainty/statistics.hpp>
)
target_include_directories(
  ${LIB_TARGET}
  INTERFACE
    $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${${PROJECT_NAME}_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(${LIB_TARGET} INTERFACE)

target_compile_features(${LIB_TARGET} INTERFACE cxx_std_17)

# Install
install(
  TARGETS ${LIB_TARGET}
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
install(
  DIRECTORY src/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h")
install(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h")
install(
  EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION cmake)
file(
  WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  "include(CMakeFindDependencyMacro)
include(\${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}Targets.cmake)
")
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPATIBILITY SameMajorVersion)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION cmake)

