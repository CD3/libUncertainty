cmake_minimum_required( VERSION 3.11 )

if( TARGET libUncertainty )
  return()
endif( TARGET libUncertainty )

OPTION( BUILD_UNIT_TESTS "Build unit tests for the library" ON )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
set(INSTALL_PATH ${INSTALL_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/externals/install/")

include(function-git_version)
GIT_VERSION(libUncertainty)
message( STATUS "libUncertainty version: ${libUncertainty_VERSION}")
project(libUncertainty VERSION ${libUncertainty_VERSION})
# need to run configure_file twice
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/src/${PROJECT_NAME}/version.h.in"
                "${CMAKE_CURRENT_BINARY_DIR}/include/${PROJECT_NAME}/version.h")
configure_file( "${CMAKE_CURRENT_BINARY_DIR}/include/${PROJECT_NAME}/version.h"
                "${CMAKE_CURRENT_BINARY_DIR}/include/${PROJECT_NAME}/version.h")

# common dependencies
# find_package( Boost REQUIRED )
# find_package( Eigen3 3.3 REQUIRED )
# load conan build information if it exists
#if( EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake" )
#  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
#  conan_basic_setup(TARGETS)
#endif()

# if( NOT TARGET CONAN_PKG::boost )
# find_package( Boost COMPONENTS filesystem REQUIRED )
# endif()


set( LIB_TARGET libUncertainty )
add_library( ${LIB_TARGET} INTERFACE )
add_library( ${PROJECT_NAME}::${LIB_TARGET} ALIAS ${LIB_TARGET} )
target_sources(
  ${LIB_TARGET}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libUncertainty/uncertain.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libUncertainty/propagate.hpp>
)
target_include_directories( ${LIB_TARGET} INTERFACE
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${${PROJECT_NAME}_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(${LIB_TARGET}
                        INTERFACE
                        # $<$<TARGET_EXISTS:CONAN_PKG::boost>:CONAN_PKG::boost>
                          #Boost::boost
                          #Eigen3::Eigen
                     )

             target_compile_features(${LIB_TARGET} INTERFACE cxx_std_17)
if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} AND BUILD_UNIT_TESTS)
add_subdirectory( tests )
endif()






# Install and Packaging
install( TARGETS ${LIB_TARGET} EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)
install(DIRECTORY src/
  DESTINATION include
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/
  DESTINATION include
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
install(EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION cmake
)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
"include(CMakeFindDependencyMacro)
# find_dependency(Boost)
include(\${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}Targets.cmake)
"
  )
include(CMakePackageConfigHelpers)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPATIBILITY SameMajorVersion
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION cmake
  )




INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "")
SET(CPACK_PACKAGE_VENDOR "")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.md")
if( ${PROJECT_NAME}_VERSION_MAJOR)
  SET(CPACK_PACKAGE_VERSION_MAJOR ${${PROJECT_NAME}_VERSION_MAJOR})
endif()
if( ${PROJECT_NAME}_VERSION_MINOR)
  SET(CPACK_PACKAGE_VERSION_MINOR ${${PROJECT_NAME}_VERSION_MINOR})
endif()
if( ${PROJECT_NAME}_VERSION_PATCH)
  SET(CPACK_PACKAGE_VERSION_PATCH ${${PROJECT_NAME}_VERSION_PATCH})
endif()
INCLUDE(CPack)



