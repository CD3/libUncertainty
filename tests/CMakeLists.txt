cmake_minimum_required( VERSION 3.1 )

project(libUncertainty-Tests)
enable_testing()

set( binDir ${CMAKE_BINARY_DIR}/testBin )
set( buildDir ${CMAKE_BINARY_DIR} )

add_subdirectory(.. libUncertainty)

find_package(Boost REQUIRED)
find_package(BoostUnitDefinitions REQUIRED)
find_package(uncertainties-cpp REQUIRED)

if( NOT TARGET libUncertainty )
  add_subdirectory(.. libUncertainty)
endif()

# Find all sources.
file( GLOB_RECURSE unitTest_SOURCES
      RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "./CatchTests/*.cpp" )
message(STATUS "Detected Catch-based Unit Test Sources:")
foreach(testSrc ${unitTest_SOURCES})
  message(STATUS "  ${testSrc}" )
endforeach()

if( EXISTS "${CMAKE_BINARY_DIR}/conanbuildinfo.cmake" )
  include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
  conan_basic_setup(TARGETS)
endif()

set( testName libUncertainty_CatchTests )
add_executable(${testName} ${unitTest_SOURCES})
target_compile_definitions( ${testName} PRIVATE -DTESTING )
target_link_libraries(${testName} PUBLIC libUncertainty
  Boost::headers
  BoostUnitDefinitions::BoostUnitDefinitions
  uncertainties-cpp::uncertainties-cpp
  )
target_include_directories( ${testName} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include" )
set_target_properties(${testName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${binDir})
set_target_properties(${testName} PROPERTIES CXX_STANDARD 11)

add_test(NAME ${testName} 
         WORKING_DIRECTORY ${buildDir}
         COMMAND ${binDir}/${testName} )

