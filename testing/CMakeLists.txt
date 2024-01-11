option(BUILD_TESTS "Build the Unit Tests" ON )
if( NOT BUILD_TESTS )
  return ()
endif()

find_package(Boost REQUIRED)
find_package(BoostUnitDefinitions REQUIRED)
find_package(Catch2 REQUIRED)

file( GLOB_RECURSE unitTest_SOURCES
      RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "./CatchTests/*.cpp" )
message(STATUS "Detected Catch-based Unit Test Sources:")
foreach(testSrc ${unitTest_SOURCES})
  message(STATUS "  ${testSrc}" )
endforeach()

add_executable(libUncertainty_CatchTests ${unitTest_SOURCES})
target_compile_definitions( libUncertainty_CatchTests PRIVATE -DTESTING )
target_link_libraries(libUncertainty_CatchTests PUBLIC libUncertainty
  Boost::headers
  BoostUnitDefinitions::BoostUnitDefinitions
  Catch2::Catch2WithMain
  )
target_include_directories( libUncertainty_CatchTests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include" )

