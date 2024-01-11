cmake_minimum_required(VERSION 3.20)

project(libUncertainty)

set(libUncertainty_VERSION 0.2 CACHE STRING "Build version number")
set(libUncertainty_VERSION $CACHE{libUncertainty_VERSION})

message(STATUS "Version Number: ${libUncertainty_VERSION}")

add_subdirectory(src)
add_subdirectory(testing)
