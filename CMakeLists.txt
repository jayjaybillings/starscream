# CMake root configuration for Starscream
# Author contact: jayjaybillings@gmail.com

# Set CMake required version and project name
cmake_minimum_required (VERSION 2.6)
project (Starscream)

# Set the version number
set (Starscream_VERSION_MAJOR 2)
set (Starscream_VERSION_MINOR 0)

# Setup a header file to hold configuration variables for Starscream
configure_file (
  "${PROJECT_SOURCE_DIR}/src/starscreamConfig.h.in"
  "${PROJECT_BINARY_DIR}/src/starscreamConfig.h"
)
# Add the binary directory to the search path so that the Starscream
# configuration file will be found
include_directories ("${PROJECT_BINARY_DIR}")

# Add the Starscream source directory to the build configuration, include files
# first and then the source files
include_directories ("${PROJECT_SOURCE_DIR}/src")
add_subdirectory (src)

# Add an Eclipse project name for the Eclipse CDT4 CMake generator
project (StarscreamSrc)
