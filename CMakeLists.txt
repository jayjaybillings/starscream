# CMake root configuration for Starscream
# Author contact: jayjaybillings@gmail.com

# Set CMake required version and project name
cmake_minimum_required (VERSION 2.6)
project (Starscream)

# Set the version number
set (Starscream_VERSION_MAJOR 2)
set (Starscream_VERSION_MINOR 0)
set (Starscream_VERSION 
     "${Starscream_VERSION_MAJOR}.${Starscream_VERSION_MINOR}")

# Say hello
MESSAGE(STATUS 
        "Configuring CMake build system for Starscream ${Starscream_VERSION}")

# Set the master build flag. This flag is used to signal if the build should be
# performed or if it should be aborted.
set (MASTER_BUILD_FLAG true)

# Check the OS type
include(CMakeMacros/CheckOS.cmake)

# Setup the build options
include(CMakeMacros/SetBuildOptions.cmake)

# Instruct CMake to create shared libraries
#set (BUILD_SHARED_LIBS ON)

# Add the binary directory to the search path so that the Starscream
# configuration file will be found
include_directories ("${PROJECT_BINARY_DIR}/src")

# Add the Starscream source directory to the build configuration, include files
# first and then the source files
include_directories ("${PROJECT_SOURCE_DIR}/src")
add_subdirectory (src)

# Add an Eclipse project name for the Eclipse CDT4 CMake generator
project (StarscreamSrc)
