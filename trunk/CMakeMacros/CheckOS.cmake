
# Check the OS and if it is Linux, set up the options
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Print a status message
  message(STATUS "Operating system ${CMAKE_SYSTEM_NAME} is supported.")
  # Setup the build options
  include(CMakeMacros/SetupLinuxOptions.cmake)
else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Print a fatal error message
  message(FATAL_ERROR "${CMAKE_SYSTEM_NAME} is not supported!")
  # Disable the build
  set (MASTER_BUILD_FLAG false)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
