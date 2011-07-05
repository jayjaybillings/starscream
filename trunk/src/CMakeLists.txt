# CMake source (/src) configuration for Starscream
# Author contact: jayjaybillings@gmail.com

# Setup a header file to hold configuration variables for Starscream
configure_file (
  "${PROJECT_SOURCE_DIR}/src/starscreamConfig.h.in"
  "${PROJECT_BINARY_DIR}/src/starscreamConfig.h"
)

# Set the targets
set (LIB starscream)
set (PRG starscream_test)

# Set the C files for the library
set (LIB_C_SRC
     starscream_init.c
     starscream_io.c
     starscream_pf.c
     starscream_structure.c
     starscream_vel.c)

# Grab the headers
FILE(GLOB LIB_HEADERS *.h)

# Add the library
add_library (${LIB} ${LIB_C_SRC})

# Link the GSL and FFTW libraries
target_link_libraries(${LIB} gsl)
target_link_libraries(${LIB} gslcblas)
if(${ENABLE_PARALLELISM} MATCHES "ON")
  target_link_libraries(${LIB} fftw3_threads)
endif(${ENABLE_PARALLELISM} MATCHES "ON")
target_link_libraries(${LIB} fftw3)

# Build the Starscream test executable
add_executable(${PRG} starscream.c)
target_link_libraries(${PRG} ${LIB} m)

# Set the install targets
install (TARGETS
  ${LIB} DESTINATION ${LIB_INSTALL_PREFIX}
)
install (FILES ${LIB_HEADERS} DESTINATION ${INCLUDE_INSTALL_PREFIX})
