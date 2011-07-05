# Set the default installation directories
set (LIB_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/lib)
set (INCLUDE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/include)

# ----- GSL ----- #

# Set the default directories to search for GSL include files
set (GSL_INCLUDES_SEARCH_DIRS 
     /usr/include/local
     /usr/include
     [CACHE INTERNAL "Search directories for GSL"])

# Locate GSL include directory
set(GSL_INCLUDE_PATH "" CACHE PATH "GSL Header File Directory")
# Include the directory in the build if it exists
if (EXISTS ${GSL_INCLUDE_PATH})
   include_directories(${GSL_INCLUDE_PATH})
endif (EXISTS ${GSL_INCLUDE_PATH})

# Set the GSL library search path
set(GSL_LIBS_PATH "" CACHE PATH "GSL Library Path")

# ----- FFTW ----- #

# Set the default directories to search for FFTW include files
set (FFTW_INCLUDES_SEARCH_DIRS 
     /usr/include/local
     /usr/include
     [CACHE INTERNAL "Search directories for FFTW"])

# Locate FFTW include directory
set(FFTW_INCLUDE_PATH "" CACHE PATH "FFTW Header File Directory")
# Include the directory in the build if it exists
if (EXISTS ${FFTW_INCLUDE_PATH})
   include_directories(${FFTW_INCLUDE_PATH})
endif (EXISTS ${FFTW_INCLUDE_PATH})

# Set the FFTW library search path
set(FFTW_LIBS_PATH "" CACHE PATH "FFTW Library Path")

# Link in the library directories
link_directories(${GSL_LIBS_PATH} ${FFTW_LIBS_PATH})