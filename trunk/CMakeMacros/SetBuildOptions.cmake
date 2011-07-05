
# ----- Define all build options ----- #

# Shared library option
option (ENABLE_SHARED "Build shared libraries")

# OpenMP threading
option (ENABLE_PARALLELISM "Enable parallel execution" OFF)
mark_as_advanced(ENABLE_PARALLELISM)

# ----- Setup all build options ----- #

# Shared library
if (${ENABLE_SHARED} MATCHES "ON")
   set (BUILD_SHARED_LIBS ON)
   message(STATUS "Shared libraries will be built")
endif (${ENABLE_SHARED} MATCHES "ON")

# OpenMP threading, default to no
set (USE_FFTW_THREADS 0)
if (${ENABLE_PARALLELISM} MATCHES "ON")
   set (USE_FFTW_THREADS 1)
   set (STATUS "Parallel execution enabled")
endif (${ENABLE_PARALLELISM} MATCHES "ON")