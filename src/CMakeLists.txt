# CMake source (/src) configuration for Starscream
# Author contact: jayjaybillings@gmail.com

# Set the target
set (LIB Starscream)

# Set the C++ src files
file (GLOB CXX_SRC *.cpp)

# Add the library
add_library (${LIB} ${CXX_SRC})

# Set the install targets
install (TARGETS
  ${LIB} DESTINATION lib
)
