
#set(ENABLE_FORTRAN OFF CACHE BOOL "Enables Fortran" FORCE)
#set(ENABLE_CXX11 OFF CACHE BOOL "Enables C++11 language support" FORCE)
#set(ENABLE_CXX14 ON CACHE BOOL "Enables C++14 language support" FORCE)

if(ENABLE_CXX11)
    MESSAGE(FATAL_ERROR "c++11 is enabled. GEOSX requires c++14")
endif(ENABLE_CXX11)

if(NOT ENABLE_CXX14)
    MESSAGE(FATAL_ERROR "c++14 is NOT enabled. GEOSX requires c++14")
endif(NOT ENABLE_CXX14)

blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -pedantic)
blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -Wabi)
blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -Wno-variadic-macros)