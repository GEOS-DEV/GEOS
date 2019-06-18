include("${CMAKE_CURRENT_LIST_DIR}/../src/coreComponents/cxx-utilities/host-configs/lassen-clang@upstream.cmake")

# asmjit doesn't work on PowerPC
set(ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE)

# Silo configure script doesn't recognize systype
set(SILO_BUILD_TYPE "powerpc64-unknown-linux-gnu" CACHE STRING "" FORCE)

set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)
