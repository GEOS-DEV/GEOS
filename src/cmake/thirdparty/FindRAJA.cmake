# first Check for RAJA_DIR

if(NOT RAJA_DIR)
    MESSAGE(FATAL_ERROR "Could not find RAJA. RAJA requires explicit RAJA_DIR.")
endif()

include(${RAJA_DIR}/share/cmake/raja/raja-config.cmake)
set(RAJA_FOUND TRUE)
set(RAJA_INCLUDE_DIRS ${RAJA_DIR}/include)
