# first Check for RAJA_DIR

if(NOT RAJA_DIR)
    MESSAGE(FATAL_ERROR "Could not find RAJA. RAJA requires explicit RAJA_DIR.")
endif()
MESSAGE("Finding RAJA in ${RAJA_DIR}")
include(${RAJA_DIR}/share/raja/cmake/raja-config.cmake OPTIONAL)
set (RAJA_FOUND FALSE)
if (RAJA_INCLUDE_PATH)
   set(RAJA_FOUND TRUE)
   set(RAJA_INCLUDE_DIRS ${RAJA_INCLUDE_PATH})
endif()
