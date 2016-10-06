####################################
# 3rd Party Dependencies
####################################

################################
# Conduit
################################
#message( "CONDUIT_DIR=${CONDUIT_DIR}")
if (CONDUIT_DIR)

  include(cmake/thirdparty/FindConduit.cmake)
  blt_register_library( NAME conduit
                        INCLUDES ${CONDUIT_INCLUDE_DIRS} 
                        LIBRARIES  conduit)
  blt_register_library( NAME conduit_io
                        INCLUDES ${CONDUIT_INCLUDE_DIRS}
                        LIBRARIES  conduit_io)
endif()


################################
# HDF5
################################
if (HDF5_DIR)
  include(cmake/thirdparty/FindHDF5.cmake)
  blt_register_library(NAME hdf5
                       INCLUDES ${HDF5_INCLUDE_DIRS}
                       LIBRARIES ${HDF5_LIBRARY} )
endif()


if (ATK_DIR)
  include(cmake/thirdparty/FindATK.cmake)
  blt_register_library( NAME sidre
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  sidre)

  blt_register_library( NAME slic
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  slic)
endif()


if (RAJA_DIR)
  include(cmake/thirdparty/FindRAJA.cmake)
 
   if (NOT RAJA_FOUND)
	MESSAGE(FATAL_ERROR "RAJA not found. Do not define RAJA_DIR to attempt internal build of RAJA.")
   endif()
elseif(EXISTS ${PROJECT_SOURCE_DIR}/thirdparty/raja/CMakeLists.txt)
  #set(RAJA_DIR ${PROJECT_SOURCE_DIR}/thirdparty/raja-install/share/raja/cmake)
  set(RAJA_DIR ${PROJECT_SOURCE_DIR}/thirdparty/raja-install)
  #find_package(RAJA)
 include(${PROJECT_SOURCE_DIR}/cmake/thirdparty/FindRAJA.cmake) 
  if (NOT RAJA_FOUND)
     MESSAGE(FATAL_ERROR "RAJA not found locally in ${RAJA_DIR}. Maybe you need to run chairajabuild in src/thirdparty?")
  endif()
else()
   MESSAGE(FATAL_ERROR "RAJA_DIR not defined and no local raja found in ${CMAKE_CURRENT_LIST_DIR}/raja/CMakeLists.txt.")
endif(RAJA_DIR)

 blt_register_library( NAME RAJA
                       INCLUDES ${RAJA_INCLUDE_DIRS} 
                       LIBRARIES  RAJA)

#if (UNCRUSTIFY_EXECUTABLE)
  include(cmake/blt/cmake/thirdparty/FindUncrustify.cmake)
#endif()

if (CHAI_DIR)
   MESSAGE("REGISTERING CHAI WITH BLT "  ${CHAI_DIR}/src)
#   set(CHAI_HEADERS ${CHAI_DIR}/src/ManagedArray.hpp ${CHAI_DIR}/src/ManagedClass.hpp ${CHAI_DIR}/src/resource_manager.hpp ${CHAI_DIR}/src/YALLPolicies.h
#                    ${CHAI_DIR}/src/execution_policies.hpp )
   set(CHAI_HEADERS ${PROJECT_SOURCE_DIR}/${CHAI_DIR}/src/)
   blt_register_library(NAME chai 
                     INCLUDES ${CHAI_HEADERS}
                     LIBRARIES ${PROJECT_SOURCE_DIR}/${CHAI_DIR}/src/libchai.a
		     DEFINES -DCHAI_DISABLE_RM=1)
endif()
