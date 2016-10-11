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







################################
# RAJA
################################
set(RAJA_LOCAL_DIR ${PROJECT_SOURCE_DIR}/thirdparty/raja-install)
if( NOT EXISTS ${RAJA_LOCAL_DIR} AND NOT EXISTS ${RAJA_DIR})
    MESSAGE(FATAL_ERROR "RAJA_DIR not defined and no locally built raja found in ${PROJECT_SOURCE_DIR}/thirdparty/raja-install. Maybe you need to run chairajabuild in src/thirdparty?")
else()
    if(EXISTS ${RAJA_LOCAL_DIR})
        MESSAGE( "Using local RAJA found at ${RAJA_LOCAL_DIR}")
        set(RAJA_DIR ${RAJA_LOCAL_DIR})
        include(${PROJECT_SOURCE_DIR}/cmake/thirdparty/FindRAJA.cmake) 
        if (NOT RAJA_FOUND)
            MESSAGE(FATAL_ERROR "RAJA not found locally in ${RAJA_DIR}. Maybe you need to run chairajabuild in src/thirdparty?")
        endif()
    elseif(EXISTS ${RAJA_DIR})
        MESSAGE( "Using system RAJA found at ${RAJA_DIR}")
        include(${PROJECT_SOURCE_DIR}/cmake/thirdparty/FindRAJA.cmake)
        if (NOT RAJA_FOUND)
            MESSAGE(FATAL_ERROR "RAJA not found in ${RAJA_DIR}. Maybe you need to build it")
        endif()
    endif()
endif()

blt_register_library( NAME RAJA
                      INCLUDES ${RAJA_INCLUDE_DIRS} 
                      LIBRARIES  RAJA)





################################
# CHAI
################################
set(CHAI_LOCAL_DIR ${PROJECT_SOURCE_DIR}/thirdparty/chai)
if( NOT EXISTS ${CHAI_LOCAL_DIR} AND NOT EXISTS ${CHAI_DIR})
    MESSAGE(FATAL_ERROR "CHAI_DIR not defined and no locally built chai found in ${PROJECT_SOURCE_DIR}/thirdparty/chai. Maybe you need to run chairajabuild in src/thirdparty?")
else()
    if(EXISTS ${CHAI_LOCAL_DIR})
        MESSAGE( "Using local CHAI found at ${CHAI_LOCAL_DIR}")
        set(CHAI_DIR ${CHAI_LOCAL_DIR})
        include(${PROJECT_SOURCE_DIR}/cmake/thirdparty/FindCHAI.cmake) 
        if (NOT CHAI_FOUND)
            MESSAGE(FATAL_ERROR "CHAI not found locally in ${CHAI_DIR}. Maybe you need to run chairajabuild in src/thirdparty?")
        endif()
    elseif(EXISTS ${CHAI_DIR})
        MESSAGE( "Using system CHAI found at ${CHAI_DIR}")
        include(${PROJECT_SOURCE_DIR}/cmake/thirdparty/FindCHAI.cmake)
        if (NOT CHAI_FOUND)
            MESSAGE(FATAL_ERROR "CHAI not found in ${CHAI_DIR}. Maybe you need to build it")
        endif()
    endif()
endif()

blt_register_library(NAME chai 
                     INCLUDES ${CHAI_INCLUDE_DIRS}
                     LIBRARIES ${CHAI_DIR}/src/libchai.a
                     DEFINES -DCHAI_DISABLE_RM=1 )
                     



#if (UNCRUSTIFY_EXECUTABLE)
  include(cmake/blt/cmake/thirdparty/FindUncrustify.cmake)
#endif()



