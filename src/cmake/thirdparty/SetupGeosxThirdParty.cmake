message("\nProcessing SetupGeosxThirdParty.cmake")


####################################
# 3rd Party Dependencies
####################################
include(ExternalProject)


################################
# Conduit
################################
if (CONDUIT_DIR)
  include(cmake/thirdparty/FindConduit.cmake)
  blt_register_library( NAME conduit
                        INCLUDES ${CONDUIT_INCLUDE_DIRS} 
                        LIBRARIES  conduit
                        TREAT_INCLUDES_AS_SYSTEM ON )
  blt_register_library( NAME conduit_io
                        INCLUDES ${CONDUIT_INCLUDE_DIRS}
                        LIBRARIES  conduit_io
                        TREAT_INCLUDES_AS_SYSTEM ON )
endif()


################################
# HDF5
################################
if (HDF5_DIR)
  include(cmake/thirdparty/FindHDF5.cmake)
  blt_register_library(NAME hdf5
                       INCLUDES ${HDF5_INCLUDE_DIRS}
                       LIBRARIES ${HDF5_LIBRARIES} 
                       TREAT_INCLUDES_AS_SYSTEM ON )
endif()


if (ATK_DIR)
  include(cmake/thirdparty/FindATK.cmake)
  blt_register_library( NAME sidre
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  sidre
                        TREAT_INCLUDES_AS_SYSTEM ON )

  blt_register_library( NAME spio
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  spio
                        TREAT_INCLUDES_AS_SYSTEM ON)

  blt_register_library( NAME slic
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  slic
                        TREAT_INCLUDES_AS_SYSTEM ON)
endif()






get_filename_component( TEMP_DIR "${CMAKE_INSTALL_PREFIX}" NAME)
set( GEOSX_TPL_DIR "../thirdPartyLibs/${TEMP_DIR}" )
################################
# SILO
################################
if( EXISTS ${SILO_DIR})
    message("Using system SILO found at ${SILO_DIR}")
else()
    message(INFO ": Using SILO from thirdPartyLibs")
    set(SILO_DIR ${GEOSX_TPL_DIR}/silo)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindSilo.cmake)
if (NOT SILO_FOUND)
    message(FATAL_ERROR ": SILO not found in ${SILO_DIR}. Maybe you need to build it")
endif()
blt_register_library( NAME silo
                      INCLUDES ${SILO_INCLUDE_DIRS}
                      LIBRARIES ${SILO_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON
                      DEPENDS_ON hdf5 )





################################
# RAJA
################################
if( EXISTS ${RAJA_DIR})
    message("Using system RAJA found at ${RAJA_DIR}")
else()
    message(INFO ": Using SILO from thirdPartyLibs")
    set(RAJA_DIR ${GEOSX_TPL_DIR}/raja)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindRAJA.cmake)
if (NOT RAJA_FOUND)
    message(FATAL_ERROR ": RAJA not found in ${RAJA_DIR}. Maybe you need to build it")
endif()    
blt_register_library( NAME raja
                      INCLUDES ${RAJA_INCLUDE_DIRS}
                      LIBRARIES ${RAJA_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON )




################################
# CHAI
################################
if( EXISTS ${CHAI_DIR})
    message("Using system CHAI found at ${CHAI_DIR}")
else()
    message(INFO ": Using SILO from thirdPartyLibs")
    set(CHAI_DIR ${GEOSX_TPL_DIR}/chai)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindCHAI.cmake)
if (NOT CHAI_FOUND)
    message(FATAL_ERROR ": CHAI not found in ${CHAI_DIR}. Maybe you need to build it")
endif()    
blt_register_library( NAME chai
                      INCLUDES ${CHAI_INCLUDE_DIRS}
                      LIBRARIES ${CHAI_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON )


################################
# FPARSER
################################
if( USE_FPARSER )

message(INFO ": Using SILO from thirdPartyLibs")
set(FPARSER_INSTALL_DIR ${GEOSX_TPL_DIR}/fparser)

find_path( FPARSER_INCLUDE_DIRS fparser.h
           PATHS  ${FPARSER_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( FPARSER_LIBRARY NAMES fparser
              PATHS ${FPARSER_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FPARSER  DEFAULT_MSG
                                  FPARSER_INCLUDE_DIRS
                                  FPARSER_LIBRARY )
if (NOT FPARSER_FOUND)
    message(FATAL_ERROR ": FPARSER not found in ${FPARSER_DIR}. Maybe you need to build it")
endif()

blt_register_library( NAME fparser
                      INCLUDES ${FPARSER_INCUDE_DIRS} 
                      LIBRARIES ${FPARSER_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON )


endif()






################################
# CALIPER
################################
if( ENABLE_CALIPER )
message( INFO ": setting up caliper" )

set(CALIPER_DIR ${GEOSX_TPL_DIR}/caliper)
find_path( CALIPER_INCLUDE_DIRS caliper.h
           PATHS  ${CALIPER_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( CALIPER_LIBRARY NAMES caliper
              PATHS ${CALIPER_DIR}/lib64
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CALIPER  DEFAULT_MSG
                                  CALIPER_INCLUDE_DIRS
                                  CALIPER_LIBRARY )
                                  
if (NOT CALIPER_FOUND)
    message(FATAL_ERROR ": CALIPER not found in ${CALIPER_DIR}. Maybe you need to build it")
endif()
blt_register_library( NAME caliper
                      INCLUDES ${CALIPER_INSTALL_DIR}/include 
                      LIBRARIES ${CALIPER_INSTALL_DIR}/lib/libcaliper.a 
                      LIBRARIES ${CALIPER_INSTALL_DIR}/lib/libcaliper-common.a
                      LIBRARIES ${CALIPER_INSTALL_DIR}/lib/libcaliper-reader.a
                      TREAT_INCLUDES_AS_SYSTEM ON )

endif()



################################
# ASMJIT / MATHPRESSO
################################
set(ENABLE_MATHPRESSO ON CACHE BOOL  "Enables mathpresso Plugin")
if( ENABLE_MATHPRESSO )

message( INFO ": setting up asmjit" )
set(ASMJIT_LOCAL_DIR ${PROJECT_BINARY_DIR}/thirdparty/asmjit/src/asmjit)
set(ASMJIT_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/asmjit)

message( INFO ": setting up MathPresso" )
set(MATHPRESSO_DIR ${GEOSX_TPL_DIR}/mathpresso)

find_path( MATHPRESSO_INCLUDE_DIRS mathpresso/mathpresso.h
           PATHS  ${MATHPRESSO_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( MATHPRESSO_LIBRARY NAMES mathpresso
              PATHS ${MATHPRESSO_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MATHPRESSO  DEFAULT_MSG
                                  MATHPRESSO_INCLUDE_DIRS
                                  MATHPRESSO_LIBRARY )
                                  
if (NOT MATHPRESSO_FOUND)
    message(FATAL_ERROR ": MATHPRESSO not found in ${MATHPRESSO_DIR}. Maybe you need to build it")
endif()

blt_register_library( NAME mathpresso
                      INCLUDES ${MATHPRESSO_INCLUDE_DIRS}
                      LIBRARIES ${MATHPRESSO_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON )

endif()



################################
# PUGIXML
################################
message( INFO ": setting up pugixml" )
set(PUGIXML_DIR ${GEOSX_TPL_DIR}/pugixml)


find_path( PUGIXML_INCLUDE_DIRS NAMES pugixml.hpp
           PATHS  ${PUGIXML_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    find_library( PUGIXML_LIBRARY NAMES pugixml
                  PATHS ${PUGIXML_DIR}/lib
                  NO_DEFAULT_PATH
                  NO_CMAKE_ENVIRONMENT_PATH
                  NO_CMAKE_PATH
                  NO_SYSTEM_ENVIRONMENT_PATH
                  NO_CMAKE_SYSTEM_PATH)

else()
    find_library( PUGIXML_LIBRARY NAMES pugixml
                  PATHS ${PUGIXML_DIR}/lib64
                  NO_DEFAULT_PATH
                  NO_CMAKE_ENVIRONMENT_PATH
                  NO_CMAKE_PATH
                  NO_SYSTEM_ENVIRONMENT_PATH
                  NO_CMAKE_SYSTEM_PATH)
endif()

find_package_handle_standard_args(PUGIXML  DEFAULT_MSG
                                  PUGIXML_INCLUDE_DIRS
                                  PUGIXML_LIBRARY )


if (NOT PUGIXML_FOUND)
    message(FATAL_ERROR ": PUGIXML not found in ${PUGIXML_DIR}. Maybe you need to build it")
endif()

blt_register_library( NAME pugixml
                      INCLUDES ${PUGIXML_INCLUDE_DIRS}
                      LIBRARIES ${PUGIXML_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON )



################################
# TRILINOS
################################
#if( ENABLE_TRILINOS )
message( INFO ": setting up TRILINOS" )

set(TRILINOS_DIR ${GEOSX_TPL_DIR}/trilinos)
include(${TRILINOS_DIR}/lib/cmake/Trilinos/TrilinosConfig.cmake)

blt_register_library( NAME trilinos
                      INCLUDES ${Trilinos_INCLUDE_DIRS} 
                      LIBRARIES ${Trilinos_LIBRARIES}
                      TREAT_INCLUDES_AS_SYSTEM ON )

#endif()



if (UNCRUSTIFY_EXECUTABLE)
  include(cmake/blt/cmake/thirdparty/FindUncrustify.cmake)
endif()

message("Leaving SetupGeosxThirdParty.cmake\n")


