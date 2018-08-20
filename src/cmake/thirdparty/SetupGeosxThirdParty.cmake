message("\nProcessing SetupGeosxThirdParty.cmake")


####################################
# 3rd Party Dependencies
####################################
include(ExternalProject)


message("GEOSX_TPL_DIR=${GEOSX_TPL_DIR}")
if( NOT GEOSX_TPL_DIR )
    message("GEOSX_TPL_ROOT_DIR=${GEOSX_TPL_ROOT_DIR}")
    get_filename_component( TEMP_DIR "${CMAKE_INSTALL_PREFIX}" NAME)
    string(REGEX REPLACE "debug" "release" TEMP_DIR2 ${TEMP_DIR})
    set( GEOSX_TPL_DIR "${GEOSX_TPL_ROOT_DIR}/${TEMP_DIR2}" )
endif()
message("GEOSX_TPL_DIR=${GEOSX_TPL_DIR}")

set(ATK_DIR "${GEOSX_TPL_DIR}/axom" CACHE PATH "")
set(CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit" CACHE PATH "")

set( thirdPartyLibs "")


################################
# Conduit
################################
if (EXISTS ${CONDUIT_DIR})
  message( "CONDUIT_DIR = ${CONDUIT_DIR}" )
  include(cmake/thirdparty/FindConduit.cmake)
  blt_register_library( NAME conduit
                        INCLUDES ${CONDUIT_INCLUDE_DIRS} 
                        LIBRARIES  conduit
                        TREAT_INCLUDES_AS_SYSTEM ON )

  blt_register_library( NAME conduit_blueprint
                        INCLUDES ${CONDUIT_INCLUDE_DIRS}
                        LIBRARIES conduit_blueprint
                        TREAT_INCLUDES_AS_SYSTEM ON )

  blt_register_library( NAME conduit_relay
                        INCLUDES ${CONDUIT_INCLUDE_DIRS}
                        LIBRARIES conduit_relay
                        TREAT_INCLUDES_AS_SYSTEM ON )
                        
  set( thirdPartyLibs ${thirdPartyLibs} conduit conduit_blueprint conduit_relay )

else()
  message( "Not using conduit" )
endif()


################################
# AXOM
################################
if (EXISTS ${ATK_DIR})
  message( "ATK_DIR = ${ATK_DIR}" )
  include(cmake/thirdparty/FindATK.cmake)
  blt_register_library( NAME sidre
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  sidre
                        TREAT_INCLUDES_AS_SYSTEM ON )

  blt_register_library( NAME slic
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  slic
                        TREAT_INCLUDES_AS_SYSTEM ON)
                        
  set( thirdPartyLibs ${thirdPartyLibs} sidre slic )
else()
  message( "Not using axom" )
endif()


set(UNCRUSTIFY_EXECUTABLE "${GEOSX_TPL_DIR}/uncrustify/bin/uncrustify" CACHE PATH "" FORCE )


################################
# HDF5
################################
if( EXISTS ${HDF5_DIR})
    message("Using system HDF5 found at ${HDF5_DIR}")
else()
    message(INFO ": Using HDF5 from thirdPartyLibs")
    set(HDF5_DIR ${GEOSX_TPL_DIR}/hdf5)
endif()

if (HDF5_DIR)
  include(cmake/thirdparty/FindHDF5.cmake)
  blt_register_library(NAME hdf5
                       INCLUDES ${HDF5_INCLUDE_DIRS}
                       LIBRARIES ${HDF5_LIBRARIES} dl
                       TREAT_INCLUDES_AS_SYSTEM ON )
                       
  set( thirdPartyLibs ${thirdPartyLibs} hdf5 )

endif()

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

set( thirdPartyLibs ${thirdPartyLibs} silo )  




################################
# RAJA
################################
if( EXISTS ${RAJA_DIR})
    message("Using system RAJA found at ${RAJA_DIR}")
else()
    message(INFO ": Using RAJA from thirdPartyLibs")
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

set( thirdPartyLibs ${thirdPartyLibs} raja )  



################################
# CHAI
################################
if( EXISTS ${CHAI_DIR})
    message("Using system CHAI found at ${CHAI_DIR}")
    set(CHAI_FOUND TRUE)
else()
    message(INFO ": Using CHAI from thirdPartyLibs")
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

set( thirdPartyLibs ${thirdPartyLibs} chai )  


################################
# FPARSER
################################
if( USE_FPARSER )

message(INFO ": Using FPARSER from thirdPartyLibs")
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

set( thirdPartyLibs ${thirdPartyLibs} fparser )  

endif()






################################
# CALIPER
################################
if( ENABLE_CALIPER )
message( INFO ": setting up caliper" )

if( EXISTS ${CALIPER_DIR} )
    message( INFO "Found system caliper" )
    message("Using system CALIPER found at ${CALIPER_DIR}")
    set(CALIPER_FOUND TRUE)    
else()
    message(INFO ": Using CALIPER from thirdPartyLibs")
    set(CALIPER_DIR ${GEOSX_TPL_DIR}/caliper)
    set(CALIPER_FOUND TRUE)
endif()

find_path( CALIPER_INCLUDE_DIRS caliper/Caliper.h
           PATHS  ${CALIPER_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

set( caliper_lib_list caliper caliper-reader caliper-common  gotcha )
                       
message(INFO "looking for libs in ${CALIPER_DIR}")
blt_find_libraries( FOUND_LIBS CALIPER_LIBRARIES
                    NAMES ${caliper_lib_list}
                    PATHS ${CALIPER_DIR}/lib ${CALIPER_DIR}/lib64
                   )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CALIPER  DEFAULT_MSG
                                  CALIPER_INCLUDE_DIRS
                                  CALIPER_LIBRARIES )


if (NOT CALIPER_FOUND)
    message(FATAL_ERROR ": CALIPER not found in ${CALIPER_DIR}. Maybe you need to build it")
else() 
    message("CALIPER_INCLUDE_DIRS = ${CALIPER_INCLUDE_DIRS}")
    message("CALIPER_LIBRARIES = ${CALIPER_LIBRARIES}")    
endif()
blt_register_library( NAME caliper
                      INCLUDES ${CALIPER_INCLUDE_DIRS}
                      LIBRARIES ${CALIPER_LIBRARIES}
                      TREAT_INCLUDES_AS_SYSTEM ON )

set( thirdPartyLibs ${thirdPartyLibs} caliper )  
                      
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

set( thirdPartyLibs ${thirdPartyLibs} mathpresso )  

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
                  PATHS ${PUGIXML_DIR}/lib64 ${PUGIXML_DIR}/lib
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

set( thirdPartyLibs ${thirdPartyLibs} pugixml )  



################################
# TRILINOS
################################
if( ENABLE_TRILINOS )

  if(EXISTS ${TRILINOS_DIR})
  
  else()
      message( INFO ": setting up TRILINOS" )
      set(TRILINOS_DIR ${GEOSX_TPL_DIR}/trilinos)
  endif()
  
  include(${TRILINOS_DIR}/lib/cmake/Trilinos/TrilinosConfig.cmake)
  
  
  blt_register_library( NAME trilinos
                        INCLUDES ${Trilinos_INCLUDE_DIRS} 
                        LIBRARIES ${Trilinos_LIBRARIES}
                        TREAT_INCLUDES_AS_SYSTEM ON )
  set( thirdPartyLibs ${thirdPartyLibs} trilinos )  

endif()

################################
# HYPRE
################################
if( ENABLE_HYPRE )
message( INFO ": setting up HYPRE" )

set(HYPRE_DIR ${GEOSX_TPL_DIR}/hypre)

find_path( HYPRE_INCLUDE_DIRS HYPRE.h
           PATHS  ${HYPRE_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( HYPRE_LIBRARY NAMES hypre
              PATHS ${HYPRE_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE  DEFAULT_MSG
                                  HYPRE_INCLUDE_DIRS
                                  HYPRE_LIBRARY )
if (NOT HYPRE_FOUND)
    message(FATAL_ERROR ": HYPRE not found in ${HYPRE_DIR}. Maybe you need to build it")
endif()

blt_register_library( NAME hypre
                      INCLUDES ${HYPRE_INCLUDE_DIRS} 
                      INCLUDES ${HYPRE_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON )

set( thirdPartyLibs ${thirdPartyLibs} hypre )  

endif()



#if (UNCRUSTIFY_EXECUTABLE)
#  include(cmake/blt/cmake/thirdparty/FindUncrustify.cmake)
#endif()
message("UNCRUSTIFY_FOUND = ${UNCRUSTIFY_FOUND}")
if(UNCRUSTIFY_FOUND)
    # targets for verifying formatting
    if( NOT TARGET uncrustify_check )
        add_custom_target(uncrustify_check)
        add_dependencies(check uncrustify_check)
    endif()

    # targets for modifying formatting
    if( NOT TARGET uncrustify_style )
        add_custom_target(uncrustify_style)
        add_dependencies(style uncrustify_style)
    endif()
    
endif()

message("Leaving SetupGeosxThirdParty.cmake\n")

