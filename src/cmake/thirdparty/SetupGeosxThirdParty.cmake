message("\nProcessing SetupGeosxThirdParty.cmake")


####################################
# 3rd Party Dependencies
####################################
include(ExternalProject)


message("GEOSX_TPL_DIR=${GEOSX_TPL_DIR}")
if( NOT GEOSX_TPL_DIR )
    message("GEOSX_TPL_ROOT_DIR=${GEOSX_TPL_ROOT_DIR}")
    get_filename_component( TEMP_DIR "${CMAKE_INSTALL_PREFIX}" NAME)
    if( CMAKE_BUILD_TYPE MATCHES "Debug" )
        string(REGEX REPLACE "debug" "release" TEMP_DIR2 ${TEMP_DIR})
    endif()
    if( CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo" )
        string(REGEX REPLACE "relwithdebinfo" "release" TEMP_DIR2 ${TEMP_DIR})
    endif()
    if( CMAKE_BUILD_TYPE MATCHES "Release" )
        set(TEMP_DIR2 ${TEMP_DIR})
    endif()
    set( GEOSX_TPL_DIR "${GEOSX_TPL_ROOT_DIR}/${TEMP_DIR2}" )
endif()
message("GEOSX_TPL_DIR=${GEOSX_TPL_DIR}")

set( thirdPartyLibs "")


################################
# Conduit
################################
if (NOT EXISTS ${CONDUIT_DIR})
    message(INFO ": Using conduit from thirdPartyLibs")
    set(CONDUIT_DIR ${GEOSX_TPL_DIR}/conduit CACHE PATH "")
endif()

if (EXISTS ${CONDUIT_DIR})
  message( "CONDUIT_DIR = ${CONDUIT_DIR}" )
  include( cmake/thirdparty/FindConduit.cmake )
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
                        
  set( CONDUIT_FOUND ON CACHE BOOL "" )
  set( thirdPartyLibs ${thirdPartyLibs} conduit conduit_blueprint conduit_relay )
else()
  set( CONDUIT_FOUND OFF CACHE BOOL "" )
  message( "Not using conduit" )
endif()

################################
# AXOM
################################
include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindATK.cmake)

################################
# uncrustify
################################
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
include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindCHAI.cmake)

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

    set( caliper_lib_list caliper caliper-reader caliper-common gotcha caliper-mpi )
                        
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
if (NOT EXISTS ${MATHPRESSO_DIR})
    message(INFO ": Using mathpresso from thirdPartyLibs")
    set(MATHPRESSO_DIR ${GEOSX_TPL_DIR}/mathpresso CACHE PATH "")
endif()

if (EXISTS ${MATHPRESSO_DIR})
    message( INFO ": setting up asmjit" )
    set(ASMJIT_LOCAL_DIR ${PROJECT_BINARY_DIR}/thirdparty/asmjit/src/asmjit)
    set(ASMJIT_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/asmjit)

    message( INFO ": setting up MathPresso" )
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
    
    set(MATHPRESSO_FOUND ON CACHE BOOL "")
    set( thirdPartyLibs ${thirdPartyLibs} mathpresso )
else()
    set(MATHPRESSO_FOUND OFF CACHE BOOL "")
    message(INFO ": Not using mathpresso")
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
  
  list(REMOVE_ITEM Trilinos_LIBRARIES "gtest")
  list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
  message(STATUS "Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
  message(STATUS "Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
  
  blt_register_library( NAME trilinos
                        INCLUDES ${Trilinos_INCLUDE_DIRS} 
                        LIBRARIES ${Trilinos_LIBRARIES}
                        TREAT_INCLUDES_AS_SYSTEM ON )
  set( thirdPartyLibs ${thirdPartyLibs} trilinos )  

endif()

################################
# METIS
################################
if( ENABLE_METIS )
    message( INFO ": setting up METIS" )

    set(METIS_DIR ${GEOSX_TPL_DIR}/metis)

    find_path( METIS_INCLUDE_DIRS metis.h
       PATHS  ${METIS_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

   find_library( METIS_LIBRARY NAMES metis
          PATHS ${METIS_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    message( "METIS_INCLUDE_DIRS = ${METIS_INCLUDE_DIRS}" )
    message( "METIS_LIBRARY = ${METIS_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(METIS DEFAULT_MSG METIS_INCLUDE_DIRS METIS_LIBRARY )
    if (NOT METIS_FOUND)
        message(FATAL_ERROR ": METIS not found in ${METIS_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME metis
                            INCLUDES ${METIS_INCLUDE_DIRS} 
                LIBRARIES ${METIS_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON )

    set( thirdPartyLibs ${thirdPartyLibs} metis )
endif()

################################
# PARMETIS
################################
if( ENABLE_PARMETIS )
    message( INFO ": setting up PARMETIS" )

    set(PARMETIS_DIR ${GEOSX_TPL_DIR}/parmetis)

    find_path( PARMETIS_INCLUDE_DIRS parmetis.h
           PATHS  ${PARMETIS_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

   find_library( PARMETIS_LIBRARY NAMES parmetis
              PATHS ${PARMETIS_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    message( "PARMETIS_INCLUDE_DIRS = ${PARMETIS_INCLUDE_DIRS}" )
    message( "PARMETIS_LIBRARY = ${PARMETIS_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(PARMETIS  DEFAULT_MSG
                                    PARMETIS_INCLUDE_DIRS
                                    PARMETIS_LIBRARY )
    if (NOT PARMETIS_FOUND)
        message(FATAL_ERROR ": PARMETIS not found in ${PARMETIS_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME parmetis
                        INCLUDES ${PARMETIS_INCLUDE_DIRS} 
                        LIBRARIES ${PARMETIS_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON )

    set( thirdPartyLibs ${thirdPartyLibs} parmetis )
endif()

################################
# SUPERLU_DIST
################################
if( ENABLE_SUPERLU_DIST)
    message( INFO ": setting up SUPERLU_DIST" )

    set(SUPERLU_DIST_DIR ${GEOSX_TPL_DIR}/superlu_dist)

    find_path( SUPERLU_DIST_INCLUDE_DIRS superlu_defs.h
        PATHS  ${SUPERLU_DIST_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

   find_library( SUPERLU_DIST_LIBRARY NAMES superlu_dist
          PATHS ${SUPERLU_DIST_DIR}/lib PATHS ${SUPERLU_DIST_DIR}/lib64
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    message( "SUPERLU_DIST_INCLUDE_DIRS = ${SUPERLU_DIST_INCLUDE_DIRS}" )
    message( "SUPERLU_DIST_LIBRARY = ${SUPERLU_DIST_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(SUPERLU_DIST  DEFAULT_MSG
                                SUPERLU_DIST_INCLUDE_DIRS
                    SUPERLU_DIST_LIBRARY )
    if (NOT SUPERLU_DIST_FOUND)
        message(FATAL_ERROR ": SUPERLU_DIST not found in ${SUPERLU_DIST_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME superlu_dist
                            INCLUDES ${SUPERLU_DIST_INCLUDE_DIRS} 
                    LIBRARIES ${SUPERLU_DIST_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON )

    set( thirdPartyLibs ${thirdPartyLibs} superlu_dist )
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

    find_library( HYPRE_LIBRARY NAMES HYPRE
                PATHS ${HYPRE_DIR}/lib
                NO_DEFAULT_PATH
                NO_CMAKE_ENVIRONMENT_PATH
                NO_CMAKE_PATH
                NO_SYSTEM_ENVIRONMENT_PATH
                NO_CMAKE_SYSTEM_PATH)

    message( "HYPRE_INCLUDE_DIRS = ${HYPRE_INCLUDE_DIRS}" )
    message( "HYPRE_LIBRARY = ${HYPRE_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(HYPRE  DEFAULT_MSG
                                    HYPRE_INCLUDE_DIRS
                                    HYPRE_LIBRARY )
    if (NOT HYPRE_FOUND)
        message(FATAL_ERROR ": HYPRE not found in ${HYPRE_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME hypre
                        INCLUDES ${HYPRE_INCLUDE_DIRS} 
                        LIBRARIES ${HYPRE_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON )

    set( thirdPartyLibs ${thirdPartyLibs} hypre )
endif()


################################
# unsrustify
################################
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
