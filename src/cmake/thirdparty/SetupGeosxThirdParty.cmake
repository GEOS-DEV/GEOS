####################################
# 3rd Party Dependencies
####################################

set( thirdPartyLibs "")

################################
# BLAS/LAPACK
################################

include(cmake/thirdparty/FindMathLibraries.cmake)

blt_register_library(NAME blas
                      TREAT_INCLUDES_AS_SYSTEM ON
                      LIBRARIES ${BLAS_LIBRARIES})

blt_register_library(NAME lapack
                     DEPENDS_ON blas
                     TREAT_INCLUDES_AS_SYSTEM ON
                     LIBRARIES ${LAPACK_LIBRARIES})

################################
# Intel MKL
################################
if(ENABLE_MKL)
    message(STATUS "Using Intel MKL")

    blt_register_library(NAME mkl
                         INCLUDES ${MKL_INCLUDE_DIRS}
                         LIBRARIES ${MKL_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(TRILINOS_DEPENDS mkl)
    set(thirdPartyLibs ${thirdPartyLibs} mkl)
endif()

################################
# IBM ESSL
################################
if(ENABLE_ESSL)
    message(STATUS "Using up IBM ESSL")

    blt_register_library(NAME essl
                         INCLUDES ${ESSL_INCLUDE_DIRS}
                         LIBRARIES ${ESSL_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(TRILINOS_DEPENDS essl)
    set(thirdPartyLibs ${thirdPartyLibs} essl)
else()
    set(TRILINOS_DEPENDS blas lapack)
    set(thirdPartyLibs ${thirdPartyLibs} blas lapack)
endif()

################################
# HDF5
################################
if(DEFINED HDF5_DIR)
    message(STATUS "HDF5_DIR = ${HDF5_DIR}")

    include(cmake/thirdparty/FindHDF5.cmake)
    blt_register_library(NAME hdf5
                        INCLUDES ${HDF5_INCLUDE_DIRS}
                        LIBRARIES ${HDF5_LIBRARIES} dl
                        TREAT_INCLUDES_AS_SYSTEM ON )

    set(ENABLE_HDF5 ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} hdf5)
else()
    message(FATAL_ERROR "GEOSX requires hdf5, set HDF5_DIR to the hdf5 installation directory.")
endif()

################################
# Conduit
################################
if(DEFINED CONDUIT_DIR)
    message(STATUS "CONDUIT_DIR = ${CONDUIT_DIR}" )
    include(cmake/thirdparty/FindConduit.cmake )
    blt_register_library(NAME conduit
                         INCLUDES ${CONDUIT_INCLUDE_DIRS} 
                         LIBRARIES conduit
                         TREAT_INCLUDES_AS_SYSTEM ON)

    blt_register_library(NAME conduit_blueprint
                         INCLUDES ${CONDUIT_INCLUDE_DIRS}
                         LIBRARIES conduit_blueprint
                         TREAT_INCLUDES_AS_SYSTEM ON)

    blt_register_library(NAME conduit_relay
                         INCLUDES ${CONDUIT_INCLUDE_DIRS}
                         LIBRARIES conduit_relay
                         TREAT_INCLUDES_AS_SYSTEM ON)

    # Removes link to Thread library (removes "-pthread" when linking to relay)
    get_target_property(_relay_libraries conduit_relay INTERFACE_LINK_LIBRARIES)
    set(_new_relay_libraries)
    string(REPLACE "Threads::Threads;" "" _new_relay_libraries "${_relay_libraries}")
    set_target_properties(conduit_relay
                          PROPERTIES INTERFACE_LINK_LIBRARIES 
                          "${_new_relay_libraries}")

    set(ENABLE_CONDUIT ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} conduit conduit_blueprint conduit_relay)
else()
    message(FATAL_ERROR "GEOSX requires conduit, set CONDUIT_DIR to the conduit installation directory.")
endif()

################################
# SILO
################################
if(DEFINED SILO_DIR)
    message(STATUS "SILO_DIR = ${SILO_DIR}")
    include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindSilo.cmake)
    if (NOT SILO_FOUND)
        message(FATAL_ERROR "SILO not found in ${SILO_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME silo
                        INCLUDES ${SILO_INCLUDE_DIRS}
                        LIBRARIES ${SILO_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON
                        DEPENDS_ON hdf5 )

    set(ENABLE_SILO ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} silo)
else()
    message(FATAL_ERROR "GEOSX requires Silo, set SILO_DIR to the Silo installation directory.")
endif()

################################
# PUGIXML
################################
if(DEFINED PUGIXML_DIR)
    message(STATUS "PUGIXML_DIR = ${PUGIXML_DIR}")

    find_path(PUGIXML_INCLUDE_DIRS NAMES pugixml.hpp
              PATHS ${PUGIXML_DIR}/include
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        find_library(PUGIXML_LIBRARY NAMES pugixml
                     PATHS ${PUGIXML_DIR}/lib
                     NO_DEFAULT_PATH
                     NO_CMAKE_ENVIRONMENT_PATH
                     NO_CMAKE_PATH
                     NO_SYSTEM_ENVIRONMENT_PATH
                     NO_CMAKE_SYSTEM_PATH)

    else()
        find_library(PUGIXML_LIBRARY NAMES pugixml
                     PATHS ${PUGIXML_DIR}/lib64 ${PUGIXML_DIR}/lib
                     NO_DEFAULT_PATH
                     NO_CMAKE_ENVIRONMENT_PATH
                     NO_CMAKE_PATH
                     NO_SYSTEM_ENVIRONMENT_PATH
                     NO_CMAKE_SYSTEM_PATH)
    endif()

    find_package_handle_standard_args(PUGIXML DEFAULT_MSG
                                      PUGIXML_INCLUDE_DIRS
                                      PUGIXML_LIBRARY )


    if (NOT PUGIXML_FOUND)
        message(FATAL_ERROR "PUGIXML not found in ${PUGIXML_DIR}. Maybe you need to build it")
    endif()

    blt_register_library(NAME pugixml
                         INCLUDES ${PUGIXML_INCLUDE_DIRS}
                         LIBRARIES ${PUGIXML_LIBRARY}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(thirdPartyLibs ${thirdPartyLibs} pugixml)
else()
    message(FATAL_ERROR "GEOSX requires pugixml, set PUGIXML_DIR to the pugixml installation directory.")
endif()

################################
# RAJA
################################
if(DEFINED RAJA_DIR)
    message(STATUS "RAJA_DIR = ${RAJA_DIR}")

    find_package(RAJA REQUIRED PATHS ${RAJA_DIR})

    get_target_property(RAJA_INCLUDE_DIRS RAJA INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(RAJA
                          PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${RAJA_INCLUDE_DIRS}")

    set(ENABLE_RAJA ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} RAJA )
else()
    message(FATAL_ERROR "GEOSX requires RAJA, set RAJA_DIR to the RAJA installation directory.")
endif()

################################
# Umpire
################################
if(DEFINED UMPIRE_DIR)
    message(STATUS "UMPIRE_DIR = ${UMPIRE_DIR}")

    find_package(umpire REQUIRED
                 PATHS ${UMPIRE_DIR})

    set(ENABLE_UMPIRE ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} umpire)
else()
    message(FATAL_ERROR "GEOSX requires Umpire, set UMPIRE_DIR to the Umpire installation directory.")
endif()

################################
# CHAI
################################
if(DEFINED CHAI_DIR)
    message(STATUS "CHAI_DIR = ${CHAI_DIR}")

    find_package(chai REQUIRED
                 PATHS ${CHAI_DIR})

    get_target_property(CHAI_INCLUDE_DIRS chai INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(chai
                        PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CHAI_INCLUDE_DIRS}")

    set(ENABLE_CHAI ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} chai)
else()
    message(FATAL_ERROR "GEOSX requires CHAI, set CHAI_DIR to the CHAI installation directory.")
endif()

################################
# Adiak
################################
if(DEFINED ADIAK_DIR)
    message(STATUS "ADIAK_DIR = ${ADIAK_DIR}")

    find_package(adiak REQUIRED
                 PATHS ${ADIAK_DIR}/lib/cmake/adiak)

    blt_register_library(NAME adiak
                         INCLUDES ${adiak_INCLUDE_DIRS}
                         LIBRARIES ${adiak_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(ENABLE_ADIAK ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} adiak)
else()
    if( ENABLE_ADIAK )
        message(WARNING "ENABLE_ADIAK is ON but ADIAK_DIR isn't defined.")
    endif()

    set(ENABLE_ADIAK OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Adiak")
endif()

################################
# Caliper
################################
if(DEFINED CALIPER_DIR)
    message(STATUS "CALIPER_DIR = ${CALIPER_DIR}")

    find_package(caliper REQUIRED
                 PATHS ${CALIPER_DIR}/share/cmake/caliper)

    blt_register_library(NAME caliper
                         INCLUDES ${caliper_INCLUDE_PATH}
                         LIBRARIES ${caliper_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(ENABLE_CALIPER ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} caliper)
else()
    if( ENABLE_CALIPER )
        message(WARNING "ENABLE_CALIPER is ON but CALIPER_DIR isn't defined.")
    endif()

    set(ENABLE_CALIPER OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Caliper")
endif()

################################
# ASMJIT / MATHPRESSO
################################
# if (NOT EXISTS ${MATHPRESSO_DIR})
#     message(STATUS "Using mathpresso from thirdPartyLibs")
#     set(MATHPRESSO_DIR ${GEOSX_TPL_DIR}/mathpresso CACHE PATH "")
# endif()

# if (EXISTS ${MATHPRESSO_DIR})
#     message( STATUS "setting up asmjit" )
#     set(ASMJIT_LOCAL_DIR ${PROJECT_BINARY_DIR}/thirdparty/asmjit/src/asmjit)
#     set(ASMJIT_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/asmjit)

#     message( STATUS "setting up MathPresso" )
#     find_path( MATHPRESSO_INCLUDE_DIRS mathpresso/mathpresso.h
#             PATHS  ${MATHPRESSO_DIR}/include
#             NO_DEFAULT_PATH
#             NO_CMAKE_ENVIRONMENT_PATH
#             NO_CMAKE_PATH
#             NO_SYSTEM_ENVIRONMENT_PATH
#             NO_CMAKE_SYSTEM_PATH)

#     find_library( MATHPRESSO_LIBRARY NAMES mathpresso
#                 PATHS ${MATHPRESSO_DIR}/lib
#                 NO_DEFAULT_PATH
#                 NO_CMAKE_ENVIRONMENT_PATH
#                 NO_CMAKE_PATH
#                 NO_SYSTEM_ENVIRONMENT_PATH
#                 NO_CMAKE_SYSTEM_PATH)


#     include(FindPackageHandleStandardArgs)
#     find_package_handle_standard_args(MATHPRESSO  DEFAULT_MSG
#                                     MATHPRESSO_INCLUDE_DIRS
#                                     MATHPRESSO_LIBRARY )
                                    
#     if (NOT MATHPRESSO_FOUND)
#         message(FATAL_ERROR "MATHPRESSO not found in ${MATHPRESSO_DIR}. Maybe you need to build it")
#     endif()

#     blt_register_library( NAME mathpresso
#                         INCLUDES ${MATHPRESSO_INCLUDE_DIRS}
#                         LIBRARIES ${MATHPRESSO_LIBRARY}
#                         TREAT_INCLUDES_AS_SYSTEM ON )
    
#     set(MATHPRESSO_FOUND ON CACHE BOOL "")
#     set( thirdPartyLibs ${thirdPartyLibs} mathpresso )
# else()
    set(MATHPRESSO_FOUND OFF CACHE BOOL "")
    set(ENABLE_MATHPRESSO OFF CACHE BOOL "")
    message(STATUS "Not using mathpresso")
# endif()

################################
# METIS
################################
if(DEFINED METIS_DIR)
    message(STATUS "METIS_DIR = ${METIS_DIR}")

    find_path(METIS_INCLUDE_DIRS metis.h
              PATHS ${METIS_DIR}/include
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    find_library(METIS_LIBRARY NAMES metis
                 PATHS ${METIS_DIR}/lib
                 NO_DEFAULT_PATH
                 NO_CMAKE_ENVIRONMENT_PATH
                 NO_CMAKE_PATH
                 NO_SYSTEM_ENVIRONMENT_PATH
                 NO_CMAKE_SYSTEM_PATH)

    message(STATUS "METIS_INCLUDE_DIRS = ${METIS_INCLUDE_DIRS}")
    message(STATUS "METIS_LIBRARY = ${METIS_LIBRARY}")
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(METIS DEFAULT_MSG METIS_INCLUDE_DIRS METIS_LIBRARY )
    if (NOT METIS_FOUND)
        message(FATAL_ERROR "METIS not found in ${METIS_DIR}. Maybe you need to build it")
    endif()

    blt_register_library(NAME metis
                         INCLUDES ${METIS_INCLUDE_DIRS}
                         LIBRARIES ${METIS_LIBRARY}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(ENABLE_METIS ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} metis)
else()
    if( ENABLE_METIS )
        message(WARNING "ENABLE_METIS is ON but METIS_DIR isn't defined.")
    endif()

    set(ENABLE_METIS OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Metis")
endif()

################################
# PARMETIS
################################
if(DEFINED PARMETIS_DIR)
    message(STATUS "PARMETIS_DIR = ${PARMETIS_DIR}")

    find_path(PARMETIS_INCLUDE_DIRS parmetis.h
              PATHS ${PARMETIS_DIR}/include
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

   find_library(PARMETIS_LIBRARY NAMES parmetis
                PATHS ${PARMETIS_DIR}/lib
                NO_DEFAULT_PATH
                NO_CMAKE_ENVIRONMENT_PATH
                NO_CMAKE_PATH
                NO_SYSTEM_ENVIRONMENT_PATH
                NO_CMAKE_SYSTEM_PATH)

    message(STATUS "PARMETIS_INCLUDE_DIRS = ${PARMETIS_INCLUDE_DIRS}" )
    message(STATUS "PARMETIS_LIBRARY = ${PARMETIS_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(PARMETIS DEFAULT_MSG PARMETIS_INCLUDE_DIRS PARMETIS_LIBRARY )
    if (NOT PARMETIS_FOUND)
        message(FATAL_ERROR "PARMETIS not found in ${PARMETIS_DIR}. Maybe you need to build it")
    endif()

    blt_register_library(NAME parmetis
                         INCLUDES ${PARMETIS_INCLUDE_DIRS} 
                         LIBRARIES ${PARMETIS_LIBRARY}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(ENABLE_PARMETIS ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} parmetis)
else()
    if( ENABLE_PARMETIS )
        message(WARNING "ENABLE_PARMETIS is ON but PARMETIS_DIR isn't defined.")
    endif()

    set(ENABLE_PARMETIS OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Metis")
endif()

################################
# SUPERLU_DIST
################################
if(DEFINED SUPERLU_DIST_DIR)
    message(STATUS "SUPERLU_DIST_DIR = ${SUPERLU_DIST_DIR}")

    find_path(SUPERLU_DIST_INCLUDE_DIRS superlu_defs.h
              PATHS ${SUPERLU_DIST_DIR}/include
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

   find_library(SUPERLU_DIST_LIBRARY NAMES superlu_dist
                PATHS ${SUPERLU_DIST_DIR}/lib PATHS ${SUPERLU_DIST_DIR}/lib64
                NO_DEFAULT_PATH
                NO_CMAKE_ENVIRONMENT_PATH
                NO_CMAKE_PATH
                NO_SYSTEM_ENVIRONMENT_PATH
                NO_CMAKE_SYSTEM_PATH)

    message(STATUS "SUPERLU_DIST_INCLUDE_DIRS = ${SUPERLU_DIST_INCLUDE_DIRS}" )
    message(STATUS "SUPERLU_DIST_LIBRARY = ${SUPERLU_DIST_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(SUPERLU_DIST DEFAULT_MSG SUPERLU_DIST_INCLUDE_DIRS SUPERLU_DIST_LIBRARY)
    if(NOT SUPERLU_DIST_FOUND)
        message(FATAL_ERROR "SUPERLU_DIST not found in ${SUPERLU_DIST_DIR}. Maybe you need to build it")
    endif()

    blt_register_library(NAME superlu_dist
                         DEPENDS_ON parmetis metis lapack blas
                         INCLUDES ${SUPERLU_DIST_INCLUDE_DIRS} 
                         LIBRARIES ${SUPERLU_DIST_LIBRARY}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(ENABLE_SUPERLU_DIST ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} superlu_dist)
else()
    if( ENABLE_SUPERLU_DIST )
        message(WARNING "ENABLE_SUPERLU_DIST is ON but SUPERLU_DIST_DIR isn't defined.")
    endif()

    set(ENABLE_SUPERLU_DIST OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using superlu_dist")
endif()

################################
# SUITESPARSE
################################
# if( ENABLE_SUITESPARSE )
#     message( STATUS "setting up SUITESPARSE" )

#     set(SUITESPARSE_DIR ${GEOSX_TPL_DIR}/suitesparse)

#     find_path( SUITESPARSE_INCLUDE_DIRS umfpack.h
#                PATHS  ${SUITESPARSE_DIR}/include
#                NO_DEFAULT_PATH
#                NO_CMAKE_ENVIRONMENT_PATH
#                NO_CMAKE_PATH
#                NO_SYSTEM_ENVIRONMENT_PATH
#                NO_CMAKE_SYSTEM_PATH)

#     find_library( SUITESPARSE_LIBRARY NAMES umfpack
#                   PATHS ${SUITESPARSE_DIR}/lib
#                   NO_DEFAULT_PATH
#                   NO_CMAKE_ENVIRONMENT_PATH
#                   NO_CMAKE_PATH
#                   NO_SYSTEM_ENVIRONMENT_PATH
#                   NO_CMAKE_SYSTEM_PATH)

#     message(STATUS "SUITESPARSE_INCLUDE_DIRS = ${SUITESPARSE_INCLUDE_DIRS}" )
#     message(STATUS "SUITESPARSE_LIBRARY = ${SUITESPARSE_LIBRARY}" )
#     include(FindPackageHandleStandardArgs)
#     find_package_handle_standard_args(SUITESPARSE DEFAULT_MSG
#                                       SUITESPARSE_INCLUDE_DIRS
#                                       SUITESPARSE_LIBRARY )
#     if (NOT SUITESPARSE_FOUND)
#         message(FATAL_ERROR "SUITESPARSE not found in ${SUITESPARSE_DIR}. Maybe you need to build it")
#     endif()
    
#     set( SUITESPARSE_DEPENDS "blas;lapack" )

#     blt_register_library( NAME suitesparse
#                           DEPENDS_ON ${SUITESPARSE_DEPENDS}
#                           INCLUDES ${SUITESPARSE_INCLUDE_DIRS}
#                           LIBRARIES ${SUITESPARSE_LIBRARY}
#                           TREAT_INCLUDES_AS_SYSTEM ON )

#     set( thirdPartyLibs ${thirdPartyLibs} suitesparse )
# endif()

################################
# HYPRE
################################
# if( ENABLE_HYPRE )
#     message( STATUS "setting up HYPRE" )

#     set(HYPRE_DIR ${GEOSX_TPL_DIR}/hypre)

#     find_path( HYPRE_INCLUDE_DIRS HYPRE.h
#             PATHS  ${HYPRE_DIR}/include
#             NO_DEFAULT_PATH
#             NO_CMAKE_ENVIRONMENT_PATH
#             NO_CMAKE_PATH
#             NO_SYSTEM_ENVIRONMENT_PATH
#             NO_CMAKE_SYSTEM_PATH)

#     find_library( HYPRE_LIBRARY NAMES HYPRE
#                 PATHS ${HYPRE_DIR}/lib
#                 NO_DEFAULT_PATH
#                 NO_CMAKE_ENVIRONMENT_PATH
#                 NO_CMAKE_PATH
#                 NO_SYSTEM_ENVIRONMENT_PATH
#                 NO_CMAKE_SYSTEM_PATH)

#     message(STATUS "HYPRE_INCLUDE_DIRS = ${HYPRE_INCLUDE_DIRS}" )
#     message(STATUS "HYPRE_LIBRARY = ${HYPRE_LIBRARY}" )
#     include(FindPackageHandleStandardArgs)
#     find_package_handle_standard_args(HYPRE  DEFAULT_MSG
#                                     HYPRE_INCLUDE_DIRS
#                                     HYPRE_LIBRARY )
#     if (NOT HYPRE_FOUND)
#         message(FATAL_ERROR "HYPRE not found in ${HYPRE_DIR}. Maybe you need to build it")
#     endif()
    
#     set( HYPRE_DEPENDS "blas;lapack" )
#     if( ENABLE_SUPERLU_DIST )
#         list( APPEND HYPRE_DEPENDS "superlu_dist" )
#     endif()

#     blt_register_library( NAME hypre
#                           DEPENDS_ON ${HYPRE_DEPENDS}
#                           INCLUDES ${HYPRE_INCLUDE_DIRS}
#                           LIBRARIES ${HYPRE_LIBRARY}
#                           TREAT_INCLUDES_AS_SYSTEM ON )

#     set( thirdPartyLibs ${thirdPartyLibs} hypre )
# endif()

################################
# TRILINOS
################################
if(DEFINED TRILINOS_DIR)
    message(STATUS "TRILINOS_DIR = ${TRILINOS_DIR}")
  
    include(${TRILINOS_DIR}/lib/cmake/Trilinos/TrilinosConfig.cmake)
  
    list(REMOVE_ITEM Trilinos_LIBRARIES "gtest")
    list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
    message(STATUS "Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
    message(STATUS "Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")

    blt_register_library(NAME trilinos
                         DEPENDS_ON ${TRILINOS_DEPENDS}
                         INCLUDES ${Trilinos_INCLUDE_DIRS} 
                         LIBRARIES ${Trilinos_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON )

    set(ENABLE_TRILINOS ON CACHE BOOL "" FORCE)
    set(thirdPartyLibs ${thirdPartyLibs} trilinos)
else()
    if( ENABLE_TRILINOS )
        message(WARNING "ENABLE_TRILINOS is ON but TRILINOS_DIR isn't defined.")
    endif()

    set(ENABLE_TRILINOS OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Trilinos")
endif()

################################
# PETSC
################################
# if( ENABLE_PETSC )
#     message( STATUS "setting up PETSC" )

#     if( EXISTS ${PETSC_DIR} )
  
#     else()
#         set(PETSC_DIR ${GEOSX_TPL_DIR}/petsc)
#     endif()

#     find_path( Petsc_INCLUDE_DIRS petscvec.h
#                PATHS  ${PETSC_DIR}/include
#                NO_DEFAULT_PATH
#                NO_CMAKE_ENVIRONMENT_PATH
#                NO_CMAKE_PATH
#                NO_SYSTEM_ENVIRONMENT_PATH
#                NO_CMAKE_SYSTEM_PATH)

#     find_library( Petsc_LIBRARIES NAMES petsc
#                   PATHS ${PETSC_DIR}/lib
#                   NO_DEFAULT_PATH
#                   NO_CMAKE_ENVIRONMENT_PATH
#                   NO_CMAKE_PATH
#                   NO_SYSTEM_ENVIRONMENT_PATH
#                   NO_CMAKE_SYSTEM_PATH)

#     message( STATUS "Petsc_INCLUDE_DIRS = ${Petsc_INCLUDE_DIRS}" )
#     message( STATUS "Petsc_LIBRARIES = ${Petsc_LIBRARIES}" )
  
  
#     blt_register_library( NAME petsc
#                           INCLUDES ${Petsc_INCLUDE_DIRS} 
#                           LIBRARIES ${Petsc_LIBRARIES}
#                           TREAT_INCLUDES_AS_SYSTEM ON )
#     set( thirdPartyLibs ${thirdPartyLibs} petsc )  

# endif()


################################
# VTK
################################
# if( ENABLE_VTK )
#   find_package(VTK REQUIRED PATHS ${GEOSX_TPL_DIR}/vtk NO_DEFAULT_PATH
#   )

#   message( STATUS "VTK_INCLUDE_DIRS = ${VTK_INCLUDE_DIRS}" )
#   message( STATUS "VTK_LIBRARIES = ${VTK_LIBRARIES}" )

#   blt_register_library( NAME vtk
#                         LIBRARIES ${VTK_LIBRARIES}
#                         TREAT_INCLUDES_AS_SYSTEM ON )
#   set( thirdPartyLibs ${thirdPartyLibs} vtk )  
# endif()

################################
# uncrustify (This is here instead of SetupGeosxThirdParty.cmake so that BLT will see it)
################################
if(DEFINED UNCRUSTIFY_EXECUTABLE)
    message( STATUS "Using system uncrustify found at ${UNCRUSTIFY_EXECUTABLE}" )
else()
    set( UNCRUSTIFY_EXECUTABLE "${GEOSX_TPL_DIR}/uncrustify/bin/uncrustify" CACHE PATH "" )
    message( STATUS "Using uncrustify from thirdPartyLibs at ${UNCRUSTIFY_EXECUTABLE}" )
endif()

################################
# unsrustify
################################
if(UNCRUSTIFY_FOUND)
    message(STATUS "UNCRUSTIFY_EXECUTABLE = ${UNCRUSTIFY_EXECUTABLE}")

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
else()
    message(STATUS "Not using uncrustify.")
endif()

string( TOUPPER "${GEOSX_LA_INTERFACE}" upper_LAI )
if( NOT ENABLE_${upper_LAI} )
  message( FATAL_ERROR "${GEOSX_LA_INTERFACE} LA interface is selected, but ENABLE_${upper_LAI} is OFF" )
endif()
option( GEOSX_LA_INTERFACE_${upper_LAI} "${upper_LAI} LA interface is selected" ON )

message(STATUS "thirdPartyLibs = ${thirdPartyLibs}")
