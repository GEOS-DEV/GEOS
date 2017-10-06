cmake_minimum_required(VERSION 3.1)

################################
# GEOS_TPL
################################
project(geosx_tpl LANGUAGES C CXX)

################################
# Include blt
################################
include(cmake/blt/SetupBLT.cmake)

include(ExternalProject)


################################
# install targets from build
################################s
#install(DIRECTORY ${PROJECT_BINARY_DIR}/bin DESTINATION . OPTIONAL)
#install(DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION . OPTIONAL)



set( dependency_list "" )
################################
# HDF5
################################
if (HDF5_DIR)
  include(../src/cmake/thirdparty/FindHDF5.cmake)
  blt_register_library(NAME hdf5
                       INCLUDES ${HDF5_INCLUDE_DIRS}
                       LIBRARIES ${HDF5_LIBRARY} 
                       TREAT_INCLUDES_AS_SYSTEM ON )
endif()


################################
# SILO
################################
if( EXISTS ${SILO_DIR})
    message("Using system SILO found at ${SILO_DIR}")
    include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindSILO.cmake)
    if (NOT SILO_FOUND)
        message(FATAL_ERROR ": SILO not found in ${SILO_DIR}. Maybe you need to build it")
    endif()
else()
    message(INFO ": Using SILO found at https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz")
    message("${CMAKE_SOURCE_DIR}")
    set(silo_install_dir ${CMAKE_INSTALL_PREFIX}/silo)
    ExternalProject_Add( silo
#                         URL https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz
#                         URL_HASH MD5=60fef9ce373daf1e9cc8320cfa509bc5
#                         DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
                          URL ${CMAKE_SOURCE_DIR}/tpl_mirror/silo-4.10.2-bsd.tar.gz
                          PREFIX ${PROJECT_BINARY_DIR}/silo
                          INSTALL_DIR ${silo_install_dir}
                          CONFIGURE_COMMAND ../silo/configure CC=${CMAKE_C_COMPILER}
                                               CXX=${CMAKE_CXX_COMPILER}
                                               --prefix=${silo_install_dir}
                                               --disable-fortran
                                               --enable-optimization
                                               --with-hdf5=${HDF5_DIR}/include,${HDF5_DIR}/lib
                                               LDFLAGS=-ldl
                                               --disable-silex 
                          BUILD_COMMAND make -j8
                          INSTALL_COMMAND make install )

    set( dependency_list "${dependency_list};silo" )
endif()
        





################################
# RAJA
################################
if( EXISTS ${RAJA_DIR})
message("${RAJA_DIR}")
    message("Using system RAJA found at ${RAJA_DIR}")
    include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindRAJA.cmake)
    if (NOT RAJA_FOUND)
        message(FATAL_ERROR ": RAJA not found in ${RAJA_DIR}. Maybe you need to build it")
    endif()    
else()
    message(INFO ": Using RAJA found at https://github.com/LLNL/RAJA/archive/develop.zip")
    set(raja_install_dir ${CMAKE_INSTALL_PREFIX}/raja)
    message("RAJA_CMAKE_ARGS=${RAJA_CMAKE_ARGS}")
    ExternalProject_Add( raja
                         GIT_REPOSITORY https://github.com/LLNL/RAJA.git
                         GIT_TAG develop
                         PATCH_COMMAND git submodule init && git submodule update                         
                         PREFIX ${PROJECT_BINARY_DIR}/raja
                         BUILD_COMMAND make -j8
                         INSTALL_COMMAND make install
                         INSTALL_DIR ${raja_install_dir}
                         CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                    -DRAJA_ENABLE_CUDA=${CUDA_ENABLED}
                                    -DRAJA_ENABLE_TESTS=${RAJA_ENABLE_TESTS}
                                    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                                    -DRAJA_ENABLE_TBB=${RAJA_ENABLE_TBB} )

    set( dependency_list "${dependency_list};raja" )
endif()




################################
# CHAI
################################
if( EXISTS ${CHAI_DIR})
    message(INFO ": Using system CHAI found at ${CHAI_DIR}")
    include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindCHAI.cmake)
    if (NOT CHAI_FOUND)
        message(FATAL_ERROR ": CHAI not found in ${CHAI_DIR}. Maybe you need to build it")
    endif()
else()
    set(CHAI_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/chai)
    message(INFO ": Using CHAI found at git@github.com:LLNL/CHAI.git")
    ExternalProject_Add( chai
                         PREFIX ${PROJECT_BINARY_DIR}/chai
                         GIT_REPOSITORY git@github.com:LLNL/CHAI.git
                         GIT_TAG develop
                         PATCH_COMMAND git submodule init && git submodule update                         
                         BUILD_COMMAND make -j8
                         INSTALL_DIR ${chai_install_dir}
                         INSTALL_COMMAND make install
                         CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                    -DENABLE_CUDA=${CUDA_ENABLED}
                                    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

    set( dependency_list "${dependency_list};chai" )

endif()


################################
# FPARSER
################################
if( USE_FPARSER )
message( INFO ": setting up fparser" )
set(FPARSER_LOCAL_DIR ${CMAKE_SOURCE_DIR}/fparser)
set(FPARSER_DIR ${FPARSER_LOCAL_DIR})
set(FPARSER_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/fparser)
set(FPARSER_INCLUDE_DIR ${fparser_install_dir}/include)

message( INFO ": FPARSER_DIR = ${FPARSER_DIR}" )
message( INFO ": FPARSER_LOCAL_DIR = ${FPARSER_LOCAL_DIR}" )
message( INFO ": FPARSER_INSTALL_DIR = ${FPARSER_INSTALL_DIR}" )

ExternalProject_Add( fparser 
                     #URL http://warp.povusers.org/FunctionParser/fparser4.5.2.zip
                     #URL_HASH MD5=60fef9ce373daf1e9cc8320cfa509bc5
                     #DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/fparser4.5.2.zip
                     PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fparser
                     INSTALL_DIR ${FPARSER_INSTALL_DIR}
                     CONFIGURE_COMMAND ""
                     BUILD_COMMAND ${CMAKE_CXX_COMPILER} -c -DFP_NO_SUPPORT_OPTIMIZER -I../fparser ../fparser/fparser.cc ../fparser/fpoptimizer.cc &&
                                   ar rcs libfparser.a fparser.o fpoptimizer.o
                     INSTALL_COMMAND mkdir -p ${FPARSER_INSTALL_DIR}/lib &&
                                     cp libfparser.a ${FPARSER_INSTALL_DIR}/lib &&
                                     cd ../fparser &&
                                     mkdir -p ${FPARSER_INSTALL_DIR}/include && 
                                     ls  &&
                                     cp fparser.hh fparser_gmpint.hh fparser_mpfr.hh fpconfig.hh ${FPARSER_INSTALL_DIR}/include;
                     )

    set( dependency_list "${dependency_list};fparser" )

endif()






################################
# CALIPER
################################
message( "ENABLE_CALIPER = ${ENABLE_CALIPER}" )
if( ENABLE_CALIPER )
message( INFO ": setting up caliper" )
set(CALIPER_LOCAL_DIR ${PROJECT_BINARY_DIR}/caliper)
set(CALIPER_DIR ${CALIPER_LOCAL_DIR})
set(CALIPER_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/caliper)

message( INFO ": CALIPER_DIR = ${CALIPER_DIR}" )
message( INFO ": CALIPER_LOCAL_DIR = ${CALIPER_LOCAL_DIR}" )
message( INFO ": CALIPER_INSTALL_DIR = ${CALIPER_INSTALL_DIR}" )

ExternalProject_Add( caliper
                     GIT_REPOSITORY git@github.com:LLNL/Caliper.git
                     GIT_TAG master
                     PATCH_COMMAND git submodule init && git submodule update                         
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/Caliper-master.zip
                     BUILD_COMMAND make -j8
                     INSTALL_DIR ${CALIPER_INSTALL_DIR}
                     INSTALL_COMMAND make install
                     CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                -DBUILD_SHARED_LIBS:BOOL=OFF
                                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

    set( dependency_list "${dependency_list};caliper" )

endif()



################################
# ASMJIT / MATHPRESSO
################################
set(ENABLE_MATHPRESSO ON CACHE BOOL  "Enables mathpresso Plugin")
if( ENABLE_MATHPRESSO )
message( INFO ": setting up asmjit" )
set(ASMJIT_LOCAL_DIR ${PROJECT_BINARY_DIR}/asmjit/src/asmjit)
set(ASMJIT_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/asmjit)

ExternalProject_Add( asmjit
                     PREFIX ${PROJECT_BINARY_DIR}/asmjit
#                     URL https://github.com/asmjit/asmjit/archive/master.zip
#                     URL_HASH MD5=3c0b3190d422240b075dfc667a081a3a
#                     DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
#                     DOWNLOAD_NAME asmjit.zip
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/asmjit-master.zip
                     BUILD_COMMAND make -j8
                     INSTALL_DIR ${ASMJIT_INSTALL_DIR}
                     INSTALL_COMMAND make install
                     CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                -DBUILD_SHARED_LIBS:BOOL=OFF
                                -DCMAKE_BUILD_TYPE=Release
                                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

message( INFO ": setting up MathPresso" )
set(MATHPRESSO_LOCAL_DIR ${PROJECT_BINARY_DIR}/mathpresso)
set(MATHPRESSO_DIR ${MATHPRESSO_LOCAL_DIR})
set(MATHPRESSO_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/mathpresso)

message( INFO ": MATHPRESSO_DIR = ${MATHPRESSO_DIR}" )
message( INFO ": MATHPRESSO_LOCAL_DIR = ${MATHPRESSO_LOCAL_DIR}" )
message( INFO ": MATHPRESSO_INSTALL_DIR = ${MATHPRESSO_INSTALL_DIR}" )

ExternalProject_Add( mathpresso
                     PREFIX ${PROJECT_BINARY_DIR}/mathpresso
#                     URL https://github.com/kobalicek/mathpresso/archive/master.zip
#                     URL_HASH MD5=b43212cafeab5e0e2ef5b87c29c15df1
#                     DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
#                     DOWNLOAD_NAME mathpresso.zip
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/mathpresso-master.zip
                     DEPENDS asmjit 
                     BUILD_COMMAND make -j8
                     INSTALL_DIR ${MATHPRESSO_INSTALL_DIR}
                     INSTALL_COMMAND mkdir -p <INSTALL_DIR>/include &&
                                     mkdir -p <INSTALL_DIR>/lib &&
                                     make INSTALL_DIR=<INSTALL_DIR> install &&
                                     cp libmathpresso.a <INSTALL_DIR>/lib/
                     CMAKE_ARGS -DMATHPRESSO_STATIC=TRUE
                                -DASMJIT_DIR=${ASMJIT_LOCAL_DIR}
                                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

set( dependency_list "${dependency_list};mathpresso" )

endif()



################################
# PUGIXML
################################
message( INFO ": setting up pugixml" )
set(PUGIXML_LOCAL_DIR ${PROJECT_BINARY_DIR}/pugixml)
set(PUGIXML_DIR ${PUGIXML_LOCAL_DIR})
set(PUGIXML_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/pugixml)

message( INFO ": PUGIXML_DIR = ${PUGIXML_DIR}" )
message( INFO ": PUGIXML_LOCAL_DIR = ${PUGIXML_LOCAL_DIR}" )
message( INFO ": PUGIXML_INSTALL_DIR = ${PUGIXML_INSTALL_DIR}" )

ExternalProject_Add( pugixml
                     PREFIX ${PROJECT_BINARY_DIR}/pugixml
#                     URL https://github.com/zeux/pugixml/archive/master.zip
#                     DOWNLOAD_NAME pugixml.zip
#                     URL_HASH MD5=0099563cb3f466fe03b10f9666c73993
#                     DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/pugixml.tar
                     BUILD_COMMAND make -j8
                     INSTALL_DIR ${PUGIXML_INSTALL_DIR}
                     INSTALL_COMMAND make install
                     CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

set( dependency_list "${dependency_list};pugixml" )




################################
# TRILINOS
################################
message( INFO ": setting up TRILINOS" )
set(TRILINOS_LOCAL_DIR ${PROJECT_BINARY_DIR}/trilinos)
set(TRILINOS_DIR ${TRILINOS_LOCAL_DIR})
set(TRILINOS_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/trilinos)

ExternalProject_Add( trilinos
                     PREFIX ${PROJECT_BINARY_DIR}/trilinos
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/trilinos-12.10.1-Source.tar.gz
                     BUILD_COMMAND make -j8
                     INSTALL_DIR ${TRILINOS_INSTALL_DIR}
                     INSTALL_COMMAND make install
                     CMAKE_ARGS -D CMAKE_C_COMPILER=${MPI_C_COMPILER}
                                -D CMAKE_CXX_COMPILER=${MPI_CXX_COMPILER}
                                -D BUILD_SHARED_LIBS:BOOL=OFF
                                -D CMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
                                -D TPL_ENABLE_MPI:BOOL=ON 
                                -D BUILD_SHARED_LIBS:BOOL=OFF
                                -D CMAKE_BUILD_TYPE:STRING=RELEASE
                                -D Trilinos_ENABLE_Fortran:BOOL=OFF 
                                -D Trilinos_WARNINGS_AS_ERRORS_FLAGS:STRING="" 
                                -D CMAKE_VERBOSE_MAKEFILE:BOOL=TRUE 
                                -D Trilinos_ENABLE_TESTS:BOOL=OFF 
                                -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=FALSE 
                                -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON 
                                -D Trilinos_ENABLE_Epetra:BOOL=ON 
                                -D Trilinos_ENABLE_EpetraExt:BOOL=ON 
                                -D Trilinos_ENABLE_Tpetra:BOOL=ON 
                                -D Trilinos_ENABLE_Jpetra:BOOL=ON 
                                -D Trilinos_ENABLE_Kokkos:BOOL=ON 
                                -D Trilinos_ENABLE_Mesquite:BOOL=ON 
                                -D Trilinos_ENABLE_Sacado:BOOL=ON 
                                -D Trilinos_ENABLE_Stratimikos:BOOL=ON 
                                -D Trilinos_ENABLE_Amesos:BOOL=ON 
                                -D Trilinos_ENABLE_AztecOO:BOOL=ON 
                                -D Trilinos_ENABLE_Ifpack:BOOL=ON 
                                -D Trilinos_ENABLE_Teuchos:BOOL=ON 
                                -D Trilinos_ENABLE_ML:BOOL=ON 
                                -D Trilinos_ENABLE_Intrepid:BOOL=ON 
                                -D Trilinos_ENABLE_Shards:BOOL=ON 
                                -D Trilinos_ENABLE_Pamgen:BOOL=ON 
                                -D Trilinos_ENABLE_Thyra:BOOL=ON 
                                -D Trilinos_ENABLE_Boost=OFF
                                -D Trilinos_ENABLE_STK=OFF
                                -D Trilinos_ENABLE_TrilinosCouplings:BOOL=ON  )

set( dependency_list "${dependency_list};trilinos" )



blt_add_executable( NAME             tpl
                    SOURCES          tpl.cpp )

add_dependencies( tpl ${dependency_list} )

