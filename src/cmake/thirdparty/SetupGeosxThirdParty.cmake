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
                       LIBRARIES ${HDF5_LIBRARY} 
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



################################
# SILO
################################
if( EXISTS ${SILO_DIR})
    if( NOT BUILD_LOCAL_SILO )
        message("Using system SILO found at ${SILO_DIR}")
        include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindSILO.cmake)
        if (NOT SILO_FOUND)
            message(FATAL_ERROR ": SILO not found in ${SILO_DIR}. Maybe you need to build it")
        endif()
    
        blt_register_library( NAME silo
                              INCLUDES ${SILO_INCLUDE_DIRS}
                              LIBRARIES ${SILO_LIBRARY}
                              TREAT_INCLUDES_AS_SYSTEM ON
                              DEPENDS_ON hdf5 )
    else()
        message(INFO ": Build SILO from source found at ${SILO_DIR}")
        set(silo_install_dir ${CMAKE_INSTALL_PREFIX}/thirdparty/silo)
        ExternalProject_Add( SILO
                             PREFIX ${PROJECT_BINARY_DIR}/thirdparty/silo
                             SOURCE_DIR ${SILO_DIR}
                             BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/silo/src
                             DEPENDS hdf5
                             INSTALL_COMMAND make install
                             INSTALL_DIR ${silo_install_dir}
		                     CONFIGURE_COMMAND ./configure CC=${CMAKE_C_COMPILER}
                   			 				   CXX=${CMAKE_CXX_COMPILER}
			                                   --prefix=$(CWD) 
			                                   --disable-fortran
                                               --with-hdf5=${HDF5_DIR}/include,${HDF5_DIR}/lib
                                               LDFLAGS=-ldl
			                                   --disable-silex 
        		             BUILD_COMMAND make
                      		 INSTALL_COMMAND make install )
                                 
        blt_register_library( NAME SILO
                              INCLUDES ${silo_install_dir}/include 
                              LIBRARIES ${silo_install_dir}/lib/libsilo.a
                              TREAT_INCLUDES_AS_SYSTEM ON
                              DEPENDS_ON hdf5 )        
    endif()
else()
    message(INFO ": Using SILO found at https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz")
    message("${CMAKE_SOURCE_DIR}")
    set(silo_install_dir ${CMAKE_INSTALL_PREFIX}/thirdparty/silo)
    ExternalProject_Add( silo
#                         URL https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2-bsd.tar.gz
#                         URL_HASH MD5=60fef9ce373daf1e9cc8320cfa509bc5
#                         DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
                         URL ${CMAKE_SOURCE_DIR}/tpl_mirror/silo-4.10.2-bsd.tar.gz
                         PREFIX ${PROJECT_BINARY_DIR}/thirdparty/silo
                         INSTALL_DIR ${silo_install_dir}
                         CONFIGURE_COMMAND ../silo/configure CC=${CMAKE_C_COMPILER}
                                               CXX=${CMAKE_CXX_COMPILER}
                                               --prefix=${silo_install_dir}
                                               --disable-fortran
                                               --enable-optimization
                                               --with-hdf5=${HDF5_DIR}/include,${HDF5_DIR}/lib
                                               LDFLAGS=-ldl
                                               --disable-silex 
                             BUILD_COMMAND make
                             INSTALL_COMMAND make install )

    blt_register_library( NAME silo
                          INCLUDES ${silo_install_dir}/include 
                          LIBRARIES ${silo_install_dir}/lib/libsiloh5.a
                          TREAT_INCLUDES_AS_SYSTEM ON
                          DEPENDS_ON hdf5  )

endif()






################################
# RAJA
################################
if( EXISTS ${RAJA_DIR})
message("${RAJA_DIR}")
    if( NOT BUILD_LOCAL_RAJA )
        message("Using system RAJA found at ${RAJA_DIR}")
        include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindRAJA.cmake)
        if (NOT RAJA_FOUND)
            message(FATAL_ERROR ": RAJA not found in ${RAJA_DIR}. Maybe you need to build it")
        endif()
    
        blt_register_library( NAME raja
                              INCLUDES ${RAJA_INCLUDE_DIRS}
                              LIBRARIES ${RAJA_LIBRARY}
                              TREAT_INCLUDES_AS_SYSTEM ON )
    else()
        message(INFO ": Build RAJA from source found at ${RAJA_DIR}")
        set(raja_install_dir ${CMAKE_INSTALL_PREFIX}/thirdparty/raja)
        ExternalProject_Add( RAJA
                             PREFIX ${PROJECT_BINARY_DIR}/thirdparty/raja
                             SOURCE_DIR ${RAJA_DIR}
                             BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/raja/src
                             INSTALL_COMMAND make install
                             INSTALL_DIR ${raja_install_dir}
                             CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                        -DRAJA_ENABLE_CUDA=${CUDA_ENABLED}
                                        -DRAJA_ENABLE_TESTS=${RAJA_ENABLE_TESTS}
                                        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>  )
                                 
        blt_register_library( NAME RAJA
                              INCLUDES ${raja_install_dir}/include 
                              LIBRARIES ${raja_install_dir}/lib/libRAJA.a
                              TREAT_INCLUDES_AS_SYSTEM ON )        
    endif()
else()
    message(INFO ": Using RAJA found at https://github.com/LLNL/RAJA/archive/develop.zip")
    set(raja_install_dir ${CMAKE_INSTALL_PREFIX}/thirdparty/raja)
    ExternalProject_Add( raja
#                         URL https://github.com/LLNL/RAJA/archive/develop.zip
#                         URL_HASH MD5=d60ad60c1c9d662893862e3b35d8b2dc
#                         DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
#                         DOWNLOAD_NAME raja.zip
                         URL ${CMAKE_SOURCE_DIR}/tpl_mirror/RAJA-develop.zip
                         PREFIX ${PROJECT_BINARY_DIR}/thirdparty/raja
                         INSTALL_COMMAND make install
                         INSTALL_DIR ${raja_install_dir}
                         CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                    -DRAJA_ENABLE_CUDA=${CUDA_ENABLED}
                                    -DRAJA_ENABLE_TESTS=${RAJA_ENABLE_TESTS}
                                    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

    blt_register_library( NAME raja
                          INCLUDES ${raja_install_dir}/include 
                          LIBRARIES ${raja_install_dir}/lib/libRAJA.a
                          TREAT_INCLUDES_AS_SYSTEM ON )

endif()




################################
# CHAI
################################
if( EXISTS ${CHAI_DIR})
    if( NOT BUILD_LOCAL_CHAI )
        message(INFO ": Using system CHAI found at ${CHAI_DIR}")
        include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindCHAI.cmake)
        if (NOT CHAI_FOUND)
            message(FATAL_ERROR ": CHAI not found in ${CHAI_DIR}. Maybe you need to build it")
        endif()
    
        blt_register_library( NAME chai
                              INCLUDES ${CHAI_INCLUDE_DIRS}
                              LIBRARIES ${CHAI_LIBRARY}
                              TREAT_INCLUDES_AS_SYSTEM ON  )
    else()
        message(INFO ": Build CHAI from source found at ${CHAI_DIR}")
        set(chai_install_dir ${CMAKE_INSTALL_PREFIX}/thirdparty/chai)
        ExternalProject_Add( chai
                             PREFIX ${PROJECT_BINARY_DIR}/thirdparty/chai
                             SOURCE_DIR ${CHAI_DIR}
                             BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/chai/src
                             INSTALL_COMMAND make install
                             INSTALL_DIR ${chai_install_dir}
                             CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                        -DENABLE_CUDA=${CUDA_ENABLED}
                                        -DENABLE_OPENMP=${ENABLE_OPENMP}
                                        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )
                                 
        blt_register_library( NAME chai
                              INCLUDES ${chai_install_dir}/include 
                              LIBRARIES ${chai_install_dir}/lib/libchai.a
                              TREAT_INCLUDES_AS_SYSTEM ON )        
    endif()
else()
    set(CHAI_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/chai)
    message( INFO ": CHAI_INSTALL_DIR = ${CHAI_INSTALL_DIR}" )

    message(INFO ": Using CHAI found at ssh://git@cz-bitbucket.llnl.gov:7999/um/chai.git")
    ExternalProject_Add( chai
                         PREFIX ${PROJECT_BINARY_DIR}/thirdparty/chai
                         GIT_REPOSITORY ssh://git@cz-bitbucket.llnl.gov:7999/um/chai.git
                         GIT_TAG develop
                         PATCH_COMMAND git submodule init && git submodule update                         
                         INSTALL_DIR ${chai_install_dir}
                         INSTALL_COMMAND make install
                         CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                    -DENABLE_CUDA=${CUDA_ENABLED}
                                    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

    blt_register_library( NAME chai
                          INCLUDES ${chai_install_dir}/include 
                          LIBRARIES ${chai_install_dir}/lib/libchai.a 
                          TREAT_INCLUDES_AS_SYSTEM ON )
endif()



################################
# FPARSER
################################
if( USE_FPARSER )
message( INFO ": setting up fparser" )
set(FPARSER_LOCAL_DIR ${CMAKE_SOURCE_DIR}/thirdparty/fparser)
set(FPARSER_DIR ${FPARSER_LOCAL_DIR})
set(FPARSER_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/fparser)
set(FPARSER_INCLUDE_DIR ${fparser_install_dir}/include)

message( INFO ": FPARSER_DIR = ${FPARSER_DIR}" )
message( INFO ": FPARSER_LOCAL_DIR = ${FPARSER_LOCAL_DIR}" )
message( INFO ": FPARSER_INSTALL_DIR = ${FPARSER_INSTALL_DIR}" )

ExternalProject_Add( fparser 
                     #URL http://warp.povusers.org/FunctionParser/fparser4.5.2.zip
                     #URL_HASH MD5=60fef9ce373daf1e9cc8320cfa509bc5
                     #DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/fparser4.5.2.zip
                     PREFIX ${CMAKE_CURRENT_BINARY_DIR}/thirdparty/fparser
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

blt_register_library( NAME fparser
                      INCLUDES ${FPARSER_INSTALL_DIR}/include 
                      LIBRARIES ${FPARSER_INSTALL_DIR}/lib/libfparser.a
                      TREAT_INCLUDES_AS_SYSTEM ON )

endif()






################################
# CALIPER
################################
if( ENABLE_CALIPER )
message( INFO ": setting up caliper" )
set(CALIPER_LOCAL_DIR ${PROJECT_BINARY_DIR}/thirdparty/caliper)
set(CALIPER_DIR ${CALIPER_LOCAL_DIR})
set(CALIPER_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/caliper)

message( INFO ": CALIPER_DIR = ${CALIPER_DIR}" )
message( INFO ": CALIPER_LOCAL_DIR = ${CALIPER_LOCAL_DIR}" )
message( INFO ": CALIPER_INSTALL_DIR = ${CALIPER_INSTALL_DIR}" )

#https://github.com/LLNL/Caliper.git
#git@github.com:LLNL/Caliper.git

ExternalProject_Add( caliper
#                     PREFIX ${PROJECT_BINARY_DIR}/thirdparty/caliper
#                     URL https://github.com/LLNL/Caliper/archive/master.zip
#                     URL_HASH MD5=60fef9ce373daf1e9cc8320cfa509bc5
#                     DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
#                     DOWNLOAD_NAME caliper.zip
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/Caliper-master.zip
                     INSTALL_DIR ${CALIPER_INSTALL_DIR}
                     INSTALL_COMMAND make install
                     CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                -DBUILD_SHARED_LIBS:BOOL=OFF
                                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

#include(${CALIPER_DIR}/share/cmake/caliper/caliper-config.cmake OPTIONAL)

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

ExternalProject_Add( asmjit
                     PREFIX ${PROJECT_BINARY_DIR}/thirdparty/asmjit
#                     URL https://github.com/asmjit/asmjit/archive/master.zip
#                     URL_HASH MD5=3c0b3190d422240b075dfc667a081a3a
#                     DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
#                     DOWNLOAD_NAME asmjit.zip
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/asmjit-master.zip
                     INSTALL_DIR ${ASMJIT_INSTALL_DIR}
                     INSTALL_COMMAND make install
                     CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                -DBUILD_SHARED_LIBS:BOOL=OFF
                                -DCMAKE_BUILD_TYPE=Release
                                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

message( INFO ": setting up MathPresso" )
set(MATHPRESSO_LOCAL_DIR ${PROJECT_BINARY_DIR}/thirdparty/mathpresso)
set(MATHPRESSO_DIR ${MATHPRESSO_LOCAL_DIR})
set(MATHPRESSO_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/mathpresso)

message( INFO ": MATHPRESSO_DIR = ${MATHPRESSO_DIR}" )
message( INFO ": MATHPRESSO_LOCAL_DIR = ${MATHPRESSO_LOCAL_DIR}" )
message( INFO ": MATHPRESSO_INSTALL_DIR = ${MATHPRESSO_INSTALL_DIR}" )

ExternalProject_Add( mathpresso
                     PREFIX ${PROJECT_BINARY_DIR}/thirdparty/mathpresso
#                     URL https://github.com/kobalicek/mathpresso/archive/master.zip
#                     URL_HASH MD5=b43212cafeab5e0e2ef5b87c29c15df1
#                     DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
#                     DOWNLOAD_NAME mathpresso.zip
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/mathpresso-master.zip
                     DEPENDS asmjit 
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

blt_register_library( NAME mathpresso
                      INCLUDES ${MATHPRESSO_INSTALL_DIR}/include
                      LIBRARIES ${MATHPRESSO_INSTALL_DIR}/lib/libmathpresso.a
                      TREAT_INCLUDES_AS_SYSTEM ON )

endif()



################################
# PUGIXML
################################
message( INFO ": setting up pugixml" )
set(PUGIXML_LOCAL_DIR ${PROJECT_BINARY_DIR}/thirdparty/pugixml)
set(PUGIXML_DIR ${PUGIXML_LOCAL_DIR})
set(PUGIXML_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/pugixml)

message( INFO ": PUGIXML_DIR = ${PUGIXML_DIR}" )
message( INFO ": PUGIXML_LOCAL_DIR = ${PUGIXML_LOCAL_DIR}" )
message( INFO ": PUGIXML_INSTALL_DIR = ${PUGIXML_INSTALL_DIR}" )

ExternalProject_Add( pugixml
                     PREFIX ${PROJECT_BINARY_DIR}/thirdparty/pugixml
#                     URL https://github.com/zeux/pugixml/archive/master.zip
#                     DOWNLOAD_NAME pugixml.zip
#                     URL_HASH MD5=0099563cb3f466fe03b10f9666c73993
#                     DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/mirror
                     URL ${CMAKE_SOURCE_DIR}/tpl_mirror/pugixml.tar
                     INSTALL_DIR ${PUGIXML_INSTALL_DIR}
                     INSTALL_COMMAND make install
                     CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	blt_register_library( NAME pugixml
    	                  INCLUDES ${PUGIXML_INSTALL_DIR}/include
        	              LIBRARIES ${PUGIXML_INSTALL_DIR}/lib/libpugixml.a
        	              TREAT_INCLUDES_AS_SYSTEM ON )
else()
	blt_register_library( NAME pugixml
    	                  INCLUDES ${PUGIXML_INSTALL_DIR}/include
        	              LIBRARIES ${PUGIXML_INSTALL_DIR}/lib64/libpugixml.a
        	              TREAT_INCLUDES_AS_SYSTEM ON )
endif()


if (UNCRUSTIFY_EXECUTABLE)
  include(cmake/blt/cmake/thirdparty/FindUncrustify.cmake)
endif()

message("Leaving SetupGeosxThirdParty.cmake\n")


