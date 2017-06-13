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
set(RAJA_LOCAL_DIR ${CMAKE_SOURCE_DIR}/thirdparty/raja)
if(EXISTS ${RAJA_DIR})
    message(INFO ": Using system RAJA found at ${RAJA_DIR}")
    include(${PROJECT_SOURCE_DIR}/cmake/thirdparty/FindRAJA.cmake)
    if (NOT RAJA_FOUND)
        message(FATAL_ERROR "RAJA not found in ${RAJA_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME raja
                          INCLUDES ${RAJA_INCLUDE_DIRS} 
                          LIBRARIES RAJA )

else()
    message(INFO ": Using RAJA found at https://github.com/LLNL/RAJA/archive/develop.zip")
    set(raja_install_dir ${CMAKE_INSTALL_PREFIX}/thirdparty/raja2)
    ExternalProject_Add( raja
                         URL https://github.com/LLNL/RAJA/archive/develop.zip
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
                          LIBRARIES ${raja_install_dir}/lib/libRAJA.a )

endif()


################################
# CHAI
################################
message( "CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}")
message( "CHAI_DIR=${CHAI_DIR}")
if( EXISTS ${CHAI_DIR})

	message( INFO ": CHAI_DIR = ${CHAI_DIR}" )

	if( NOT BUILD_LOCAL_CHAI )
	    message(INFO ": Using system CHAI found at ${CHAI_DIR}")
	    include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindCHAI.cmake)
	    if (NOT CHAI_FOUND)
	        message(FATAL_ERROR ": CHAI not found in ${CHAI_DIR}. Maybe you need to build it")
	    endif()
	
	    blt_register_library( NAME chai
	                          INCLUDES ${CHAI_INCLUDE_DIRS}
	                          LIBRARIES ${CHAI_LIBRARY}  )
	else()
		message(INFO ": Build CHAI from source found at {CHAI_DIR}")
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
	                                    -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )
	                             
	    blt_register_library( NAME chai
	                          INCLUDES ${chai_install_dir}/include 
	                          LIBRARIES ${chai_install_dir}/lib/libchai.a )		
	endif()
else()
	set(CHAI_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/chai)
	message( INFO ": CHAI_INSTALL_DIR = ${CHAI_INSTALL_DIR}" )

    message(INFO ": Using CHAI found at ssh://git@cz-bitbucket.llnl.gov:7999/um/chai.git")
    ExternalProject_Add( chai
                         PREFIX ${PROJECT_BINARY_DIR}/thirdparty/chai
       					 GIT_REPOSITORY ssh://git@cz-bitbucket.llnl.gov:7999/um/chai.git
        				 GIT_TAG master
                         INSTALL_DIR ${chai_install_dir}
                         INSTALL_COMMAND make install
	                     CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        	                        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        	                        -DENABLE_CUDA=${CUDA_ENABLED}
        	                        -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR> )

    blt_register_library( NAME chai
                          INCLUDES ${chai_install_dir}/include 
                          LIBRARIES ${chai_install_dir}/lib/libchai.a  )
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
                     URL http://warp.povusers.org/FunctionParser/fparser4.5.2.zip
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
                      LIBRARIES ${FPARSER_INSTALL_DIR}/lib/libfparser.a )

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
                     PREFIX ${PROJECT_BINARY_DIR}/thirdparty/caliper
                     GIT_REPOSITORY https://github.com/LLNL/Caliper.git
                     GIT_TAG master
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
                      LIBRARIES ${CALIPER_INSTALL_DIR}/lib/libcaliper-reader.a )

endif()

if (UNCRUSTIFY_EXECUTABLE)
  include(cmake/blt/cmake/thirdparty/FindUncrustify.cmake)
endif()

message("Leaving SetupGeosxThirdParty.cmake\n")


