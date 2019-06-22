if( NOT LAPACK_DIR )

    find_package(LAPACK REQUIRED)
    
    list( GET LAPACK_LIBRARIES 0 LAPACK_LIB0 )
    get_filename_component( LAPACK_LIBRARY_DIR ${LAPACK_LIB0} DIRECTORY CACHE )
    get_filename_component( LAPACK_DIR ${LAPACK_LIBRARY_DIR} DIRECTORY CACHE )
    set( LAPACK_INCLUDE_DIR ${LAPACK_DIR}/include  CACHE PATH "")
    
else()
    find_path( LAPACK_INCLUDE_DIR 
               NAMES  lapack.h  mkl_lapack.h
               PATHS  ${LAPACK_DIR}/include
               NO_DEFAULT_PATH
               NO_CMAKE_ENVIRONMENT_PATH
               NO_CMAKE_PATH
               NO_SYSTEM_ENVIRONMENT_PATH
               NO_CMAKE_SYSTEM_PATH)
     
     find_library( LAPACK_LIBRARIES 
              NAMES lapack libmkl_rt.so 
              PATHS ${LAPACK_DIR}/lib ${LAPACK_DIR}/lib64
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)
              
    list( GET LAPACK_LIBRARIES 0 LAPACK_LIB0 )
    get_filename_component( LAPACK_LIBRARY_DIR ${LAPACK_LIB0} DIRECTORY CACHE )

endif()

find_package_handle_standard_args( LAPACK 
                                   DEFAULT_MSG
                                   LAPACK_DIR
                                   LAPACK_INCLUDE_DIR
                                   LAPACK_LIBRARY_DIR
                                   LAPACK_LIBRARIES
                                   )

set(LAPACK_LIBRARY_NAMES "" CACHE PATH "")
foreach( lib ${LAPACK_LIBRARIES} )
    get_filename_component( libname ${lib} NAME )
    list( APPEND LAPACK_LIBRARY_NAMES ${libname} )
endforeach()

message( "    LAPACK_FOUND         = ${LAPACK_FOUND}" )
message( "    LAPACK_DIR           = ${LAPACK_DIR}" )
message( "    LAPACK_INCLUDE_DIR   = ${LAPACK_INCLUDE_DIR}" )
message( "    LAPACK_LIBRARY_DIR   = ${LAPACK_LIBRARY_DIR}" )
message( "    LAPACK_LIBRARIES     = ${LAPACK_LIBRARIES}" )
message( "    LAPACK_LIBRARY_NAMES = ${LAPACK_LIBRARY_NAMES}" )
message( "    LAPACK_LINKER_FLAGS  = ${LAPACK_LINKER_FLAGS}" )
