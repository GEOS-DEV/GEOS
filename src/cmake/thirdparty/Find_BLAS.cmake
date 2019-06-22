
if( NOT BLAS_DIR )

    find_package(BLAS REQUIRED)
    
    list( GET BLAS_LIBRARIES 0 BLAS_LIB0 )
    get_filename_component( BLAS_LIBRARY_DIR ${BLAS_LIB0} DIRECTORY CACHE )
    get_filename_component( BLAS_DIR ${BLAS_LIBRARY_DIR} DIRECTORY CACHE )
    set( BLAS_INCLUDE_DIR ${BLAS_DIR}/include CACHE PATH "" )
    
else()
    find_path( BLAS_INCLUDE_DIR 
               NAMES  blas.h cblas.h mkl_blas.h mkl_cblas.h
               PATHS  ${BLAS_DIR}/include
               NO_DEFAULT_PATH
               NO_CMAKE_ENVIRONMENT_PATH
               NO_CMAKE_PATH
               NO_SYSTEM_ENVIRONMENT_PATH
               NO_CMAKE_SYSTEM_PATH)
     
    if( NOT DEFINED BLAS_LIBRARY_NAMES )
        set( BLAS_LIBRARY_NAMES "blas;cblas;libmkl_rt.so" CACHE STRING "")
    endif()
    
    find_library( BLAS_LIBRARIES 
                  NAMES ${BLAS_LIBRARY_NAMES}
                  PATHS ${BLAS_DIR}/lib ${BLAS_DIR}/lib64
                  NO_DEFAULT_PATH
                  NO_CMAKE_ENVIRONMENT_PATH
                  NO_CMAKE_PATH
                  NO_SYSTEM_ENVIRONMENT_PATH
                  NO_CMAKE_SYSTEM_PATH)
              
    list( GET BLAS_LIBRARIES 0 BLAS_LIB0 )
    get_filename_component( BLAS_LIBRARY_DIR ${BLAS_LIB0} DIRECTORY CACHE )

endif()

find_package_handle_standard_args( BLAS 
                                   DEFAULT_MSG
                                   BLAS_DIR
                                   BLAS_INCLUDE_DIR
                                   BLAS_LIBRARY_DIR
                                   BLAS_LIBRARIES
                                   )

set(BLAS_LIBRARY_NAMES "" CACHE STRING "" FORCE)
foreach( lib ${BLAS_LIBRARIES} )
    get_filename_component( libname ${lib} NAME )
    list( APPEND BLAS_LIBRARY_NAMES ${libname} )
endforeach()

message( "    BLAS_FOUND         = ${BLAS_FOUND}" )
message( "    BLAS_DIR           = ${BLAS_DIR}" )
message( "    BLAS_INCLUDE_DIR   = ${BLAS_INCLUDE_DIR}" )
message( "    BLAS_LIBRARY_DIR   = ${BLAS_LIBRARY_DIR}" )
message( "    BLAS_LIBRARIES     = ${BLAS_LIBRARIES}" )
message( "    BLAS_LIBRARY_NAMES = ${BLAS_LIBRARY_NAMES}" )
message( "    BLAS_LINKER_FLAGS  = ${BLAS_LINKER_FLAGS}" )
