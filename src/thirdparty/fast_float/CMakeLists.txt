message("adding fast_float library")
    
blt_add_library( NAME fast_float
                 HEADERS include/fast_float.h
                 )
                 
target_include_directories( fast_float SYSTEM INTERFACE ${CMAKE_SOURCE_DIR}/thirdparty/fast_float/include )
                 