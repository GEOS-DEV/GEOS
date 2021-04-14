message("adding optionparser library")
#blt_import_library( 
#    NAME optionparser
#    INCLUDES optionparser/optionparser.h
#    TREAT_INCLUDES_AS_SYSTEM ON
#    )
    
blt_add_library( NAME     optionparser
                 HEADERS  src/optionparser.h
                 )
                 
target_include_directories( optionparser SYSTEM INTERFACE ${CMAKE_SOURCE_DIR}/thirdparty/optionparser/src)
                 