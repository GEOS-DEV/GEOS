set(pugixml_headers
    src/pugixml.hpp
    src/pugiconfig.hpp )

set(pugixml_sources
    src/pugixml.cpp )


message("adding pugixml library")
blt_add_library( NAME             pugixml
                 SOURCES          ${pugixml_sources}
                 HEADERS          ${pugixml_headers}
                 HEADERS_OUTPUT_SUBDIR pugixml
                 DEPENDS_ON        
                 SYSTEM_DEPENDS_ON 
                 )


if( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
set_target_properties(pugixml PROPERTIES COMPILE_FLAGS "-Wno-everything -Wno-switch-enum -Wno-float-equal -Wno-deprecated -Wno-reserved-id-macro -Wno-non-virtual-dtor")
endif()

#target_include_directories( pugixml PUBLIC ${CMAKE_CURRENT_LIST_DIR}/src )

