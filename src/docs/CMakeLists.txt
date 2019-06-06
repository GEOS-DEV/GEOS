#add_subdirectory(sphinx)

if(DOXYGEN_FOUND)
    add_subdirectory(doxygen)
endif()


  add_custom_target( geosx_update_rst_tables
                     ALL
                     COMMAND bash ../scripts/updateRstTables.bash > /dev/null 2>&1
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     DEPENDS geosx
                     )