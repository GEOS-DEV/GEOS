# Specify list of tests
set( gtest_geosx_tests
     testFEMConsistency.cpp )

set( tplDependencyList ${parallelDeps} gtest )

set( dependencyList mainInterface testingUtilities )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

# Add gtest C++ based tests
foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${decoratedDependencies} ${tplDependencyList} )

  geos_add_test( NAME ${test_name}
                 COMMAND ${test_name} )

  set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "TEST_BINARY_DIR=$<TARGET_FILE_DIR:${test_name}>")

  # Add TEST_BINARY_DIR as a compile definition
  target_compile_definitions(${test_name} PRIVATE TEST_BINARY_DIR=\"$<TARGET_FILE_DIR:${test_name}>\")

  # --- Copy mesh files to test binary dir ---
  set( mesh_files
       fractured_mesh_hex_DFN_1.vtu
       fractured_mesh_hex_DFN_12.vtu
       fractured_mesh_hex_DFN_123.vtu
       fractured_mesh_hex_DFN_13.vtu
       fractured_mesh_hex_DFN_2.vtu
       fractured_mesh_hex_DFN_23.vtu
       fractured_mesh_hex_DFN_3.vtu
       fractured_mesh_tet_DFN_1.vtu
       fractured_mesh_tet_DFN_12.vtu
       fractured_mesh_tet_DFN_123.vtu
       fractured_mesh_tet_DFN_13.vtu
       fractured_mesh_tet_DFN_2.vtu
       fractured_mesh_tet_DFN_23.vtu
       fractured_mesh_tet_DFN_3.vtu
       fractured_wavy_mesh_hex_DFN_1.vtu
       fractured_wavy_mesh_hex_DFN_12.vtu
       fractured_wavy_mesh_hex_DFN_123.vtu
       fractured_wavy_mesh_hex_DFN_13.vtu
       fractured_wavy_mesh_hex_DFN_2.vtu
       fractured_wavy_mesh_hex_DFN_23.vtu
       fractured_wavy_mesh_hex_DFN_3.vtu
       fractured_wavy_mesh_tet_DFN_1.vtu
       fractured_wavy_mesh_tet_DFN_12.vtu
       fractured_wavy_mesh_tet_DFN_123.vtu
       fractured_wavy_mesh_tet_DFN_13.vtu
       fractured_wavy_mesh_tet_DFN_2.vtu
       fractured_wavy_mesh_tet_DFN_23.vtu
       fractured_wavy_mesh_tet_DFN_3.vtu
  )

  foreach( mesh_file ${mesh_files} )
    add_custom_command(
        TARGET ${test_name} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_CURRENT_SOURCE_DIR}/${mesh_file}
                $<TARGET_FILE_DIR:${test_name}>
    )
  endforeach()
endforeach()
