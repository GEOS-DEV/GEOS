#
# Mesh benchmarks
#

set( dependencyList mesh ${parallelDeps} ${tplDependencyList} )

if( ENABLE_MPI )
  blt_add_executable(
    NAME benchmarkScatterMethods
    SOURCES benchmarkScatterMethods.cpp
    OUTPUT_DIR ${BENCHMARK_OUTPUT_DIRECTORY}
    DEPENDS_ON ${dependencyList}
    FOLDER coreComponents/mesh/benchmarks
  )
  set_source_files_properties( benchmarkScatterMethods.cpp PROPERTIES LANGUAGE CXX )
  install(TARGETS benchmarkScatterMethods DESTINATION ${CMAKE_INSTALL_PREFIX}/benchmarks)
endif()

blt_add_executable(
  NAME benchmarkUniformRefinement
  SOURCES benchmarkUniformRefinement.cpp
  OUTPUT_DIR ${BENCHMARK_OUTPUT_DIRECTORY}
  DEPENDS_ON ${dependencyList}
  FOLDER coreComponents/mesh/benchmarks
)
set_source_files_properties( benchmarkUniformRefinement.cpp PROPERTIES LANGUAGE CXX )
install(TARGETS benchmarkUniformRefinement DESTINATION ${CMAKE_INSTALL_PREFIX}/benchmarks)

# Dynamic-call audit for Linux shared builds; it forwards all calls unchanged.
if( CMAKE_SYSTEM_NAME STREQUAL "Linux" AND GEOS_BUILD_SHARED_LIBS )
  blt_add_library( NAME auditUniformRefinementCalls
                   SOURCES auditUniformRefinementCalls.cpp
                   DEPENDS_ON ${dependencyList} ${CMAKE_DL_LIBS}
                   SHARED TRUE
                   FOLDER coreComponents/mesh/benchmarks )
  set_source_files_properties( auditUniformRefinementCalls.cpp PROPERTIES LANGUAGE CXX )
  # Interposed constructors forward with RTLD_NEXT. Keep their libraries loaded
  # even when the executable itself does not use those classes (the audit probe
  # deliberately calls them). GCC's --as-needed otherwise drops these targets.
  target_link_options( auditUniformRefinementCalls PRIVATE "LINKER:--no-as-needed" )
  if( GEOS_ENABLE_TESTS )
    find_package( Python3 REQUIRED COMPONENTS Interpreter )
    set( audit_runner ${CMAKE_SOURCE_DIR}/../benchmarks/auditUniformRefinementCalls.py )
    set( audit_filter --gtest_filter=VTKUniformRefinement.ConnectedHexesMarkersFieldsAndLineageThroughTwoLevels )
    add_test( NAME testUniformRefinementAudit
              COMMAND ${Python3_EXECUTABLE} ${audit_runner}
                      --library $<TARGET_FILE:auditUniformRefinementCalls> --
                      $<TARGET_FILE:testVTKUniformRefinement> ${audit_filter} )
    add_test( NAME testUniformRefinementAuditProbe
              COMMAND ${Python3_EXECUTABLE} ${audit_runner}
                      --library $<TARGET_FILE:auditUniformRefinementCalls> --expect-probe --
                      $<TARGET_FILE:testVTKUniformRefinement> ${audit_filter} )
    set_tests_properties( testUniformRefinementAudit testUniformRefinementAuditProbe PROPERTIES TIMEOUT 60 )
    if( ENABLE_MPI AND NOT ENABLE_CUDA AND NOT ENABLE_HIP )
      add_test( NAME testUniformRefinementAudit_4ranks
                COMMAND ${Python3_EXECUTABLE} ${audit_runner}
                        --library $<TARGET_FILE:auditUniformRefinementCalls> --minimum-calls 4 --
                        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:testVTKUniformRefinement> -x 4 ${audit_filter} ${MPIEXEC_POSTFLAGS} )
      set_tests_properties( testUniformRefinementAudit_4ranks PROPERTIES TIMEOUT 60 )
      add_test( NAME testUniformRefinementAuditFractures_4ranks
                COMMAND ${Python3_EXECUTABLE} ${audit_runner}
                        --library $<TARGET_FILE:auditUniformRefinementCalls> --minimum-calls 8 --
                        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                        $<TARGET_FILE:testVTKRefinedFractureImport_mpi> -x 4
                        --gtest_filter=VTKRefinedFractureImport.FullInitializationOfNamedJunctionRegions
                        ${MPIEXEC_POSTFLAGS} )
      set_tests_properties( testUniformRefinementAuditFractures_4ranks PROPERTIES TIMEOUT 60 )
    endif()
  endif()
endif()
