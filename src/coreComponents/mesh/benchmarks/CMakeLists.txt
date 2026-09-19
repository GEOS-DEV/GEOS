#
# Mesh benchmarks
#

set( dependencyList mesh ${parallelDeps} ${tplDependencyList} )

blt_add_executable(
  NAME benchmarkScatterMethods
  SOURCES benchmarkScatterMethods.cpp
  OUTPUT_DIR ${BENCHMARK_OUTPUT_DIRECTORY}
  DEPENDS_ON ${dependencyList}
  FOLDER coreComponents/mesh/benchmarks
)

install(TARGETS benchmarkScatterMethods DESTINATION ${CMAKE_INSTALL_PREFIX}/benchmarks)