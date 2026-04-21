#
# Mesh benchmarks
#

blt_add_executable(
  NAME benchmarkScatterMethods
  SOURCES benchmarkScatterMethods.cpp
  OUTPUT_DIR ${BENCHMARK_OUTPUT_DIRECTORY}
  DEPENDS_ON mesh
  FOLDER coreComponents/mesh/benchmarks
)

install(TARGETS benchmarkScatterMethods DESTINATION ${CMAKE_INSTALL_PREFIX}/benchmarks)