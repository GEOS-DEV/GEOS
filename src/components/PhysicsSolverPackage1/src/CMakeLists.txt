
#
# Specify all headers
#
set(componentHeaders "SolidMechanicsLagrangianFEM.hpp" )

#
# Specify all sources
#

set(componentSources "SolidMechanicsLagrangianFEM.cpp" )

blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)

message("adding physicsSolverPackage1 library")
blt_add_library( NAME             PhysicsSolverPackage1
                 SOURCES          ${componentSources}
                 HEADERS          ${componentHeaders}
                 DEPENDS_ON        geosx_core cxx-utilities
                 SYSTEM_DEPENDS_ON slic sidre conduit hdf5 RAJA chai
                 )

target_include_directories( PhysicsSolverPackage1 PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
target_include_directories( PhysicsSolverPackage1 PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src )
