#--------------------------------------------------------------------------------------------------
#[[
Package: mixedMimetic

Contains the mixed mimetic finite difference discretization infrastructure:
the discretization object registered under the NumericalMethods block, the
residual-based Global Adaptation kernels, and the dispatch over the mixed-form
inner products.
#]]
#--------------------------------------------------------------------------------------------------

#
# Specify all headers
#
set( mixedMimetic_headers
     MixedMimeticDiscretization.hpp
     MixedMimeticDiscretizationManager.hpp
     MixedMimeticDispatch.hpp
     MixedMimeticFields.hpp
     adaptivity/GlobalAdaptationKernels.hpp )

#
# Specify all sources
#
set( mixedMimetic_sources
     MixedMimeticDiscretization.cpp
     MixedMimeticDiscretizationManager.cpp )

set( dependencyList ${parallelDeps} finiteVolume )

geos_decorate_link_dependencies( LIST decoratedDependencies DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       mixedMimetic
                 SOURCES    ${mixedMimetic_sources}
                 HEADERS    ${mixedMimetic_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS} )

target_include_directories( mixedMimetic PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS mixedMimetic LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
