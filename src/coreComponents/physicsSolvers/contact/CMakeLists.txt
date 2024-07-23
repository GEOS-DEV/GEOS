# Specify solver headers
set( physicsSolversContact_headers
     contact/ContactSolverBase.hpp
     contact/ContactFields.hpp
     contact/SolidMechanicsEFEMKernelsBase.hpp
     contact/SolidMechanicsEFEMKernels.hpp
     contact/SolidMechanicsEFEMStaticCondensationKernels.hpp
     contact/SolidMechanicsEFEMKernelsHelper.hpp
     contact/SolidMechanicsEmbeddedFractures.hpp
     contact/SolidMechanicsLagrangeContact.hpp
     contact/SolidMechanicsAugmentedLagrangianContact.hpp
     contact/SolidMechanicsALMKernelsBase.hpp
     contact/SolidMechanicsALMKernels.hpp
     contact/SolidMechanicsALMKernelsHelper.hpp
     contact/SolidMechanicsALMJumpUpdateKernels.hpp
     contact/SolidMechanicsALMBubbleKernels.hpp
 )

# Specify solver sources
set( physicsSolversContact_sources
     contact/ContactSolverBase.cpp
     contact/SolidMechanicsEmbeddedFractures.cpp
     contact/SolidMechanicsLagrangeContact.cpp
     contact/SolidMechanicsAugmentedLagrangianContact.cpp
      )

set( dependencyList ${parallelDeps} constitutive mesh linearAlgebra discretizationMethods events )

blt_add_library( NAME       physicsSolversContact
                 SOURCES    ${physicsSolversContact_sources}
                 HEADERS    ${physicsSolversContact_headers}
                 DEPENDS_ON ${dependencyList} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS} )

target_include_directories( physicsSolversContact PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if( externalComponentDeps )
  target_include_directories( physicsSolversContact PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()


