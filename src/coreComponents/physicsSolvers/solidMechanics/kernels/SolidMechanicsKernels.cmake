



set( kernelPath "coreComponents/physicsSolvers/solidMechanics/kernels" )

set( ExplicitSmallStrainPolicy "geos::parallelDevicePolicy< ${GEOSX_BLOCK_SIZE} >" )
set( ExplicitFiniteStrainPolicy "geos::parallelDevicePolicy< ${GEOSX_BLOCK_SIZE} >" )
set( FixedStressThermoPoromechanicsPolicy "geos::parallelDevicePolicy< ${GEOSX_BLOCK_SIZE} >" )
set( ImplicitSmallStrainNewmarkPolicy "geos::parallelDevicePolicy< ${GEOSX_BLOCK_SIZE} >" )
set( ImplicitSmallStrainQuasiStaticPolicy "geos::parallelDevicePolicy< ${GEOSX_BLOCK_SIZE} >" )
set( SmallStrainResidualPolicy "geos::parallelDevicePolicy< ${GEOSX_BLOCK_SIZE} >" )


configure_file( ${CMAKE_SOURCE_DIR}/${kernelPath}/policies.hpp.in
                ${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/policies.hpp )


set( kernelNames SolidMechanicsKernels )
set( subregionList CellElementSubRegion )
set( solidBaseDispatch DamageSpectral<ElasticIsotropic>
                       DamageVolDev<ElasticIsotropic>
                       Damage<ElasticIsotropic>
                       DruckerPragerExtended
                       ModifiedCamClay
                       DelftEgg
                       DruckerPrager
                       ElasticIsotropic
                       ElasticTransverseIsotropic
                       ElasticIsotropicPressureDependent
                       ElasticOrthotropic )

set( finiteElementDispatch H1_Hexahedron_Lagrange1_GaussLegendre2
                           H1_Wedge_Lagrange1_Gauss6
                           H1_Tetrahedron_Lagrange1_Gauss1
                           H1_Pyramid_Lagrange1_Gauss5
                           H1_Tetrahedron_VEM_Gauss1
                           H1_Prism5_VEM_Gauss1
                           H1_Prism6_VEM_Gauss1
                           H1_Prism7_VEM_Gauss1
                           H1_Prism8_VEM_Gauss1
                           H1_Prism9_VEM_Gauss1
                           H1_Prism10_VEM_Gauss1 )

if ( NOT ${ENABLE_HIP} )
  list(APPEND finiteElementDispatch
              # Q1_Hexahedron_Lagrange_GaussLobatto
              # Q2_Hexahedron_Lagrange_GaussLobatto
              # Q3_Hexahedron_Lagrange_GaussLobatto
              # Q4_Hexahedron_Lagrange_GaussLobatto
              # Q5_Hexahedron_Lagrange_GaussLobatto
              H1_Hexahedron_VEM_Gauss1
              H1_Wedge_VEM_Gauss1
              H1_Prism11_VEM_Gauss1 )
endif( )

  foreach( KERNELNAME ${kernelNames} )
    foreach( SUBREGION_TYPE  ${subregionList} )
      foreach( FE_TYPE ${finiteElementDispatch} )
        foreach( CONSTITUTIVE_TYPE ${solidBaseDispatch} )

        set( filename "${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/${KERNELNAME}_${SUBREGION_TYPE}_${CONSTITUTIVE_TYPE}_${FE_TYPE}.cpp" )
        string(REPLACE "<" "-" filename ${filename})
        string(REPLACE ">" "-" filename ${filename})
        string(REPLACE "," "-" filename ${filename})
        string(REPLACE " " "" filename ${filename})
        message( " -- Generating file: ${filename}")
        configure_file( ${CMAKE_SOURCE_DIR}/${kernelPath}/SolidMechanicsKernels.cpp.template
                        ${filename} )

          list( APPEND physicsSolvers_sources ${filename} )
        endforeach()
      endforeach()
    endforeach()
  endforeach()


  configure_file( ${CMAKE_SOURCE_DIR}/coreComponents/physicsSolvers/solidMechanics/kernels/SolidMechanicsSmallStrainResidualKernels.cpp
                  ${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/SolidMechanicsSmallStrainResidualKernels.cpp )

  list( APPEND physicsSolvers_sources ${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/SolidMechanicsSmallStrainResidualKernels.cpp )
  
  set( porousSolidDispatch PorousSolid<ElasticIsotropic> )

  set( kernelNames SolidMechanicsFixedStressThermoPoroElasticKernels )


  foreach( KERNELNAME ${kernelNames} )
    foreach( SUBREGION_TYPE  ${subregionList} )
      foreach( FE_TYPE ${finiteElementDispatch} )
        foreach( CONSTITUTIVE_TYPE ${porousSolidDispatch} )

        set( filename "${CMAKE_BINARY_DIR}/generatedSrc/${kernelPath}/${KERNELNAME}_${SUBREGION_TYPE}_${CONSTITUTIVE_TYPE}_${FE_TYPE}.cpp" )
        string(REPLACE "<" "-" filename ${filename})
        string(REPLACE ">" "-" filename ${filename})
        string(REPLACE "," "-" filename ${filename})
        string(REPLACE " " "" filename ${filename})
        message( " -- Generating file: ${filename}")
        configure_file( ${CMAKE_SOURCE_DIR}/${kernelPath}/SolidMechanicsFixedStressThermoPoromechanicsKernels.cpp.template
                        ${filename} )
          list( APPEND physicsSolvers_sources ${filename} )

        endforeach()
      endforeach()
    endforeach()
  endforeach()
