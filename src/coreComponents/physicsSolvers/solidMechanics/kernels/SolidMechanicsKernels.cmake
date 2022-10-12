



set( kernelPath "coreComponents/physicsSolvers/solidMechanics/kernels" )

set( ExplicitSmallStrainPolicy "geosx::parallelDevicePolicy<32>" )
set( ExplicitFiniteStrainPolicy "geosx::parallelDevicePolicy<32>" )
set( ImplicitSmallStrainQuasiStaticPolicy "geosx::parallelDevicePolicy<32>" )


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
                           Q3_Hexahedron_Lagrange_GaussLobatto
                           H1_Tetrahedron_VEM_Gauss1
                           H1_Wedge_VEM_Gauss1
                           H1_Hexahedron_VEM_Gauss1
                           H1_Prism5_VEM_Gauss1
                           H1_Prism6_VEM_Gauss1
                           H1_Prism7_VEM_Gauss1
                           H1_Prism8_VEM_Gauss1
                           H1_Prism9_VEM_Gauss1
                           H1_Prism10_VEM_Gauss1
                           H1_Prism11_VEM_Gauss1 )

  foreach( KERNELNAME ${kernelNames} )
    foreach( SUBREGION_TYPE  ${subregionList} )
      foreach( CONSTITUTIVE_TYPE ${solidBaseDispatch} )
        foreach( FE_TYPE ${finiteElementDispatch} )

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