set( dispatchImplPath "coreComponents/finiteElement/dispatchImpl" )

set( finiteElementDispatch3D H1_Hexahedron_Lagrange1_GaussLegendre2
                             H1_Hexahedron_Lagrange1_GaussLegendre2_DEBUG1
                             H1_Hexahedron_Lagrange1_GaussLegendre2_DEBUG2
                             H1_Hexahedron_Lagrange1_GaussLegendre2_DEBUG3
                             H1_Hexahedron_Lagrange1_GaussLegendre2_DEBUG4
                             H1_Wedge_Lagrange1_Gauss6
                             H1_Tetrahedron_Lagrange1_Gauss1
                             H1_Pyramid_Lagrange1_Gauss5
                             H1_Tetrahedron_VEM_Gauss1
                             H1_Wedge_VEM_Gauss1
                             H1_Hexahedron_VEM_Gauss1
                             H1_Prism5_VEM_Gauss1
                             H1_Prism6_VEM_Gauss1
                             H1_Prism7_VEM_Gauss1
                             H1_Prism8_VEM_Gauss1
                             H1_Prism9_VEM_Gauss1
                             H1_Prism10_VEM_Gauss1
                             H1_Prism11_VEM_Gauss1
                             Q3_Hexahedron_Lagrange_GaussLobatto )
set( finiteElementDispatch2D H1_QuadrilateralFace_Lagrange1_GaussLegendre2
                             H1_TriangleFace_Lagrange1_Gauss1 )

  foreach( FE_TYPE ${finiteElementDispatch3D} )

    set( filename "${CMAKE_BINARY_DIR}/generatedSrc/${dispatchImplPath}/FiniteElementDispatch3D_${FE_TYPE}.cpp" )
    string(REPLACE "<" "-" filename ${filename})
    string(REPLACE ">" "-" filename ${filename})
    string(REPLACE "," "-" filename ${filename})
    string(REPLACE " " "" filename ${filename})
    message( " -- Generating file: ${filename}")
    configure_file( ${CMAKE_SOURCE_DIR}/${dispatchImplPath}/FiniteElementDispatch3D.cpp.template
                    ${filename} )

    list( APPEND finiteElement_sources ${filename} )
  endforeach()

  foreach( FE_TYPE ${finiteElementDispatch2D} )

    set( filename "${CMAKE_BINARY_DIR}/generatedSrc/${dispatchImplPath}/FiniteElementDispatch2D_${FE_TYPE}.cpp" )
    string(REPLACE "<" "-" filename ${filename})
    string(REPLACE ">" "-" filename ${filename})
    string(REPLACE "," "-" filename ${filename})
    string(REPLACE " " "" filename ${filename})
    message( " -- Generating file: ${filename}")
    configure_file( ${CMAKE_SOURCE_DIR}/${dispatchImplPath}/FiniteElementDispatch2D.cpp.template
                    ${filename} )

    list( APPEND finiteElement_sources ${filename} )
  endforeach()
