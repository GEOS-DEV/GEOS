/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FiniteElementDispatch
 */

#ifndef GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_
#define GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_


#include "elementFormulations/BiLinearQuadrilateralFaceShapeFunctionKernel.hpp"
#include "elementFormulations/BiLinearWedgeShapeFunctionKernel.hpp"
#include "elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp"
#include "elementFormulations/H1_Tetrahedron_Lagrange1_Gauss1.hpp"
#include "elementFormulations/LinearTriangleFaceShapeFunctionKernel.hpp"
#include "elementFormulations/PyramidShapeFunctionKernel.hpp"



namespace geosx
{
namespace finiteElement
{

struct ParentElementTypeStrings
{
  static constexpr auto Tetrahedon    = "C3D4";
  static constexpr auto Pyramid       = "C3D5";
  static constexpr auto Prism         = "C3D6";
  static constexpr auto Hexahedron    = "C3D8";
  static constexpr auto Quadralateral = "C2D4";
  static constexpr auto Triangle      = "C2D3";
  static constexpr auto Polyhedral    = "POLYHEDRAL";
  static constexpr auto Polytope      = "POLYTOPE";
};

template< typename LAMBDA >
void
dispatch3D( FiniteElementBase const & input,
            LAMBDA && lambda )
{
  if( dynamic_cast<H1_Hexahedron_Lagrange1_GaussLegendre2 const *>(&input) )
  {
    lambda( static_cast<H1_Hexahedron_Lagrange1_GaussLegendre2 const &>(input) );
  }
  else if( dynamic_cast<BiLinearWedgeShapeFunctionKernel const *>(&input) )
  {
    lambda( static_cast<BiLinearWedgeShapeFunctionKernel const &>(input) );
  }
  else if( dynamic_cast<H1_Tetrahedron_Lagrange1_Gauss1 const *>(&input) )
  {
    lambda( static_cast<H1_Tetrahedron_Lagrange1_Gauss1 const &>(input) );
  }
  else if( dynamic_cast<PyramidShapeFunctionKernel const *>(&input) )
  {
    lambda( static_cast<PyramidShapeFunctionKernel const &>(input) );
  }
  else
  {
    GEOSX_ERROR( "finiteElement::dispatch3D() is not implemented for input of "<<typeid(input).name() );
  }
}

template< typename LAMBDA >
void
dispatch2D( FiniteElementBase const & input,
            LAMBDA && lambda )
{
  if( dynamic_cast<BiLinearQuadrilateralFaceShapeFunctionKernel const *>(&input) )
  {
    lambda( static_cast<BiLinearQuadrilateralFaceShapeFunctionKernel const &>(input) );
  }
  else if( dynamic_cast<LinearTriangleFaceShapeFunctionKernel const *>(&input) )
  {
    lambda( static_cast<LinearTriangleFaceShapeFunctionKernel const &>(input) );
  }
  else
  {
    GEOSX_ERROR( "finiteElement::dispatch2D() is not implemented for input of: "<<typeid(input).name() );
  }
}

}
}



#endif /* GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_ */
