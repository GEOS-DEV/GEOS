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


#include "elementFormulations/BiLinearWedgeShapeFunctionKernel.hpp"
#include "elementFormulations/LinearTetrahedronShapeFunctionKernel.hpp"
#include "elementFormulations/TrilinearHexahedronShapeFunctionKernel.hpp"
#include "elementFormulations/PyramidShapeFunctionKernel.hpp"
#include "elementFormulations/BiLinearQuadrilateralFaceShapeFunctionKernel.hpp"
#include "elementFormulations/LinearTriangleFaceShapeFunctionKernel.hpp"



namespace geosx
{
namespace finiteElement
{

struct ParentElementTypeStrings
{
  static constexpr auto Tetrahedal    = "C3D4";
  static constexpr auto Pyramid       = "C3D5";
  static constexpr auto Prism         = "C3D6";
  static constexpr auto Hexahedral    = "C3D8";
  static constexpr auto Quadralateral = "C2D4";
  static constexpr auto Triangle      = "C2D3";
  static constexpr auto Polyhedral    = "POLYHEDRAL";
  static constexpr auto Polytope      = "POLYTOPE";
};


//template< typename LAMBDA >
//void
//dispatch( string const & input,
//          LAMBDA && lambda )
//{
//  if( input ==  ParentElementTypeStrings::Hexahedral )
//  {
//    lambda( TrilinearHexahedronShapeFunctionKernel() );
//  }
//  else if( input == ParentElementTypeStrings::Tetrahedal )
//  {
//    lambda( LinearTetrahedronShapeFunctionKernel() );
//  }
//  else if( input == ParentElementTypeStrings::Prism )
//  {
//    lambda( BiLinearWedgeShapeFunctionKernel() );
//  }
//  else
//  {
//    GEOSX_ERROR( "integralTypeDispatch() is not implemented for value of: "<<input );
//  }
//}

template< typename LAMBDA >
void
dispatch3D( FiniteElementShapeFunctionKernelBase const & input,
            LAMBDA && lambda )
{
  if( dynamic_cast<TrilinearHexahedronShapeFunctionKernel const *>(&input) )
  {
    lambda( static_cast<TrilinearHexahedronShapeFunctionKernel const &>(input) );
  }
  else if( dynamic_cast<BiLinearWedgeShapeFunctionKernel const *>(&input) )
  {
    lambda( static_cast<BiLinearWedgeShapeFunctionKernel const &>(input) );
  }
  else if( dynamic_cast<LinearTetrahedronShapeFunctionKernel const *>(&input) )
  {
    lambda( static_cast<LinearTetrahedronShapeFunctionKernel const &>(input) );
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
dispatch2D( FiniteElementShapeFunctionKernelBase const & input,
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
