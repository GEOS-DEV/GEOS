/*
 * FiniteElementDispatch.hpp
 *
 *  Created on: Jul 6, 2020
 *      Author: settgast
 */

#ifndef GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_
#define GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_


#include "elementFormulations/BiLinearWedgeShapeFunctionKernel.hpp"
#include "elementFormulations/LinearTetrahedronShapeFunctionKernel.hpp"
#include "elementFormulations/TrilinearHexahedronShapeFunctionKernel.hpp"



namespace geosx
{
namespace finiteElement
{

enum class ParentElementType
{
  Tetrahedal,
  Pyramid,
  Prism,
  Hexahedral,
  Polyhedral,
  Polytope,
  INVALID
};

struct ParentElementTypeStrings
{
  static constexpr auto Tetrahedal  = "C3D4";
  static constexpr auto Pyramid     = "C3D5";
  static constexpr auto Prism       = "C3D6";
  static constexpr auto Hexahedral  = "C3D8";
  static constexpr auto Polyhedral  = "POLYHEDRAL";
  static constexpr auto Polytope    = "POLYTOPE";
};



template< typename LAMBDA >
void
dispatch( ParentElementType const input,
          LAMBDA && lambda )
{
  switch( input )
  {
    case ParentElementType::Hexahedral:
    {
      lambda( TrilinearHexahedronShapeFunctionKernel() );
      break;
    }
    case ParentElementType::Tetrahedal:
    {
      lambda( LinearTetrahedronShapeFunctionKernel() );
      break;
    }
    case ParentElementType::Prism:
    {
      lambda( BiLinearWedgeShapeFunctionKernel() );
      break;
    }
    default:
      GEOSX_ERROR( "integralTypeDispatch() is not implemented for value of: " );
  }
}


template< typename LAMBDA >
void
dispatch( string const & input,
          LAMBDA && lambda )
{
  if( input ==  ParentElementTypeStrings::Hexahedral )
  {
    lambda( TrilinearHexahedronShapeFunctionKernel() );
  }
  else if( input == ParentElementTypeStrings::Tetrahedal )
  {
    lambda( LinearTetrahedronShapeFunctionKernel() );
  }
  else if( input == ParentElementTypeStrings::Prism )
  {
    lambda( BiLinearWedgeShapeFunctionKernel() );
  }
  else
  {
    GEOSX_ERROR( "integralTypeDispatch() is not implemented for value of: "<<input );
  }
}

}
}



#endif /* GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_ */
