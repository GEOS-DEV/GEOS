/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FiniteElementDispatch.hpp
 */

#ifndef GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_
#define GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_


#include "elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp"
#include "elementFormulations/H1_Pyramid_Lagrange1_Gauss5.hpp"
#include "elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp"
#include "elementFormulations/H1_Tetrahedron_Lagrange1_Gauss1.hpp"
#include "elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp"
#include "elementFormulations/H1_Wedge_Lagrange1_Gauss6.hpp"
#include "elementFormulations/Q3_Hexahedron_Lagrange_GaussLobatto.hpp"

#include "LvArray/src/system.hpp"



namespace geosx
{
namespace finiteElement
{

template< typename LAMBDA >
void
dispatch3D( FiniteElementBase const & input,
            LAMBDA && lambda )
{
  if( auto const * const ptr1 = dynamic_cast< H1_Hexahedron_Lagrange1_GaussLegendre2 const * >(&input) )
  {
    lambda( *ptr1 );
  }
  else if( auto const * const ptr2 = dynamic_cast< H1_Wedge_Lagrange1_Gauss6 const * >(&input) )
  {
    lambda( *ptr2 );
  }
  else if( auto const * const ptr3 = dynamic_cast< H1_Tetrahedron_Lagrange1_Gauss1 const * >(&input) )
  {
    lambda( *ptr3 );
  }
  else if( auto const * const ptr4 = dynamic_cast< H1_Pyramid_Lagrange1_Gauss5 const * >(&input) )
  {
    lambda( *ptr4 );
  }
  else
  {
    GEOSX_ERROR( "finiteElement::dispatch3D() is not implemented for input of "<<typeid(input).name() );
  }
}


template< typename LAMBDA >
void
dispatch3D( FiniteElementBase & input,
            LAMBDA && lambda )
{
  if( auto * const ptr1 = dynamic_cast< H1_Hexahedron_Lagrange1_GaussLegendre2 * >(&input) )
  {
    lambda( *ptr1 );
  }
  else if( auto * const ptr2 = dynamic_cast< H1_Wedge_Lagrange1_Gauss6 * >(&input) )
  {
    lambda( *ptr2 );
  }
  else if( auto * const ptr3 = dynamic_cast< H1_Tetrahedron_Lagrange1_Gauss1 * >(&input) )
  {
    lambda( *ptr3 );
  }
  else if( auto * const ptr4 = dynamic_cast< H1_Pyramid_Lagrange1_Gauss5 * >(&input) )
  {
    lambda( *ptr4 );
  }
  else if( auto * const ptr5 = dynamic_cast< Q3_Hexahedron_Lagrange_GaussLobatto * >(&input) )
  {
    lambda( *ptr5 );
  }
  else
  {
    GEOSX_ERROR( "finiteElement::dispatch3D() is not implemented for input of "<<LvArray::system::demangleType( &input ) );
  }
}

template< typename LAMBDA >
void
dispatch2D( FiniteElementBase const & input,
            LAMBDA && lambda )
{
  if( auto const * const ptr1 = dynamic_cast< H1_QuadrilateralFace_Lagrange1_GaussLegendre2 const * >(&input) )
  {
    lambda( *ptr1 );
  }
  else if( auto const * const ptr2 = dynamic_cast< H1_TriangleFace_Lagrange1_Gauss1 const * >(&input) )
  {
    lambda( *ptr2 );
  }
  else
  {
    GEOSX_ERROR( "finiteElement::dispatch2D() is not implemented for input of: "<<LvArray::system::demangleType( &input ) );
  }
}

}
}



#endif /* GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_ */
