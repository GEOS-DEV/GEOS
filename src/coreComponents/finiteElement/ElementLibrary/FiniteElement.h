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

#ifndef FINITE_ELEMENT_H
#define FINITE_ELEMENT_H

/**
 * @file FiniteElement.h
 */

#include "finiteElement/basis/BasisBase.hpp"
#include "finiteElement/quadrature/QuadratureBase.hpp"
#include "FiniteElementBase.h"
#include "common/TimingMacros.hpp"
#include "LvArray/src/tensorOps.hpp"

/**
 * Class representing a generic finite element.  Its constructor
 * takes a specific interpolation basis and quadrature rule in
 * order to define a complete element.
 *
 * The class assume that the
 * mapping from parent coordinates to real coordinates is
 * iso-parametric, and therefore the same basis is used for both
 * interpolation and mapping.
 *
 * The class also defines a generic interface for accessing finite
 * element data.  In the future, more sophisticated element
 * definitions that do not fit within the current class can be
 * defined through derived classes.
 */
namespace geosx
{
template< int dim >
class FiniteElement : public FiniteElementBase
{
public:

  /**
   * Constructor.  Takes an interpolation basis and quadrature rule,
   * and pre-computes all static finite element data.  Any data
   * that depends on the mapped configuration of the element, however,
   * is left uninitialized until reinit(...) is called.
   */
  FiniteElement( BasisBase const & basis,
                 QuadratureBase const & quadrature,
                 const int num_zero_energy_modes = 0 ):
    FiniteElementBase( dim, quadrature.size(), basis.size(), num_zero_energy_modes )
  {
    for( int q=0; q<n_q_points; ++q )
    {
      data[q].parent_q_point = quadrature.integration_point( q );
      data[q].parent_q_weight = quadrature.integration_weight( q );

      for( int i=0; i<n_dofs; ++i )
      {
        data[q].parent_values[i] = basis.value( i, data[q].parent_q_point );
        R1Tensor const gradient = basis.gradient( i, data[q].parent_q_point );
        for( int j = 0; j < 3; ++j )
        {
          data[q].parent_gradients[i][j] = gradient[j];
        }
      }
    }
  }

  ~FiniteElement() override
  {}

  static string CatalogName() { return "C3D8"; }

  void reinit( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
               arraySlice1d< localIndex const, -1 > const & mapped_support_points ) override
  { return reinitPrivate( X, mapped_support_points ); }

  void reinit( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
               arraySlice1d< localIndex const, 0 > const & mapped_support_points ) override
  { return reinitPrivate( X, mapped_support_points ); }

private:
  /**
   * Reinitialize the finite element basis on a particular element.
   * We use the coordinates of the support points in real space to
   * construct the forward mapping from the parent coordinate system.  The
   * support points are assumed to follow a lexicographic ordering:
   * On the parent element, we loop over the x-coordinate fastest,
   * the y, then z (depending on the desired spatial dimension of the
   * element).
   */
  template< int USD >
  void reinitPrivate( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                      arraySlice1d< localIndex const, USD > const & mapped_support_points )
  {
    real64 jacobian[ 3 ][ 3 ];

    for( int q=0; q<n_q_points; ++q )
    {
      arrayView2d< real64 const > const & parentGradients = data[ q ].parent_gradients;

      LvArray::tensorOps::AiBj< 3, 3 >( jacobian, X[ mapped_support_points[ 0 ] ], parentGradients[ 0 ] );

      for( int i=1; i<n_dofs; ++i )
      {
        LvArray::tensorOps::plusAiBj< 3, 3 >( jacobian, X[ mapped_support_points[ i ] ], parentGradients[ i ] );
      }

      if( dim == 2 )
      {
        jacobian[ 2 ][ 2 ] = 1;
      }

      data[ q ].jacobian_determinant = LvArray::tensorOps::invert< 3 >( jacobian );

      for( int i=0; i<n_dofs; ++i )
      {
        LvArray::tensorOps::AjiBj< 3, 3 >( data[ q ].mapped_gradients[ i ], jacobian, parentGradients[ i ] );
      }
    }
  }
};

template< typename KERNELWRAPPER, typename ... PARAMS >
inline real64
finiteElementLaunchDispatch( localIndex NUM_NODES_PER_ELEM,
                             localIndex NUM_QUADRATURE_POINTS,
                             PARAMS && ... params )
{
  if( NUM_NODES_PER_ELEM == 8 && NUM_QUADRATURE_POINTS == 8 )
  {
    return KERNELWRAPPER::template Launch< 8, 8 >( std::forward< PARAMS >( params )... );
  }

  GEOSX_ERROR( "Not implemented!" );

  return 0;
}

}

#endif
