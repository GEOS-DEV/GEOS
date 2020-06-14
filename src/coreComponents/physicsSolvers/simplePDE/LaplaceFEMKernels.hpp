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
 * @file SolidMechanicsLagrangianFEMKernels.hpp
 */

#pragma once

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "finiteElement/FiniteElementShapeFunctionKernel.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
namespace LaplaceFEMKernels
{

struct ImplicitKernel
{

  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
  static real64 Launch( arrayView4d< real64 const > const & dNdX,
                        arrayView2d< real64 const > const & detJ,
                        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemNodes,
                        arrayView1d< globalIndex const > const & dofIndex,
                        globalIndex const dofRankOffset,
                        CRSMatrixView< real64, globalIndex const > const & matrix )
  {
    localIndex const numElems = dNdX.size( 0 );
    GEOSX_ERROR_IF_NE( dNdX.size( 0 ), numElems );
    GEOSX_ERROR_IF_NE( dNdX.size( 1 ), NUM_QUADRATURE_POINTS );
    GEOSX_ERROR_IF_NE( dNdX.size( 2 ), NUM_NODES_PER_ELEM );

    GEOSX_ERROR_IF_NE( detJ.size( 0 ), numElems );
    GEOSX_ERROR_IF_NE( detJ.size( 1 ), NUM_QUADRATURE_POINTS );

    GEOSX_ERROR_IF_NE( elemNodes.size( 0 ), numElems );
    GEOSX_ERROR_IF_NE( elemNodes.size( 1 ), NUM_NODES_PER_ELEM );

    // begin element loop, skipping ghost elements
    forAll< parallelDevicePolicy< 32 > >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      globalIndex dofIndices[ NUM_NODES_PER_ELEM ];
      real64 elementMatrix[ NUM_NODES_PER_ELEM ][ NUM_NODES_PER_ELEM ] = { { 0 } };

      for( localIndex q = 0; q < NUM_QUADRATURE_POINTS; ++q )
      {
        for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
        {
          dofIndices[ a ] = dofIndex[ elemNodes( k, a ) ];

          real64 diffusion = 1.0;
          for( localIndex b = 0; b < NUM_NODES_PER_ELEM; ++b )
          {
            elementMatrix[ a ][ b ] += detJ( k, q ) *
                                       diffusion *
                                       + LvArray::tensorOps::AiBi< 3 >( dNdX[ k ][ q ][ a ], dNdX[ k ][ q ][ b ] );
          }

        }
      }

      for( localIndex i = 0; i < NUM_NODES_PER_ELEM; ++i )
      {
        globalIndex const dof = dofIndices[ i ] - dofRankOffset;
        if( dof < 0 || dof >= matrix.numRows() ) continue;
        matrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                     dofIndices,
                                                                     elementMatrix[ i ],
                                                                     NUM_NODES_PER_ELEM );
      }
    } );

    return 0;
  }
};

} // namespace LaplaceFEMKernels
} // namespace geosx
