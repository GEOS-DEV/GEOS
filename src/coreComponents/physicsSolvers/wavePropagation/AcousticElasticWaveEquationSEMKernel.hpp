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
 * @file AcousticElasticWaveEquationSEMKernel.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEMKERNEL_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#if !defined( GEOS_USE_HIP )
#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#endif

#include <utility>

namespace geos
{

namespace acousticElasticWaveEquationSEMKernels
{

template< typename FE_TYPE >
struct CouplingKernel
{
  static constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;

  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const,
                       nodes::REFERENCE_POSITION_USD > const nodeCoords,
          localIndex const regionIndex,
          localIndex const subRegionIndex,
          arrayView2d< localIndex const > const faceToSubRegion,
          arrayView2d< localIndex const > const faceToRegion,
          arrayView2d< localIndex const > const faceToElement,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView2d< real64 const > const faceNormals,
          arrayView1d< real32 > const couplingVectorx,
          arrayView1d< real32 > const couplingVectory,
          arrayView1d< real32 > const couplingVectorz )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const f )
    {
      localIndex e0 = faceToElement( f, 0 ), e1 = faceToElement( f, 1 );
      localIndex er0 = faceToRegion( f, 0 ), er1 = faceToRegion( f, 1 );
      localIndex esr0 = faceToSubRegion( f, 0 ), esr1 = faceToSubRegion( f, 1 );

      if( e0 != -1 && e1 != -1 && er0 != er1 )  // an interface is defined as a transition between regions
      {
        // check that one of the region is the fluid subregion for the fluid -> solid coupling term
        if((er0 == regionIndex && esr0 == subRegionIndex) || (er1 == regionIndex && esr1 == subRegionIndex))
        {
          real64 xLocal[ 4 ][ 3 ];
          for( localIndex a = 0; a < 4; ++a )
          {
            localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
            for( localIndex i = 0; i < 3; ++i )
            {
              xLocal[a][i] = nodeCoords( nodeIndex, i );
            }
          }

          // determine normal sign for fluid -> solid coupling
          localIndex sgn = er0 == regionIndex ? 1 : (er1 == regionIndex ? -1 : 0);

          for( localIndex q = 0; q < numNodesPerFace; ++q )
          {
            real64 const aux = FE_TYPE::computeDampingTerm( q, xLocal );

            real32 const localIncrementx = aux * (sgn * faceNormals( f, 0 ));
            real32 const localIncrementy = aux * (sgn * faceNormals( f, 1 ));
            real32 const localIncrementz = aux * (sgn * faceNormals( f, 2 ));

            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVectorx[facesToNodes( f, q )], localIncrementx );
            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVectory[facesToNodes( f, q )], localIncrementy );
            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVectorz[facesToNodes( f, q )], localIncrementz );
          }
        }
      }
    } );

  }
};

} /* namespace acousticElasticWaveEquationSEMKernels */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEMKERNEL_HPP_ */
