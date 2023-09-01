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
                       nodes::REFERENCE_POSITION_USD > const X,
          arrayView1d< localIndex const > fluid_indices,
          arrayView2d< localIndex const > const faceToRegion,
          arrayView2d< localIndex const > const faceToElement,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView2d< real64 const > const faceNormals,
          arrayView1d< real32 > const couplingVectorx,
          arrayView1d< real32 > const couplingVectory,
          arrayView1d< real32 > const couplingVectorz )
  {
    array1d< localIndex > count( 3 );
    count.zero();
    arrayView1d< localIndex > const count_view = count.toView();

    bool const dump = helpers::ienv( "DUMP" ) > 0;

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const f )
    {
      localIndex e0 = faceToElement( f, 0 ), e1 = faceToElement( f, 1 );
      localIndex er0 = faceToRegion( f, 0 ), er1 = faceToRegion( f, 1 );
      // localIndex esr0 = faceToSubRegion( f, 0 ), esr1 = faceToSubRegion( f, 1 );

      if((e0 != -1 && e1 == -1) || (e0 == -1 && e1 != -1))
      {
        // printf("\t[CouplingKernel::launch] debug\n");
        RAJA::atomicInc< ATOMIC_POLICY >( &count_view[0] );
      }

      if( e0 != -1 && e1 != -1 )
      {
        RAJA::atomicInc< ATOMIC_POLICY >( &count_view[1] );
        // printf("\t[CouplingKernel::launch] f=%i -> (e0=%i, e1=%i)\n", f, e0, e1);
        // NOTE: subregion check doesn't work: esr0 != esr1
        if( er0 != er1 )  // should define an interface
        {
          // determine normal sign for fluid -> solid coupling
          localIndex sgn = -1;
          for( auto const & idx : fluid_indices )
          {
            if( er0 == idx )
            {
              sgn = +1;
              break;
            }
          }
          RAJA::atomicInc< ATOMIC_POLICY >( &count_view[2] );
          // printf("\t[CouplingKernel::launch] interface found for f=%i\n", f);
          real64 xLocal[ numNodesPerFace ][ 3 ];
          for( localIndex a = 0; a < numNodesPerFace; ++a )
          {
            for( localIndex i = 0; i < 3; ++i )
            {
              xLocal[a][i] = X( facesToNodes( f, a ), i );
            }
          }

          for( localIndex q = 0; q < numNodesPerFace; ++q )
          {
            real64 const aux = -FE_TYPE::computeDampingTerm( q, xLocal );

            // real64 const mag = sqrt(pow(faceNormals[f][0], 2) + pow(faceNormals[f][1], 2) + pow(faceNormals[f][2], 2));

            real32 const localIncrementx = aux * (sgn * faceNormals[f][0]);
            real32 const localIncrementy = aux * (sgn * faceNormals[f][1]);
            real32 const localIncrementz = aux * (sgn * faceNormals[f][2]);

            if( dump && q == 0 )
              printf(
                "\t[CouplingKernel::launch] nx=%g ny=%g nz=%g\n",
                faceNormals[f][0], faceNormals[f][1], faceNormals[f][2]
                );
            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVectorx[facesToNodes[f][q]], localIncrementx );
            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVectory[facesToNodes[f][q]], localIncrementy );
            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVectorz[facesToNodes[f][q]], localIncrementz );
          }
        }
      }
    } );

    printf(
      "\t[CouplingKernel::launch] n_faces=%i n_boundary_faces=%i n_internal_faces=%i n_interface_faces=%i\n",
      size, count[0], count[1], count[2]
      );


  }
};

} /* namespace acousticElasticWaveEquationSEMKernels */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEMKERNEL_HPP_ */
