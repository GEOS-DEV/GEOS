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

  CouplingKernel( FE_TYPE const & finiteElement ) : m_finiteElement( finiteElement )
  {}

  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const > const faceToRegion,
          arrayView2d< localIndex const > const faceToSubRegion,
          arrayView2d< localIndex const > const faceToElement,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView1d< real32 > const couplingVector )
  {
    real32 rhof = 1020;  // hardcoded until github.com/GEOS-DEV/GEOS/pull/2548 is merged

    array1d< localIndex > count( 3 );
    count.zero();
    arrayView1d< localIndex > const count_view = count.toView();

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const f )
    {
      localIndex e0 = faceToElement( f, 0 ), e1 = faceToElement( f, 1 );
      localIndex er0 = faceToRegion( f, 0 ), er1 = faceToRegion( f, 1 );
      localIndex esr0 = faceToSubRegion( f, 0 ), esr1 = faceToSubRegion( f, 1 );

      if ((e0 != -1 && e1 == -1) || (e0 == -1 && e1 != -1))
      {
        // printf( "\t[CouplingKernel::launch] debug\n" );
        RAJA::atomicInc< ATOMIC_POLICY >( &count_view[0] );
      }

      if (e0 != -1 && e1 != -1)
      {
        RAJA::atomicInc< ATOMIC_POLICY >( &count_view[1] );
        printf( "\t[CouplingKernel::launch] f=%i -> (e0=%i, e1=%i)\n", f, e0, e1 );
        if (er0 != er1 && esr0 != esr1)  /* should define an interface */
        {
          RAJA::atomicInc< ATOMIC_POLICY >( &count_view[2] );
          printf( "\t[CouplingKernel::launch] interface found for f=%i\n", f );
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
            real32 const localIncrement = -rhof * m_finiteElement.computeDampingTerm( q, xLocal );
            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVector[facesToNodes[f][q]], localIncrement );
          }
        }
      }
    } );

    printf(
      "\t[CouplingKernel::launch] n_faces=%i n_boundary_faces=%i n_internal_faces=%i n_interface_faces=%i\n",
      size, count_view[0], count_view[1], count_view[2]
    );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

} /* namespace acousticElasticWaveEquationSEMKernels */

namespace acousticElasticWaveEquationSEMKernels_old
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class CouplingSEMKernel : public finiteElement::KernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, 1, 1 >
{
public:
  using Base = finiteElement::KernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, 1, 1 >;

  static constexpr localIndex numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  static constexpr localIndex num1dNodes = FE_TYPE::num1dNodes;

  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables(): xLocal() {}
    real64 xLocal[ numNodesPerElem ][ 3 ];
  };

  CouplingSEMKernel( NodeManager & nodeManager,
                     EdgeManager const & edgeManager,
                     FaceManager const & faceManager,
                     localIndex const targetRegionIndex,
                     SUBREGION_TYPE const & elementSubRegion,
                     FE_TYPE const & finiteElementSpace,
                     CONSTITUTIVE_TYPE & inputConstitutiveType,
                     SortedArrayView< localIndex const > const & interfaceNodesSet,
                     real64 const dt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.getField< fields::referencePosition32 >() ),
    m_p_n( nodeManager.getField< fields::Pressure_n >() ),
    m_ux_nm1( nodeManager.getField< fields::Displacementx_nm1 >() ),
    m_uy_nm1( nodeManager.getField< fields::Displacementy_nm1 >() ),
    m_uz_nm1( nodeManager.getField< fields::Displacementz_nm1 >() ),
    m_ux_n( nodeManager.getField< fields::Displacementx_n >() ),
    m_uy_n( nodeManager.getField< fields::Displacementy_n >() ),
    m_uz_n( nodeManager.getField< fields::Displacementz_n >() ),
    m_p_np1( nodeManager.getField< fields::Pressure_np1 >() ),
    m_ux_np1( nodeManager.getField< fields::Displacementx_np1 >() ),
    m_uy_np1( nodeManager.getField< fields::Displacementy_np1 >() ),
    m_uz_np1( nodeManager.getField< fields::Displacementz_np1 >() ),
    m_interfaceNodesSet( interfaceNodesSet ),
    m_dt( dt )
  {
    GEOS_UNUSED_VAR( edgeManager );
    GEOS_UNUSED_VAR( faceManager );
    GEOS_UNUSED_VAR( targetRegionIndex );
  }

  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    /// numDofPerTrialSupportPoint = 1
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = Base::m_elemsToNodes( k, a );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X[ nodeIndex ][ i ];
      }
    }
  }

protected:
  arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const m_X;

  arrayView1d< real32 const > const m_p_n;
  arrayView1d< real32 const > const m_ux_nm1;
  arrayView1d< real32 const > const m_uy_nm1;
  arrayView1d< real32 const > const m_uz_nm1;
  arrayView1d< real32 const > const m_ux_n;
  arrayView1d< real32 const > const m_uy_n;
  arrayView1d< real32 const > const m_uz_n;

  // writing coupling terms in these arrays
  arrayView1d< real32 > const m_p_np1;
  arrayView1d< real32 > const m_ux_np1;
  arrayView1d< real32 > const m_uy_np1;
  arrayView1d< real32 > const m_uz_np1;

  SortedArrayView< localIndex const > const & m_interfaceNodesSet;

  real64 const m_dt;
};

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class AcousticElasticSEM : public CouplingSEMKernel< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >
{
public:
  using Base = CouplingSEMKernel< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >;
  using Base::Base::m_elemsToNodes;
  using Base::Base::m_finiteElementSpace;
  using typename Base::StackVariables;

  template< typename ... Args >
  AcousticElasticSEM( Args &&... args ): Base( std::forward< Args >( args )... ) {}

  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // printf( "\t[AcousticElasticSEM::quadraturePointKernel]\n" );
    // TODO: restrict to interface nodeset
    real64 const val = m_finiteElementSpace.computeMassTerm( q, stack.xLocal );
    real32 const rhof = 1020.0;  // hardcoded until github.com/GEOS-DEV/GEOS/pull/2548 is merged

    for( int i = 0; i < Base::num1dNodes; ++i )
    {
      for( int j = 0; j < Base::num1dNodes; ++j )
      {
        localIndex const ai = m_elemsToNodes[k][i];
        localIndex const aj = m_elemsToNodes[k][j];

        real32 const localIncrementx = val * (this->m_p_n[aj] / rhof);
        real32 const localIncrementy = val * (this->m_p_n[aj] / rhof);
        real32 const localIncrementz = val * (this->m_p_n[aj] / rhof);

        RAJA::atomicAdd< parallelDeviceAtomic >( &this->m_ux_np1[ai], localIncrementx );
        RAJA::atomicAdd< parallelDeviceAtomic >( &this->m_uy_np1[ai], localIncrementy );
        RAJA::atomicAdd< parallelDeviceAtomic >( &this->m_ux_np1[ai], localIncrementz );
      }
    }
  }

};

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ElasticAcousticSEM : public CouplingSEMKernel< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >
{
public:
  using Base = CouplingSEMKernel< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >;
  using Base::Base::m_elemsToNodes;
  using Base::Base::m_finiteElementSpace;
  using typename Base::StackVariables;

  template< typename ... Args >
  ElasticAcousticSEM( Args &&... args ): Base( std::forward< Args >( args )... ) {}

  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // printf( "\t[ElasticAcousticSEM::quadraturePointKernel]\n" );
    // TODO: restrict to interface nodeset

    real64 const val = m_finiteElementSpace.computeMassTerm( q, stack.xLocal );
    real32 const rhof = 1020.0;  // hardcoded until github.com/GEOS-DEV/GEOS/pull/2548 is merged
    real32 const dt2 = pow( this->m_dt, 2 );

    for( int i = 0; i < Base::num1dNodes; ++i )
    {
      for( int j = 0; j < Base::num1dNodes; ++j )
      {
        localIndex const ai = m_elemsToNodes[k][i];
        localIndex const aj = m_elemsToNodes[k][j];

        real32 const ux_dt2 = ( this->m_ux_np1[aj] - 2.0 * this->m_ux_n[aj] + this->m_ux_nm1[aj] ) / dt2;
        real32 const uy_dt2 = ( this->m_uy_np1[aj] - 2.0 * this->m_uy_n[aj] + this->m_uy_nm1[aj] ) / dt2;
        real32 const uz_dt2 = ( this->m_uz_np1[aj] - 2.0 * this->m_uz_n[aj] + this->m_uz_nm1[aj] ) / dt2;

        real32 const localIncrement = val * rhof * ( ux_dt2 + uy_dt2 + uz_dt2 );
        RAJA::atomicAdd< parallelDeviceAtomic >( &this->m_p_np1[ai], localIncrement );
      }
    }
  }

};

using AcousticElasticSEMFactory = finiteElement::KernelFactory< AcousticElasticSEM, SortedArrayView< localIndex const > const &, real64 >;
using ElasticAcousticSEMFactory = finiteElement::KernelFactory< ElasticAcousticSEM, SortedArrayView< localIndex const > const &, real64 >;

} // namespace acousticElasticWaveEquationSEMKernels

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEMKERNEL_HPP_ */
