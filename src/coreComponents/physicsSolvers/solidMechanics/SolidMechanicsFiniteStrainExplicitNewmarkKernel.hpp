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
 * @file SolidMechanicsSmallStrainExplicitNewmarkKernels.hpp
 */

#pragma once

#include "SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp"
#include "finiteElement/Kinematics.h"


namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

class ExplicitFiniteStrain : public ExplicitSmallStrain
{
#if defined(GEOSX_USE_CUDA)
  #define CALCFEMSHAPE
#endif
  // If UPDATE_STRESS is undef, then stress is not updated at all.
//  #define UPDATE_STRESS 1 // uses total displacement to and adds material stress state to integral for nodalforces.
#define UPDATE_STRESS 2 // uses velocity*dt and updates material stress state.

public:
  using KernelBase = ExplicitSmallStrain;
//  using Base::numTestDofPerSP;
//  using Base::numTrialDofPerSP;

  template< int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
            int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >
  struct StackVariables : public KernelBase::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                                                             NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >
  {
public:
    using StackVariablesBase = KernelBase::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM,
                                                           NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >;

    using StackVariablesBase::numNodes;
    using StackVariablesBase::fLocal;
    using StackVariablesBase::varLocal;

#if defined(CALCFEMSHAPE)
    using StackVariablesBase::xLocal;
    using StackVariablesBase::dNdX;
    using StackVariablesBase::detJ;
#endif


    GEOSX_HOST_DEVICE
    StackVariables():
      StackVariablesBase(),
      uLocal{ {0.0} }
    {}

    real64 uLocal[ numNodes ][ 3 ];
  };



  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  class Components : public KernelBase::Components< SUBREGION_TYPE,
                                                    CONSTITUTIVE_TYPE,
                                                    NUM_NODES_PER_ELEM,
                                                    NUM_NODES_PER_ELEM >
  {
public:
    using ComponentsBase = KernelBase::Components< SUBREGION_TYPE,
                                                   CONSTITUTIVE_TYPE,
                                                   NUM_NODES_PER_ELEM,
                                                   NUM_NODES_PER_ELEM >;
    using ComponentsBase::numNodesPerElem;
    using ComponentsBase::m_dt;
    using ComponentsBase::elemsToNodes;
    using ComponentsBase::u;
    using ComponentsBase::vel;
    using ComponentsBase::X;
    using ComponentsBase::acc;
    using ComponentsBase::constitutiveUpdate;
    using ComponentsBase::dNdX;
    using ComponentsBase::detJ;

    using StackVars = StackVariables< numNodesPerElem,
                                      numNodesPerElem >;

    Components( arrayView1d< globalIndex const > const & inputDofNumber,
                ParallelMatrix & inputMatrix,
                ParallelVector & inputRhs,
                NodeManager & nodeManager,
                SUBREGION_TYPE const & elementSubRegion,
                FiniteElementBase const * const finiteElementSpace,
                CONSTITUTIVE_TYPE * const inputConstitutiveType,
                real64 const dt,
                string const & elementListName ):
      ComponentsBase( inputDofNumber,
                      inputMatrix,
                      inputRhs,
                      nodeManager,
                      elementSubRegion,
                      finiteElementSpace,
                      inputConstitutiveType,
                      dt,
                      elementListName )
    {}

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void setup( localIndex const k,
                STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a< NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        for( int i=0; i<3; ++i )
        {
#if defined(CALCFEMSHAPE)
          stack.xLocal[ a ][ i ] = X[ nodeIndex ][ i ];
#endif
          stack.uLocal[ a ][ i ] = u[ nodeIndex ][ i ];
          stack.varLocal[ a ][ i ] = vel[ nodeIndex ][ i ];
        }
      }
    }

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void quadraturePointStateUpdate( localIndex const k,
                                     localIndex const q,
                                     STACK_VARIABLE_TYPE & stack ) const
    {
#if defined(CALCFEMSHAPE)
      real64 dNdX[ 8 ][ 3 ];
      real64 const detJ = FiniteElementShapeKernel::shapeFunctionDerivatives( q, X_local, dNdX );
#define DNDX dNdX
#define DETJ detJ
#else
#define DNDX dNdX[k][q]
#define DETJ detJ( k, q )
#endif
      R2Tensor dUhatdX, dUdX;
      real64 * const GEOSX_RESTRICT g0 = dUhatdX.Data();
      real64 * const GEOSX_RESTRICT g1 = dUdX.Data();

      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        g0[0] += stack.varLocal[a][0]*DNDX[a][0];
        g0[1] += stack.varLocal[a][0]*DNDX[a][1];
        g0[2] += stack.varLocal[a][0]*DNDX[a][2];
        g0[3] += stack.varLocal[a][1]*DNDX[a][0];
        g0[4] += stack.varLocal[a][1]*DNDX[a][1];
        g0[5] += stack.varLocal[a][1]*DNDX[a][2];
        g0[6] += stack.varLocal[a][2]*DNDX[a][0];
        g0[7] += stack.varLocal[a][2]*DNDX[a][1];
        g0[8] += stack.varLocal[a][2]*DNDX[a][2];


        g1[0] += stack.uLocal[a][0]*DNDX[a][0];
        g1[1] += stack.uLocal[a][0]*DNDX[a][1];
        g1[2] += stack.uLocal[a][0]*DNDX[a][2];
        g1[3] += stack.uLocal[a][1]*DNDX[a][0];
        g1[4] += stack.uLocal[a][1]*DNDX[a][1];
        g1[5] += stack.uLocal[a][1]*DNDX[a][2];
        g1[6] += stack.uLocal[a][2]*DNDX[a][0];
        g1[7] += stack.uLocal[a][2]*DNDX[a][1];
        g1[8] += stack.uLocal[a][2]*DNDX[a][2];
      }


      dUhatdX *= m_dt;

      R2Tensor F, Ldt, fInv;

      // calculate du/dX
      F = dUhatdX;
      F *= 0.5;
      F += dUdX;
      F.PlusIdentity( 1.0 );
      fInv.Inverse( F );

      // chain rule: calculate dv/du = dv/dX * dX/du
      Ldt.AijBjk( dUhatdX, fInv );

      // calculate gradient (end of step)
      F = dUhatdX;
      F += dUdX;
      F.PlusIdentity( 1.0 );
      real64 const detF = F.Det();
      fInv.Inverse( F );


      R2Tensor Rot;
      R2SymTensor Dadt;
      HughesWinget( Rot, Dadt, Ldt );

      constitutiveUpdate.HypoElastic( k, q, Dadt.Data(), Rot );

      real64 const integrationFactor = -DETJ * detF;

      real64 const * const stress = constitutiveUpdate.m_stress[k][q];

      real64 P[ 3 ][ 3 ];
      P[ 0 ][ 0 ] = ( stress[ 0 ] * fInv( 0, 0 ) + stress[ 5 ] * fInv( 0, 1 ) + stress[ 4 ] * fInv( 0, 2 ) ) * integrationFactor;
      P[ 0 ][ 1 ] = ( stress[ 0 ] * fInv( 1, 0 ) + stress[ 5 ] * fInv( 1, 1 ) + stress[ 4 ] * fInv( 1, 2 ) ) * integrationFactor;
      P[ 0 ][ 2 ] = ( stress[ 0 ] * fInv( 2, 0 ) + stress[ 5 ] * fInv( 2, 1 ) + stress[ 4 ] * fInv( 2, 2 ) ) * integrationFactor;

      P[ 1 ][ 0 ] = ( stress[ 5 ] * fInv( 0, 0 ) + stress[ 1 ] * fInv( 0, 1 ) + stress[ 3 ] * fInv( 0, 2 ) ) * integrationFactor;
      P[ 1 ][ 1 ] = ( stress[ 5 ] * fInv( 1, 0 ) + stress[ 1 ] * fInv( 1, 1 ) + stress[ 3 ] * fInv( 1, 2 ) ) * integrationFactor;
      P[ 1 ][ 2 ] = ( stress[ 5 ] * fInv( 2, 0 ) + stress[ 1 ] * fInv( 2, 1 ) + stress[ 3 ] * fInv( 2, 2 ) ) * integrationFactor;

      P[ 2 ][ 0 ] = ( stress[ 4 ] * fInv( 0, 0 ) + stress[ 3 ] * fInv( 0, 1 ) + stress[ 2 ] * fInv( 0, 2 ) ) * integrationFactor;
      P[ 2 ][ 1 ] = ( stress[ 4 ] * fInv( 1, 0 ) + stress[ 3 ] * fInv( 1, 1 ) + stress[ 2 ] * fInv( 1, 2 ) ) * integrationFactor;
      P[ 2 ][ 2 ] = ( stress[ 4 ] * fInv( 2, 0 ) + stress[ 3 ] * fInv( 2, 1 ) + stress[ 2 ] * fInv( 2, 2 ) ) * integrationFactor;

      for( int a=0; a<NUM_NODES_PER_ELEM; ++a )      // loop through all shape functions in element
      {
        stack.fLocal[a][0] = stack.fLocal[a][0] + P[ 0 ][ 0 ] * DNDX[ a ][ 0 ] + P[ 0 ][ 1 ] * DNDX[ a ][ 1 ] + P[ 0 ][ 2 ] * DNDX[ a ][ 2 ];
        stack.fLocal[a][1] = stack.fLocal[a][1] + P[ 1 ][ 0 ] * DNDX[ a ][ 0 ] + P[ 1 ][ 1 ] * DNDX[ a ][ 1 ] + P[ 1 ][ 2 ] * DNDX[ a ][ 2 ];
        stack.fLocal[a][2] = stack.fLocal[a][2] + P[ 2 ][ 0 ] * DNDX[ a ][ 0 ] + P[ 2 ][ 1 ] * DNDX[ a ][ 1 ] + P[ 2 ][ 2 ] * DNDX[ a ][ 2 ];
      }

    }

    using ComponentsBase::quadraturePointJacobianContribution;
    using ComponentsBase::quadraturePointResidualContribution;
    using ComponentsBase::complete;
    using ComponentsBase::Launch;

  };
#undef CALCFEMSHAPE
#undef DNDX
#undef DETJ
#undef UPDATE_STRESS

};

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
