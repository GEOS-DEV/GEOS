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

#if defined(GEOSX_USE_CUDA)
  #define CALCFEMSHAPE
#endif
// If UPDATE_STRESS is undef, then stress is not updated at all.
//  #define UPDATE_STRESS 1 // uses total displacement to and adds material stress state to integral for nodalforces.
#define UPDATE_STRESS 2 // uses velocity*dt and updates material stress state.

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          int NUM_NODES_PER_ELEM,
          int >
class ExplicitFiniteStrain : public ExplicitSmallStrain< SUBREGION_TYPE,
                                                         CONSTITUTIVE_TYPE,
                                                         NUM_NODES_PER_ELEM,
                                                         NUM_NODES_PER_ELEM >
{
public:
  using Base = ExplicitSmallStrain< SUBREGION_TYPE,
                                    CONSTITUTIVE_TYPE,
                                    NUM_NODES_PER_ELEM,
                                    NUM_NODES_PER_ELEM >;
  using Base::numTestDofPerSP;
  using Base::numTrialDofPerSP;
  using Base::numNodesPerElem;

  using Base::m_dt;
  using Base::elemsToNodes;
  using Base::u;
  using Base::vel;
  using Base::X;
  using Base::acc;
  using Base::constitutiveUpdate;
#if !defined(CALCFEMSHAPE)
  using Base::dNdX;
  using Base::detJ;
#endif
  using ConstructorParams = ExplicitSmallStrainConstructorParams;

  ExplicitFiniteStrain( NodeManager & nodeManager,
                        EdgeManager const & edgeManager,
                        FaceManager const & faceManager,
                        SUBREGION_TYPE const & elementSubRegion,
                        FiniteElementBase const * const finiteElementSpace,
                        CONSTITUTIVE_TYPE * const inputConstitutiveType,
                        real64 const dt,
                        string const & elementListName ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          dt,
          elementListName )
  {}

  ExplicitFiniteStrain( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       SUBREGION_TYPE const & elementSubRegion,
                       FiniteElementBase const * const finiteElementSpace,
                       CONSTITUTIVE_TYPE * const inputConstitutiveType,
                       ConstructorParams & params ):
    ExplicitFiniteStrain( nodeManager,
                         edgeManager,
                         faceManager,
                         elementSubRegion,
                         finiteElementSpace,
                         inputConstitutiveType,
                         params.m_dt,
                         params.m_elementListName )
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
    real64 const detJ = FiniteElementShapeKernel::shapeFunctionDerivatives( q, stack.xLocal, dNdX );
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

    auto const & stress = constitutiveUpdate.m_stress[k][q];

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

  using Base::quadraturePointJacobianContribution;
  using Base::quadraturePointResidualContribution;
  using Base::complete;
  using Base::Launch;


  //*****************************************************************************
    struct StackVariables : public Base::StackVariables
    {
  public:
      using Base::StackVariables::fLocal;
      using Base::StackVariables::varLocal;

  #if defined(CALCFEMSHAPE)
      using Base::StackVariables::xLocal;
      using Base::StackVariables::dNdX;
      using Base::StackVariables::detJ;
  #endif


      GEOSX_HOST_DEVICE
      StackVariables():
        Base::StackVariables(),
        uLocal{ {0.0} }
      {}

      real64 uLocal[ numNodesPerElem ][ numTrialDofPerSP ];
    };
  //*****************************************************************************





};
#undef CALCFEMSHAPE
#undef DNDX
#undef DETJ
#undef UPDATE_STRESS


} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
