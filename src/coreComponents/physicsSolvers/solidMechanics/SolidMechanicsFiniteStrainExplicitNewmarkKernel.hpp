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
 * @file SolidMechanicsFiniteStrainExplicitNewmarkKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSFINITESTRAINEXPLICITNEWMARK_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSFINITESTRAINEXPLICITNEWMARK_HPP_

#include "SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp"
#include "finiteElement/Kinematics.h"

namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

/// If UPDATE_STRESS is undef, uses total displacement and stress is not
/// updated at all.
/// If UPDATE_STRESS 1, uses total displacement to and adds material stress
/// state to integral for nodalforces.
/// If UPDATE_STRESS 2 then velocity*dt is used to update material stress state
#define UPDATE_STRESS 2


/**
 * @brief Implements kernels for solving the equations of motion using the
 *   explicit Newmark method under the finite strain assumption.
 * @copydoc ExplicitSmallStrain
 *
 * ### Explicit Small Strain Description
 * Finite strain implementation.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitFiniteStrain : public ExplicitSmallStrain< SUBREGION_TYPE,
                                                         CONSTITUTIVE_TYPE,
                                                         FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = ExplicitSmallStrain< SUBREGION_TYPE,
                                    CONSTITUTIVE_TYPE,
                                    FE_TYPE >;

  using Base::numNodesPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

  using Base::m_dt;
  using Base::m_u;
  using Base::m_vel;
  using Base::m_acc;
#if !defined(CALCFEMSHAPE)
  using Base::m_X;
#endif

  /**
   * @copydoc ExplicitSmallStrain
   */
  ExplicitFiniteStrain( NodeManager & nodeManager,
                        EdgeManager const & edgeManager,
                        FaceManager const & faceManager,
                        localIndex const targetRegionIndex,
                        SUBREGION_TYPE const & elementSubRegion,
                        FE_TYPE const & finiteElementSpace,
                        CONSTITUTIVE_TYPE & inputConstitutiveType,
                        real64 const dt,
                        string const & elementListName ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          dt,
          elementListName )
  {}


  //*****************************************************************************
  /**
   * @copydoc ExplicitSmallStrain::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
    using Base::StackVariables::fLocal;
    using Base::StackVariables::varLocal;
    using Base::StackVariables::xLocal;


    /**
     * @brief constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            uLocal{ {0.0} }
    {}

    /// Local stack storage for nodal displacements.
    real64 uLocal[ numNodesPerElem ][ numDofPerTrialSupportPoint ];
  };
  //*****************************************************************************


  /**
   * @copydoc ExplicitSmallStrain::setup
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int i=0; i<numDofPerTrialSupportPoint; ++i )
      {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
        stack.xLocal[ a ][ i ] = m_X[ nodeIndex ][ i ];
#endif
        stack.uLocal[ a ][ i ] = m_u[ nodeIndex ][ i ];
        stack.varLocal[ a ][ i ] = m_vel[ nodeIndex ][ i ];
      }
    }
  }

  /**
   * @copydoc ExplicitSmallStrain::quadraturePointKernel
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    real64 dUhatdX[3][3] = { {0} };
    real64 dUdX[3][3] = { {0} };
    real64 F[3][3] = { {0} };
    real64 Ldt[3][3] = { {0} };
    real64 fInv[3][3] = { {0} };

    FE_TYPE::gradient( dNdX, stack.varLocal, dUhatdX );
    FE_TYPE::gradient( dNdX, stack.uLocal, dUdX );

    LvArray::tensorOps::scale< 3, 3 >( dUhatdX, m_dt );

    // calculate du/dX
    LvArray::tensorOps::scaledCopy< 3, 3 >( F, dUhatdX, 0.5 );
    LvArray::tensorOps::add< 3, 3 >( F, dUdX );
    LvArray::tensorOps::addIdentity< 3 >( F, 1.0 );
    LvArray::tensorOps::invert< 3 >( fInv, F );

    // chain rule: calculate dv/dx^(n+1/2) = dv/dX * dX/dx^(n+1/2)
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( Ldt, dUhatdX, fInv );

    // calculate gradient (end of step)
    LvArray::tensorOps::copy< 3, 3 >( F, dUhatdX );
    LvArray::tensorOps::add< 3, 3 >( F, dUdX );
    LvArray::tensorOps::addIdentity< 3 >( F, 1.0 );
    real64 const detF = LvArray::tensorOps::invert< 3 >( fInv, F );

    real64 Rot[ 3 ][ 3 ];
    real64 Dadt[ 6 ];
    HughesWinget( Rot, Dadt, Ldt );

    real64 stress[ 6 ] = { };
    m_constitutiveUpdate.hypoUpdate_StressOnly( k, q, Dadt, Rot, stress );

    real64 P[ 3 ][ 3 ];
    LvArray::tensorOps::Rij_eq_symAikBjk< 3 >( P, stress, fInv );
    LvArray::tensorOps::scale< 3, 3 >( P, -detJ * detF );

    FE_TYPE::plusGradNajAij( dNdX, P, stack.fLocal );
  }
};
#undef UPDATE_STRESS

/// The factory used to construct a ExplicitFiniteStrain kernel.
using ExplicitFiniteStrainFactory = finiteElement::KernelFactory< ExplicitFiniteStrain,
                                                                  real64,
                                                                  string const & >;

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSFINITESTRAINEXPLICITNEWMARK_HPP_
