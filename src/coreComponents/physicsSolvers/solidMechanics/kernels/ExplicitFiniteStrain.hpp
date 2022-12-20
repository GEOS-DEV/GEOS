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
 * @file ExplicitFiniteStrain.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITFINITESTRAIN_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITFINITESTRAIN_HPP_

#include "ExplicitSmallStrain.hpp"
#include "finiteElement/Kinematics.h"

namespace geosx
{

namespace solidMechanicsLagrangianFEMKernels
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
                        string const elementListName );


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
  void setup( localIndex const k,
              StackVariables & stack ) const;

  /**
   * @copydoc ExplicitSmallStrain::quadraturePointKernel
   */
  GEOSX_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;


<<<<<<< HEAD:src/coreComponents/physicsSolvers/solidMechanics/SolidMechanicsFiniteStrainExplicitNewmarkKernel.hpp
    // chain rule: calculate dv/dx^(n+1/2) = dv/dX * dX/dx^(n+1/2)
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( Ldt, dUhatdX, fInv );

    // calculate gradient (end of step)
    LvArray::tensorOps::copy< 3, 3 >( F, dUhatdX );
    LvArray::tensorOps::add< 3, 3 >( F, dUdX );
    LvArray::tensorOps::addIdentity< 3 >( F, 1.0 );
    real64 const detF = LvArray::tensorOps::invert< 3 >( fInv, F );

    real64 Rot[ 3 ][ 3 ];
    real64 Dadt[ 6 ];
    real64 timeIncrement;
    HughesWinget( Rot, Dadt, Ldt );

    real64 stress[ 6 ] = { };
    m_constitutiveUpdate.hypoUpdate_StressOnly( k, q, timeIncrement, Dadt, Rot, stress );
=======
  /**
   * @copydoc ExplicitSmallStrain::kernelLaunch
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );

>>>>>>> develop:src/coreComponents/physicsSolvers/solidMechanics/kernels/ExplicitFiniteStrain.hpp


};
/// The factory used to construct a ExplicitFiniteStrain kernel.
using ExplicitFiniteStrainFactory = finiteElement::KernelFactory< ExplicitFiniteStrain,
                                                                  real64,
                                                                  string const >;

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITFINITESTRAIN_HPP_
