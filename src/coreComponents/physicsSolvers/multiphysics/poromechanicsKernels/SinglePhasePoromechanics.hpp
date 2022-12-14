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
 * @file SinglePhasePoromechanics.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{

namespace poromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static single-phase poromechanics.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### SinglePhasePoroelastic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static single-phase poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhasePoromechanics :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;


  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  SinglePhasePoromechanics( NodeManager const & nodeManager,
                            EdgeManager const & edgeManager,
                            FaceManager const & faceManager,
                            localIndex const targetRegionIndex,
                            SUBREGION_TYPE const & elementSubRegion,
                            FE_TYPE const & finiteElementSpace,
                            CONSTITUTIVE_TYPE & inputConstitutiveType,
                            arrayView1d< globalIndex const > const inputDispDofNumber,
                            globalIndex const rankOffset,
                            CRSMatrixView< real64, globalIndex const > const inputMatrix,
                            arrayView1d< real64 > const inputRhs,
                            real64 const (&inputGravityVector)[3],
                            string const inputFlowDofKey,
                            string const fluidModelKey );

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    static constexpr int numDispDofPerElem =  Base::StackVariables::maxNumRows;

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            u_local(),
            uhat_local(),
            localResidualMomentum( Base::StackVariables::localResidual ),
      dLocalResidualMomentum_dDisplacement( Base::StackVariables::localJacobian ),
      dLocalResidualMomentum_dPressure{ {0.0} },
      localResidualMass{ 0.0 },
      dLocalResidualMass_dDisplacement{ {0.0} },
      dLocalResidualMass_dPressure{ {0.0} },
      localFlowDofIndex{ 0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][3];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    real64 ( &localResidualMomentum )[numDispDofPerElem];
    real64 ( &dLocalResidualMomentum_dDisplacement )[numDispDofPerElem][numDispDofPerElem];
    real64 dLocalResidualMomentum_dPressure[numDispDofPerElem][1];

    real64 localResidualMass[1];
    real64 dLocalResidualMass_dDisplacement[1][numDispDofPerElem];
    real64 dLocalResidualMass_dPressure[1][1];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localFlowDofIndex[1];

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the SinglePhase implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  GEOSX_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::kernelLaunch
   *
   * ### SinglePhasePoromechancis Description
   * Copy of the KernelBase::kernelLaunch function
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];
  real64 const m_gravityAcceleration;

  /// The rank global densities
  arrayView2d< real64 const > const m_solidDensity;
  arrayView2d< real64 const > const m_fluidDensity;
  arrayView2d< real64 const > const m_fluidDensity_n;
  arrayView2d< real64 const > const m_dFluidDensity_dPressure;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_flowDofNumber;

  /// The rank-global fluid pressure arrays.
  arrayView1d< real64 const > const m_fluidPressure_n;
  arrayView1d< real64 const > const m_fluidPressure;

};

using SinglePhaseKernelFactory = finiteElement::KernelFactory< SinglePhasePoromechanics,
                                                               arrayView1d< globalIndex const > const,
                                                               globalIndex const,
                                                               CRSMatrixView< real64, globalIndex const > const,
                                                               arrayView1d< real64 > const,
                                                               real64 const (&)[3],
                                                               string const,
                                                               string const >;

} // namespace poromechanicsKernels

} // namespace geosx


#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_HPP_
