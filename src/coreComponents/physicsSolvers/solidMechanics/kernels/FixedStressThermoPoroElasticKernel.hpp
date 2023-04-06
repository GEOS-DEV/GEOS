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
 * @file SolidMechanicsFixedStressThermoPoroElasticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_FIXEDSTRESSTHERMOPOROELASTIC_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_FIXEDSTRESSTHERMOPOROELASTIC_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geosx
{

namespace solidMechanicsLagrangianFEMKernels
{

/**
 * @brief Implements kernels for solving the solid part of the fixed-stress thermoporomechanics problem.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### FixedStressThermoPoroElastic Description
 * Implements the KernelBase interface functions required for using the
 * effective stress for the integration of the stress divergence. This is
 * templated on one of the "finite element kernel application" functions
 * such as geosx::finiteElement::RegionBasedKernelApplication.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class FixedStressThermoPoroElastic :
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
  using Base::m_meshData;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  FixedStressThermoPoroElastic( NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager const & faceManager,
                                localIndex const targetRegionIndex,
                                SUBREGION_TYPE const & elementSubRegion,
                                FE_TYPE const & finiteElementSpace,
                                CONSTITUTIVE_TYPE & inputConstitutiveType,
                                arrayView1d< globalIndex const > const inputDofNumber,
                                globalIndex const rankOffset,
                                CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                arrayView1d< real64 > const inputRhs,
                                real64 const (&inputGravityVector)[3] );

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

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
                                       xLocal(),
                                       u_local(),
                                       uhat_local(),
                                       constitutiveStiffness()
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[ 6 ][ 6 ];
  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the FixedStressThermoPoroElastic implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   * For solid mechanics kernels, the strain increment is calculated, and the
   * constitutive update is called. In addition, the constitutive stiffness
   * stack variable is filled by the constitutive model.
   */
  GEOSX_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::kernelLaunch
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

  /// The rank global densities
  arrayView2d< real64 const > const m_density;

  /// The rank-global fluid pressure arrays.
  arrayView1d< real64 const > const m_pressure_n;
  arrayView1d< real64 const > const m_pressure;

  /// The rank-global initial temperature array
  arrayView1d< real64 const > const m_initialTemperature;

  /// The rank-global temperature arrays.
  arrayView1d< real64 const > const m_temperature_n;
  arrayView1d< real64 const > const m_temperature;

  /**
   * @brief Get a parameter representative of the stiffness, used as physical scaling for the
   * stabilization matrix.
   * @param[in] k Element index.
   * @return A parameter representative of the stiffness matrix dstress/dstrain
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 computeStabilizationScaling( localIndex const k ) const
  {
    // TODO: generalize this to other constitutive models (currently we assume linear elasticity).
    return 2.0 * m_constitutiveUpdate.getShearModulus( k );
  }
};

/// The factory used to construct a FixedStressThermoPoroElastic kernel.
using FixedStressThermoPoroElasticFactory = finiteElement::KernelFactory< FixedStressThermoPoroElastic,
                                                                          arrayView1d< globalIndex const > const,
                                                                          globalIndex,
                                                                          CRSMatrixView< real64, globalIndex const > const,
                                                                          arrayView1d< real64 > const,
                                                                          real64 const (&)[3] >;

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_FIXEDSTRESSTHERMOPOROELASTIC_HPP_
