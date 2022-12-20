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
 * @file MultiphasePoromechanics.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{

namespace poromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static multiphase poromechanics.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### MultiphasePoroelastic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static multiphase poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class MultiphasePoromechanics :
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
  static constexpr int numMaxComponents = 3;
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
  MultiphasePoromechanics( NodeManager const & nodeManager,
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
                           localIndex const numComponents,
                           localIndex const numPhases,
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
      dLocalResidualMomentum_dComponents{ {0.0} },
      localResidualMass{ 0.0 },
      dLocalResidualMass_dDisplacement{ {0.0} },
      dLocalResidualMass_dPressure{ {0.0} },
      dLocalResidualMass_dComponents{ {0.0} },
      localPoreVolumeConstraint{ 0.0 },
      dLocalPoreVolumeConstraint_dPressure{ {0.0} },
      dLocalPoreVolumeConstraint_dComponents{ {0.0} },
      localPressureDofIndex{ 0 },
      localComponentDofIndices{ 0 }
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
    real64 dLocalResidualMomentum_dComponents[numDispDofPerElem][numMaxComponents];

    real64 localResidualMass[numMaxComponents];
    real64 dLocalResidualMass_dDisplacement[numMaxComponents][numDispDofPerElem];
    real64 dLocalResidualMass_dPressure[numMaxComponents][1];
    real64 dLocalResidualMass_dComponents[numMaxComponents][numMaxComponents];

    real64 localPoreVolumeConstraint[1];
    real64 dLocalPoreVolumeConstraint_dPressure[1][1];
    real64 dLocalPoreVolumeConstraint_dComponents[1][numMaxComponents];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localPressureDofIndex[1];
    globalIndex localComponentDofIndices[numMaxComponents];

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the MultiphasePoromechanics implementation, global values from the displacement,
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
   * @brief Get a parameter representative of the stiffness, used as physical scaling for the
   * stabilization matrix.
   * @param[in] k Element index.
   * @return A parameter representative of the stiffness matrix dstress/dstrain
   */
  GEOSX_HOST_DEVICE
  real64 computeStabilizationScaling( localIndex const k ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::kernelLaunch
   *
   * ### MultiphasePoromechanics Description
   * Copy of the KernelBase::keranelLaunch function
   * elements.
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];
  real64 const m_gravityAcceleration;

  /// The rank global density
  arrayView2d< real64 const > m_solidDensity;

  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseDensity;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseDensity_n;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dFluidPhaseDensity;

  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > m_fluidPhaseCompFrac;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > m_fluidPhaseCompFrac_n;
  arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > m_dFluidPhaseCompFrac;

  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseMassDensity;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dFluidPhaseMassDensity;

  arrayView2d< real64 const, compflow::USD_PHASE > m_fluidPhaseSaturation;
  arrayView2d< real64 const, compflow::USD_PHASE > m_fluidPhaseSaturation_n;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dFluidPhaseSaturation;

  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dGlobalCompFraction_dGlobalCompDensity;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > m_flowDofNumber;

  /// The rank-global fluid pressure arrays.
  arrayView1d< real64 const > m_fluidPressure_n;
  arrayView1d< real64 const > m_fluidPressure;

  /// Number of components
  localIndex const m_numComponents;

  /// Number of phases
  localIndex const m_numPhases;

};

using MultiphaseKernelFactory = finiteElement::KernelFactory< MultiphasePoromechanics,
                                                              arrayView1d< globalIndex const > const,
                                                              globalIndex const,
                                                              CRSMatrixView< real64, globalIndex const > const,
                                                              arrayView1d< real64 > const,
                                                              real64 const (&)[3],
                                                              string const,
                                                              localIndex const,
                                                              localIndex const,
                                                              string const >;

} // namespace poromechanicsKernels

} // namespace geosx

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_HPP_
