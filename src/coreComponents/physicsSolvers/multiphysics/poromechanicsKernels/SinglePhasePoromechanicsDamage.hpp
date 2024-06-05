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
 * @file SinglePhasePoromechanicsDamageKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSDAMAGE_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSDAMAGE_HPP_

#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geos
{

namespace poromechanicsDamageKernels
{

/**
 * @brief Implements kernels for solving quasi-static single-phase poromechanics with phase-field damage.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### SinglePhasePoroelasticDamage Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static single-phase phase-field poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhasePoromechanicsDamage :
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
  using Base::m_dt;


  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  SinglePhasePoromechanicsDamage( NodeManager const & nodeManager,
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
                                  real64 const inputDt,
                                  real64 const (&gravityVector)[3],
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
    GEOS_HOST_DEVICE
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
      localPressureDofIndex{ 0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][3];
#endif

    // Storage for displacements

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint]{};
    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint]{};

    // Storage for helper variables used in the quadrature point kernel

    /// Total stress
    real64 totalStress[6]{};
    /// Derivative of total stress wrt displacement
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;
    /// Derivative of total stress wrt pressure
    real64 dTotalStress_dPressure[6]{};
    /// Derivative of total stress wrt temperature
    real64 dTotalStress_dTemperature[6]{};

    /// Body force
    real64 bodyForce[3]{};
    /// Derivative of body force wrt volumetric strain increment
    real64 dBodyForce_dVolStrainIncrement[3]{};
    /// Derivative of body force wrt pressure
    real64 dBodyForce_dPressure[3]{};

    /// Delta temperature since the beginning of the simulation
    real64 deltaTemperatureFromInit{}; // for stress computation
    /// Delta temperature since last time step
    real64 deltaTemperatureFromLastStep{}; // for porosity update

    // Storage for residual and degrees of freedom

    /// Linear momentum balance residual
    real64 ( &localResidualMomentum )[numDispDofPerElem];
    /// Derivative of linear momentum balance residual wrt displacement
    real64 ( &dLocalResidualMomentum_dDisplacement )[numDispDofPerElem][numDispDofPerElem];
    /// Derivative of linear momentum balance residual wrt pressure
    real64 dLocalResidualMomentum_dPressure[numDispDofPerElem][1]{};

    /// Mass accumulation
    real64 fluidMassIncrement{};
    /// Derivative of mass accumulation wrt volumetric strain increment
    real64 dFluidMassIncrement_dVolStrainIncrement{};
    /// Derivative of mass accumulation wrt pressure
    real64 dFluidMassIncrement_dPressure{};

    // Storage for mass residual and degrees of freedom

    /// Mass balance residual
    real64 localResidualMass[1];
    /// Derivative of mass balance residual wrt displacement
    real64 dLocalResidualMass_dDisplacement[1][numDispDofPerElem]{};
    /// Derivative of mass balance residual wrt pressure
    real64 dLocalResidualMass_dPressure[1][1]{};

    /// Fracture flow term
    real64 fractureFlowTerm[3]{};
    /// Derivative of the fracture flow term wrt pressure
    real64 dFractureFlowTerm_dPressure[3]{};

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localPressureDofIndex{};

  };
  //*****************************************************************************

  /**
   * @brief Helper function to compute 1) the total stress, 2) the body force term, and 3) the fluidMassIncrement
   * using quantities returned by the PorousDamageSolid constitutive model.
   * This function also computes the derivatives of these three quantities wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] strainIncrement the strain increment used in total stress and porosity computation
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const ( &strainIncrement )[6],
                          StackVariables & stack ) const;

  /**
   * @brief Helper function to compute the body force term (\rho g) and its derivatives wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[in] dSolidDensity_dPressure the derivative of solid density wrt pressure
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeBodyForce( localIndex const k,
                         localIndex const q,
                         real64 const & porosity,
                         real64 const & dPorosity_dVolStrain,
                         real64 const & dPorosity_dPressure,
                         real64 const & dPorosity_dTemperature,
                         real64 const & dSolidDensity_dPressure,
                         StackVariables & stack ) const;

  /**
   * @brief Helper function to compute the fluid mass increment and its derivatives wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] porosity_n the element porosity at the previous converged time step
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFluidIncrement( localIndex const k,
                              localIndex const q,
                              real64 const & porosity,
                              real64 const & porosity_n,
                              real64 const & dPorosity_dVolStrain,
                              real64 const & dPorosity_dPressure,
                              real64 const & dPorosity_dTemperature,
                              StackVariables & stack ) const;

  /**
   * @brief Assemble the local linear momentum balance residual and derivatives using total stress and body force terms
   * @param[in] N displacement finite element basis functions
   * @param[in] dNdX basis function derivatives
   * @param[in] detJxW determinant of the Jacobian transformation matrix times the quadrature weight
   * @param[inout] stack the stack variables
   * @detail This function assembles the discretized form of the following equation
   *   divergence( totalStress ) + bodyForce = 0
   * with the following dependencies on the strainIncrement tensor and pressure
   *   totalStress = totalStress( strainIncrement, pressure)
   *   bodyForce   = bodyForce( strainIncrement, pressure)
   */
  GEOS_HOST_DEVICE
  void assembleMomentumBalanceTerms( real64 const ( &N )[numNodesPerElem],
                                     real64 const ( &dNdX )[numNodesPerElem][3],
                                     real64 const & detJxW,
                                     StackVariables & stack ) const;

  /**
   * @brief Assemble the local mass balance residual and derivatives using fluid mass/energy increment
   * @param[in] dNdX basis function derivatives
   * @param[in] detJxW determinant of the Jacobian transformation matrix times the quadrature weight
   * @param[inout] stack the stack variables
   * @detail This function assembles the discretized form of the temporal derivative in the following equation
   *   dFluidMass_dTime + divergence( fluidMassFlux ) = source  (fluid phase mass balance)
   * with the following dependencies on the strainIncrement tensor and pressure
   *   fluidMass = fluidMass( strainIncrement, pressure)
   *   fluidMassFlux = fluidMassFlux( pressure)
   */
  GEOS_HOST_DEVICE
  void assembleElementBasedFlowTerms( real64 const ( &dNdX )[numNodesPerElem][3],
                                      real64 const & detJxW,
                                      StackVariables & stack ) const;

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the SinglePhasePoromechanicsDamage implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOS_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::kernelLaunch
   *
   * ### SinglePhasePoromechancisDamage Description
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
  real64 const m_gravityVector[3]{};
  /// The L2-norm of the gravity vector
  real64 const m_gravityAcceleration;

  /// The rank global density
  arrayView2d< real64 const > m_solidDensity;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_flowDofNumber;

  /// The rank-global fluid pressure at the previous converged time step
  arrayView1d< real64 const > const m_pressure_n;
  /// The rank-global fluid pressure
  arrayView1d< real64 const > const m_pressure;

  /// Fluid density
  arrayView2d< real64 const > const m_fluidDensity;
  /// Fluid density at the previous converged time step
  arrayView2d< real64 const > const m_fluidDensity_n;
  /// Derivative of fluid density wrt pressure
  arrayView2d< real64 const > const m_dFluidDensity_dPressure;

  arrayView2d< real64 const > const m_fluidPressureGradient;

};

using SinglePhasePoromechanicsDamageKernelFactory =
  finiteElement::KernelFactory< SinglePhasePoromechanicsDamage,
                                arrayView1d< globalIndex const > const,
                                globalIndex const,
                                CRSMatrixView< real64, globalIndex const > const,
                                arrayView1d< real64 > const,
                                real64 const,
                                real64 const (&)[3],
                                string const,
                                string const >;

} // namespace poromechanicsDamageKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSDAMAGE_HPP_
