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

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICS_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICS_HPP_

#include "codingUtilities/Utilities.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/PoromechanicsBase.hpp"

namespace geosx
{

namespace poromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static multiphase poromechanics.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
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
  public PoromechanicsBase< SUBREGION_TYPE,
                            CONSTITUTIVE_TYPE,
                            FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = PoromechanicsBase< SUBREGION_TYPE,
                                  CONSTITUTIVE_TYPE,
                                  FE_TYPE >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  static constexpr int maxNumComponents = 3;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_gravityAcceleration;
  using Base::m_gravityVector;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::m_solidDensity;
  using Base::m_pressure_n;
  using Base::m_pressure;
  using Base::m_flowDofNumber;
  using Base::m_meshData;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param gravityVector The gravity vector.
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
                           real64 const (&gravityVector)[3],
                           string const inputFlowDofKey,
                           localIndex const numComponents,
                           localIndex const numPhases,
                           string const fluidModelKey,
                           string const capPressureModelKey );

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
      Base::StackVariables()
    {}

    /// Mass accumulation
    real64 compMassIncrement[maxNumComponents]{};
    /// Derivative of mass accumulation wrt volumetric strain increment
    real64 dCompMassIncrement_dVolStrainIncrement[maxNumComponents]{};
    /// Derivative of mass accumulation wrt pressure
    real64 dCompMassIncrement_dPressure[maxNumComponents]{};
    /// Derivative of mass accumulation wrt comp density
    real64 dCompMassIncrement_dComponents[maxNumComponents][maxNumComponents]{};

    /// Pore volume constraint
    real64 poreVolConstraint{};
    /// Derivative of pore volume constraint wrt pressure
    real64 dPoreVolConstraint_dPressure{};
    /// Derivative of pore volume constraint wrt comp density
    real64 dPoreVolConstraint_dComponents[1][maxNumComponents]{};

    /// Derivative of body force wrt comp density
    real64 dBodyForce_dComponents[3][maxNumComponents]{};

    // Storage for mass residual and degrees of freedom

    /// Derivative of momemtum balance residual wrt comp density
    real64 dLocalResidualMomentum_dComponents[Base::StackVariables::numDispDofPerElem][maxNumComponents]{};

    /// Mass balance residual
    real64 localResidualMass[maxNumComponents]{};
    /// Derivative of mass balance residual wrt volumetric strain increment
    real64 dLocalResidualMass_dDisplacement[maxNumComponents][Base::StackVariables::numDispDofPerElem]{};
    /// Derivative of mass balance residual wrt pressure
    real64 dLocalResidualMass_dPressure[maxNumComponents][1]{};
    /// Derivative of mass balance residual wrt components
    real64 dLocalResidualMass_dComponents[maxNumComponents][maxNumComponents]{};

    /// Pore volume constraint residual (sum of saturations = 1)
    real64 localResidualPoreVolConstraint[1]{};
    /// Derivative of pore volume constraint residual wrt pressure
    real64 dLocalResidualPoreVolConstraint_dPressure[1][1]{};
    /// Derivative of pore volume constraint residual wrt components
    real64 dLocalResidualPoreVolConstraint_dComponents[1][maxNumComponents]{};

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localComponentDofIndices[maxNumComponents]{};

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


  /**
   * @brief Helper function to compute 1) the total stress, 2) the body force term, and 3) the fluidMassIncrement
   * using quantities returned by the PorousSolid constitutive model.
   * This function also computes the derivatives of these three quantities wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          StackVariables & stack ) const;

  /**
   * @brief Helper function to compute the body force term (\rho g) and its derivatives wrt primary variables
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[in] dSolidDensity_dPressure the derivative of solid density wrt pressure
   * @param[inout] stack the stack variables
   * @param[inout] bodyForceKernelOp the lambda function to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeBodyForce( localIndex const k,
                         localIndex const q,
                         real64 const & porosity,
                         real64 const & dPorosity_dVolStrain,
                         real64 const & dPorosity_dPressure,
                         real64 const & dPorosity_dTemperature,
                         real64 const & dSolidDensity_dPressure,
                         StackVariables & stack,
                         FUNC && bodyForceKernelOp = NoOpFunc{} ) const;

  /**
   * @brief Helper function to compute the fluid mass increment and its derivatives wrt primary variables
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] porosity_n the element porosity at the previous converged time step
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[inout] stack the stack variables
   * @param[inout] fluidIncrementKernelOp the lambda function to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeFluidIncrement( localIndex const k,
                              localIndex const q,
                              real64 const & porosity,
                              real64 const & porosity_n,
                              real64 const & dPorosity_dVolStrain,
                              real64 const & dPorosity_dPressure,
                              real64 const & dPorosity_dTemperature,
                              StackVariables & stack,
                              FUNC && fluidIncrementKernelOp = NoOpFunc{} ) const;

  /**
   * @brief Helper function to compute the pore-volume constraint and its derivatives wrt primary variables
   * @param[in] k the element index
   * @param[in] porosity the element porosity
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void computePoreVolumeConstraint( localIndex const k,
                                    real64 const & porosity,
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
  GEOSX_HOST_DEVICE
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
  GEOSX_HOST_DEVICE
  void assembleElementBasedFlowTerms( real64 const ( &dNdX )[numNodesPerElem][3],
                                      real64 const & detJxW,
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

  /// Views on phase densities and derivatives
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseDensity;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseDensity_n;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dFluidPhaseDensity;

  /// Views on phase component fractions and derivatives
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > m_fluidPhaseCompFrac;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > m_fluidPhaseCompFrac_n;
  arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > m_dFluidPhaseCompFrac;

  /// Views on phase mass densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseMassDensity;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dFluidPhaseMassDensity;

  /// Views on phase volume fraction (saturation)
  arrayView2d< real64 const, compflow::USD_PHASE > const m_fluidPhaseVolFrac;
  arrayView2d< real64 const, compflow::USD_PHASE > const m_fluidPhaseVolFrac_n;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > const m_dFluidPhaseVolFrac;

  /// Views on capillary pressure
  arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > m_fluidPhaseCapPressure;
  arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > m_dFluidPhaseCapPressure_dPhaseVolFrac;

  /// Views on derivatives of global comp fraction wrt global comp density
  arrayView3d< real64 const, compflow::USD_COMP_DC > const m_dGlobalCompFraction_dGlobalCompDensity;

  /// Number of components
  localIndex const m_numComponents;

  /// Number of phases
  localIndex const m_numPhases;

};

using MultiphasePoromechanicsKernelFactory =
  finiteElement::KernelFactory< MultiphasePoromechanics,
                                arrayView1d< globalIndex const > const,
                                globalIndex const,
                                CRSMatrixView< real64, globalIndex const > const,
                                arrayView1d< real64 > const,
                                real64 const (&)[3],
                                string const,
                                localIndex const,
                                localIndex const,
                                string const,
                                string const >;

} // namespace poromechanicsKernels

} // namespace geosx

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICS_HPP_
