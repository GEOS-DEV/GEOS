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
 * @file IsothermalCompositionalMultiphaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureExtrinsicData.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/MultiFluidExtrinsicData.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityExtrinsicData.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geosx
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

using namespace constitutive;

/******************************** PhaseMobilityKernel ********************************/

/**
 * @class PhaseMobilityKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase mobilities
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseMobilityKernel : public isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  PhaseMobilityKernel( ObjectManagerBase & subRegion,
                       MultiFluidBase const & fluid,
                       RelativePermeabilityBase const & relperm )
    : Base(),
    m_phaseVolFrac( subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac_dPres( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure >() ),
    m_dPhaseVolFrac_dComp( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity >() ),
    m_dCompFrac_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseVisc( fluid.phaseViscosity() ),
    m_dPhaseVisc( fluid.dPhaseViscosity() ),
    m_phaseRelPerm( relperm.phaseRelPerm() ),
    m_dPhaseRelPerm_dPhaseVolFrac( relperm.dPhaseRelPerm_dPhaseVolFraction() ),
    m_phaseMob( subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseMobility >() ),
    m_dPhaseMob_dPres( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseMobility_dPressure >() ),
    m_dPhaseMob_dComp( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseMobility_dGlobalCompDensity >() )
  {}

  /**
   * @brief Compute the phase mobilities in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseMobilityKernelOp the function used to customize the kernel
   */
  template< typename FUNC = isothermalCompositionalMultiphaseBaseKernels::NoOpFunc >
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && phaseMobilityKernelOp = isothermalCompositionalMultiphaseBaseKernels::NoOpFunc{} ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const phaseVisc = m_phaseVisc[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const dPhaseVisc = m_dPhaseVisc[ei][0];
    arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > const phaseRelPerm = m_phaseRelPerm[ei][0];
    arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > const dPhaseRelPerm_dPhaseVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[ei][0];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac_dComp = m_dPhaseVolFrac_dComp[ei];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const phaseMob = m_phaseMob[ei];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const dPhaseMob_dPres = m_dPhaseMob_dPres[ei];
    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseMob_dComp = m_dPhaseMob_dComp[ei];

    real64 dRelPerm_dC[numComp]{};
    real64 dDens_dC[numComp]{};
    real64 dVisc_dC[numComp]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {

      // compute the phase mobility only if the phase is present
      bool const phaseExists = (phaseVolFrac[ip] > 0);
      if( !phaseExists )
      {
        phaseMob[ip] = 0.0;
        dPhaseMob_dPres[ip] = 0.0;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseMob_dComp[ip][jc] = 0.0;
        }
        continue;
      }

      real64 const density = phaseDens[ip];
      real64 const dDens_dP = dPhaseDens[ip][Deriv::dP];
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseDens[ip], dDens_dC, Deriv::dC );

      real64 const viscosity = phaseVisc[ip];
      real64 const dVisc_dP = dPhaseVisc[ip][Deriv::dP];
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseVisc[ip], dVisc_dC, Deriv::dC );

      real64 const relPerm = phaseRelPerm[ip];
      real64 dRelPerm_dP = 0.0;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dRelPerm_dC[ic] = 0.0;
      }

      for( integer jp = 0; jp < numPhase; ++jp )
      {
        real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
        dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_dPres[jp];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_dComp[jp][jc];
        }
      }

      real64 const mobility = relPerm * density / viscosity;

      phaseMob[ip] = mobility;
      dPhaseMob_dPres[ip] = dRelPerm_dP * density / viscosity
                            + mobility * (dDens_dP / density - dVisc_dP / viscosity);

      // compositional derivatives
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] * density / viscosity
                                  + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
      }

      // call the lambda in the phase loop to allow the reuse of the relperm, density, viscosity, and mobility
      // possible use: assemble the derivatives wrt temperature
      phaseMobilityKernelOp( ip, phaseMob[ip], dPhaseMob_dPres[ip], dPhaseMob_dComp[ip] );
    }
  }

protected:

  // inputs

  /// Views on the phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView2d< real64 const, compflow::USD_PHASE > m_dPhaseVolFrac_dPres;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dPhaseVolFrac_dComp;
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

  /// Views on the phase densities
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseDens;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseDens;

  /// Views on the phase viscosities
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseVisc;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseVisc;

  /// Views on the phase relative permeabilities
  arrayView3d< real64 const, relperm::USD_RELPERM > m_phaseRelPerm;
  arrayView4d< real64 const, relperm::USD_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;

  // outputs

  /// Views on the phase mobilities
  arrayView2d< real64, compflow::USD_PHASE > m_phaseMob;
  arrayView2d< real64, compflow::USD_PHASE > m_dPhaseMob_dPres;
  arrayView3d< real64, compflow::USD_PHASE_DC > m_dPhaseMob_dComp;

};

/**
 * @class PhaseMobilityKernelFactory
 */
class PhaseMobilityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   MultiFluidBase const & fluid,
                   RelativePermeabilityBase const & relperm )
  {
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 2 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 3 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};


/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @brief Base class for FaceBasedAssemblyKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of components/dofs).
 */
class FaceBasedAssemblyKernelBase
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using CompFlowAccessors =
    StencilAccessors< extrinsicMeshData::ghostRank,
                      extrinsicMeshData::flow::gravityCoefficient,
                      extrinsicMeshData::flow::pressure,
                      extrinsicMeshData::flow::deltaPressure,
                      extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity,
                      extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure,
                      extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity,
                      extrinsicMeshData::flow::phaseMobility,
                      extrinsicMeshData::flow::dPhaseMobility_dPressure,
                      extrinsicMeshData::flow::dPhaseMobility_dGlobalCompDensity >;
  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              extrinsicMeshData::multifluid::phaseMassDensity,
                              extrinsicMeshData::multifluid::dPhaseMassDensity,
                              extrinsicMeshData::multifluid::phaseCompFraction,
                              extrinsicMeshData::multifluid::dPhaseCompFraction >;

  using CapPressureAccessors =
    StencilMaterialAccessors< CapillaryPressureBase,
                              extrinsicMeshData::cappres::phaseCapPressure,
                              extrinsicMeshData::cappres::dPhaseCapPressure_dPhaseVolFraction >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              extrinsicMeshData::permeability::permeability,
                              extrinsicMeshData::permeability::dPerm_dPressure >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] capPressureFlag flag specifying whether capillary pressure is used or not
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernelBase( integer const numPhases,
                               globalIndex const rankOffset,
                               integer const capPressureFlag,
                               DofNumberAccessor const & dofNumberAccessor,
                               CompFlowAccessors const & compFlowAccessors,
                               MultiFluidAccessors const & multiFluidAccessors,
                               CapPressureAccessors const & capPressureAccessors,
                               PermeabilityAccessors const & permeabilityAccessors,
                               real64 const & dt,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs );

protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Flag to specify whether capillary pressure is used or not
  integer const m_capPressureFlag;

  /// Time step size
  real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > m_permeability;
  ElementViewConst< arrayView3d< real64 const > > m_dPerm_dPres;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  // Primary and secondary variables

  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres;
  ElementViewConst< arrayView1d< real64 const > > const m_dPres;

  /// Views on derivatives of phase volume fractions and comp fractions
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const m_dCompFrac_dCompDens;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_dPhaseVolFrac_dPres;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseVolFrac_dCompDens;

  /// Views on phase mobilities
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseMob;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_dPhaseMob_dPres;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseMob_dCompDens;

  /// Views on phase mass densities
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseMassDens;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const m_dPhaseMassDens;

  /// Views on phase component fractions
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const m_phaseCompFrac;
  ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const m_dPhaseCompFrac;

  /// Views on phase capillary pressure
  ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const m_phaseCapPressure;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac;

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;
};

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public FaceBasedAssemblyKernelBase
{
public:

  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations (all of them, except the volume balance equation)
  static constexpr integer numEqn = NUM_DOF-1;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::maxNumPointsInFlux;

  /// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::maxNumConnections;

  /// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::maxStencilSize;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] capPressureFlag flag specifying whether capillary pressure is used or not
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( integer const numPhases,
                           globalIndex const rankOffset,
                           integer const capPressureFlag,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           CompFlowAccessors const & compFlowAccessors,
                           MultiFluidAccessors const & multiFluidAccessors,
                           CapPressureAccessors const & capPressureAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    : FaceBasedAssemblyKernelBase( numPhases,
                                   rankOffset,
                                   capPressureFlag,
                                   dofNumberAccessor,
                                   compFlowAccessors,
                                   multiFluidAccessors,
                                   capPressureAccessors,
                                   permeabilityAccessors,
                                   dt,
                                   localMatrix,
                                   localRhs ),
    m_stencilWrapper( stencilWrapper ),
    m_seri( stencilWrapper.getElementRegionIndices() ),
    m_sesri( stencilWrapper.getElementSubRegionIndices() ),
    m_sei( stencilWrapper.getElementIndices() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOSX_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : stencilSize( size ),
      numFluxElems( numElems ),
      compFlux( numComp ),
      dCompFlux_dP( size, numComp ),
      dCompFlux_dC( size, numComp, numComp ),
      dofColIndices( size * numDof ),
      localFlux( numElems * numEqn ),
      localFluxJacobian( numElems * numEqn, size * numDof )
    {}

    // Stencil information

    /// Stencil size for a given connection
    localIndex const stencilSize;

    /// Number of elements for a given connection
    localIndex const numFluxElems;

    // Transmissibility and derivatives

    /// Transmissibility
    real64 transmissibility[maxNumConns][2]{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dPres[maxNumConns][2]{};

    // Component fluxes and derivatives

    /// Component fluxes
    stackArray1d< real64, numComp > compFlux;
    /// Derivatives of component fluxes wrt pressure
    stackArray2d< real64, maxStencilSize * numComp > dCompFlux_dP;
    /// Derivatives of component fluxes wrt component densities
    stackArray3d< real64, maxStencilSize * numComp * numComp > dCompFlux_dC;

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    stackArray1d< globalIndex, maxNumElems * numDof > dofColIndices;

    /// Storage for the face local residual vector (all equations except volume balance)
    stackArray1d< real64, maxNumElems * numEqn > localFlux;
    /// Storage for the face local Jacobian matrix
    stackArray2d< real64, maxNumElems * numEqn * maxStencilSize * numDof > localFluxJacobian;

  };


  /**
   * @brief Getter for the stencil size at this connection
   * @param[in] iconn the connection index
   * @return the size of the stencil at this connection
   */
  GEOSX_HOST_DEVICE
  localIndex stencilSize( localIndex const iconn ) const
  { return meshMapUtilities::size1( m_sei, iconn ); }

  /**
   * @brief Getter for the number of elements at this connection
   * @param[in] iconn the connection index
   * @return the number of elements at this connection
   */
  GEOSX_HOST_DEVICE
  localIndex numPointsInFlux( localIndex const iconn ) const
  { return m_stencilWrapper.numPointsInFlux( iconn ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // set degrees of freedom indices for this face
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      globalIndex const offset = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];

      for( integer jdof = 0; jdof < numDof; ++jdof )
      {
        stack.dofColIndices[i * numDof + jdof] = offset + jdof;
      }
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = isothermalCompositionalMultiphaseBaseKernels::NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && compFluxKernelOp = isothermalCompositionalMultiphaseBaseKernels::NoOpFunc{} ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    // first, compute the transmissibilities at this face
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      // clear working arrays
      real64 densMean{};
      stackArray1d< real64, maxNumElems > dDensMean_dP( stack.numFluxElems );
      stackArray2d< real64, maxNumElems * numComp > dDensMean_dC( stack.numFluxElems, numComp );

      // create local work arrays
      real64 phaseFlux{};
      real64 dPhaseFlux_dP[maxStencilSize]{};
      real64 dPhaseFlux_dC[maxStencilSize][numComp]{};

      real64 presGrad{};
      stackArray1d< real64, maxStencilSize > dPresGrad_dP( stack.stencilSize );
      stackArray2d< real64, maxStencilSize *numComp > dPresGrad_dC( stack.stencilSize, numComp );

      real64 gravHead{};
      stackArray1d< real64, maxNumElems > dGravHead_dP( stack.numFluxElems );
      stackArray2d< real64, maxNumElems * numComp > dGravHead_dC( stack.numFluxElems, numComp );

      real64 dCapPressure_dC[numComp]{};

      // Working array
      real64 dProp_dC[numComp]{};

      // calculate quantities on primary connected cells
      for( integer i = 0; i < stack.numFluxElems; ++i )
      {
        localIndex const er  = m_seri( iconn, i );
        localIndex const esr = m_sesri( iconn, i );
        localIndex const ei  = m_sei( iconn, i );

        // density
        real64 const density  = m_phaseMassDens[er][esr][ei][0][ip];
        real64 const dDens_dP = m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];

        applyChainRule( numComp,
                        m_dCompFrac_dCompDens[er][esr][ei],
                        m_dPhaseMassDens[er][esr][ei][0][ip],
                        dProp_dC,
                        Deriv::dC );

        // average density and derivatives
        densMean += 0.5 * density;
        dDensMean_dP[i] = 0.5 * dDens_dP;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dDensMean_dC[i][jc] = 0.5 * dProp_dC[jc];
        }
      }

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      for( integer i = 0; i < stack.stencilSize; ++i )
      {
        localIndex const er  = m_seri( iconn, i );
        localIndex const esr = m_sesri( iconn, i );
        localIndex const ei  = m_sei( iconn, i );

        // capillary pressure
        real64 capPressure     = 0.0;
        real64 dCapPressure_dP = 0.0;

        for( integer ic = 0; ic < numComp; ++ic )
        {
          dCapPressure_dC[ic] = 0.0;
        }

        if( m_capPressureFlag )
        {
          capPressure = m_phaseCapPressure[er][esr][ei][0][ip];

          for( integer jp = 0; jp < m_numPhases; ++jp )
          {
            real64 const dCapPressure_dS = m_dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
            dCapPressure_dP += dCapPressure_dS * m_dPhaseVolFrac_dPres[er][esr][ei][jp];

            for( integer jc = 0; jc < numComp; ++jc )
            {
              dCapPressure_dC[jc] += dCapPressure_dS * m_dPhaseVolFrac_dCompDens[er][esr][ei][jp][jc];
            }
          }
        }

        presGrad += stack.transmissibility[0][i] * (m_pres[er][esr][ei] + m_dPres[er][esr][ei] - capPressure);
        dPresGrad_dP[i] += stack.transmissibility[0][i] * (1 - dCapPressure_dP)
                           + stack.dTrans_dPres[0][i] * (m_pres[er][esr][ei] + m_dPres[er][esr][ei] - capPressure);
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPresGrad_dC[i][jc] += -stack.transmissibility[0][i] * dCapPressure_dC[jc];
        }

        real64 const gravD     = stack.transmissibility[0][i] * m_gravCoef[er][esr][ei];
        real64 const dGravD_dP = stack.dTrans_dPres[0][i] * m_gravCoef[er][esr][ei];

        // the density used in the potential difference is always a mass density
        // unlike the density used in the phase mobility, which is a mass density
        // if useMass == 1 and a molar density otherwise
        gravHead += densMean * gravD;

        // need to add contributions from both cells the mean density depends on
        for( integer j = 0; j < stack.numFluxElems; ++j )
        {
          dGravHead_dP[j] += dDensMean_dP[j] * gravD + dGravD_dP * densMean;
          for( integer jc = 0; jc < numComp; ++jc )
          {
            dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
          }
        }
      }

      // *** upwinding ***

      // compute phase potential gradient
      real64 const potGrad = presGrad - gravHead;

      // choose upstream cell
      localIndex const k_up = (potGrad >= 0) ? 0 : 1;

      localIndex const er_up  = m_seri( iconn, k_up );
      localIndex const esr_up = m_sesri( iconn, k_up );
      localIndex const ei_up  = m_sei( iconn, k_up );

      real64 const mobility = m_phaseMob[er_up][esr_up][ei_up][ip];

      // skip the phase flux if phase not present or immobile upstream
      if( LvArray::math::abs( mobility ) < 1e-20 ) // TODO better constant
      {
        continue;
      }

      // pressure gradient depends on all points in the stencil
      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc];
        }
      }

      // gravitational head depends only on the two cells connected (same as mean density)
      for( integer ke = 0; ke < stack.numFluxElems; ++ke )
      {
        dPhaseFlux_dP[ke] -= dGravHead_dP[ke];
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] -= dGravHead_dC[ke][jc];
        }
      }

      // compute the phase flux and derivatives using upstream cell mobility
      phaseFlux = mobility * potGrad;
      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        dPhaseFlux_dP[ke] *= mobility;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] *= mobility;
        }
      }

      real64 const dMob_dP  = m_dPhaseMob_dPres[er_up][esr_up][ei_up][ip];
      arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 > dPhaseMob_dCompSub =
        m_dPhaseMob_dCompDens[er_up][esr_up][ei_up][ip];

      // add contribution from upstream cell mobility derivatives
      dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[k_up][jc] += dPhaseMob_dCompSub[jc] * potGrad;
      }

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d< real64 const, multifluid::USD_PHASE_COMP-3 > phaseCompFracSub =
        m_phaseCompFrac[er_up][esr_up][ei_up][0][ip];
      arraySlice2d< real64 const, multifluid::USD_PHASE_COMP_DC-3 > dPhaseCompFracSub =
        m_dPhaseCompFrac[er_up][esr_up][ei_up][0][ip];

      // compute component fluxes and derivatives using upstream cell composition
      for( integer ic = 0; ic < numComp; ++ic )
      {
        real64 const ycp = phaseCompFracSub[ic];
        stack.compFlux[ic] += phaseFlux * ycp;

        // derivatives stemming from phase flux
        for( integer ke = 0; ke < stack.stencilSize; ++ke )
        {
          stack.dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
          }
        }

        // additional derivatives stemming from upstream cell phase composition
        stack.dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFracSub[ic][Deriv::dP];

        // convert derivatives of comp fraction w.r.t. comp fractions to derivatives w.r.t. comp densities
        applyChainRule( numComp,
                        m_dCompFrac_dCompDens[er_up][esr_up][ei_up],
                        dPhaseCompFracSub[ic],
                        dProp_dC,
                        Deriv::dC );
        for( integer jc = 0; jc < numComp; ++jc )
        {
          stack.dCompFlux_dC[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
        }
      }

      // call the lambda in the phase loop to allow the reuse of the phase fluxes and their derivatives
      // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase
      compFluxKernelOp( ip, k_up, er_up, esr_up, ei_up, potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

    }

    // *** end of upwinding

    // populate local flux vector and derivatives
    for( integer ic = 0; ic < numComp; ++ic )
    {
      stack.localFlux[ic]           =  m_dt * stack.compFlux[ic];
      stack.localFlux[numComp + ic] = -m_dt * stack.compFlux[ic];

      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        localIndex const localDofIndexPres = ke * numDof;
        stack.localFluxJacobian[ic][localDofIndexPres]           =  m_dt * stack.dCompFlux_dP[ke][ic];
        stack.localFluxJacobian[numComp + ic][localDofIndexPres] = -m_dt * stack.dCompFlux_dP[ke][ic];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
          stack.localFluxJacobian[ic][localDofIndexComp]           =  m_dt * stack.dCompFlux_dC[ke][ic][jc];
          stack.localFluxJacobian[numComp + ic][localDofIndexComp] = -m_dt * stack.dCompFlux_dC[ke][ic][jc];
        }
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    using namespace compositionalMultiphaseUtilities;

    // Apply equation/variable change transformation(s)
    stackArray1d< real64, maxStencilSize * numDof > work( stack.stencilSize * numDof );
    shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numDof*stack.stencilSize, stack.numFluxElems,
                                                             stack.localFluxJacobian, work );
    shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( numComp, stack.numFluxElems,
                                                               stack.localFlux );

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numComp-1)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    for( integer i = 0; i < stack.numFluxElems; ++i )
    {
      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( m_localMatrix.numRows(), localRow + numComp );

        for( integer ic = 0; ic < numComp; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + ic], stack.localFlux[i * numComp + ic] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
            ( localRow + ic,
            stack.dofColIndices.data(),
            stack.localFluxJacobian[i * numComp + ic].dataIfContiguous(),
            stack.stencilSize * numDof );
        }
      }
    }
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numConnections the number of connections
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numConnections,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    forAll< POLICY >( numConnections, [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent.stencilSize( iconn ),
                                                  kernelComponent.numPointsInFlux( iconn ) );

      kernelComponent.setup( iconn, stack );
      kernelComponent.computeFlux( iconn, stack );
      kernelComponent.complete( iconn, stack );
    } );
  }

protected:

  // Stencil information

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
};

/**
 * @class FaceBasedAssemblyKernelFactory
 */
class FaceBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] capPressureFlag flag specifying whether capillary pressure is used or not
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const capPressureFlag,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC()+1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      using KERNEL_TYPE = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

      KERNEL_TYPE kernel( numPhases, rankOffset, capPressureFlag, stencilWrapper, dofNumberAccessor,
                          compFlowAccessors, multiFluidAccessors, capPressureAccessors, permeabilityAccessors,
                          dt, localMatrix, localRhs );
      KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

/******************************** CFLFluxKernel ********************************/

/**
 * @brief Functions to compute the (outflux) total volumetric flux needed in the calculation of CFL numbers
 */
struct CFLFluxKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename VIEWTYPE >
  using ElementView = ElementRegionManager::ElementView< VIEWTYPE >;

  using CompFlowAccessors =
    StencilAccessors< extrinsicMeshData::flow::pressure,
                      extrinsicMeshData::flow::gravityCoefficient,
                      extrinsicMeshData::flow::phaseVolumeFraction,
                      extrinsicMeshData::flow::phaseOutflux,
                      extrinsicMeshData::flow::componentOutflux >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              extrinsicMeshData::multifluid::phaseViscosity,
                              extrinsicMeshData::multifluid::phaseDensity,
                              extrinsicMeshData::multifluid::phaseMassDensity,
                              extrinsicMeshData::multifluid::phaseCompFraction >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              extrinsicMeshData::permeability::permeability,
                              extrinsicMeshData::permeability::dPerm_dPressure >;


  using RelPermAccessors =
    StencilMaterialAccessors< RelativePermeabilityBase, extrinsicMeshData::relperm::phaseRelPerm >;

  template< integer NC, localIndex NUM_ELEMS, localIndex maxStencilSize >
  GEOSX_HOST_DEVICE
  static void
  compute( integer const numPhases,
           localIndex const stencilSize,
           real64 const & dt,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           real64 const (&transmissibility)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
           ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & phaseRelPerm,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseVisc,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux,
           ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux );

  template< integer NC, typename STENCILWRAPPER_TYPE >
  static void
  launch( integer const numPhases,
          real64 const & dt,
          STENCILWRAPPER_TYPE const & stencil,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & phaseRelPerm,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseVisc,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux,
          ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux );
};

/******************************** CFLKernel ********************************/

/**
 * @brief Functions to compute the CFL number using the phase volumetric outflux and the component mass outflux in each cell
 */
struct CFLKernel
{

  static constexpr real64 minPhaseMobility = 1e-12;
  static constexpr real64 minComponentFraction = 1e-12;

  template< integer NP >
  GEOSX_HOST_DEVICE
  static void
  computePhaseCFL( real64 const & poreVol,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
                   arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > phaseRelPerm,
                   arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > dPhaseRelPerm_dPhaseVolFrac,
                   arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseVisc,
                   arraySlice1d< real64 const, compflow::USD_PHASE- 1 > phaseOutflux,
                   real64 & phaseCFLNumber );

  template< integer NC >
  GEOSX_HOST_DEVICE
  static void
  computeCompCFL( real64 const & poreVol,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compDens,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compFrac,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compOutflux,
                  real64 & compCFLNumber );

  template< integer NC, integer NP >
  static void
  launch( localIndex const size,
          arrayView1d< real64 const > const & volume,
          arrayView2d< real64 const > const & porosity,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux,
          arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux,
          arrayView1d< real64 > const & phaseCFLNumber,
          arrayView1d< real64 > const & compCFLNumber,
          real64 & maxPhaseCFLNumber,
          real64 & maxCompCFLNumber );

};

/******************************** AquiferBCKernel ********************************/

/**
 * @brief Functions to assemble aquifer boundary condition contributions to residual and Jacobian
 */
struct AquiferBCKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using CompFlowAccessors =
    StencilAccessors< extrinsicMeshData::ghostRank,
                      extrinsicMeshData::flow::pressure,
                      extrinsicMeshData::flow::deltaPressure,
                      extrinsicMeshData::flow::gravityCoefficient,
                      extrinsicMeshData::flow::phaseVolumeFraction,
                      extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure,
                      extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity,
                      extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              extrinsicMeshData::multifluid::phaseDensity,
                              extrinsicMeshData::multifluid::dPhaseDensity,
                              extrinsicMeshData::multifluid::phaseCompFraction,
                              extrinsicMeshData::multifluid::dPhaseCompFraction >;

  template< integer NC >
  GEOSX_HOST_DEVICE
  static void
    compute( integer const numPhases,
             integer const ipWater,
             bool const allowAllPhasesIntoAquifer,
             real64 const & aquiferVolFlux,
             real64 const & dAquiferVolFlux_dPres,
             real64 const & aquiferWaterPhaseDens,
             arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
             arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens,
             arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dPres,
             arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac_dCompDens,
             arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac,
             arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac,
             arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens,
             real64 const & dt,
             real64 ( &localFlux )[NC],
             real64 ( &localFluxJacobian )[NC][NC+1] );

  template< integer NC >
  static void
  launch( integer const numPhases,
          integer const ipWater,
          bool const allowAllPhasesIntoAquifer,
          BoundaryStencil const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
          real64 const & aquiferWaterPhaseDens,
          arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseVolFrac_dPres,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
          real64 const & timeAtBeginningOfStep,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
