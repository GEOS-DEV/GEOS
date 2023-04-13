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

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geos
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
    m_phaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::flow::dPhaseVolumeFraction >() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseVisc( fluid.phaseViscosity() ),
    m_dPhaseVisc( fluid.dPhaseViscosity() ),
    m_phaseRelPerm( relperm.phaseRelPerm() ),
    m_dPhaseRelPerm_dPhaseVolFrac( relperm.dPhaseRelPerm_dPhaseVolFraction() ),
    m_phaseMob( subRegion.getField< fields::flow::phaseMobility >() ),
    m_dPhaseMob( subRegion.getField< fields::flow::dPhaseMobility >() )
  {}

  /**
   * @brief Compute the phase mobilities in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseMobilityKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && phaseMobilityKernelOp = NoOpFunc{} ) const
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
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dPhaseVolFrac[ei];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const phaseMob = m_phaseMob[ei];
    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseMob = m_dPhaseMob[ei];

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
        for( integer jc = 0; jc < numComp + 2; ++jc )
        {
          dPhaseMob[ip][jc] = 0.0;
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
        dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac[jp][Deriv::dP];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac[jp][Deriv::dC+jc];
        }
      }

      real64 const mobility = relPerm * density / viscosity;

      phaseMob[ip] = mobility;
      dPhaseMob[ip][Deriv::dP] = dRelPerm_dP * density / viscosity
                                 + mobility * (dDens_dP / density - dVisc_dP / viscosity);

      // compositional derivatives
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseMob[ip][Deriv::dC+jc] = dRelPerm_dC[jc] * density / viscosity
                                      + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
      }

      // call the lambda in the phase loop to allow the reuse of the relperm, density, viscosity, and mobility
      // possible use: assemble the derivatives wrt temperature
      phaseMobilityKernelOp( ip, phaseMob[ip], dPhaseMob[ip] );
    }
  }

protected:

  // inputs

  /// Views on the phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dPhaseVolFrac;
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
  arrayView3d< real64, compflow::USD_PHASE_DC > m_dPhaseMob;

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
    StencilAccessors< fields::ghostRank,
                      fields::flow::gravityCoefficient,
                      fields::flow::pressure,
                      fields::flow::dGlobalCompFraction_dGlobalCompDensity,
                      fields::flow::dPhaseVolumeFraction,
                      fields::flow::phaseMobility,
                      fields::flow::dPhaseMobility >;
  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseMassDensity,
                              fields::multifluid::dPhaseMassDensity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

  using CapPressureAccessors =
    StencilMaterialAccessors< CapillaryPressureBase,
                              fields::cappres::phaseCapPressure,
                              fields::cappres::dPhaseCapPressure_dPhaseVolFraction >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] dofNumberAccessor accessor for the dof numbers
   * @param[in] compFlowAccessors accessor for wrappers registered by the solver
   * @param[in] multiFluidAccessors accessor for wrappers registered by the multifluid model
   * @param[in] capPressureAccessors accessor for wrappers registered by the cap pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernelBase( integer const numPhases,
                               globalIndex const rankOffset,
                               integer const hasCapPressure,
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
  integer const m_hasCapPressure;

  /// Time step size
  real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > const m_permeability;
  ElementViewConst< arrayView3d< real64 const > > const m_dPerm_dPres;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  // Primary and secondary variables

  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres;

  /// Views on derivatives of phase volume fractions and comp fractions
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const m_dCompFrac_dCompDens;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseVolFrac;

  /// Views on phase mobilities
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseMob;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseMob;

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

  /// Number of flux support points (hard-coded for TFPA)
  static constexpr integer numFluxSupportPoints = 2;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
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
                           integer const hasCapPressure,
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
                                   hasCapPressure,
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
      numConnectedElems( numElems ),
      dofColIndices( size * numDof ),
      localFlux( numElems * numEqn ),
      localFluxJacobian( numElems * numEqn, size * numDof )
    {}

    // Stencil information

    /// Stencil size for a given connection
    localIndex const stencilSize;
    /// Number of elements connected at a given connection
    localIndex const numConnectedElems;

    // Transmissibility and derivatives

    /// Transmissibility
    real64 transmissibility[maxNumConns][numFluxSupportPoints]{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dPres[maxNumConns][numFluxSupportPoints]{};

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
  { return m_sei[iconn].size(); }

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
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    // first, compute the transmissibilities at this face
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );


    localIndex k[numFluxSupportPoints];
    localIndex connectionIndex = 0;
    for( k[0] = 0; k[0] < stack.numConnectedElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numConnectedElems; ++k[1] )
      {
        /// cell indices
        localIndex const seri[numFluxSupportPoints]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[numFluxSupportPoints] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[numFluxSupportPoints]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // clear working arrays
        real64 compFlux[numComp]{};
        real64 dCompFlux_dP[numFluxSupportPoints][numComp]{};
        real64 dCompFlux_dC[numFluxSupportPoints][numComp][numComp]{};

        real64 const trans[numFluxSupportPoints] = { stack.transmissibility[connectionIndex][0],
                                                     stack.transmissibility[connectionIndex][1] };

        real64 const dTrans_dPres[numFluxSupportPoints] = { stack.dTrans_dPres[connectionIndex][0],
                                                            stack.dTrans_dPres[connectionIndex][1] };

        //***** calculation of flux *****
        // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          // create local work arrays
          real64 densMean = 0.0;
          real64 dDensMean_dP[numFluxSupportPoints]{};
          real64 dDensMean_dC[numFluxSupportPoints][numComp]{};

          real64 phaseFlux = 0.0;
          real64 dPhaseFlux_dP[numFluxSupportPoints]{};
          real64 dPhaseFlux_dC[numFluxSupportPoints][numComp]{};

          real64 presGrad = 0.0;
          real64 dPresGrad_dP[numFluxSupportPoints]{};
          real64 dPresGrad_dC[numFluxSupportPoints][numComp]{};

          real64 gravHead = 0.0;
          real64 dGravHead_dP[numFluxSupportPoints]{};
          real64 dGravHead_dC[numFluxSupportPoints][numComp]{};

          real64 dCapPressure_dC[numComp]{};

          // Working array
          real64 dProp_dC[numComp]{};

          // calculate quantities on primary connected cells
          for( integer i = 0; i < numFluxSupportPoints; ++i )
          {
            localIndex const er  = seri[i];
            localIndex const esr = sesri[i];
            localIndex const ei  = sei[i];

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

          /// compute the TPFA potential difference
          for( integer i = 0; i < numFluxSupportPoints; i++ )
          {
            localIndex const er  = seri[i];
            localIndex const esr = sesri[i];
            localIndex const ei  = sei[i];

            // capillary pressure
            real64 capPressure     = 0.0;
            real64 dCapPressure_dP = 0.0;

            for( integer ic = 0; ic < numComp; ++ic )
            {
              dCapPressure_dC[ic] = 0.0;
            }

            if( m_hasCapPressure )
            {
              capPressure = m_phaseCapPressure[er][esr][ei][0][ip];

              for( integer jp = 0; jp < m_numPhases; ++jp )
              {
                real64 const dCapPressure_dS = m_dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
                dCapPressure_dP += dCapPressure_dS * m_dPhaseVolFrac[er][esr][ei][jp][Deriv::dP];

                for( integer jc = 0; jc < numComp; ++jc )
                {
                  dCapPressure_dC[jc] += dCapPressure_dS * m_dPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
                }
              }
            }

            presGrad += trans[i] * (m_pres[er][esr][ei] - capPressure);
            dPresGrad_dP[i] += trans[i] * (1 - dCapPressure_dP)
                               + dTrans_dPres[i] * (m_pres[er][esr][ei] - capPressure);
            for( integer jc = 0; jc < numComp; ++jc )
            {
              dPresGrad_dC[i][jc] += -trans[i] * dCapPressure_dC[jc];
            }

            real64 const gravD     = trans[i] * m_gravCoef[er][esr][ei];
            real64 const dGravD_dP = dTrans_dPres[i] * m_gravCoef[er][esr][ei];

            // the density used in the potential difference is always a mass density
            // unlike the density used in the phase mobility, which is a mass density
            // if useMass == 1 and a molar density otherwise
            gravHead += densMean * gravD;

            // need to add contributions from both cells the mean density depends on
            for( integer j = 0; j < numFluxSupportPoints; ++j )
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

          localIndex const er_up  = seri[k_up];
          localIndex const esr_up = sesri[k_up];
          localIndex const ei_up  = sei[k_up];

          real64 const mobility = m_phaseMob[er_up][esr_up][ei_up][ip];

          // skip the phase flux if phase not present or immobile upstream
          if( LvArray::math::abs( mobility ) < 1e-20 ) // TODO better constant
          {
            continue;
          }

          // pressure gradient depends on all points in the stencil
          for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
          {
            dPhaseFlux_dP[ke] += dPresGrad_dP[ke] - dGravHead_dP[ke];
            dPhaseFlux_dP[ke] *= mobility;
            for( integer jc = 0; jc < numComp; ++jc )
            {
              dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc] - dGravHead_dC[ke][jc];
              dPhaseFlux_dC[ke][jc] *= mobility;
            }
          }
          // compute phase flux using upwind mobility.
          phaseFlux = mobility * potGrad;

          real64 const dMob_dP  = m_dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
          arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 > dPhaseMobSub =
            m_dPhaseMob[er_up][esr_up][ei_up][ip];

          // add contribution from upstream cell mobility derivatives
          dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
          for( integer jc = 0; jc < numComp; ++jc )
          {
            dPhaseFlux_dC[k_up][jc] += dPhaseMobSub[Deriv::dC+jc] * potGrad;
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
            compFlux[ic] += phaseFlux * ycp;

            // derivatives stemming from phase flux
            for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
            {
              dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
              for( integer jc = 0; jc < numComp; ++jc )
              {
                dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
              }
            }

            // additional derivatives stemming from upstream cell phase composition
            dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFracSub[ic][Deriv::dP];

            // convert derivatives of comp fraction w.r.t. comp fractions to derivatives w.r.t. comp densities
            applyChainRule( numComp,
                            m_dCompFrac_dCompDens[er_up][esr_up][ei_up],
                            dPhaseCompFracSub[ic],
                            dProp_dC,
                            Deriv::dC );
            for( integer jc = 0; jc < numComp; ++jc )
            {
              dCompFlux_dC[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
            }
          } // *** end of upwinding

          // call the lambda in the phase loop to allow the reuse of the phase fluxes and their derivatives
          // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase
          compFluxKernelOp( ip, k, seri, sesri, sei, connectionIndex,
                            k_up, er_up, esr_up, ei_up, potGrad,
                            phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

        } // loop over phases

        // populate local flux vector and derivatives
        for( integer ic = 0; ic < numComp; ++ic )
        {
          integer const eqIndex0 = k[0] * numEqn + ic;
          integer const eqIndex1 = k[1] * numEqn + ic;

          stack.localFlux[eqIndex0]  +=  m_dt * compFlux[ic];
          stack.localFlux[eqIndex1]  -=  m_dt * compFlux[ic];

          for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
          {
            localIndex const localDofIndexPres = k[ke] * numDof;
            stack.localFluxJacobian[eqIndex0][localDofIndexPres] += m_dt * dCompFlux_dP[ke][ic];
            stack.localFluxJacobian[eqIndex1][localDofIndexPres] -= m_dt * dCompFlux_dP[ke][ic];

            for( integer jc = 0; jc < numComp; ++jc )
            {
              localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
              stack.localFluxJacobian[eqIndex0][localDofIndexComp] += m_dt * dCompFlux_dC[ke][ic][jc];
              stack.localFluxJacobian[eqIndex1][localDofIndexComp] -= m_dt * dCompFlux_dC[ke][ic][jc];
            }
          }
        }
        connectionIndex++;
      }   // loop over k[1]
    }   // loop over k[0]

  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && assemblyKernelOp = NoOpFunc{} ) const
  {
    using namespace compositionalMultiphaseUtilities;

    // Apply equation/variable change transformation(s)
    stackArray1d< real64, maxStencilSize * numDof > work( stack.stencilSize * numDof );
    shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numEqn, numDof*stack.stencilSize, stack.numConnectedElems,
                                                             stack.localFluxJacobian, work );
    shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( numComp, numEqn, stack.numConnectedElems,
                                                               stack.localFlux );

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numComp-1)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    for( integer i = 0; i < stack.numConnectedElems; ++i )
    {
      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( m_localMatrix.numRows(), localRow + numComp );

        for( integer ic = 0; ic < numComp; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + ic], stack.localFlux[i * numEqn + ic] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
            ( localRow + ic,
            stack.dofColIndices.data(),
            stack.localFluxJacobian[i * numEqn + ic].dataIfContiguous(),
            stack.stencilSize * numDof );
        }

        // call the lambda to assemble additional terms, such as thermal terms
        assemblyKernelOp( i, localRow );
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
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
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
                   integer const hasCapPressure,
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

      using kernelType = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

      kernelType kernel( numPhases, rankOffset, hasCapPressure, stencilWrapper, dofNumberAccessor,
                         compFlowAccessors, multiFluidAccessors, capPressureAccessors, permeabilityAccessors,
                         dt, localMatrix, localRhs );
      kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

/******************************** DirichletFaceBasedAssemblyKernel ********************************/

/**
 * @class DirichFaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename FLUIDWRAPPER >
class DirichletFaceBasedAssemblyKernel : public FaceBasedAssemblyKernel< NUM_COMP,
                                                                         NUM_DOF,
                                                                         BoundaryStencilWrapper >
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

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_numPhases;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;
  using AbstractBase::m_phaseMob;
  using AbstractBase::m_dPhaseMob;
  using AbstractBase::m_phaseMassDens;
  using AbstractBase::m_dPhaseMassDens;
  using AbstractBase::m_phaseCompFrac;
  using AbstractBase::m_dPhaseCompFrac;
  using AbstractBase::m_dCompFrac_dCompDens;
  using AbstractBase::m_permeability;
  using AbstractBase::m_dPerm_dPres;
  using AbstractBase::m_localMatrix;
  using AbstractBase::m_localRhs;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, BoundaryStencilWrapper >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] faceManager the face manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  DirichletFaceBasedAssemblyKernel( integer const numPhases,
                                    globalIndex const rankOffset,
                                    integer const hasCapPressure,
                                    FaceManager const & faceManager,
                                    BoundaryStencilWrapper const & stencilWrapper,
                                    FLUIDWRAPPER const & fluidWrapper,
                                    DofNumberAccessor const & dofNumberAccessor,
                                    CompFlowAccessors const & compFlowAccessors,
                                    MultiFluidAccessors const & multiFluidAccessors,
                                    CapPressureAccessors const & capPressureAccessors,
                                    PermeabilityAccessors const & permeabilityAccessors,
                                    real64 const & dt,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )
    : Base( numPhases,
            rankOffset,
            hasCapPressure,
            stencilWrapper,
            dofNumberAccessor,
            compFlowAccessors,
            multiFluidAccessors,
            capPressureAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_facePres( faceManager.getField< fields::flow::facePressure >() ),
    m_faceTemp( faceManager.getField< fields::flow::faceTemperature >() ),
    m_faceCompFrac( faceManager.getField< fields::flow::faceGlobalCompFraction >() ),
    m_faceGravCoef( faceManager.getField< fields::flow::gravityCoefficient >() ),
    m_fluidWrapper( fluidWrapper )
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
    StackVariables( localIndex const GEOSX_UNUSED_PARAM( size ),
                    localIndex GEOSX_UNUSED_PARAM( numElems ) )
    {}

    // Transmissibility
    real64 transmissibility = 0.0;

    // Component fluxes and derivatives

    /// Component fluxes
    real64 compFlux[numComp]{};
    /// Derivatives of component fluxes wrt pressure
    real64 dCompFlux_dP[numComp]{};
    /// Derivatives of component fluxes wrt component densities
    real64 dCompFlux_dC[numComp][numComp]{};

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    globalIndex dofColIndices[numDof]{};

    /// Storage for the face local residual vector
    real64 localFlux[numEqn]{};
    /// Storage for the face local Jacobian matrix
    real64 localFluxJacobian[numEqn][numDof]{};

  };


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    globalIndex const offset =
      m_dofNumber[m_seri( iconn, BoundaryStencil::Order::ELEM )][m_sesri( iconn, BoundaryStencil::Order::ELEM )][m_sei( iconn, BoundaryStencil::Order::ELEM )];

    for( integer jdof = 0; jdof < numDof; ++jdof )
    {
      stack.dofColIndices[jdof] = offset + jdof;
    }
  }


  /**
   * @brief Compute the local Dirichlet face flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {
    using Deriv = multifluid::DerivativeOffset;
    using Order = BoundaryStencil::Order;

    localIndex const er  = m_seri( iconn, Order::ELEM );
    localIndex const esr = m_sesri( iconn, Order::ELEM );
    localIndex const ei  = m_sei( iconn, Order::ELEM );
    localIndex const kf  = m_sei( iconn, Order::FACE );

    // Step 1: compute the transmissibility at the boundary face

    real64 dTrans_dPerm[3]{};
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     stack.transmissibility,
                                     dTrans_dPerm );
    real64 const dTrans_dPres = LvArray::tensorOps::AiBi< 3 >( dTrans_dPerm, m_dPerm_dPres[er][esr][ei][0] );

    // Step 2: compute the fluid properties on the face
    // This is needed to get the phase mass density and the phase comp fraction at the face
    // Because we approximate the face mobility using the total element mobility

    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > facePhaseFrac( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > facePhaseDens( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > facePhaseMassDens( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > facePhaseVisc( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > facePhaseEnthalpy( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > facePhaseInternalEnergy( 1, 1, m_numPhases );
    StackArray< real64, 4, constitutive::MultiFluidBase::MAX_NUM_PHASES *NUM_COMP,
                multifluid::LAYOUT_PHASE_COMP > facePhaseCompFrac( 1, 1, m_numPhases, NUM_COMP );
    real64 faceTotalDens = 0.0;

    m_fluidWrapper.compute( m_facePres[kf],
                            m_faceTemp[kf],
                            m_faceCompFrac[kf],
                            facePhaseFrac[0][0],
                            facePhaseDens[0][0],
                            facePhaseMassDens[0][0],
                            facePhaseVisc[0][0],
                            facePhaseEnthalpy[0][0],
                            facePhaseInternalEnergy[0][0],
                            facePhaseCompFrac[0][0],
                            faceTotalDens );

    // Step 3: loop over phases, compute and upwind phase flux and sum contributions to each component's flux

    for( integer ip = 0; ip < m_numPhases; ++ip )
    {

      // working variables
      real64 dDensMean_dC[numComp]{};
      real64 dF_dC[numComp]{};
      real64 dProp_dC[numComp]{};

      real64 phaseFlux = 0.0; // for the lambda
      real64 dPhaseFlux_dP = 0.0;
      real64 dPhaseFlux_dC[numComp]{};


      // Step 3.1: compute the average phase mass density at the face

      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[er][esr][ei],
                      m_dPhaseMassDens[er][esr][ei][0][ip],
                      dProp_dC,
                      Deriv::dC );

      // average density and derivatives
      real64 const densMean = 0.5 * ( m_phaseMassDens[er][esr][ei][0][ip] + facePhaseMassDens[0][0][ip] );
      real64 const dDensMean_dP = 0.5 * m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dC[jc] = 0.5 * dProp_dC[jc];
      }


      // Step 3.2: compute the (TPFA) potential difference at the face

      real64 const gravTimesDz = m_gravCoef[er][esr][ei] - m_faceGravCoef[kf];
      real64 const potDif = m_pres[er][esr][ei] - m_facePres[kf] - densMean * gravTimesDz;
      real64 const f = stack.transmissibility * potDif;
      real64 const dF_dP = stack.transmissibility * ( 1.0 - dDensMean_dP * gravTimesDz ) + dTrans_dPres * potDif;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dF_dC[jc] = -stack.transmissibility * dDensMean_dC[jc] * gravTimesDz;
      }

      // Step 3.3: computation of the mobility
      // We do that before the if/else statement to be able to pass it to the compFluxOpKernel

      // recomputing the exact mobility at the face would be quite complex, as it would require:
      //   1) computing the saturation
      //   2) computing the relperm
      //   3) computing the mobility as \lambda_p = \rho_p kr_p( S_p ) / \mu_p
      // the second step in particular would require yet another dispatch to get the relperm model
      // so, for simplicity, we approximate the face mobility as
      //    \lambda^approx_p = \rho_p S_p / \mu_p
      //                     = \rho_p ( (nu_p / rho_p) * rho_t ) / \mu_p (plugging the expression of saturation)
      //                     = \nu_p * rho_t / \mu_p
      // fortunately, we don't need the derivatives
      real64 const facePhaseMob = ( facePhaseFrac[0][0][ip] > 0.0 )
  ? facePhaseFrac[0][0][ip] * faceTotalDens / facePhaseVisc[0][0][ip]
  : 0.0;

      // *** upwinding ***
      // Step 3.4: upwinding based on the sign of the phase potential gradient
      // It is easier to hard-code the if/else because it is difficult to address elem and face variables in a uniform way

      if( potDif >= 0 ) // the element is upstream
      {

        // compute the phase flux and derivatives using the element mobility
        phaseFlux = m_phaseMob[er][esr][ei][ip] * f;
        dPhaseFlux_dP = m_phaseMob[er][esr][ei][ip] * dF_dP + m_dPhaseMob[er][esr][ei][ip][Deriv::dP] * f;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[jc] =
            m_phaseMob[er][esr][ei][ip] * dF_dC[jc] + m_dPhaseMob[er][esr][ei][ip][Deriv::dC+jc] * f;
        }

        // slice some constitutive arrays to avoid too much indexing in component loop
        arraySlice1d< real64 const, multifluid::USD_PHASE_COMP-3 > phaseCompFracSub =
          m_phaseCompFrac[er][esr][ei][0][ip];
        arraySlice2d< real64 const, multifluid::USD_PHASE_COMP_DC-3 > dPhaseCompFracSub =
          m_dPhaseCompFrac[er][esr][ei][0][ip];

        // compute component fluxes and derivatives using element composition
        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const ycp = phaseCompFracSub[ic];
          stack.compFlux[ic] += phaseFlux * ycp;
          stack.dCompFlux_dP[ic] += dPhaseFlux_dP * ycp + phaseFlux * dPhaseCompFracSub[ic][Deriv::dP];

          applyChainRule( numComp,
                          m_dCompFrac_dCompDens[er][esr][ei],
                          dPhaseCompFracSub[ic],
                          dProp_dC,
                          Deriv::dC );
          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dCompFlux_dC[ic][jc] += dPhaseFlux_dC[jc] * ycp + phaseFlux * dProp_dC[jc];
          }
        }

      }
      else // the face is upstream
      {

        // compute the phase flux and derivatives using the approximated face mobility
        // we only have to take derivatives of the potential gradient in this case
        phaseFlux = facePhaseMob * f;
        dPhaseFlux_dP = facePhaseMob * dF_dP;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[jc] = facePhaseMob * dF_dC[jc];
        }

        // compute component fluxes and derivatives using the face composition
        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const ycp = facePhaseCompFrac[0][0][ip][ic];
          stack.compFlux[ic] += phaseFlux * ycp;
          stack.dCompFlux_dP[ic] += dPhaseFlux_dP * ycp;
          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dCompFlux_dC[ic][jc] += dPhaseFlux_dC[jc] * ycp;
          }
        }
      }

      // call the lambda in the phase loop to allow the reuse of the phase fluxes and their derivatives
      // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase
      compFluxKernelOp( ip, er, esr, ei, kf, f,
                        facePhaseMob, facePhaseEnthalpy[0][0], facePhaseCompFrac[0][0],
                        phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

    }

    // *** end of upwinding

    // Step 4: populate local flux vector and derivatives
    for( integer ic = 0; ic < numComp; ++ic )
    {
      stack.localFlux[ic]            = m_dt * stack.compFlux[ic];
      stack.localFluxJacobian[ic][0] = m_dt * stack.dCompFlux_dP[ic];
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localFluxJacobian[ic][jc+1] = m_dt * stack.dCompFlux_dC[ic][jc];
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && assemblyKernelOp = NoOpFunc{} ) const
  {
    using namespace compositionalMultiphaseUtilities;
    using Order = BoundaryStencil::Order;

    // Apply equation/variable change transformation(s)
    real64 work[numDof]{};
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numDof, stack.localFluxJacobian, work );
    shiftElementsAheadByOneAndReplaceFirstElementWithSum( numComp, stack.localFlux );

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numComp-1)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    if( m_ghostRank[m_seri( iconn, Order::ELEM )][m_sesri( iconn, Order::ELEM )][m_sei( iconn, Order::ELEM )] < 0 )
    {
      globalIndex const globalRow = m_dofNumber[m_seri( iconn, Order::ELEM )][m_sesri( iconn, Order::ELEM )][m_sei( iconn, Order::ELEM )];
      localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
      GEOSX_ASSERT_GE( localRow, 0 );
      GEOSX_ASSERT_GT( AbstractBase::m_localMatrix.numRows(), localRow + numComp );

      for( integer ic = 0; ic < numComp; ++ic )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + ic], stack.localFlux[ic] );
        AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
          ( localRow + ic,
          stack.dofColIndices,
          stack.localFluxJacobian[ic],
          numDof );
      }

      // call the lambda to assemble additional terms, such as thermal terms
      assemblyKernelOp( localRow );
    }
  }

protected:

  /// Views on face pressure, temperature, and composition
  arrayView1d< real64 const > const m_facePres;
  arrayView1d< real64 const > const m_faceTemp;
  arrayView2d< real64 const, compflow::USD_COMP > const m_faceCompFrac;

  /// View on the face gravity coefficient
  arrayView1d< real64 const > const m_faceGravCoef;

  /// Reference to the fluid wrapper
  FLUIDWRAPPER const m_fluidWrapper;

};


/**
 * @class DirichletFaceBasedAssemblyKernelFactory
 */
class DirichletFaceBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] faceManager reference to the face manager
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the boundary stencil wrapper
   * @param[in] fluidBase the multifluid constitutive model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   BoundaryStencilWrapper const & stencilWrapper,
                   MultiFluidBase & fluidBase,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    constitutive::constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
    {
      using FluidType = TYPEOFREF( fluid );
      typename FluidType::KernelWrapper const fluidWrapper = fluid.createKernelWrapper();

      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        integer constexpr NUM_DOF = NC()+1;

        ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
          elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
        dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

        using kernelType = DirichletFaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, typename FluidType::KernelWrapper >;
        typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
        typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
        typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
        typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

        // for now, we neglect capillary pressure in the kernel
        bool const hasCapPressure = false;

        kernelType kernel( numPhases, rankOffset, hasCapPressure, faceManager, stencilWrapper, fluidWrapper,
                           dofNumberAccessor, compFlowAccessors, multiFluidAccessors, capPressureAccessors, permeabilityAccessors,
                           dt, localMatrix, localRhs );
        kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
      } );
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
    StencilAccessors< fields::flow::pressure,
                      fields::flow::gravityCoefficient,
                      fields::flow::phaseVolumeFraction,
                      fields::flow::phaseOutflux,
                      fields::flow::componentOutflux >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseViscosity,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::phaseMassDensity,
                              fields::multifluid::phaseCompFraction >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;


  using RelPermAccessors =
    StencilMaterialAccessors< RelativePermeabilityBase, fields::relperm::phaseRelPerm >;

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
    StencilAccessors< fields::ghostRank,
                      fields::flow::pressure,
                      fields::flow::pressure_n,
                      fields::flow::gravityCoefficient,
                      fields::flow::phaseVolumeFraction,
                      fields::flow::dPhaseVolumeFraction,
                      fields::flow::dGlobalCompFraction_dGlobalCompDensity >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

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
             arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac,
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
          ElementViewConst< arrayView1d< real64 const > > const & pres_n,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
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

} // namespace geos


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
