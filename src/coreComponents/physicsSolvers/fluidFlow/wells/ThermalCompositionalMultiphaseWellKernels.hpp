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
 * @file ThermalCompositionalMultiphaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALCOMPOSITIONALMULTIPHASEWELLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALCOMPOSITIONALMULTIPHASEWELLKERNELS_HPP

#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseWellKernels
{

using namespace constitutive;

/******************************** TotalMassDensityKernel ****************************/

/**
 * @class TotalMassDensityKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the total mass density
 */
template< integer NUM_COMP, integer NUM_PHASE >
class TotalMassDensityKernel : public compositionalMultiphaseWellKernels::TotalMassDensityKernel< NUM_COMP, NUM_PHASE >
{
public:
  using Base = compositionalMultiphaseWellKernels::TotalMassDensityKernel< NUM_COMP, NUM_PHASE >;
  using Base::m_dCompFrac_dCompDens;
  using Base::m_dPhaseMassDens;
  using Base::m_dPhaseVolFrac;
  using Base::m_dTotalMassDens;
  using Base::m_dTotalMassDens_dCompDens;
  using Base::m_dTotalMassDens_dPres;
  using Base::m_phaseMassDens;
  using Base::m_phaseVolFrac;
  using Base::m_totalMassDens;
  using Base::numComp;
  using Base::numPhase;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  TotalMassDensityKernel( ObjectManagerBase & subRegion,
                          MultiFluidBase const & fluid )
    : Base( subRegion, fluid ),
    m_dTotalMassDens_dTemp( subRegion.getField< fields::well::dTotalMassDensity_dTemperature >())
  {}

  /**
   * @brief Compute the total mass density in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] totalMassDensityKernelOp the function used to customize the kernel
   */
  GEOS_HOST_DEVICE inline void compute( localIndex const ei ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseMassDens = m_phaseMassDens[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseMassDens = m_dPhaseMassDens[ei][0];

    real64 & dTotalMassDens_dTemp = m_dTotalMassDens_dTemp[ei];
    real64 & dTotalMassDens_dT = m_dTotalMassDens[ei][Deriv::dT];

    // Call the base compute the compute the total mass density and derivatives
    return Base::compute( ei, [&]( localIndex const ip )
    {
      dTotalMassDens_dTemp += dPhaseVolFrac[ip][Deriv::dT] * phaseMassDens[ip] + phaseVolFrac[ip] * dPhaseMassDens[ip][Deriv::dT];
      dTotalMassDens_dT += dPhaseVolFrac[ip][Deriv::dT] * phaseMassDens[ip] + phaseVolFrac[ip] * dPhaseMassDens[ip][Deriv::dT];
    } );
  }

protected:
  // outputs
  arrayView1d< real64 > m_dTotalMassDens_dTemp;
};

/**
 * @class TotalMassDensityKernelFactory
 */
class TotalMassDensityKernelFactory
{
public:
  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   MultiFluidBase const & fluid )
  {
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&]( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&]( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of thermal accumulation and volume balance
 */
template< localIndex NUM_COMP, localIndex NUM_DOF >
class ElementBasedAssemblyKernel : public compositionalMultiphaseWellKernels::ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
{
public:
  using Base = compositionalMultiphaseWellKernels::ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >;
  using Base::m_dCompFrac_dCompDens;
  using Base::m_dofNumber;
  using Base::m_dPhaseCompFrac;
  using Base::m_dPhaseDens;
  using Base::m_dPhaseVolFrac;
  using Base::m_dPoro_dPres;
  using Base::m_elemGhostRank;
  using Base::m_localMatrix;
  using Base::m_localRhs;
  using Base::m_numPhases;
  using Base::m_phaseCompFrac;
  using Base::m_phaseCompFrac_n;
  using Base::m_phaseDens;
  using Base::m_phaseDens_n;
  using Base::m_phaseVolFrac;
  using Base::m_phaseVolFrac_n;
  using Base::m_porosity;
  using Base::m_porosity_n;
  using Base::m_rankOffset;
  using Base::m_volume;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( localIndex const numPhases,
                              globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              MultiFluidBase const & fluid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > const kernelFlags )
    : Base( numPhases, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs, kernelFlags ),
    m_phaseInternalEnergy_n( fluid.phaseInternalEnergy_n()),
    m_phaseInternalEnergy( fluid.phaseInternalEnergy()),
    m_dPhaseInternalEnergy( fluid.dPhaseInternalEnergy())
  {}

  struct StackVariables : public Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables()
      : Base::StackVariables()
    {}
    using Base::StackVariables::eqnRowIndices;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localJacobian;
    using Base::StackVariables::localResidual;
    using Base::StackVariables::localRow;
    using Base::StackVariables::volume;



  };
  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    Base::computeAccumulation( ei, stack, [&]( integer const ip, real64 const & phaseAmount, real64 const & phaseAmount_n, real64 const & dPhaseAmount_dP, real64 const (&dPhaseAmount_dC)[numComp] )
    {
      // We are in the loop over phases, ip provides the current phase index.
      // We have to do two things:
      //   1- Assemble the derivatives of the component mass balance equations with respect to temperature
      //   2- Assemble the phase-dependent part of the accumulation term of the energy equation

      real64 dPhaseInternalEnergy_dC[numComp]{};

      // construct the slices
      arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens = m_dPhaseDens[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac = m_phaseCompFrac[ei][0];
      arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac = m_dPhaseCompFrac[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy_n = m_phaseInternalEnergy_n[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy = m_phaseInternalEnergy[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseInternalEnergy = m_dPhaseInternalEnergy[ei][0];

      // Step 1: assemble the derivatives of the component mass balance equations with respect to temperature

      real64 const dPhaseAmount_dT = stack.volume * (dPhaseVolFrac[ip][Deriv::dT] * phaseDens[ip] + phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dT] );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        stack.localJacobian[ic][numDof-1] += dPhaseAmount_dT * phaseCompFrac[ip][ic]
                                             + phaseAmount * dPhaseCompFrac[ip][ic][Deriv::dT];
      }

      // Step 2: assemble the phase-dependent part of the accumulation term of the energy equation

      real64 const phaseEnergy = phaseAmount * phaseInternalEnergy[ip];
      real64 const phaseEnergy_n = phaseAmount_n * phaseInternalEnergy_n[ip];
      real64 const dPhaseEnergy_dP = dPhaseAmount_dP * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dP];
      real64 const dPhaseEnergy_dT = dPhaseAmount_dT * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dT];

      // local accumulation
      stack.localResidual[numEqn-1] += phaseEnergy - phaseEnergy_n;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[numEqn-1][0]        += dPhaseEnergy_dP;
      stack.localJacobian[numEqn-1][numDof-1] += dPhaseEnergy_dT;

      // derivatives w.r.t. component densities
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseInternalEnergy[ip], dPhaseInternalEnergy_dC, Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localJacobian[numEqn-1][jc + 1] += phaseInternalEnergy[ip] * dPhaseAmount_dC[jc]
                                                 + dPhaseInternalEnergy_dC[jc] * phaseAmount;
      }
    } );


  }

  /**
   * @brief Compute the local volume balance contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeVolumeBalance( localIndex const ei,
                             StackVariables & stack ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    Base::computeVolumeBalance( ei, stack, [&]( real64 const & oneMinusPhaseVolFraction )
    {
      GEOS_UNUSED_VAR( oneMinusPhaseVolFraction );

      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];

      for( integer ip = 0; ip < m_numPhases; ++ip )
      {
        stack.localJacobian[numEqn-2][numDof-1] -= dPhaseVolFrac[ip][Deriv::dT];
      }
    } );
  }

  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // Step 1: assemble the component mass balance equations and volume balance equations
    Base::complete( ei, stack );
    /* fix me tjb
       // Step 2: assemble the energy equation
       m_localRhs[stack.localRow + numEqn - 1] += stack.localResidual[numEqn - 1];
       m_localMatrix.template addToRow<serialAtomic>(stack.localRow + numEqn - 1,
                                                  stack.dofIndices,
                                                  stack.localJacobian[numEqn - 1],
                                                  numDof);
     */
  }

protected:
  /// View on derivative of porosity w.r.t temperature
  arrayView2d< real64 const > const m_dPoro_dTemp;

  /// Views on phase internal energy
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseInternalEnergy_n;
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseInternalEnergy;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseInternalEnergy;

};

/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:
  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( localIndex const numComps,
                   localIndex const numPhases,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   MultiFluidBase const & fluid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      localIndex constexpr NUM_COMP = NC();
      localIndex constexpr NUM_DOF = NC()+2;

      BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags::TotalMassEquation );

      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
      kernel( numPhases, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs, kernelFlags );
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template
      launch< POLICY, ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF > >( subRegion.size(), kernel );
    } );
  }
};
/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NC, integer NUM_DOF >
class FaceBasedAssemblyKernel : public compositionalMultiphaseWellKernels::FaceBasedAssemblyKernel< NC, NUM_DOF >
{
public:

  using Base  = compositionalMultiphaseWellKernels::FaceBasedAssemblyKernel< NC, NUM_DOF >;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;
  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using TAG = compositionalMultiphaseWellKernels::ElemTag;


  /// Compile time value for the number of components
  static constexpr integer numComp = NC;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;


  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  FaceBasedAssemblyKernel( real64 const dt,
                           globalIndex const rankOffset,
                           string const wellDofKey,
                           WellControls const & wellControls,
                           ElementSubRegionBase const & subRegion,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs,
                           BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > kernelFlags )
    : Base( dt , rankOffset , subRegion.getReference< array1d< globalIndex > >( wellDofKey )  
      , subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString()) 
      , subRegion.getField< fields::well::mixtureConnectionRate >() 
      , subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >() 
      , localMatrix
      , localRhs
      , kernelFlags.isSet( isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags::TotalMassEquation )
      , wellControls.isProducer()
      , wellControls.getInjectionStream() )
  { }


  GEOS_HOST_DEVICE
  inline
  void
  computeExit( real64 const & dt,
               real64 const ( &compFlux )[NC],
               real64 const ( &dCompFlux_dRate )[NC],
               real64 const ( &dCompFlux_dPresUp )[NC],
               real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
               real64 ( & oneSidedFlux )[NC],
               real64 ( & oneSidedFluxJacobian_dRate )[NC][1],
               real64 ( & oneSidedFluxJacobian_dPresCompUp )[NC][NC + 1] ) const
  {
    for( integer ic = 0; ic < NC; ++ic )
    {
      oneSidedFlux[ic] = -dt * compFlux[ic];

      // derivative with respect to rate
      oneSidedFluxJacobian_dRate[ic][0] = -dt * dCompFlux_dRate[ic];

      // derivative with respect to upstream pressure
      oneSidedFluxJacobian_dPresCompUp[ic][0] = -dt * dCompFlux_dPresUp[ic];

      // derivatives with respect to upstream component densities
      for( integer jdof = 0; jdof < NC; ++jdof )
      {
        oneSidedFluxJacobian_dPresCompUp[ic][jdof+1] = -dt * dCompFlux_dCompDensUp[ic][jdof];
      }
    }
  }

  GEOS_HOST_DEVICE
  inline
  void
  compute( real64 const & dt,
           real64 const ( &compFlux )[NC],
           real64 const ( &dCompFlux_dRate )[NC],
           real64 const ( &dCompFlux_dPresUp )[NC],
           real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
           real64 ( & localFlux )[2*NC],
           real64 ( & localFluxJacobian_dRate )[2*NC][1],
           real64 ( & localFluxJacobian_dPresCompUp )[2*NC][NC + 1] ) const
  {
    // flux terms
    for( integer ic = 0; ic < NC; ++ic )
    {
      localFlux[TAG::NEXT *NC+ic]    = dt * compFlux[ic];
      localFlux[TAG::CURRENT *NC+ic] = -dt * compFlux[ic];

      // derivative with respect to rate
      localFluxJacobian_dRate[TAG::NEXT *NC+ic][0]    = dt * dCompFlux_dRate[ic];
      localFluxJacobian_dRate[TAG::CURRENT *NC+ic][0] = -dt * dCompFlux_dRate[ic];

      // derivative with respect to upstream pressure
      localFluxJacobian_dPresCompUp[TAG::NEXT *NC+ic][0]    = dt * dCompFlux_dPresUp[ic];
      localFluxJacobian_dPresCompUp[TAG::CURRENT *NC+ic][0] = -dt * dCompFlux_dPresUp[ic];

      // derivatives with respect to upstream component densities
      for( integer jdof = 0; jdof < NC; ++jdof )
      {
        localFluxJacobian_dPresCompUp[TAG::NEXT *NC+ic][jdof+1]    =  dt * dCompFlux_dCompDensUp[ic][jdof];
        localFluxJacobian_dPresCompUp[TAG::CURRENT *NC+ic][jdof+1] = -dt * dCompFlux_dCompDensUp[ic][jdof];
      }
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] ie the element index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void computeFlux( localIndex const iwelem,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {

    using namespace compositionalMultiphaseUtilities;

    // create local work arrays
    real64 compFracUp[NC]{};
    real64 dCompFrac_dCompDensUp[NC][NC]{};

    real64 compFlux[NC]{};
    real64 dCompFlux_dRate[NC]{};
    real64 dCompFlux_dPresUp[NC]{};
    real64 dCompFlux_dCompDensUp[NC][NC]{};

    // Step 1) decide the upwind well element

    /*  currentConnRate < 0 flow from iwelem to iwelemNext
     *  currentConnRate > 0 flow from iwelemNext to iwelem
     *  With this convention, currentConnRate < 0 at the last connection for a producer
     *                        currentConnRate > 0 at the last connection for a injector
     */

    localIndex const iwelemNext = m_nextWellElemIndex[iwelem];
    real64 const currentConnRate = m_connRate[iwelem];
    localIndex iwelemUp = -1;

    if( iwelemNext < 0 && !m_isProducer )  // exit connection, injector
    {
      // we still need to define iwelemUp for Jacobian assembly
      iwelemUp = iwelem;

      // just copy the injection stream into compFrac
      for( integer ic = 0; ic < NC; ++ic )
      {
        compFracUp[ic] = m_injection[ic];
        for( integer jc = 0; jc < NC; ++jc )
        {
          dCompFrac_dCompDensUp[ic][jc] = 0.0;
        }
      }
    }
    else
    {
      // first set iwelemUp to the upstream cell
      if( ( iwelemNext < 0 && m_isProducer )  // exit connection, producer
          || currentConnRate < 0 ) // not an exit connection, iwelem is upstream
      {
        iwelemUp = iwelem;
      }
      else // not an exit connection, iwelemNext is upstream
      {
        iwelemUp = iwelemNext;
      }

      // copy the vars of iwelemUp into compFrac
      for( integer ic = 0; ic < NC; ++ic )
      {
        compFracUp[ic] = m_wellElemCompFrac[iwelemUp][ic];
        for( integer jc = 0; jc < NC; ++jc )
        {
          dCompFrac_dCompDensUp[ic][jc] = m_dWellElemCompFrac_dCompDens[iwelemUp][ic][jc];
        }
      }
    }

    // Step 2) compute upstream transport coefficient

    for( integer ic = 0; ic < NC; ++ic )
    {
      compFlux[ic]          = compFracUp[ic] * currentConnRate;
      dCompFlux_dRate[ic]   = compFracUp[ic];
      dCompFlux_dPresUp[ic] = 0.0; // none of these quantities depend on pressure
      for( integer jc = 0; jc < NC; ++jc )
      {
        dCompFlux_dCompDensUp[ic][jc] = dCompFrac_dCompDensUp[ic][jc] * currentConnRate;
      }
    }

    globalIndex const offsetUp = m_wellElemDofNumber[iwelemUp];
    globalIndex const offsetCurrent = m_wellElemDofNumber[iwelem];

    if( iwelemNext < 0 )  // exit connection
    {
      // for this case, we only need NC mass conservation equations
      // so we do not use the arrays initialized before the loop
      real64 oneSidedFlux[NC]{};
      real64 oneSidedFluxJacobian_dRate[NC][1]{};
      real64 oneSidedFluxJacobian_dPresCompUp[NC][NC+1]{};

      computeExit ( m_dt,
                    compFlux,
                    dCompFlux_dRate,
                    dCompFlux_dPresUp,
                    dCompFlux_dCompDensUp,
                    oneSidedFlux,
                    oneSidedFluxJacobian_dRate,
                    oneSidedFluxJacobian_dPresCompUp );


      globalIndex oneSidedEqnRowIndices[NC]{};
      globalIndex oneSidedDofColIndices_dPresCompUp[NC+1]{};
      globalIndex oneSidedDofColIndices_dRate = 0;

      // jacobian indices
      for( integer ic = 0; ic < NC; ++ic )
      {
        // mass balance equations for all components
        oneSidedEqnRowIndices[ic] = offsetUp + ROFFSET::MASSBAL + ic - m_rankOffset;
      }

      // in the dof ordering used in this class, there are 1 pressure dofs
      // and NC compDens dofs before the rate dof in this block
      localIndex const dRateColOffset = COFFSET::DCOMP + NC;
      oneSidedDofColIndices_dRate = offsetCurrent + dRateColOffset;

      for( integer jdof = 0; jdof < NC+1; ++jdof )
      {
        // dofs are the **upstream** pressure and component densities
        oneSidedDofColIndices_dPresCompUp[jdof] = offsetUp + COFFSET::DPRES + jdof;
      }

      if( m_useTotalMassEquation > 0 )
      {
        // Apply equation/variable change transformation(s)
        real64 work[NC + 1]{};
        shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, 1, oneSidedFluxJacobian_dRate, work );
        shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC + 1, oneSidedFluxJacobian_dPresCompUp, work );
        shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, oneSidedFlux );
      }

      for( integer i = 0; i < NC; ++i )
      {
        if( oneSidedEqnRowIndices[i] >= 0 && oneSidedEqnRowIndices[i] < m_localMatrix.numRows() )
        {
          m_localMatrix.addToRow< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                          &oneSidedDofColIndices_dRate,
                                                          oneSidedFluxJacobian_dRate[i],
                                                          1 );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                                              oneSidedDofColIndices_dPresCompUp,
                                                                              oneSidedFluxJacobian_dPresCompUp[i],
                                                                              NC+1 );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[oneSidedEqnRowIndices[i]], oneSidedFlux[i] );
        }
      }
    }
    else // not an exit connection
    {
      real64 localFlux[2*NC]{};
      real64 localFluxJacobian_dRate[2*NC][1]{};
      real64 localFluxJacobian_dPresCompUp[2*NC][NC+1]{};

      compute( m_dt,
               compFlux,
               dCompFlux_dRate,
               dCompFlux_dPresUp,
               dCompFlux_dCompDensUp,
               localFlux,
               localFluxJacobian_dRate,
               localFluxJacobian_dPresCompUp );


      globalIndex eqnRowIndices[2*NC]{};
      globalIndex dofColIndices_dPresCompUp[NC+1]{};
      globalIndex dofColIndices_dRate = 0;

      globalIndex const offsetNext = m_wellElemDofNumber[iwelemNext];

      // jacobian indices
      for( integer ic = 0; ic < NC; ++ic )
      {
        // mass balance equations for all components
        eqnRowIndices[TAG::NEXT *NC+ic]    = offsetNext + ROFFSET::MASSBAL + ic - m_rankOffset;
        eqnRowIndices[TAG::CURRENT *NC+ic] = offsetCurrent + ROFFSET::MASSBAL + ic - m_rankOffset;
      }

      // in the dof ordering used in this class, there are 1 pressure dofs
      // and NC compDens dofs before the rate dof in this block
      localIndex const dRateColOffset = COFFSET::DCOMP + NC;
      dofColIndices_dRate = offsetCurrent + dRateColOffset;

      for( integer jdof = 0; jdof < NC+1; ++jdof )
      {
        // dofs are the **upstream** pressure and component densities
        dofColIndices_dPresCompUp[jdof] = offsetUp + COFFSET::DPRES + jdof;
      }

      if( m_useTotalMassEquation > 0 )
      {
        // Apply equation/variable change transformation(s)
        real64 work[NC + 1]{};
        shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC, 1, 2, localFluxJacobian_dRate, work );
        shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC, NC + 1, 2, localFluxJacobian_dPresCompUp, work );
        shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( NC, NC, 2, localFlux );
      }

      for( integer i = 0; i < 2*NC; ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < m_localMatrix.numRows() )
        {
          m_localMatrix.addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                          &dofColIndices_dRate,
                                                          localFluxJacobian_dRate[i],
                                                          1 );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                              dofColIndices_dPresCompUp,
                                                                              localFluxJacobian_dPresCompUp[i],
                                                                              NC+1 );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[eqnRowIndices[i]], localFlux[i] );
        }
      }
    }
  }


  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElements the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElements,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const ie )
    {
      //typename KERNEL_TYPE::StackVariables stack( kernelComponent.stencilSize( iconn ),
      //                                            kernelComponent.numPointsInFlux( iconn ) );

      //kernelComponent.setup( iconn, stack );
      kernelComponent.computeFlux( ie );
      //kernelComponent.complete( iconn, stack );
    } );
  }

protected:
  /// Time step size
  real64 const m_dt;
  /// Rank offset for calculating row/col Jacobian indices
  integer const m_rankOffset;

  /// Reference to the degree-of-freedom numbers
  arrayView1d< globalIndex const > const m_wellElemDofNumber;
  /// Next element index, needed since iterating over element nodes, not edges
  arrayView1d< localIndex const > const m_nextWellElemIndex;

  /// Connection rate
  arrayView1d< real64 const > const m_connRate;

  /// Element component fraction
  arrayView2d< real64 const, compflow::USD_COMP > const m_wellElemCompFrac;
  /// Element component fraction derivatives
  arrayView3d< real64 const, compflow::USD_COMP_DC > const m_dWellElemCompFrac_dCompDens;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  /// Kernel option flag
  integer const m_useTotalMassEquation;

  /// Well type
  bool const m_isProducer;

  /// Injection stream composition
  arrayView1d< real64 const > const m_injection;


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
   * @param[in] numComps the number of fluid components
   * @param[in] dt time step size
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] useTotalMassEquation flag specifying whether to replace one component bal eqn with total mass eqn
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] wellControls object holding well control/constraint information
   * @param[in] subregion well subregion
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComps,
                   real64 const dt,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const dofKey,
                   WellControls const & wellControls,
                   ElementSubRegionBase const & subRegion,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 2;

      BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags::TotalMassEquation );


      using kernelType = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF >;


      kernelType kernel( dt, rankOffset, dofKey, wellControls, subRegion, localMatrix, localRhs, kernelFlags );
      kernelType::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }
};

};   // end namespace thermalCompositionalMultiphaseWellKernels

} // end namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALCOMPOSITIONALMULTIPHASEWELLKERNELS_HPP
