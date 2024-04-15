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
 * @file PerforationFluxKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_PERFORATIONFLUXLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_PERFORATIONFLUXLKERNELS_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/KernelLaunchSelectors.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
//#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
//#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellTags.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"


namespace geos
{

struct NoOpStuct
{
  NoOpStuct(){}
};

namespace isothermalPerforationFluxKernels
{

using namespace constitutive;



/******************************** PerforationFluxKernel ********************************/

template< integer NC, integer NP, integer IS_THERMAL >
class PerforationFluxKernel
{
public:
  /// Compile time value for the number of components
  static constexpr integer numComp = NC;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NP;

  /// Compile time value for thermal option
  static constexpr integer isThermal = IS_THERMAL;

  using TAG = wellTags::SubRegionTag;

  using CompFlowAccessors =
    StencilAccessors< fields::flow::pressure,
                      fields::flow::phaseVolumeFraction,
                      fields::flow::dPhaseVolumeFraction,
                      fields::flow::dGlobalCompFraction_dGlobalCompDensity >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseViscosity,
                              fields::multifluid::dPhaseViscosity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

  using RelPermAccessors =
    StencilMaterialAccessors< RelativePermeabilityBase,
                              fields::relperm::phaseRelPerm,
                              fields::relperm::dPhaseRelPerm_dPhaseVolFraction >;


  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  PerforationFluxKernel ( PerforationData * const perforationData,
                          ElementSubRegionBase const & subRegion,
                          CompFlowAccessors const & compFlowAccessors,
                          MultiFluidAccessors const & multiFluidAccessors,
                          RelPermAccessors const & relPermAccessors,
                          bool const disableReservoirToWellFlow ):
    m_resPres( compFlowAccessors.get( fields::flow::pressure {} )),
    m_resPhaseVolFrac( compFlowAccessors.get( fields::flow::phaseVolumeFraction {} )),
    m_dResPhaseVolFrac(compFlowAccessors.get( fields::flow::dPhaseVolumeFraction {} )),
    m_dResCompFrac_dCompDens( compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity {} )),
    m_resPhaseDens( multiFluidAccessors.get( fields::multifluid::phaseDensity {} )),
    m_dResPhaseDens( multiFluidAccessors.get( fields::multifluid::dPhaseDensity {} )),
    m_resPhaseVisc( multiFluidAccessors.get( fields::multifluid::phaseViscosity {} )),
    m_dResPhaseVisc( multiFluidAccessors.get( fields::multifluid::dPhaseViscosity {} )),
    m_resPhaseCompFrac( multiFluidAccessors.get( fields::multifluid::phaseCompFraction {} )),
    m_dResPhaseCompFrac( multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction {} )),
    m_resPhaseRelPerm( relPermAccessors.get( fields::relperm::phaseRelPerm {} )),
    m_dResPhaseRelPerm_dPhaseVolFrac( relPermAccessors.get( fields::relperm::dPhaseRelPerm_dPhaseVolFraction {} )),
    m_wellElemGravCoef( subRegion.getField< fields::well::gravityCoefficient >()),
    m_wellElemPres( subRegion.getField< fields::well::pressure >()),
    m_wellElemCompDens( subRegion.getField< fields::well::globalCompDensity >()),
    m_wellElemTotalMassDens( subRegion.getField< fields::well::totalMassDensity >()),
    m_dWellElemTotalMassDens( subRegion.getField< fields::well::dTotalMassDensity >()),
    m_dWellElemTotalMassDens_dPres( subRegion.getField< fields::well::dTotalMassDensity_dPressure >()),
    m_dWellElemTotalMassDens_dCompDens( subRegion.getField< fields::well::dTotalMassDensity_dGlobalCompDensity >()),
    m_wellElemCompFrac( subRegion.getField< fields::well::globalCompFraction >()),
    m_dWellElemCompFrac_dCompDens( subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >()),
    m_perfGravCoef( perforationData->getField< fields::well::gravityCoefficient >()),
    m_perfWellElemIndex( perforationData->getField< fields::perforation::wellElementIndex >()),
    m_perfTrans( perforationData->getField< fields::perforation::wellTransmissibility >()),
    m_resElementRegion( perforationData->getField< fields::perforation::reservoirElementRegion >()),
    m_resElementSubRegion( perforationData->getField< fields::perforation::reservoirElementSubRegion >()),
    m_resElementIndex( perforationData->getField< fields::perforation::reservoirElementIndex >()),
    m_compPerfRate( perforationData->getField< fields::well::compPerforationRate >()),
    m_dCompPerfRate( perforationData->getField< fields::well::dCompPerforationRate >()),
    m_dCompPerfRate_dPres( perforationData->getField< fields::well::dCompPerforationRate_dPres >()),
    m_dCompPerfRate_dComp( perforationData->getField< fields::well::dCompPerforationRate_dComp >()),
    m_disableReservoirToWellFlow( disableReservoirToWellFlow )
  {} 

struct StackVariables
  {
public:
    /**
     * @brief Constructor for the stack variables
     */

    GEOS_HOST_DEVICE
    StackVariables() {}
   
  };

  template< typename STRUCT = NoOpFunc >  
  GEOS_HOST_DEVICE
  inline
  void
  computeFlux( localIndex const iperf, STRUCT && stack= NoOpFunc{}  ) const
  {
   // get the index of the reservoir elem
      localIndex const er  = m_resElementRegion[iperf];
      localIndex const esr = m_resElementSubRegion[iperf];
      localIndex const ei  = m_resElementIndex[iperf];

      // get the index of the well elem
      localIndex const iwelem = m_perfWellElemIndex[iperf];
 
    using Deriv = multifluid::DerivativeOffset;
    using CP_Deriv = multifluid::DerivativeOffsetC< NC, IS_THERMAL >;

    // local working variables and arrays
    real64 pres[2]{};
    real64 dPres_dP[2]{};
    real64 dPres_dC[2][NC]{};
    real64 dFlux_dP[2]{};
    real64 dFlux_dC[2][NC]{};
    real64 dMult_dP[2]{};
    real64 dMult_dC[2][NC]{};
    real64 dPotDiff_dP[2]{};
    real64 dPotDiff_dC[2][NC]{};
    real64 multiplier[2]{};

    real64 dResTotalMob_dC[NC]{};
    real64 dDens_dC[NC]{};
    real64 dVisc_dC[NC]{};
    real64 dRelPerm_dC[NC]{};
    real64 dMob_dC[NC]{};
    real64 dCompFrac_dCompDens[NC]{};

    // local working variables - compact
    // All derivative quantiites generated are stored in arrays using CP_Deriv offsets
    // The input well/reservoir quantites use the Deriv offsets
    // The arrays using the deriv offsets have extra column for dT in isothermal cases

    real64 dPres[2][CP_Deriv::nDer]{};
    real64 dFlux[2][CP_Deriv::nDer]{};
    real64 dMob[CP_Deriv::nDer]{};
    real64 dPotDiff[2][CP_Deriv::nDer]{};
    real64 dCompFrac[CP_Deriv::nDer]{};

    // Step 1: reset the perforation rates
    // tjb - remove when safe
    for( integer ic = 0; ic < NC; ++ic )
    {
      m_compPerfRate[iperf][ic] = 0.0;
      for( integer ke = 0; ke < 2; ++ke )
      {
        m_dCompPerfRate_dPres[iperf][ke][ic] = 0.0;
        for( integer jc = 0; jc < NC; ++jc )
        {
          m_dCompPerfRate_dComp[iperf][ke][ic][jc] = 0.0;
        }
      }
    }
    for( integer ic = 0; ic < NC; ++ic )
    {
      m_compPerfRate[iperf][ic] = 0.0;
      for( integer ke = 0; ke < 2; ++ke )
      {
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          m_dCompPerfRate[iperf][ke][ic][jc] = 0.0;
        }
      }
    }
    if constexpr (IS_THERMAL ) 
    {
      // initialize outputs
      stack.m_energyPerfFlux[iperf]=0;
      for( integer ke = 0; ke < 2; ++ke )
      {
        for ( integer i = 0; i < CP_Deriv::nDer; ++i )
        {
          stack.m_dEnergyPerfFlux[iperf][ke][i]=0;
        }
      }
    }

    // Step 2: copy the variables from the reservoir and well element

    // a) get reservoir variables

    pres[TAG::RES] = m_resPres[er][esr][ei];
    dPres_dP[TAG::RES] = 1.0;
    dPres[TAG::RES][CP_Deriv::dP] = 1.0;
    multiplier[TAG::RES] = 1.0;

    // Here in the absence of a buoyancy term we assume that the reservoir cell is perforated at its center
    // TODO: add a buoyancy term for the reservoir side here


    // b) get well variables

    pres[TAG::WELL] = m_wellElemPres[iwelem];
    dPres_dP[TAG::WELL] = 1.0;
    dPres[TAG::WELL][CP_Deriv::dP] = 1.0;
    multiplier[TAG::WELL] = -1.0;

    real64 const gravD = ( m_perfGravCoef[iperf] - m_wellElemGravCoef[iwelem] );

    pres[TAG::WELL] +=  m_wellElemTotalMassDens[iwelem] * gravD;
    // Note RHS uses CP_Deriv while LHS uses Deriv !!!
    dPres_dP[TAG::WELL] +=  m_dWellElemTotalMassDens_dPres[iwelem] * gravD;
    dPres[TAG::WELL][CP_Deriv::dP] +=  m_dWellElemTotalMassDens_dPres[iwelem] * gravD;
    if constexpr ( IS_THERMAL )
    {
      dPres[TAG::WELL][CP_Deriv::dT] += m_dWellElemTotalMassDens[iwelem][Deriv::dT] * gravD;
    }
    for( integer ic = 0; ic < NC; ++ic )
    {
      dPres_dC[TAG::WELL][ic] += m_dWellElemTotalMassDens_dCompDens[iwelem][ic] * gravD;
      dPres[TAG::WELL][CP_Deriv::dC+ic] += m_dWellElemTotalMassDens[iwelem][Deriv::dC+ic] * gravD;
    }



    // Step 3: compute potential difference

    real64 potDiff = 0.0;
    for( integer i = 0; i < 2; ++i )
    {
      potDiff += multiplier[i] * m_perfTrans[iperf] * pres[i];
      dPotDiff_dP[i] += multiplier[i] * m_perfTrans[iperf] * dPres_dP[i];

      for( integer ic = 0; ic < NC; ++ic )
      {
        dPotDiff_dC[i][ic] += multiplier[i] * m_perfTrans[iperf] * dPres_dC[i][ic];
      }
      // LHS & RHS both use CP_Deriv
      for( integer ic = 0; ic < CP_Deriv::nDer; ++ic )
      {
        dPotDiff[i][ic] += multiplier[i] * m_perfTrans[iperf] * dPres[i][ic];
      }
    }


    // Step 4: upwinding based on the flow direction

    real64 flux = 0.0;
    if( potDiff >= 0 )  // ** reservoir cell is upstream **
    {

      // loop over phases, compute and upwind phase flux
      // and sum contributions to each component's perforation rate
      for( integer ip = 0; ip < NP; ++ip )
      {

        // skip the rest of the calculation if the phase is absent
        // or if crossflow is disabled for injectors
        bool const phaseExists = (m_resPhaseVolFrac[er][esr][ei][ip] > 0);
        if( !phaseExists || m_disableReservoirToWellFlow )
        {
          continue;
        }

        // here, we have to recompute the reservoir phase mobility (not including density)

        // density
        real64 const resDens = m_resPhaseDens[er][esr][ei][0][ip];
        real64 dDens[CP_Deriv::nDer]{};

        real64 const dResDens_dP  = m_dResPhaseDens[er][esr][ei][0][ip][Deriv::dP];
        dDens[CP_Deriv::dP]  = m_dResPhaseDens[er][esr][ei][0][ip][Deriv::dP];
        if constexpr ( IS_THERMAL )
        {
          dDens[CP_Deriv::dT]  = m_dResPhaseDens[er][esr][ei][0][ip][Deriv::dT];
        }
        applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                        m_dResPhaseDens[er][esr][ei][0][ip],
                        dDens_dC,
                        Deriv::dC );
        applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                        m_dResPhaseDens[er][esr][ei][0][ip],
                        &dDens[CP_Deriv::dC],
                        Deriv::dC );
        // viscosity
        real64 const resVisc = m_resPhaseVisc[er][esr][ei][0][ip];
        real64 dVisc[CP_Deriv::nDer]{};
        real64 const dResVisc_dP  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dP];
        dVisc[CP_Deriv::dP]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dP];
        if constexpr ( IS_THERMAL )
        {
          dVisc[CP_Deriv::dT]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dT];
        }
        applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                        m_dResPhaseVisc[er][esr][ei][0][ip],
                        dVisc_dC,
                        Deriv::dC );
        applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                        m_dResPhaseVisc[er][esr][ei][0][ip],
                        &dVisc[CP_Deriv::dC],
                        Deriv::dC );

        // relative permeability
        real64 const resRelPerm = m_resPhaseRelPerm[er][esr][ei][0][ip];
        real64 dRelPerm[CP_Deriv::nDer]{};
        real64 dResRelPerm_dP = 0.0;
        for( integer jc = 0; jc < NC; ++jc )
        {
          dRelPerm_dC[jc] = 0;
        }
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          dRelPerm[jc]=0;
        }
        for( integer jp = 0; jp < NP; ++jp )
        {
          real64 const dResRelPerm_dS = m_dResPhaseRelPerm_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dResRelPerm_dP += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dP];
          dRelPerm[CP_Deriv::dP] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dP];
          if constexpr ( IS_THERMAL )
          {
            dRelPerm[CP_Deriv::dT] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dT];
          }
          for( integer jc = 0; jc < NC; ++jc )
          {
            dRelPerm_dC[jc] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
            dRelPerm[CP_Deriv::dC+jc] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
          }
        }

        // compute the reservoir phase mobility, including phase density
        real64 const resPhaseMob = resDens * resRelPerm / resVisc;
        real64 const dResPhaseMob_dPres = dResRelPerm_dP * resDens / resVisc
                                          + resPhaseMob * (dResDens_dP / resDens - dResVisc_dP / resVisc);
        for( integer jc = 0; jc < NC; ++jc )
        {
          dMob_dC[jc] = dRelPerm_dC[jc] * resDens / resVisc
                        + resPhaseMob * (dDens_dC[jc] / resDens - dVisc_dC[jc] / resVisc);
        }
        // Handles all dependencies
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          dMob[jc] = dRelPerm[jc] * resDens / resVisc
                    + resPhaseMob * (dDens[jc] / resDens - dVisc[jc] / resVisc);
        }

        // compute the phase flux and derivatives using upstream cell mobility
        flux = resPhaseMob * potDiff;
        dFlux_dP[TAG::RES]  = dResPhaseMob_dPres * potDiff + resPhaseMob * dPotDiff_dP[TAG::RES];
        dFlux_dP[TAG::WELL] = resPhaseMob *  dPotDiff_dP[TAG::WELL];

        for( integer ic = 0; ic < NC; ++ic )
        {
          dFlux_dC[TAG::RES][ic] = dMob_dC[ic] * potDiff + resPhaseMob * dPotDiff_dC[TAG::RES][ic];
          dFlux_dC[TAG::WELL][ic] = resPhaseMob * dPotDiff_dC[TAG::WELL][ic];
        }

        // Handles all dependencies
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          dFlux[TAG::RES][jc]  = dMob[jc] * potDiff + resPhaseMob * dPotDiff[TAG::RES][jc];
          dFlux[TAG::WELL][jc] = resPhaseMob * dPotDiff[TAG::WELL][jc];
        }
        // increment component fluxes
        for( integer ic = 0; ic < NC; ++ic )
        {
          m_compPerfRate[iperf][ic] += flux *  m_resPhaseCompFrac[er][esr][ei][0][ip][ic];

          m_dCompPerfRate_dPres[iperf][TAG::RES][ic]  +=  m_resPhaseCompFrac[er][esr][ei][0][ip][ic] * dFlux_dP[TAG::RES];
          m_dCompPerfRate_dPres[iperf][TAG::RES][ic]  += m_dResPhaseCompFrac[er][esr][ei][0][ip][ic][Deriv::dP] * flux;
          m_dCompPerfRate_dPres[iperf][TAG::WELL][ic] += m_resPhaseCompFrac[er][esr][ei][0][ip][ic] * dFlux_dP[TAG::WELL];

          applyChainRule( NC,
                          m_dResCompFrac_dCompDens[er][esr][ei],
                          m_dResPhaseCompFrac[er][esr][ei][0][ip][ic],
                          dCompFrac_dCompDens,
                          Deriv::dC );

          for( integer jc = 0; jc < NC; ++jc )
          {
            m_dCompPerfRate_dComp[iperf][TAG::RES][ic][jc]  += dFlux_dC[TAG::RES][jc] *  m_resPhaseCompFrac[er][esr][ei][0][ip][ic];
            m_dCompPerfRate_dComp[iperf][TAG::RES][ic][jc]  += flux * dCompFrac_dCompDens[jc];
            m_dCompPerfRate_dComp[iperf][TAG::WELL][ic][jc] += dFlux_dC[TAG::WELL][jc] *  m_resPhaseCompFrac[er][esr][ei][0][ip][ic];
          }
        }
        // increment component fluxes
        for( integer ic = 0; ic < NC; ++ic )
        {
          // Note this needs to be uncommented out
          // m_compPerfRate[iperf][ic] += flux *  m_resPhaseCompFrac[er][esr][ei][0][ip][ic];
          dCompFrac[CP_Deriv::dP] = m_dResPhaseCompFrac[er][esr][ei][0][ip][ic][Deriv::dP];
          if constexpr (IS_THERMAL)
          {
            dCompFrac[CP_Deriv::dT] = m_dResPhaseCompFrac[er][esr][ei][0][ip][ic][Deriv::dT];
          }

          applyChainRule( NC,
                          m_dResCompFrac_dCompDens[er][esr][ei],
                          m_dResPhaseCompFrac[er][esr][ei][0][ip][ic],
                          &dCompFrac[CP_Deriv::dC],
                          Deriv::dC );

          for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
          {
            m_dCompPerfRate[iperf][TAG::RES][ic][jc]  += dFlux[TAG::RES][jc] *  m_resPhaseCompFrac[er][esr][ei][0][ip][ic];
            m_dCompPerfRate[iperf][TAG::RES][ic][jc]  += flux * dCompFrac[jc];
            m_dCompPerfRate[iperf][TAG::WELL][ic][jc] += dFlux[TAG::WELL][jc] *  m_resPhaseCompFrac[er][esr][ei][0][ip][ic];
          }
        }

        if constexpr ( IS_THERMAL ) 
        {
          real64 const res_enthalpy = stack.m_resPhaseEnthalpy[er][esr][ei][0][ip];

          stack.m_energyPerfFlux[iperf] += flux * res_enthalpy;
          // energy equation derivatives WRT res P & T
          stack.m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dP] += dFlux[TAG::RES][CP_Deriv::dP] * res_enthalpy + 
                                                      flux * stack.m_dResPhaseEnthalpy[er][esr][ei][0][ip][Deriv::dP];
          stack.m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dT] += dFlux[TAG::RES][CP_Deriv::dT] * res_enthalpy + 
                                                      flux * stack.m_dResPhaseEnthalpy[er][esr][ei][0][ip][Deriv::dT];
          // energy equation derivatives WRT well P
          stack.m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dP] += dFlux[TAG::WELL][CP_Deriv::dP] * res_enthalpy ;
          // stack.m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dT] += dFlux[TAG::WELL][CP_Deriv::dT] * res_enthalpy ;   

          // energy equation derivatives WRT reservoir dens
          real64 dProp_dC[numComp]{};
          applyChainRule( NC,
                          m_dResCompFrac_dCompDens[er][esr][ei],
                          stack.m_dResPhaseEnthalpy[er][esr][ei][0][ip],
                          dProp_dC,
                          Deriv::dC );   

          for( integer jc = 0; jc < NC; ++jc )
          {
            stack.m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dC+jc] += flux * dProp_dC[jc];
          }                                      
        }
      }  // end resevoir is upstream phase loop

      // tjb- remove when safe
      for( integer ic = 0; ic < NC; ic++ )
      {
        assert( fabs(  m_dCompPerfRate[iperf][TAG::RES][ic][CP_Deriv::dP] -m_dCompPerfRate_dPres[iperf][TAG::RES][ic] ) < FLT_EPSILON );
        assert( fabs(  m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dP] -m_dCompPerfRate_dPres[iperf][TAG::WELL][ic] ) < FLT_EPSILON );
        for( integer jc = 0; jc < NC; ++jc )
        {
          assert( fabs(  m_dCompPerfRate[iperf][TAG::RES][ic][CP_Deriv::dC+jc]  -m_dCompPerfRate_dComp[iperf][TAG::RES][ic][jc] ) < FLT_EPSILON );
          assert( fabs(  m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dC+jc]  -m_dCompPerfRate_dComp[iperf][TAG::WELL][ic][jc] ) < FLT_EPSILON );
        }
      }
    }
    else // ** well is upstream **
    {

      real64 resTotalMob     = 0.0;
      real64 dResTotalMob_dP = 0.0;

      // we re-compute here the total mass (when useMass == 1) or molar (when useMass == 0) density
      real64 wellElemTotalDens = 0;
      for( integer ic = 0; ic < NC; ++ic )
      {
        wellElemTotalDens += m_wellElemCompDens[iwelem][ic];
      }

      // first, compute the reservoir total mobility (excluding phase density)
      for( integer ip = 0; ip < NP; ++ip )
      {

        // skip the rest of the calculation if the phase is absent
        bool const phaseExists = (m_resPhaseVolFrac[er][esr][ei][ip] > 0);
        if( !phaseExists )
        {
          continue;
        }

        // viscosity
        real64 const resVisc = m_resPhaseVisc[er][esr][ei][0][ip];
        real64 dVisc[CP_Deriv::nDer]{};
        real64 const dResVisc_dP  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dP];
        dVisc[CP_Deriv::dP]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dP];
        if constexpr ( IS_THERMAL )
        {
          dVisc[CP_Deriv::dT]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dT];
        }
        applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                        m_dResPhaseVisc[er][esr][ei][0][ip],
                        dVisc_dC,
                        Deriv::dC );
        applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                        m_dResPhaseVisc[er][esr][ei][0][ip],
                        &dVisc[CP_Deriv::dC],
                        Deriv::dC );


        // relative permeability
        real64 const resRelPerm = m_resPhaseRelPerm[er][esr][ei][0][ip];
        real64 dRelPerm[CP_Deriv::nDer]{};
        real64 dResRelPerm_dP = 0.0;
        for( integer jc = 0; jc < NC; ++jc )
        {
          dRelPerm_dC[jc] = 0;
        }
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          dRelPerm[jc]=0;
        }
        for( integer jp = 0; jp < NP; ++jp )
        {
          real64 const dResRelPerm_dS = m_dResPhaseRelPerm_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dResRelPerm_dP += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dP];
          dRelPerm[CP_Deriv::dP] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dP];
          if constexpr ( IS_THERMAL )
          {
            dRelPerm[CP_Deriv::dT] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dT];
          }
          for( integer jc = 0; jc < NC; ++jc )
          {
            dRelPerm_dC[jc] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
            dRelPerm[CP_Deriv::dC+jc] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
          }
        }
        // increment total mobility
        resTotalMob     += resRelPerm / resVisc;

        dResTotalMob_dP += ( dResRelPerm_dP * resVisc - resRelPerm * dResVisc_dP )
                          / ( resVisc * resVisc );
        for( integer ic = 0; ic < NC; ++ic )
        {
          dResTotalMob_dC[ic] += ( dRelPerm_dC[ic] * resVisc - resRelPerm * dVisc_dC[ic] )
                                / ( resVisc * resVisc );
        }
        // Handles all dependencies
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          dMob[jc] += (dRelPerm[jc] *resVisc -  resRelPerm * dVisc[jc] )
                      / ( resVisc * resVisc);
        }
      } // end well is upstream phase loop

      // compute a potdiff multiplier = wellElemTotalDens * resTotalMob
      // wellElemTotalDens is a mass density if useMass == 1 and a molar density otherwise
      real64 const mult   = wellElemTotalDens * resTotalMob;
      dMult_dP[TAG::RES]  = wellElemTotalDens * dResTotalMob_dP;
      dMult_dP[TAG::WELL] = 0.0; // because totalDens does not depend on pressure
      for( integer ic = 0; ic < NC; ++ic )
      {
        dMult_dC[TAG::RES][ic]  = wellElemTotalDens * dResTotalMob_dC[ic];
        dMult_dC[TAG::WELL][ic] = resTotalMob;
      }

      real64 dMult[2][CP_Deriv::nDer]{};
      dMult[TAG::WELL][CP_Deriv::dP] = 0.0;
      if constexpr ( IS_THERMAL )
      {
        dMult[TAG::WELL][CP_Deriv::dT] = 0.0;
      }
      for( integer ic = 0; ic < NC; ++ic )
      {
        dMult[TAG::WELL][CP_Deriv::dC+ic] = resTotalMob;
      }
      for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
      {
        dMult[TAG::RES][jc] = wellElemTotalDens * dMob[jc];
      }


      // compute the volumetric flux and derivatives using upstream cell mobility
      flux = mult * potDiff;
      dFlux_dP[TAG::RES]  = dMult_dP[TAG::RES] * potDiff + mult * dPotDiff_dP[TAG::RES];
      dFlux_dP[TAG::WELL] = dMult_dP[TAG::WELL] * potDiff + mult * dPotDiff_dP[TAG::WELL];

      for( integer ic = 0; ic < NC; ++ic )
      {
        dFlux_dC[TAG::RES][ic]  = dMult_dC[TAG::RES][ic] * potDiff + mult * dPotDiff_dC[TAG::RES][ic];
        dFlux_dC[TAG::WELL][ic] = dMult_dC[TAG::WELL][ic] * potDiff + mult * dPotDiff_dC[TAG::WELL][ic];
      }

      for( integer ic = 0; ic < CP_Deriv::nDer; ++ic )
      {
        dFlux[TAG::RES][ic]  = dMult[TAG::RES][ic] * potDiff + mult * dPotDiff[TAG::RES][ic];
        dFlux[TAG::WELL][ic] = dMult[TAG::WELL][ic] * potDiff + mult * dPotDiff[TAG::WELL][ic];
      }
      // compute component fluxes
      for( integer ic = 0; ic < NC; ++ic )
      {
        m_compPerfRate[iperf][ic] += m_wellElemCompFrac[iwelem][ic] * flux;
        m_dCompPerfRate_dPres[iperf][TAG::RES][ic]  = m_wellElemCompFrac[iwelem][ic] * dFlux_dP[TAG::RES];
        m_dCompPerfRate_dPres[iperf][TAG::WELL][ic] = m_wellElemCompFrac[iwelem][ic] * dFlux_dP[TAG::WELL];

        for( integer jc = 0; jc < NC; ++jc )
        {
          m_dCompPerfRate_dComp[iperf][TAG::RES][ic][jc]  += m_wellElemCompFrac[iwelem][ic] * dFlux_dC[TAG::RES][jc];
          m_dCompPerfRate_dComp[iperf][TAG::WELL][ic][jc] += m_wellElemCompFrac[iwelem][ic] * dFlux_dC[TAG::WELL][jc];
          m_dCompPerfRate_dComp[iperf][TAG::WELL][ic][jc] += m_dWellElemCompFrac_dCompDens[iwelem][ic][jc] * flux;
        }
      }
      for( integer ic = 0; ic < NC; ++ic )
      {
        // Note this needs to be include below when this code above is removed
        // m_compPerfRate[iperf][ic] += m_wellElemCompFrac[iwelem][ic] * flux;
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          m_dCompPerfRate[iperf][TAG::RES][ic][jc]  = m_wellElemCompFrac[iwelem][ic] * dFlux[TAG::RES][jc];
        }
      }
      for( integer ic = 0; ic < NC; ++ic )
      {
        m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dP] = m_wellElemCompFrac[iwelem][ic] * dFlux[TAG::WELL][CP_Deriv::dP];
        if constexpr ( IS_THERMAL )
        {
          m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dT] = m_wellElemCompFrac[iwelem][ic] * dFlux[TAG::WELL][CP_Deriv::dT];
        }
        for( integer jc = 0; jc < NC; ++jc )
        {
          m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dC+jc] += m_wellElemCompFrac[iwelem][ic] * dFlux[TAG::WELL][CP_Deriv::dC+jc];
          m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dC+jc] += m_dWellElemCompFrac_dCompDens[iwelem][ic][jc] * flux;
        }
      }
      if constexpr ( IS_THERMAL )
      {
        for( integer ip = 0; ip < NP; ++ip )
        {
          real64 const wellelem_enthalpy = stack.m_wellElemPhaseEnthalpy[iwelem][0][ip];
          stack.m_energyPerfFlux[iperf] += flux * wellelem_enthalpy;
          // energy equation derivatives WRT res P & T
          stack.m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dP] += dFlux[TAG::RES][CP_Deriv::dP] * wellelem_enthalpy ;
          stack.m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dT] += dFlux[TAG::RES][CP_Deriv::dT] * wellelem_enthalpy ;
          // energy equation derivatives WRT well P & T
          stack.m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dP] += dFlux[TAG::WELL][CP_Deriv::dP] * wellelem_enthalpy 
                          +  flux * stack.m_dWellElemPhaseEnthalpy[iwelem][0][ip][CP_Deriv::dP];
          stack.m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dT] += dFlux[TAG::WELL][CP_Deriv::dT] * wellelem_enthalpy 
                          +  flux * stack.m_dWellElemPhaseEnthalpy[iwelem][0][ip][CP_Deriv::dT];   

          // energy equation derivatives WRT reservoir dens
          real64 dProp_dC[numComp]{};
          applyChainRule( NC,
                          m_dWellElemCompFrac_dCompDens[iwelem],
                          &stack.m_dWellElemPhaseEnthalpy[iwelem][0][ip][CP_Deriv::dC],
                          dProp_dC,
                          Deriv::dC );   

          for( integer jc = 0; jc < NC; ++jc )
          {
            stack.m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dC+jc] += flux * dProp_dC[jc];
          }  
        }
      }
      // tjb- remove when safe
      for( integer ic = 0; ic < NC; ic++ )
      {
        if( fabs(  m_dCompPerfRate[iperf][TAG::RES][ic][CP_Deriv::dP] -m_dCompPerfRate_dPres[iperf][TAG::RES][ic] ) > FLT_EPSILON )
        {
          std::cout << ic << " " <<  m_dCompPerfRate[iperf][TAG::RES][ic][CP_Deriv::dP] << " " << m_dCompPerfRate_dPres[iperf][TAG::RES][ic] << std::endl;
        }
        assert( fabs(  m_dCompPerfRate[iperf][TAG::RES][ic][CP_Deriv::dP] -m_dCompPerfRate_dPres[iperf][TAG::RES][ic] ) < FLT_EPSILON );
        assert( fabs(  m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dP] -m_dCompPerfRate_dPres[iperf][TAG::WELL][ic] ) < FLT_EPSILON );
        for( integer jc = 0; jc < NC; ++jc )
        {
          assert( fabs(  m_dCompPerfRate[iperf][TAG::RES][ic][CP_Deriv::dC+jc]  -m_dCompPerfRate_dComp[iperf][TAG::RES][ic][jc] ) < FLT_EPSILON );
          assert( fabs(  m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dC+jc]  -m_dCompPerfRate_dComp[iperf][TAG::WELL][ic][jc] ) < FLT_EPSILON );
        }
      }
      //compFluxKernelOp(er, esr, ei, potDiff, potGrad, flux, dFlux);
    } // end phase loop
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
    forAll< POLICY >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
    {

      kernelComponent.computeFlux(iperf, kernelComponent.m_stackVariables);
     
    } );
  }


StackVariables m_stackVariables;

protected:
  ElementViewConst< arrayView1d< real64 const > > const  m_resPres;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const  m_resPhaseVolFrac;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const  m_dResPhaseVolFrac;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const  m_dResCompFrac_dCompDens;
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const  m_resPhaseDens;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const  m_dResPhaseDens;
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const  m_resPhaseVisc;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const  m_dResPhaseVisc;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const  m_resPhaseCompFrac;
  ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const  m_dResPhaseCompFrac;
  ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const  m_resPhaseRelPerm;
  ElementViewConst< arrayView4d< real64 const, relperm::USD_RELPERM_DS > > const  m_dResPhaseRelPerm_dPhaseVolFrac;
  arrayView1d< real64 const > const  m_wellElemGravCoef;
  arrayView1d< real64 const > const  m_wellElemPres;
  arrayView2d< real64 const, compflow::USD_COMP > const  m_wellElemCompDens;
  arrayView1d< real64 const > const  m_wellElemTotalMassDens;
  arrayView2d< real64 const, compflow::USD_FLUID_DC > const  m_dWellElemTotalMassDens;
  arrayView1d< real64 const > const  m_dWellElemTotalMassDens_dPres;
  arrayView2d< real64 const, compflow::USD_FLUID_DC > const  m_dWellElemTotalMassDens_dCompDens;
  arrayView2d< real64 const, compflow::USD_COMP > const  m_wellElemCompFrac;
  arrayView3d< real64 const, compflow::USD_COMP_DC > const  m_dWellElemCompFrac_dCompDens;
  arrayView1d< real64 const > const  m_perfGravCoef;
  arrayView1d< localIndex const > const  m_perfWellElemIndex;
  arrayView1d< real64 const > const  m_perfTrans;
  arrayView1d< localIndex const > const  m_resElementRegion;
  arrayView1d< localIndex const > const  m_resElementSubRegion;
  arrayView1d< localIndex const > const  m_resElementIndex;
  arrayView2d< real64 > const  m_compPerfRate;
  arrayView4d< real64 > const  m_dCompPerfRate;
  arrayView3d< real64 > const  m_dCompPerfRate_dPres;
  arrayView4d< real64 > const  m_dCompPerfRate_dComp;

  bool const m_disableReservoirToWellFlow;


};

/**
 * @class PerforationKernelFactory
 */
class PerforationFluxKernelFactory
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
  createAndLaunch( integer const numComp,
                   integer const numPhases,
                   string const flowSolverName,
                   PerforationData * const perforationData,
                   ElementSubRegionBase const & subRegion,
                   ElementRegionManager & elemManager,
                   integer const disableReservoirToWellFlow )
  {
    geos::internal::kernelLaunchSelectorCompPhaseSwitch( numComp, numPhases , [&]( auto NC, auto NP )
        {
          integer constexpr NUM_COMP = NC();
          integer constexpr NUM_PHASE = NP();
          integer constexpr IS_THERMAL = 0;

          using kernelType = PerforationFluxKernel< NUM_COMP, NUM_PHASE, IS_THERMAL >;
          typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, flowSolverName );
          typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, flowSolverName );
          //typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, flowSolverName );
          typename kernelType::RelPermAccessors relPermAccessors( elemManager ,flowSolverName );
          //ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
          // elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
          //dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

          kernelType kernel( perforationData, subRegion, compFlowAccessors, multiFluidAccessors, relPermAccessors, disableReservoirToWellFlow );
          kernelType::template launch< POLICY >( perforationData->size(), kernel );
        } );
  }
};

}; // end namespace isothermalPerforationFluxKernels

namespace thermalPerforationFluxKernels
{

using namespace constitutive;

/******************************** PerforationFluxKernel ********************************/

template< integer NC, integer NP, integer IS_THERMAL >
class PerforationFluxKernel : public isothermalPerforationFluxKernels::PerforationFluxKernel<NC,NP,IS_THERMAL>
{
public:

  using Base = isothermalPerforationFluxKernels::PerforationFluxKernel<NC,NP,IS_THERMAL>;
  //using AbstractBase::m_dPhaseVolFrac;
  using Base::m_resPhaseCompFrac;
  using Base::m_dResCompFrac_dCompDens;
  //using AbstractBase::m_dPhaseCompFrac;
  //using AbstractBase::m_dCompFrac_dCompDens;
  /// Compile time value for the number of components
  static constexpr integer numComp = NC;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NP;

  /// Compile time value for thermal option
  static constexpr integer isThermal = IS_THERMAL;
 
  using TAG =  typename Base::TAG;
  using CompFlowAccessors = typename Base::CompFlowAccessors ;
  using MultiFluidAccessors = typename Base::MultiFluidAccessors;
  using RelPermAccessors = typename Base::RelPermAccessors;


  using ThermalCompFlowAccessors =
    StencilAccessors< fields::flow::temperature >;

  using ThermalMultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseEnthalpy,
                              fields::multifluid::dPhaseEnthalpy >;

  //using ThermalConductivityAccessors =
  //  StencilMaterialAccessors< MultiPhaseThermalConductivityBase,
  //                            fields::thermalconductivity::effectiveConductivity >;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  PerforationFluxKernel ( PerforationData * const perforationData,
                          ElementSubRegionBase const & subRegion,
                          MultiFluidBase const & fluid,
                          CompFlowAccessors const & compFlowAccessors,
                          MultiFluidAccessors const & multiFluidAccessors,
                          RelPermAccessors const & relPermAccessors,
                          bool const disableReservoirToWellFlow,
                          ThermalCompFlowAccessors const & thermalCompFlowAccessors,
                          ThermalMultiFluidAccessors const & thermalMultiFluidAccessors)
                          : Base(perforationData,
                          subRegion,
                          compFlowAccessors,
                          multiFluidAccessors,
                          relPermAccessors,
                          disableReservoirToWellFlow ),
                          m_stack(perforationData,
                          fluid,
                          thermalCompFlowAccessors,
                          thermalMultiFluidAccessors) 
                          
    //m_wellElemPhaseEnthalpy( fluid.phaseEnthalpy()),
    //m_dWellElemPhaseEnthalpy( fluid.dPhaseEnthalpy()),
   // m_energyPerfFlux( perforationData->getField< fields::well::energyPerforationFlux >()) ,
    //m_dEnergyPerfFlux( perforationData->getField< fields::well::dEnergyPerforationFlux >()) ,
    //m_temp( thermalCompFlowAccessors.get( fields::flow::temperature {} ) ),
    //m_resPhaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::phaseEnthalpy {} ) ),
    //m_dResPhaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::dPhaseEnthalpy {} ) ) 

  {} 


  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( PerforationData * const perforationData,
                    MultiFluidBase const & fluid,
                    ThermalCompFlowAccessors const & thermalCompFlowAccessors,
                    ThermalMultiFluidAccessors const & thermalMultiFluidAccessors )
      : Base::StackVariables(),
      m_wellElemPhaseEnthalpy( fluid.phaseEnthalpy()),
      m_dWellElemPhaseEnthalpy( fluid.dPhaseEnthalpy()),
      m_energyPerfFlux( perforationData->getField< fields::well::energyPerforationFlux >()) ,
      m_dEnergyPerfFlux( perforationData->getField< fields::well::dEnergyPerforationFlux >()) ,
      m_temp( thermalCompFlowAccessors.get( fields::flow::temperature {} ) ),
      m_resPhaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::phaseEnthalpy {} ) ),
      m_dResPhaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::dPhaseEnthalpy {} ) ) 
    {}

      /// Views on phase enthalpy
  arrayView3d< real64 const, multifluid::USD_PHASE >  m_wellElemPhaseEnthalpy;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC >  m_dWellElemPhaseEnthalpy;
 
  /// Views on energy flux
  arrayView1d< real64 > const & m_energyPerfFlux;
  arrayView3d< real64 > const & m_dEnergyPerfFlux;

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const  m_temp;

  /// Views on phase enthalpies
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const  m_resPhaseEnthalpy;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const  m_dResPhaseEnthalpy;
  };

  template< typename FUNC = NoOpFunc >  
  GEOS_HOST_DEVICE
  inline
  void
  computeFlux( localIndex const iperf ) const
  {
    Base::computeFlux(iperf, m_stack );
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
    forAll< POLICY >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
    {
      //kernelComponent.setupStack(iperf);
      kernelComponent.computeFlux(iperf);
     
    } );
  }

protected:

  /// Stack continaining thermal variables
  StackVariables m_stack;
 

};

/**
 * @class PerforationKernelFactory
 */
class PerforationFluxKernelFactory
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
  createAndLaunch( integer const numComp,
                   integer const numPhases,
                   string const flowSolverName,
                   PerforationData * const perforationData,
                   ElementSubRegionBase const & subRegion,
                   MultiFluidBase const & fluid,
                   ElementRegionManager & elemManager,
                   integer const disableReservoirToWellFlow )
  {
    geos::internal::kernelLaunchSelectorCompPhaseSwitch( numComp, numPhases , [&]( auto NC, auto NP )
        {
          integer constexpr NUM_COMP = NC();
          integer constexpr NUM_PHASE = NP();
          integer constexpr IS_THERMAL = 1;

          using kernelType = PerforationFluxKernel< NUM_COMP, NUM_PHASE, IS_THERMAL >;
          typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, flowSolverName );
          typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, flowSolverName );
          //typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, flowSolverName );
          typename kernelType::RelPermAccessors relPermAccessors( elemManager ,flowSolverName );
          typename kernelType::ThermalCompFlowAccessors  thermalCompFlowAccessors( elemManager, flowSolverName );
          typename kernelType::ThermalMultiFluidAccessors  thermalMultiFluidAccessors( elemManager, flowSolverName );
          //typename kernelType::ThermalConductivityAccessors  thermalConductivityAccessors( elemManager, flowSolverName );
          //ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
          // elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
          //dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

          kernelType kernel( perforationData, subRegion, fluid, compFlowAccessors, multiFluidAccessors, 
                            relPermAccessors, disableReservoirToWellFlow, 
                            thermalCompFlowAccessors, 
                            thermalMultiFluidAccessors );
          kernelType::template launch< POLICY >( perforationData->size(), kernel );
        } );
  }
};

} // end namespace thermalPerforationFluxKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_PERFORATIONFLUXLKERNELS_HPP
