/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
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
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "physicsSolvers/KernelLaunchSelectors.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
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
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseViscosity,
                              fields::multifluid::dPhaseViscosity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

  using RelPermAccessors =
    StencilMaterialAccessors< constitutive::RelativePermeabilityBase,
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
                          bool const isInjector,
                          bool const isCrossflowEnabled ):
    m_resPres( compFlowAccessors.get( fields::flow::pressure {} )),
    m_resPhaseVolFrac( compFlowAccessors.get( fields::flow::phaseVolumeFraction {} )),
    m_dResPhaseVolFrac( compFlowAccessors.get( fields::flow::dPhaseVolumeFraction {} )),
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
    m_perfStatus( perforationData->getField< fields::perforation::perforationStatus >()),
    m_isInjector( isInjector ),
    m_isCrossflowEnabled( isCrossflowEnabled )
  {}

  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void
  computeFlux( localIndex const iperf, FUNC && fluxKernelOp= NoOpFunc {} ) const
  {
    using CP_Deriv = constitutive::multifluid::DerivativeOffsetC< NC, IS_THERMAL >;

    if( !m_perfStatus[iperf] )
    {
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
      return;
    }
    // get the index of the reservoir elem
    localIndex const er  = m_resElementRegion[iperf];
    localIndex const esr = m_resElementSubRegion[iperf];
    localIndex const ei  = m_resElementIndex[iperf];

    // get the index of the well elem
    localIndex const iwelem = m_perfWellElemIndex[iperf];

    using Deriv = constitutive::multifluid::DerivativeOffset;


    // local working variables and arrays
    real64 pres[2]{};
    real64 multiplier[2]{};

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

    // Step 2: copy the variables from the reservoir and well element

    // a) get reservoir variables

    pres[TAG::RES] = m_resPres[er][esr][ei];
    dPres[TAG::RES][CP_Deriv::dP] = 1.0;
    multiplier[TAG::RES] = 1.0;

    // Here in the absence of a buoyancy term we assume that the reservoir cell is perforated at its center
    // TODO: add a buoyancy term for the reservoir side here


    // b) get well variables

    pres[TAG::WELL] = m_wellElemPres[iwelem];
    dPres[TAG::WELL][CP_Deriv::dP] = 1.0;
    multiplier[TAG::WELL] = -1.0;

    real64 const gravD = ( m_perfGravCoef[iperf] - m_wellElemGravCoef[iwelem] );

    pres[TAG::WELL] +=  m_wellElemTotalMassDens[iwelem] * gravD;
    // Note LHS uses CP_Deriv while RHS uses Deriv !!!
    dPres[TAG::WELL][CP_Deriv::dP] +=  m_dWellElemTotalMassDens[iwelem][Deriv::dP] * gravD;
    if constexpr ( IS_THERMAL )
    {
      dPres[TAG::WELL][CP_Deriv::dT] += m_dWellElemTotalMassDens[iwelem][Deriv::dT] * gravD;
    }
    for( integer ic = 0; ic < NC; ++ic )
    {
      dPres[TAG::WELL][CP_Deriv::dC+ic] += m_dWellElemTotalMassDens[iwelem][Deriv::dC+ic] * gravD;
    }

    // Step 3: compute potential difference

    real64 potDiff = 0.0;
    for( integer i = 0; i < 2; ++i )
    {
      potDiff += multiplier[i] * m_perfTrans[iperf] * pres[i];
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
        if( !phaseExists || (m_isInjector && !m_isCrossflowEnabled) )
        {
          continue;
        }

        // here, we have to recompute the reservoir phase mobility (not including density)

        // density
        real64 const resDens = m_resPhaseDens[er][esr][ei][0][ip];
        real64 dDens[CP_Deriv::nDer]{};

        dDens[CP_Deriv::dP]  = m_dResPhaseDens[er][esr][ei][0][ip][Deriv::dP];
        if constexpr ( IS_THERMAL )
        {
          dDens[CP_Deriv::dT]  = m_dResPhaseDens[er][esr][ei][0][ip][Deriv::dT];
        }
        applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                        m_dResPhaseDens[er][esr][ei][0][ip],
                        &dDens[CP_Deriv::dC],
                        Deriv::dC );
        // viscosity
        real64 const resVisc = m_resPhaseVisc[er][esr][ei][0][ip];
        real64 dVisc[CP_Deriv::nDer]{};
        dVisc[CP_Deriv::dP]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dP];
        if constexpr ( IS_THERMAL )
        {
          dVisc[CP_Deriv::dT]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dT];
        }

        applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                        m_dResPhaseVisc[er][esr][ei][0][ip],
                        &dVisc[CP_Deriv::dC],
                        Deriv::dC );

        // relative permeability
        real64 const resRelPerm = m_resPhaseRelPerm[er][esr][ei][0][ip];
        real64 dRelPerm[CP_Deriv::nDer]{};
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          dRelPerm[jc]=0;
        }
        for( integer jp = 0; jp < NP; ++jp )
        {
          real64 const dResRelPerm_dS = m_dResPhaseRelPerm_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dRelPerm[CP_Deriv::dP] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dP];
          if constexpr ( IS_THERMAL )
          {
            dRelPerm[CP_Deriv::dT] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dT];
          }
          for( integer jc = 0; jc < NC; ++jc )
          {
            dRelPerm[CP_Deriv::dC+jc] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
          }
        }

        // compute the reservoir phase mobility, including phase density
        real64 const resPhaseMob = resDens * resRelPerm / resVisc;

        // Handles all dependencies
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          dMob[jc] = dRelPerm[jc] * resDens / resVisc
                     + resPhaseMob * (dDens[jc] / resDens - dVisc[jc] / resVisc);
        }
        // compute the phase flux and derivatives using upstream cell mobility
        flux = resPhaseMob * potDiff;
        std::cout << "tjb resPhaseMob: " << resPhaseMob << " potDiff " << potDiff << " resDens " << resDens << " resRelPerm " << resRelPerm << " resVisc " << resVisc << std::endl;
        std::cout << "tjb flux: " << flux << " potdiff " << potDiff << std::endl;
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

      }  // end resevoir is upstream phase loop

    }
    else // ** well is upstream **
    {

      real64 resTotalMob     = 0.0;

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
        // or if crossflow is disabled for non-injectors aka producers
        bool const phaseExists = (m_resPhaseVolFrac[er][esr][ei][ip] > 0);
        if( !phaseExists || (!m_isInjector && !m_isCrossflowEnabled) )
        {
          continue;
        }

        // viscosity
        real64 const resVisc = m_resPhaseVisc[er][esr][ei][0][ip];
        real64 dVisc[CP_Deriv::nDer]{};
        dVisc[CP_Deriv::dP]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dP];
        if constexpr ( IS_THERMAL )
        {
          dVisc[CP_Deriv::dT]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dT];
        }

        applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                        m_dResPhaseVisc[er][esr][ei][0][ip],
                        &dVisc[CP_Deriv::dC],
                        Deriv::dC );


        // relative permeability
        real64 const resRelPerm = m_resPhaseRelPerm[er][esr][ei][0][ip];
        real64 dRelPerm[CP_Deriv::nDer]{};
        for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
        {
          dRelPerm[jc]=0;
        }
        for( integer jp = 0; jp < NP; ++jp )
        {
          real64 const dResRelPerm_dS = m_dResPhaseRelPerm_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dRelPerm[CP_Deriv::dP] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dP];
          if constexpr ( IS_THERMAL )
          {
            dRelPerm[CP_Deriv::dT] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dT];
          }
          for( integer jc = 0; jc < NC; ++jc )
          {
            dRelPerm[CP_Deriv::dC+jc] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
          }
        }
        // increment total mobility
        resTotalMob     += resRelPerm / resVisc;
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

      for( integer ic = 0; ic < CP_Deriv::nDer; ++ic )
      {
        dFlux[TAG::RES][ic]  = dMult[TAG::RES][ic] * potDiff + mult * dPotDiff[TAG::RES][ic];
        std::cout << "tjb dFlux RES: " << ic << " " << dFlux[TAG::RES][ic] << " dMult " << dMult[TAG::RES][ic] << " potDiff " << potDiff << " mult " << mult << " dPotDiff " <<
            dPotDiff[TAG::RES][ic] << std::endl;
        dFlux[TAG::WELL][ic] = dMult[TAG::WELL][ic] * potDiff + mult * dPotDiff[TAG::WELL][ic];
      }
      // compute component fluxes
      for( integer ic = 0; ic < NC; ++ic )
      {
        m_compPerfRate[iperf][ic] += m_wellElemCompFrac[iwelem][ic] * flux;
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
          std::cout << "tjb dFlux WELL dT: " << m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dT] << "  m_wellElemCompFrac[iwelem][ic] " <<  m_wellElemCompFrac[iwelem][ic] <<
              " dFlux[TAG::WELL][CP_Deriv::dT] " << dFlux[TAG::WELL][CP_Deriv::dT] << std::endl;
        }
        for( integer jc = 0; jc < NC; ++jc )
        {
          m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dC+jc] += m_wellElemCompFrac[iwelem][ic] * dFlux[TAG::WELL][CP_Deriv::dC+jc];
          m_dCompPerfRate[iperf][TAG::WELL][ic][CP_Deriv::dC+jc] += m_dWellElemCompFrac_dCompDens[iwelem][ic][jc] * flux;
        }
      }

    } // end upstream
    if constexpr ( IS_THERMAL )
    {
      fluxKernelOp( iwelem, er, esr, ei, potDiff, dPotDiff, flux, dFlux );
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
    forAll< POLICY >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
    {
      kernelComponent.computeFlux( iperf );
    } );
  }


protected:
  ElementViewConst< arrayView1d< real64 const > > const m_resPres;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_resPhaseVolFrac;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dResPhaseVolFrac;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const m_dResCompFrac_dCompDens;
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_resPhaseDens;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dResPhaseDens;
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_resPhaseVisc;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dResPhaseVisc;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const m_resPhaseCompFrac;
  ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const m_dResPhaseCompFrac;
  ElementViewConst< arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > > const m_resPhaseRelPerm;
  ElementViewConst< arrayView4d< real64 const, constitutive::relperm::USD_RELPERM_DS > > const m_dResPhaseRelPerm_dPhaseVolFrac;
  arrayView1d< real64 const > const m_wellElemGravCoef;
  arrayView1d< real64 const > const m_wellElemPres;
  arrayView2d< real64 const, compflow::USD_COMP > const m_wellElemCompDens;
  arrayView1d< real64 const > const m_wellElemTotalMassDens;
  arrayView2d< real64 const, compflow::USD_FLUID_DC > const m_dWellElemTotalMassDens;
  arrayView2d< real64 const, compflow::USD_COMP > const m_wellElemCompFrac;
  arrayView3d< real64 const, compflow::USD_COMP_DC > const m_dWellElemCompFrac_dCompDens;
  arrayView1d< real64 const > const m_perfGravCoef;
  arrayView1d< localIndex const > const m_perfWellElemIndex;
  arrayView1d< real64 const > const m_perfTrans;
  arrayView1d< localIndex const > const m_resElementRegion;
  arrayView1d< localIndex const > const m_resElementSubRegion;
  arrayView1d< localIndex const > const m_resElementIndex;
  arrayView2d< real64 > const m_compPerfRate;
  arrayView4d< real64 > const m_dCompPerfRate;
  arrayView3d< real64 > const m_dCompPerfRate_dPres;
  arrayView4d< real64 > const m_dCompPerfRate_dComp;
  arrayView1d< integer > const m_perfStatus;

  bool const m_isInjector;
  bool const m_isCrossflowEnabled;

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
   * @param[in] numPhases the number of phases
   * @param[in] flowSolverName the name of the flow solver
   * @param[in] perforationData the perforation data
   * @param[in] subregion well subregion
   * @param[in] elemManager element region manager
   * @param[in] isInjector flag indicating whether the well is an injector
   * @param[in] isCrossflowEnabled flag indicating whether crossflow is enabled
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhases,
                   string const flowSolverName,
                   PerforationData * const perforationData,
                   ElementSubRegionBase const & subRegion,
                   ElementRegionManager const & elemManager,
                   bool const isInjector,
                   bool const isCrossflowEnabled )
  {
    geos::internal::kernelLaunchSelectorCompPhaseSwitch( numComp, numPhases, [&]( auto NC, auto NP )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_PHASE = NP();
      integer constexpr IS_THERMAL = 0;

      using kernelType = PerforationFluxKernel< NUM_COMP, NUM_PHASE, IS_THERMAL >;
      typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, flowSolverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, flowSolverName );
      typename kernelType::RelPermAccessors relPermAccessors( elemManager, flowSolverName );

      kernelType kernel( perforationData, subRegion,
                         compFlowAccessors, multiFluidAccessors, relPermAccessors,
                         isInjector, isCrossflowEnabled );
      kernelType::template launch< POLICY >( perforationData->size(), kernel );
    } );
  }
};

} // end namespace isothermalPerforationFluxKernels

namespace thermalPerforationFluxKernels
{

using namespace constitutive;

/******************************** PerforationFluxKernel ********************************/

template< integer NC, integer NP, integer IS_THERMAL >
class PerforationFluxKernel : public isothermalPerforationFluxKernels::PerforationFluxKernel< NC, NP, IS_THERMAL >
{
public:

  using Base = isothermalPerforationFluxKernels::PerforationFluxKernel< NC, NP, IS_THERMAL >;
  using Base::m_dWellElemCompFrac_dCompDens;
  using Base::m_resPres;
  using Base::m_resPhaseVolFrac;
  using Base::m_dResPhaseVolFrac;
  using Base::m_dResCompFrac_dCompDens;
  using Base::m_resPhaseDens;
  using Base::m_dResPhaseDens;
  using Base::m_resPhaseVisc;
  using Base::m_dResPhaseVisc;
  using Base::m_resPhaseCompFrac;
  using Base::m_dResPhaseCompFrac;
  using Base::m_resPhaseRelPerm;
  using Base::m_dResPhaseRelPerm_dPhaseVolFrac;
  using Base::m_isInjector;
  using Base::m_isCrossflowEnabled;

  /// Compile time value for the number of components
  static constexpr integer numComp = NC;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NP;

  /// Compile time value for thermal option
  static constexpr integer isThermal = IS_THERMAL;

  using TAG =  typename Base::TAG;
  using CompFlowAccessors = typename Base::CompFlowAccessors;
  using MultiFluidAccessors = typename Base::MultiFluidAccessors;
  using RelPermAccessors = typename Base::RelPermAccessors;


  using ThermalCompFlowAccessors =
    StencilAccessors< fields::flow::temperature >;

  using ThermalMultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseFraction,
                              fields::multifluid::dPhaseFraction,
                              fields::multifluid::phaseEnthalpy,
                              fields::multifluid::dPhaseEnthalpy >;

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
                          ElementSubRegionBase & subRegion,
                          MultiFluidBase const & wellFluid,
                          CompFlowAccessors const & compFlowAccessors,
                          MultiFluidAccessors const & multiFluidAccessors,
                          RelPermAccessors const & relPermAccessors,
                          ThermalCompFlowAccessors const & thermalCompFlowAccessors,
                          ThermalMultiFluidAccessors const & thermalMultiFluidAccessors,
                          bool const isInjector,
                          bool const isCrossflowEnabled )
    : Base( perforationData,
            subRegion,
            compFlowAccessors,
            multiFluidAccessors,
            relPermAccessors,
            isInjector,
            isCrossflowEnabled ),
    m_wellElemPhaseFrac( wellFluid.phaseFraction() ),
    m_dWellElemPhaseFrac( wellFluid.dPhaseFraction() ),
    m_wellElemPhaseEnthalpy( wellFluid.phaseEnthalpy()),
    m_dWellElemPhaseEnthalpy( wellFluid.dPhaseEnthalpy()),
    m_energyPerfFlux( perforationData->getField< fields::well::energyPerforationFlux >() ),
    m_dEnergyPerfFlux( perforationData->getField< fields::well::dEnergyPerforationFlux >() ),
    m_temp( thermalCompFlowAccessors.get( fields::flow::temperature {} ) ),
    m_resPhaseFraction( thermalMultiFluidAccessors.get( fields::multifluid::phaseFraction {} ) ),
    m_dResPhaseFraction( thermalMultiFluidAccessors.get( fields::multifluid::dPhaseFraction {} )),
    m_resPhaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::phaseEnthalpy {} ) ),
    m_dResPhaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::dPhaseEnthalpy {} ) )
  {}

  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void
  computeFlux( localIndex const iperf ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;
    using CP_Deriv =constitutive::multifluid::DerivativeOffsetC< NC, IS_THERMAL >;
    // initialize outputs
    m_energyPerfFlux[iperf]=0;
    for( integer ke = 0; ke < 2; ++ke )
    {
      for( integer i = 0; i < CP_Deriv::nDer; ++i )
      {
        m_dEnergyPerfFlux[iperf][ke][i]=0;
      }
    }

    Base::computeFlux ( iperf, [&]( localIndex const iwelem, localIndex const er, localIndex const esr, localIndex const ei
                                    , real64 const potDiff, real64 const dPotDiff[2][CP_Deriv::nDer],
                                    real64 const flux, real64 const (&dFlux)[2][CP_Deriv::nDer] )
    {
      real64 dMob[CP_Deriv::nDer]{};
      GEOS_UNUSED_VAR( iwelem );
      if( potDiff >= 0 )    // ** reservoir cell is upstream **
      {
        real64 eflux = 0.0;
        // loop over phases, compute and upwind phase flux
        // and sum contributions to each component's perforation rate
        for( integer ip = 0; ip < NP; ++ip )
        {
          // skip the rest of the calculation if the phase is absent
          // or if crossflow is disabled for injectors
          real64 const resPhaseVolFrac = m_resPhaseVolFrac[er][esr][ei][ip];
          bool const phaseExists = (resPhaseVolFrac > 0 );
          if( !phaseExists || (m_isInjector && !m_isCrossflowEnabled) )
          {
            continue;
          }

          // here, we have to recompute the reservoir phase mobility (not including density)

          // density
          real64 const resDens = m_resPhaseDens[er][esr][ei][0][ip];
          real64 dDens[CP_Deriv::nDer]{};
          dDens[CP_Deriv::dP]  = m_dResPhaseDens[er][esr][ei][0][ip][Deriv::dP];
          dDens[CP_Deriv::dT]  = m_dResPhaseDens[er][esr][ei][0][ip][Deriv::dT];
          applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                          m_dResPhaseDens[er][esr][ei][0][ip],
                          &dDens[CP_Deriv::dC],
                          Deriv::dC );
          // viscosity
          real64 const resVisc = m_resPhaseVisc[er][esr][ei][0][ip];
          real64 dVisc[CP_Deriv::nDer]{};
          dVisc[CP_Deriv::dP]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dP];
          dVisc[CP_Deriv::dT]  = m_dResPhaseVisc[er][esr][ei][0][ip][Deriv::dT];
          applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                          m_dResPhaseVisc[er][esr][ei][0][ip],
                          &dVisc[CP_Deriv::dC],
                          Deriv::dC );
          // enthalpy
          real64 const resEnthalpy =  m_resPhaseEnthalpy[er][esr][ei][0][ip];
          real64 dResEnthalpy[CP_Deriv::nDer]{};
          dResEnthalpy[CP_Deriv::dP]  = m_dResPhaseEnthalpy[er][esr][ei][0][ip][Deriv::dP];
          dResEnthalpy[CP_Deriv::dT]  = m_dResPhaseEnthalpy[er][esr][ei][0][ip][Deriv::dT];

          applyChainRule( NC, m_dResCompFrac_dCompDens[er][esr][ei],
                          m_dResPhaseEnthalpy[er][esr][ei][0][ip],
                          &dResEnthalpy[CP_Deriv::dC],
                          Deriv::dC );

          // relative permeability
          real64 const resRelPerm = m_resPhaseRelPerm[er][esr][ei][0][ip];
          real64 dRelPerm[CP_Deriv::nDer]{};
          for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
          {
            dRelPerm[jc]=0;
          }
          for( integer jp = 0; jp < NP; ++jp )
          {
            real64 const dResRelPerm_dS = m_dResPhaseRelPerm_dPhaseVolFrac[er][esr][ei][0][ip][jp];
            dRelPerm[CP_Deriv::dP] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dP];
            dRelPerm[CP_Deriv::dT] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dT];
            for( integer jc = 0; jc < NC; ++jc )
            {
              dRelPerm[CP_Deriv::dC+jc] += dResRelPerm_dS * m_dResPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
            }
          }

          // compute the reservoir phase mobility, including phase density and enthalpy
          real64 const resPhaseMob = resDens * resRelPerm / resVisc;
          real64 const resPhaseMobE = resEnthalpy * resPhaseMob;
          std::cout << "tjb resPhaseMob: " << resPhaseMob << ", resPhaseMobE: " << resPhaseMobE << " resDens " << resDens << " resRelPerm " << resRelPerm << " resVisc " << resVisc << " resEnthalpy " << resEnthalpy <<  std::endl;
          // Handles all dependencies
          for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
          {
            dMob[jc] = resPhaseMob*dResEnthalpy[jc]  +
                       resPhaseMobE*(dRelPerm[jc]/(resRelPerm+0.000000001) + dDens[jc] / (resDens +0.000000001)  - dVisc[jc]/(resVisc+0.000000001));
          }
          // compute the phase flux and derivatives using upstream cell mobility
          eflux = resPhaseMobE * potDiff;

          real64 dEFlux[2][CP_Deriv::nDer]{};
          dEFlux[TAG::WELL][CP_Deriv::dP] = resPhaseMobE * dPotDiff[TAG::WELL][CP_Deriv::dP];
          dEFlux[TAG::WELL][CP_Deriv::dT] = resPhaseMobE * dPotDiff[TAG::WELL][CP_Deriv::dT];
          // Handles all dependencies
          for( integer jc = 0; jc < CP_Deriv::nDer; ++jc )
          {
            dEFlux[TAG::RES][jc]  = dMob[jc] * potDiff + resPhaseMobE * dPotDiff[TAG::RES][jc];
            m_dEnergyPerfFlux[iperf][TAG::WELL][jc] = resPhaseMobE * dPotDiff[TAG::WELL][jc];
            std::cout << "tjb TAG::WELL " << jc <<  " m_dEnergyPerfFlux[iperf][TAG::WELL][jc] " << m_dEnergyPerfFlux[iperf][TAG::WELL][jc] <<   std::endl;
          }
          // m_energyPerfFlux[iperf] += resPhaseVolFrac * eflux;
          m_energyPerfFlux[iperf] +=   eflux;
          // energy equation derivatives WRT res P & T
          std::cout << "tjb m_dEnergyPerfFlux phase : " << ip << m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dP] << std::endl;
          m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dP] += dEFlux[TAG::RES][CP_Deriv::dP] * resPhaseVolFrac +
                                                              eflux *   m_dResPhaseFraction[er][esr][ei][0][ip][Deriv::dP];
          std::cout << "tjb m_dEnergyPerfFlux 1 phase : " << ip << " resPhaseVolFrac " << resPhaseVolFrac << " dEFlux[TAG::RES][CP_Deriv::dP] " << dEFlux[TAG::RES][CP_Deriv::dP] << std::endl;
          std::cout << "tjb m_dEnergyPerfFlux 2 phase : " << ip << " eflux " << eflux << "  m_dResPhaseFraction[er][esr][ei][0][ip][Deriv::dP] " << m_dResPhaseFraction[er][esr][ei][0][ip][Deriv::dP] << std::endl;
          m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dT] += dEFlux[TAG::RES][CP_Deriv::dT] * resPhaseVolFrac +
                                                              eflux *m_dResPhaseFraction[er][esr][ei][0][ip][Deriv::dT];
          // energy equation derivatives WRT well P - //tjb
          //m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dP] += dEFlux[TAG::WELL][CP_Deriv::dP];
          // m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dT] += dEFlux[TAG::WELL][CP_Deriv::dT];
          std::cout << "tjb TAG::WELL m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dP]  " << m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dP] <<  " m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dT]  " <<
              m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dT] << std::endl;
          real64 dProp_dC[numComp]{};
          applyChainRule( NC,
                          m_dResCompFrac_dCompDens[er][esr][ei],
                          m_dResPhaseFraction[er][esr][ei][0][ip],
                          dProp_dC,
                          Deriv::dC );
          for( integer jc = 0; jc < NC; ++jc )
          {
            //real64 const resPhaseCompFrac = m_resPhaseCompFrac[er][esr][ei][0][ip][jc];
            // m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dC+jc] +=  eflux * dProp_dC[jc] +
            //                                                       resPhaseCompFrac * dEFlux[TAG::RES][CP_Deriv::dC+jc];
            m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dC+jc] +=  dEFlux[TAG::RES][CP_Deriv::dC+jc];
          }

        } // end resevoir is upstream phase loop

      }
      else   // ** reservoir cell is downstream
      {
        for( integer iphase = 0; iphase < NP; ++iphase )
        {
          bool const phaseExists = (m_wellElemPhaseFrac[iwelem][0][iphase] > 0);
          if( !phaseExists || (!m_isInjector && !m_isCrossflowEnabled) )
          {
            continue;
          }
          real64 const phaseFlux =  flux;
          real64 const wellelem_enthalpy = m_wellElemPhaseEnthalpy[iwelem][0][iphase];
          m_energyPerfFlux[iperf] += phaseFlux * wellelem_enthalpy;
          std::cout << "tjbi phaseFlux " << phaseFlux << " wellelem_enthalpy " << wellelem_enthalpy << " m_energyPerfFlux[iperf] " << m_energyPerfFlux[iperf] << " m_wellElemPhaseFrac[iwelem][0][iphase] " <<
              m_wellElemPhaseFrac[iwelem][0][iphase] << std::endl;
          // energy equation derivatives WRT res P & T
          std::cout << "tjbi m_dEnergyPerfFlux downstream phase : " << iphase << " " << m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dP] << " dFlux[TAG::RES][CP_Deriv::dP] " <<
              dFlux[TAG::RES][CP_Deriv::dP] << std::endl;
          m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dP] +=   dFlux[TAG::RES][CP_Deriv::dP] * wellelem_enthalpy;
          std::cout << "tjbi m_dEnergyPerfFlux downstream phase : " << iphase << " " << m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dP] << std::endl;
          m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dT] +=  dFlux[TAG::RES][CP_Deriv::dT] * wellelem_enthalpy;
          std::cout << "tjbi m_dEnergyPerfFlux dt downstream phase : " << iphase << " " << m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dT] << std::endl;
          m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dP] +=   dFlux[TAG::WELL][CP_Deriv::dP] * wellelem_enthalpy
                                                               + phaseFlux * m_dWellElemPhaseEnthalpy[iwelem][0][iphase][Deriv::dP];
          // + phaseFlux * wellelem_enthalpy * m_dWellElemPhaseFrac[iwelem][0][iphase][Deriv::dP];
          m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dT] +=   dFlux[TAG::WELL][CP_Deriv::dT] * wellelem_enthalpy
                                                               + phaseFlux * m_dWellElemPhaseEnthalpy[iwelem][0][iphase][Deriv::dT];
          //+ phaseFlux * wellelem_enthalpy * m_dWellElemPhaseFrac[iwelem][0][iphase][Deriv::dT];

          //energy e
          real64 dPVF_dC[numComp]{};
          applyChainRule( NC,
                          m_dWellElemCompFrac_dCompDens[iwelem],
                          m_dWellElemPhaseFrac[iwelem][0][iphase],
                          dPVF_dC,
                          Deriv::dC );
          for( integer ic=0; ic<NC; ic++ )
          {
            m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dC+ic]  += wellelem_enthalpy * dFlux[TAG::WELL][CP_Deriv::dC+ic];// * m_wellElemPhaseFrac[iwelem][0][iphase];
            //m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dC+ic]  += wellelem_enthalpy * phaseFlux * dPVF_dC[ic];
          }
          // energy equation enthalpy derivatives WRT well dens
          real64 dProp_dC[numComp]{};
          applyChainRule( NC,
                          m_dWellElemCompFrac_dCompDens[iwelem],
                          m_dWellElemPhaseEnthalpy[iwelem][0][iphase],
                          dProp_dC,
                          Deriv::dC );

          for( integer jc = 0; jc < NC; ++jc )
          {
            m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dC+jc] += phaseFlux * dProp_dC[jc];
          }
        }

      }
    } );
  }



  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElements the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack
   * variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElements,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
    {
      kernelComponent.computeFlux( iperf );
    } );
  }

protected:

  /// Views on well element properties
  /// Element phase fraction
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_wellElemPhaseFrac;
  arrayView2d< real64 const, compflow::USD_PHASE > const m_wellElemPhaseVolFrac;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const m_dWellElemPhaseFrac;
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_wellElemPhaseEnthalpy;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const m_dWellElemPhaseEnthalpy;

  /// Views on energy flux
  arrayView1d< real64 > const m_energyPerfFlux;
  arrayView3d< real64 > const m_dEnergyPerfFlux;

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;

  /// Views on phase enthalpies
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_resPhaseFraction;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const m_dResPhaseFraction;
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_resPhaseEnthalpy;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const m_dResPhaseEnthalpy;


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
   * @param[in] numPhases the number of phases
   * @param[in] flowSolverName the name of the flow solver
   * @param[in] perforationData the perforation data
   * @param[in] subregion well subregion
   * @param[in] fluid the multifluid object
   * @param[in] elemManager element region manager
   * @param[in] isInjector flag indicating whether the well is an injector
   * @param[in] isCrossflowEnabled flag indicating whether crossflow is enabled
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhases,
                   string const flowSolverName,
                   PerforationData * const perforationData,
                   ElementSubRegionBase & subRegion,
                   MultiFluidBase const & fluid,
                   ElementRegionManager const & elemManager,
                   bool const isInjector,
                   bool const isCrossflowEnabled )
  {
    geos::internal::kernelLaunchSelectorCompPhaseSwitch( numComp, numPhases, [&]( auto NC, auto NP )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_PHASE = NP();
      integer constexpr IS_THERMAL = 1;

      using kernelType = PerforationFluxKernel< NUM_COMP, NUM_PHASE, IS_THERMAL >;
      typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, flowSolverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, flowSolverName );
      typename kernelType::RelPermAccessors relPermAccessors( elemManager, flowSolverName );
      typename kernelType::ThermalCompFlowAccessors thermalCompFlowAccessors( elemManager, flowSolverName );
      typename kernelType::ThermalMultiFluidAccessors thermalMultiFluidAccessors( elemManager, flowSolverName );

      kernelType kernel( perforationData, subRegion, fluid, compFlowAccessors, multiFluidAccessors,
                         relPermAccessors, thermalCompFlowAccessors, thermalMultiFluidAccessors,
                         isInjector, isCrossflowEnabled );
      kernelType::template launch< POLICY >( perforationData->size(), kernel );
    } );
  }
};

}   // end namespace thermalPerforationFluxKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_PERFORATIONFLUXLKERNELS_HPP
