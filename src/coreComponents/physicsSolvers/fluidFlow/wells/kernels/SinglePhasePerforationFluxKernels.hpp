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

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEPERFORATIONFLUXLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEPERFORATIONFLUXLKERNELS_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"

#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"

#include "mesh/ElementRegionManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "physicsSolvers/KernelLaunchSelectors.hpp"

#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellTags.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"


namespace geos
{

struct NoOpStuct
{
  NoOpStuct(){}
};

namespace isothermalSinglePhasePerforationFluxKernels
{

/******************************** PerforationFluxKernel ********************************/

template< integer IS_THERMAL >
class PerforationFluxKernel
{
public:

  /// Compile time value for thermal option
  static constexpr integer isThermal = IS_THERMAL;

  using TAG = wellTags::SubRegionTag;
  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::pressure >;

  using SingleFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity,
                              fields::singlefluid::viscosity,
                              fields::singlefluid::dViscosity >;


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
                          constitutive::SingleFluidBase const & fluid,
                          SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                          SingleFluidAccessors const & singleFluidAccessors
                          ):
    m_resPres( singlePhaseFlowAccessors.get( fields::flow::pressure {} )),
    m_resDens( singleFluidAccessors.get( fields::singlefluid::density {} )),
    m_dResDens( singleFluidAccessors.get( fields::singlefluid::dDensity {} )),
    m_resVisc( singleFluidAccessors.get( fields::singlefluid::viscosity {} )),
    m_dResVisc( singleFluidAccessors.get( fields::singlefluid::dViscosity {} )),
    m_wellElemGravCoef( subRegion.getField< fields::well::gravityCoefficient >()),
    m_wellElemPres( subRegion.getField< fields::well::pressure >()),
    m_wellElemDens( fluid.density()),
    m_dWellElemDens( fluid.dDensity()),
    m_wellElemVisc( fluid.viscosity()),
    m_dWellElemVisc( fluid.dViscosity()),
    m_perfGravCoef( perforationData->getField< fields::well::gravityCoefficient >()),
    m_perfWellElemIndex( perforationData->getField< fields::perforation::wellElementIndex >()),
    m_perfTrans( perforationData->getField< fields::perforation::wellTransmissibility >()),
    m_resElementRegion( perforationData->getField< fields::perforation::reservoirElementRegion >()),
    m_resElementSubRegion( perforationData->getField< fields::perforation::reservoirElementSubRegion >()),
    m_resElementIndex( perforationData->getField< fields::perforation::reservoirElementIndex >()),
    m_perfRate( perforationData->getField< fields::well::perforationRate >()),
    m_dPerfRate( perforationData->getField< fields::well::dPerforationRate >())
  {}

  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void
  computeFlux( localIndex const iperf, FUNC && fluxKernelOp= NoOpFunc {} ) const
  {

    // get the reservoir (sub)region and element indices
    localIndex const er  = m_resElementRegion[iperf];
    localIndex const esr = m_resElementSubRegion[iperf];
    localIndex const ei  = m_resElementIndex[iperf];

    // get the local index of the well element
    localIndex const iwelem = m_perfWellElemIndex[iperf];

    using Deriv = constitutive::singlefluid::DerivativeOffset;
    using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

    // local working variables and arrays
    real64 pressure[2]{};
    real64 dPressure[2][DerivOffset::nDer]{};

    real64 multiplier[2]{};

    m_perfRate[iperf] = 0.0;
    for( integer i = 0; i < 2; ++i )
    {
      for( integer j=0; j<DerivOffset::nDer; j++ )
      {
        m_dPerfRate[iperf][i][j] = 0.0;
      }
    }

    // 1) Reservoir side
    // get reservoir variables
    pressure[TAG::RES] = m_resPres[er][esr][ei];
    dPressure[TAG::RES][DerivOffset::dP] = 1.0;

    // TODO: add a buoyancy term for the reservoir side here

    // multiplier for reservoir side in the flux
    multiplier[TAG::RES] = 1;

    // 2) Well side

    // get well variables
    pressure[TAG::WELL] = m_wellElemPres[iwelem];
    dPressure[TAG::WELL][DerivOffset::dP] = 1.0;

    real64 const gravD = ( m_perfGravCoef[iperf] - m_wellElemGravCoef[iwelem] );
    pressure[TAG::WELL]     += m_wellElemDens[iwelem][0] * gravD;
    dPressure[TAG::WELL][DerivOffset::dP] +=  m_dWellElemDens[iwelem][0][Deriv::dP] * gravD;
    if constexpr (IS_THERMAL)
    {
      dPressure[TAG::WELL][DerivOffset::dT] = m_dWellElemDens[iwelem][0][DerivOffset::dT] * gravD;
    }

    // multiplier for well side in the flux
    multiplier[TAG::WELL] = -1;

    // compute potential difference
    real64 potDif = 0.0;
    for( localIndex i = 0; i < 2; ++i )
    {
      potDif += multiplier[i] * m_perfTrans[iperf] * pressure[i];
      for( integer j=0; j<DerivOffset::nDer; j++ )
      {
        m_dPerfRate[iperf][i][j] = multiplier[i] * m_perfTrans[iperf] * dPressure[i][j];
      }
    }

    // choose upstream cell based on potential difference
    localIndex const k_up = (potDif >= 0) ? TAG::RES : TAG::WELL;

    // compute upstream density, viscosity, and mobility
    real64 densityUp       = 0.0;
    real64 viscosityUp     = 0.0;
    real64 dDensityUp[DerivOffset::nDer]{};
    real64 dViscosityUp[DerivOffset::nDer]{};
    for( integer j=0; j<DerivOffset::nDer; j++ )
    {
      dDensityUp[j]=0.0;
      dViscosityUp[j]=0.0;
    }

    // upwinding the variables
    if( k_up == TAG::RES ) // use reservoir vars
    {
      densityUp     = m_resDens[er][esr][ei][0];
      viscosityUp     = m_resVisc[er][esr][ei][0];
      for( integer j=0; j<DerivOffset::nDer; j++ )
      {
        dDensityUp[j] = m_dResDens[er][esr][ei][0][j];
        dViscosityUp[j] = m_dResVisc[er][esr][ei][0][j];
      }
    }
    else // use well vars
    {
      densityUp = m_wellElemDens[iwelem][0];
      viscosityUp = m_wellElemVisc[iwelem][0];
      for( integer j=0; j<DerivOffset::nDer; j++ )
      {
        dDensityUp[j] = m_dWellElemDens[iwelem][0][j];
        dViscosityUp[j] = m_dWellElemVisc[iwelem][0][j];
      }
    }

    // compute mobility
    real64 const mobilityUp     = densityUp / viscosityUp;
    real64 dMobilityUp[DerivOffset::nDer]{};
    for( integer j=0; j<DerivOffset::nDer; j++ )
    {
      dMobilityUp[j] = dDensityUp[j] / viscosityUp
                       - mobilityUp / viscosityUp * dViscosityUp[j];
    }

    m_perfRate[iperf] = mobilityUp * potDif;
    for( localIndex ke = 0; ke < 2; ++ke )
    {
      for( integer j=0; j<DerivOffset::nDer; j++ )
      {
        m_dPerfRate[iperf][ke][j]  *= mobilityUp;
      }
    }
    for( integer j=0; j<DerivOffset::nDer; j++ )
    {
      m_dPerfRate[iperf][k_up][j]  += dMobilityUp[j]* potDif;
    }
    if constexpr ( IS_THERMAL )
    {
      fluxKernelOp( iwelem, er, esr, ei, potDif, m_perfRate[iperf], m_dPerfRate[iperf] );
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
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_resDens;
  ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const m_dResDens;
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_resVisc;
  ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const m_dResVisc;
  arrayView1d< real64 const > const m_wellElemGravCoef;
  arrayView1d< real64 const > const m_wellElemPres;
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_wellElemDens;
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const m_dWellElemDens;
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_wellElemVisc;
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const m_dWellElemVisc;

  arrayView1d< real64 const > const m_perfGravCoef;
  arrayView1d< localIndex const > const m_perfWellElemIndex;
  arrayView1d< real64 const > const m_perfTrans;
  arrayView1d< localIndex const > const m_resElementRegion;
  arrayView1d< localIndex const > const m_resElementSubRegion;
  arrayView1d< localIndex const > const m_resElementIndex;
  arrayView1d< real64 > const m_perfRate;
  arrayView3d< real64 > const m_dPerfRate;

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
   * @param[in] flowSolverName the flow solver
   * @param[in] perforationData container with perforation variables
   * @param[in] subRegion well subregion
   * @param[in] fluid  subregion fluid model
   * @param[in] subregion well subregion
   * @param[in] elemManager  region element managers
   */
  template< typename POLICY >
  static void
  createAndLaunch( string const flowSolverName,
                   PerforationData * const perforationData,
                   ElementSubRegionBase const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   ElementRegionManager & elemManager )
  {
    integer constexpr IS_THERMAL = 0;
    using kernelType = PerforationFluxKernel< IS_THERMAL >;
    typename kernelType::SinglePhaseFlowAccessors singlePhaseFlowAccessors( elemManager, flowSolverName );
    typename kernelType::SingleFluidAccessors singleFluidAccessors( elemManager, flowSolverName );
    kernelType kernel( perforationData, subRegion, fluid, singlePhaseFlowAccessors, singleFluidAccessors );
    kernelType::template launch< POLICY >( perforationData->size(), kernel );
  }
};

} // end namespace isothermalPerforationFluxKernels

namespace thermalSinglePhasePerforationFluxKernels
{

using namespace constitutive;

/******************************** PerforationFluxKernel ********************************/

template< integer IS_THERMAL >
class PerforationFluxKernel : public geos::isothermalSinglePhasePerforationFluxKernels::PerforationFluxKernel< IS_THERMAL >
{
public:

  using Base = geos::isothermalSinglePhasePerforationFluxKernels::PerforationFluxKernel< IS_THERMAL >;
  using SinglePhaseFlowAccessors = typename Base::SinglePhaseFlowAccessors;
  using SingleFluidAccessors = typename Base::SingleFluidAccessors;

  /// Compile time value for thermal option
  static constexpr integer isThermal = IS_THERMAL;

  using TAG =  typename Base::TAG;

  using ThermalSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::temperature >;

  using ThermalSingleFluidAccessors =
    StencilMaterialAccessors< SingleFluidBase,
                              fields::singlefluid::enthalpy,
                              fields::singlefluid::dEnthalpy >;


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
                          SingleFluidBase const & fluid,
                          SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                          SingleFluidAccessors const & singleFluidAccessors,
                          ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                          ThermalSingleFluidAccessors const & thermalSingleFluidAccessors )
    : Base( perforationData,
            subRegion,
            fluid,
            singlePhaseFlowAccessors,
            singleFluidAccessors ),
    m_wellElemEnthalpy( fluid.enthalpy()),
    m_dWellElemEnthalpy( fluid.dEnthalpy()),
    m_energyPerfFlux( perforationData->getField< fields::well::energyPerforationFlux >()),
    m_dEnergyPerfFlux( perforationData->getField< fields::well::dEnergyPerforationFlux >()),
    m_temp( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) ),
    m_resEnthalpy( thermalSingleFluidAccessors.get( fields::singlefluid::enthalpy {} ) ),
    m_dResEnthalpy( thermalSingleFluidAccessors.get( fields::singlefluid::dEnthalpy {} ) )

  {}

  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void
  computeFlux( localIndex const iperf ) const
  {
    using Deriv = constitutive::singlefluid::DerivativeOffset;
    using CP_Deriv =constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;
    // initialize outputs
    m_energyPerfFlux[iperf]=0;
    for( integer ke = 0; ke < 2; ++ke )
    {
      for( integer i = 0; i < CP_Deriv::nDer; ++i )
      {
        m_dEnergyPerfFlux[iperf][ke][i]=0;
      }
    }

    Base::computeFlux ( iperf, [&]( localIndex const iwelem, localIndex const er, localIndex const esr, localIndex const ei,
                                    real64 const potDiff, real64 const flux, arraySlice2d< real64 > const dFlux )
    {
      if( potDiff >= 0 )      // ** reservoir cell is upstream **
      {

        real64 const res_enthalpy =  m_resEnthalpy[er][esr][ei][0];

        m_energyPerfFlux[iperf] += flux * res_enthalpy;

        // energy equation derivatives WRT res P & T
        m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dP] += dFlux[TAG::RES][CP_Deriv::dP] * res_enthalpy +
                                                            flux *  m_dResEnthalpy[er][esr][ei][0][Deriv::dP];
        m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dT] += dFlux[TAG::RES][CP_Deriv::dT] * res_enthalpy +
                                                            flux *  m_dResEnthalpy[er][esr][ei][0][Deriv::dT];
        // energy equation derivatives WRT well P
        m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dP] += dFlux[TAG::WELL][CP_Deriv::dP] * res_enthalpy;
        m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dT] += dFlux[TAG::WELL][CP_Deriv::dT] * res_enthalpy;

      }
      else     // ** reservoir cell is downstream
      {
        real64 const wellelem_enthalpy = m_wellElemEnthalpy[iwelem][0];
        m_energyPerfFlux[iperf] += flux * wellelem_enthalpy;

        // energy equation derivatives WRT res P & T
        m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dP] += dFlux[TAG::RES][CP_Deriv::dP] * wellelem_enthalpy;
        m_dEnergyPerfFlux[iperf][TAG::RES][CP_Deriv::dT] += dFlux[TAG::RES][CP_Deriv::dT] * wellelem_enthalpy;

        m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dP] += dFlux[TAG::WELL][CP_Deriv::dP] * wellelem_enthalpy
                                                             +  flux * m_dWellElemEnthalpy[iwelem][0][Deriv::dP];
        m_dEnergyPerfFlux[iperf][TAG::WELL][CP_Deriv::dT] += dFlux[TAG::WELL][CP_Deriv::dT] * wellelem_enthalpy
                                                             +  flux * m_dWellElemEnthalpy[iwelem][0][Deriv::dT];
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
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_wellElemEnthalpy;
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const m_dWellElemEnthalpy;

  /// Views on energy flux
  arrayView1d< real64 > const m_energyPerfFlux;
  arrayView3d< real64 > const m_dEnergyPerfFlux;

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;

  /// Views on phase enthalpies
  ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const m_resEnthalpy;
  ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const m_dResEnthalpy;

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
   * @param[in] flowSolverName the flow solver
   * @param[in] perforationData container with perforation variables
   * @param[in] subRegion well subregion
   * @param[in] fluid  subregion fluid model
   * @param[in] subregion well subregion
   * @param[in] elemManager  region element managers
   */
  template< typename POLICY >
  static void
  createAndLaunch( string const flowSolverName,
                   PerforationData * const perforationData,
                   ElementSubRegionBase const & subRegion,
                   SingleFluidBase const & fluid,
                   ElementRegionManager & elemManager )
  {
    integer constexpr IS_THERMAL = 1;
    using kernelType = PerforationFluxKernel< IS_THERMAL >;
    typename kernelType::SinglePhaseFlowAccessors singlePhaseFlowAccessors( elemManager, flowSolverName );
    typename kernelType::SingleFluidAccessors singleFluidAccessors( elemManager, flowSolverName );
    typename kernelType::ThermalSinglePhaseFlowAccessors thermalSinglePhaseFlowAccessors( elemManager, flowSolverName );
    typename kernelType::ThermalSingleFluidAccessors thermalSingleFluidAccessors( elemManager, flowSolverName );
    kernelType kernel( perforationData, subRegion, fluid, singlePhaseFlowAccessors, singleFluidAccessors, thermalSinglePhaseFlowAccessors, thermalSingleFluidAccessors );
    kernelType::template launch< POLICY >( perforationData->size(), kernel );
  }
};

}   // end namespace thermalPerforationFluxKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEPERFORATIONFLUXLKERNELS_HPP
