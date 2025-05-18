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
 * @file AquiferBCKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_AQUIFERBCKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_AQUIFERBCKERNEL_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

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
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

  template< integer NC >
  GEOS_HOST_DEVICE
  inline
  static void
    compute( integer const numPhases,
             integer const ipWater,
             bool const allowAllPhasesIntoAquifer,
             real64 const aquiferVolFlux,
             real64 const dAquiferVolFlux_dPres,
             real64 const aquiferWaterPhaseDens,
             arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
             arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > phaseDens,
             arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > dPhaseDens,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
             arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac,
             arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > phaseCompFrac,
             arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac,
             arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens,
             real64 const dt,
             real64 ( &localFlux )[NC],
             real64 ( &localFluxJacobian )[NC][NC+1] );

  template< integer NC >
  static void
  launch( integer const numPhases,
          integer const ipWater,
          bool const allowAllPhasesIntoAquifer,
          integer const useTotalMassEquation,
          BoundaryStencil const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
          real64 const aquiferWaterPhaseDens,
          arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
          MeshLevel & mesh,
          real64 const timeAtBeginningOfStep,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};

/******************************** AquiferBCKernel ********************************/

template< integer NC >
GEOS_HOST_DEVICE
void
AquiferBCKernel::
  compute( integer const numPhases,
           integer const ipWater,
           bool const allowAllPhasesIntoAquifer,
           real64 const aquiferVolFlux,
           real64 const dAquiferVolFlux_dPres,
           real64 const aquiferWaterPhaseDens,
           arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
           arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > phaseDens,
           arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > dPhaseDens,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac,
           arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > phaseCompFrac,
           arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens,
           real64 const dt,
           real64 (& localFlux)[NC],
           real64 (& localFluxJacobian)[NC][NC+1] )
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  real64 dProp_dC[NC]{};
  real64 dPhaseFlux_dCompDens[NC]{};

  if( aquiferVolFlux > 0 ) // aquifer is upstream
  {
    // in this case, we assume that:
    //    - only the water phase is present in the aquifer
    //    - the aquifer water phase composition is constant

    for( integer ic = 0; ic < NC; ++ic )
    {
      real64 const phaseFlux = aquiferVolFlux * aquiferWaterPhaseDens;
      localFlux[ic] -= dt * phaseFlux * aquiferWaterPhaseCompFrac[ic];
      localFluxJacobian[ic][0] -= dt * dAquiferVolFlux_dPres * aquiferWaterPhaseDens * aquiferWaterPhaseCompFrac[ic];
    }
  }
  else // reservoir is upstream
  {
    for( integer ip = 0; ip < numPhases; ++ip )
    {

      // Why two options below:
      //   - The aquifer model assumes single-phase water flow, so ideally, we should only allow water phase flow from the reservoir to the
      // aquifer
      //   - But, if/when the CO2 plume reaches the reservoir cell connected to the aquifer and saturates it, the aquifer flux becomes zero
      //     if we don't let some CO2 go into the aquifer

      if( ip == ipWater || allowAllPhasesIntoAquifer )
      {
        real64 const phaseDensVolFrac = phaseDens[ip] * phaseVolFrac[ip];
        real64 const phaseFlux = aquiferVolFlux * phaseDensVolFrac;
        real64 const dPhaseFlux_dPres = dAquiferVolFlux_dPres * phaseDensVolFrac
                                        + aquiferVolFlux * ( dPhaseDens[ip][Deriv::dP] * phaseVolFrac[ip] + phaseDens[ip] * dPhaseVolFrac[ip][Deriv::dP] );

        applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens[ip], dProp_dC, Deriv::dC );
        for( integer ic = 0; ic < NC; ++ic )
        {
          dPhaseFlux_dCompDens[ic] = aquiferVolFlux * ( dProp_dC[ic] * phaseVolFrac[ip] + phaseDens[ip] * dPhaseVolFrac[ip][Deriv::dC+ic] );
        }

        for( integer ic = 0; ic < NC; ++ic )
        {
          localFlux[ic] -= dt * phaseFlux * phaseCompFrac[ip][ic];
          localFluxJacobian[ic][0] -= dt * ( dPhaseFlux_dPres * phaseCompFrac[ip][ic] + phaseFlux * dPhaseCompFrac[ip][ic][Deriv::dP] );

          applyChainRule( NC, dCompFrac_dCompDens, dPhaseCompFrac[ip][ic], dProp_dC, Deriv::dC );
          for( integer jc = 0; jc < NC; ++jc )
          {
            localFluxJacobian[ic][jc+1] -= dt * ( dPhaseFlux_dCompDens[jc] * phaseCompFrac[ip][ic] + phaseFlux * dProp_dC[jc] );
          }
        }
      }
    }
  }
}

template< integer NC >
void
AquiferBCKernel::
  launch( integer const numPhases,
          integer const ipWater,
          bool const allowAllPhasesIntoAquifer,
          integer const useTotalMassEquation,
          BoundaryStencil const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
          real64 const aquiferWaterPhaseDens,
          arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
          MeshLevel & mesh,
          real64 const timeAtBeginningOfStep,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{

  using namespace compositionalMultiphaseUtilities;
  using Order = BoundaryStencil::Order;

  CompFlowAccessors compFlowAccessors( mesh.getElemManager(), "compFlowAccessors" );
  MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), "multiFluidAccessors" );

  ElementViewConst< arrayView1d< integer const > > const & ghostRank = compFlowAccessors.get( fields::ghostRank{} );
  ElementViewConst< arrayView1d< real64 const > > const & pres = compFlowAccessors.get( fields::flow::pressure{} );
  ElementViewConst< arrayView1d< real64 const > > const & presOld = compFlowAccessors.get( fields::flow::pressure_n{} );
  ElementViewConst< arrayView1d< real64 const > > const & gravCoef = compFlowAccessors.get( fields::flow::gravityCoefficient{} );
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac = compFlowAccessors.get( fields::flow::phaseVolumeFraction{} );
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac = compFlowAccessors.get( fields::flow::dPhaseVolumeFraction{} );
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens = compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity{} );
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens = multiFluidAccessors.get( fields::multifluid::phaseDensity{} );
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseDens = multiFluidAccessors.get( fields::multifluid::dPhaseDensity{} );
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac = multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} );
  ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac = multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction{} );

  BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
  BoundaryStencil::WeightContainerViewConstType const & weight = stencil.getWeights();

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    constexpr integer NDOF = NC + 1;

    // working arrays
    globalIndex dofColIndices[NDOF]{};
    real64 localFlux[NC]{};
    real64 localFluxJacobian[NC][NDOF]{};

    localIndex const er  = seri( iconn, Order::ELEM );
    localIndex const esr = sesri( iconn, Order::ELEM );
    localIndex const ei  = sefi( iconn, Order::ELEM );
    real64 const areaFraction = weight( iconn, Order::ELEM );

    // compute the aquifer influx rate using the pressure influence function and the aquifer props
    real64 dAquiferVolFlux_dPres = 0.0;
    real64 const aquiferVolFlux = aquiferBCWrapper.compute( timeAtBeginningOfStep,
                                                            dt,
                                                            pres[er][esr][ei],
                                                            presOld[er][esr][ei],
                                                            gravCoef[er][esr][ei],
                                                            areaFraction,
                                                            dAquiferVolFlux_dPres );

    // compute the phase/component aquifer flux
    AquiferBCKernel::compute< NC >( numPhases,
                                    ipWater,
                                    allowAllPhasesIntoAquifer,
                                    aquiferVolFlux,
                                    dAquiferVolFlux_dPres,
                                    aquiferWaterPhaseDens,
                                    aquiferWaterPhaseCompFrac,
                                    phaseDens[er][esr][ei][0],
                                    dPhaseDens[er][esr][ei][0],
                                    phaseVolFrac[er][esr][ei],
                                    dPhaseVolFrac[er][esr][ei],
                                    phaseCompFrac[er][esr][ei][0],
                                    dPhaseCompFrac[er][esr][ei][0],
                                    dCompFrac_dCompDens[er][esr][ei],
                                    dt,
                                    localFlux,
                                    localFluxJacobian );

    // populate dof indices
    globalIndex const offset = dofNumber[er][esr][ei];
    for( integer jdof = 0; jdof < NDOF; ++jdof )
    {
      dofColIndices[jdof] = offset + jdof;
    }

    if( useTotalMassEquation > 0 )
    {
      // Apply equation/variable change transformation(s)
      real64 work[NDOF];
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NDOF, localFluxJacobian, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, localFlux );
    }

    // Add to residual/jacobian
    if( ghostRank[er][esr][ei] < 0 )
    {
      globalIndex const globalRow = dofNumber[er][esr][ei];
      localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
      GEOS_ASSERT_GE( localRow, 0 );
      GEOS_ASSERT_GT( localMatrix.numRows(), localRow + NC );

      for( integer ic = 0; ic < NC; ++ic )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow + ic], localFlux[ic] );
        localMatrix.addToRow< parallelDeviceAtomic >( localRow + ic,
                                                      dofColIndices,
                                                      localFluxJacobian[ic],
                                                      NDOF );
      }
    }
  } );
}

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_AQUIFERBCKERNEL_HPP
