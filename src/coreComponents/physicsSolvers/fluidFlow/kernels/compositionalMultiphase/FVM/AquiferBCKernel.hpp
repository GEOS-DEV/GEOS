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
 * @file AquiferBCKernel.hpp
 */

#ifndef GEOSX_AQUIFERBCKERNEL_HPP
#define GEOSX_AQUIFERBCKERNEL_HPP

#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/fields/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/fields/FlowSolverBaseFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

using namespace constitutive;

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
  using ElementViewConst = geos::ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using CompFlowAccessors =
    geos::StencilAccessors< geos::fields::ghostRank,
                            geos::fields::flow::pressure,
                            geos::fields::flow::pressure_n,
                            geos::fields::flow::gravityCoefficient,
                            geos::fields::flow::phaseVolumeFraction,
                            geos::fields::flow::dPhaseVolumeFraction,
                            geos::fields::flow::dGlobalCompFraction_dGlobalCompDensity >;

  using MultiFluidAccessors =
    geos::StencilMaterialAccessors< geos::constitutive::MultiFluidBase,
                                    geos::fields::multifluid::phaseDensity,
                                    geos::fields::multifluid::dPhaseDensity,
                                    geos::fields::multifluid::phaseCompFraction,
                                    geos::fields::multifluid::dPhaseCompFraction >;

  template< geos::integer NC >
  GEOS_HOST_DEVICE
  inline
  static void
    compute( geos::integer const numPhases,
             geos::integer const ipWater,
             bool const allowAllPhasesIntoAquifer,
             geos::real64 const & aquiferVolFlux,
             geos::real64 const & dAquiferVolFlux_dPres,
             geos::real64 const & aquiferWaterPhaseDens,
             geos::arrayView1d< geos::real64 const > const & aquiferWaterPhaseCompFrac,
             arraySlice1d< geos::real64 const, multifluid::USD_PHASE - 2 > phaseDens,
             arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
             arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac,
             arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac,
             arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac,
             arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens,
             real64 const & dt,
             real64 ( &localFlux )[NC],
             real64 ( &localFluxJacobian )[NC][NC + 1] );

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

}

}

#endif //GEOSX_AQUIFERBCKERNEL_HPP
