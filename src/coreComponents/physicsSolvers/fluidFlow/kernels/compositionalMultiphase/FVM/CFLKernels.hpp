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
 * @file CFLFluxKernels.hpp
 */

#ifndef GEOSX_CFLKERNELS_HPP
#define GEOSX_CFLKERNELS_HPP

#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/fields/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/fields/FlowSolverBaseFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

using namespace constitutive;

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
  GEOS_HOST_DEVICE
  inline
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
  GEOS_HOST_DEVICE
  inline
  static void
  computePhaseCFL( real64 const & poreVol,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
                   arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > phaseRelPerm,
                   arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > dPhaseRelPerm_dPhaseVolFrac,
                   arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseVisc,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseOutflux,
                   real64 & phaseCFLNumber );

  template< integer NC >
  GEOS_HOST_DEVICE
  inline
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

}

}

#endif //GEOSX_CFLKERNELS_HPP
