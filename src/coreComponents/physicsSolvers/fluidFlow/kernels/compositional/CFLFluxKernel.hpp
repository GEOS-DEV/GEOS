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
 * @file CFLFluxKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_CFLFLUXKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_CFLFLUXKERNEL_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

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
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::phaseViscosity,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::phaseMassDensity,
                              fields::multifluid::phaseCompFraction >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;


  using RelPermAccessors =
    StencilMaterialAccessors< constitutive::RelativePermeabilityBase, fields::relperm::phaseRelPerm >;

  template< integer NC >
  GEOS_HOST_DEVICE inline static void
  compute( integer const numPhases,
           integer const checkPhasePresenceInGravity,
           localIndex const stencilSize,
           real64 const dt,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           real64 const (&transmissibility)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
           ElementViewConst< arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > > const & phaseRelPerm,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseVisc,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux,
           ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux );

  GEOS_HOST_DEVICE inline static void
  calculateMeanDensity( integer const checkPhasePresenceInGravity,
                        integer const ip, localIndex const stencilSize,
                        arraySlice1d< localIndex const > const seri,
                        arraySlice1d< localIndex const > const sesri,
                        arraySlice1d< localIndex const > const sei,
                        ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                        ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                        real64 & densMean );

  template< integer NC, typename STENCILWRAPPER_TYPE >
  static void
  launch( integer const numPhases,
          integer const checkPhasePresenceInGravity,
          real64 const dt,
          STENCILWRAPPER_TYPE const & stencil,
          MeshLevel & mesh );
};

/******************************** CFLFluxKernel ********************************/

template< integer NC >
GEOS_HOST_DEVICE
inline
void
CFLFluxKernel::
  compute( integer const numPhases,
           integer const checkPhasePresenceInGravity,
           localIndex const stencilSize,
           real64 const dt,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           real64 const (&transmissibility)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
           ElementViewConst< arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > > const & phaseRelPerm,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseVisc,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux,
           ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux )
{
  // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
  for( integer ip = 0; ip < numPhases; ++ip )
  {

    // clear working arrays
    real64 densMean{};

    // create local work arrays
    real64 presGrad{};
    real64 gravHead{};

    // calculate quantities on primary connected cells
    calculateMeanDensity( checkPhasePresenceInGravity, ip, stencilSize, seri, sesri, sei, phaseVolFrac, phaseMassDens, densMean );

    //***** calculation of phase volumetric flux *****

    // compute potential difference MPFA-style
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      presGrad += transmissibility[i] * pres[er][esr][ei];
      gravHead += transmissibility[i] * densMean * gravCoef[er][esr][ei];
    }

    // *** upwinding ***

    // compute phase potential gradient
    real64 const potGrad = presGrad - gravHead;

    // choose upstream cell
    localIndex const k_up = (potGrad >= 0) ? 0 : 1;

    localIndex const er_up  = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up  = sei[k_up];

    // compute the phase flux only if the phase is present
    bool const phaseExists = (phaseVolFrac[er_up][esr_up][ei_up][ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    real64 const mobility = phaseRelPerm[er_up][esr_up][ei_up][0][ip] / phaseVisc[er_up][esr_up][ei_up][0][ip];

    // increment the phase (volumetric) outflux of the upstream cell
    real64 const absPhaseFlux = LvArray::math::abs( dt * mobility * potGrad );
    RAJA::atomicAdd( parallelDeviceAtomic{}, &phaseOutflux[er_up][esr_up][ei_up][ip], absPhaseFlux );

    // increment the component (mass/molar) outflux of the upstream cell
    for( integer ic = 0; ic < NC; ++ic )
    {
      real64 const absCompFlux = phaseCompFrac[er_up][esr_up][ei_up][0][ip][ic]
                                 * phaseDens[er_up][esr_up][ei_up][0][ip]
                                 * absPhaseFlux;
      RAJA::atomicAdd( parallelDeviceAtomic{}, &compOutflux[er_up][esr_up][ei_up][ic], absCompFlux );
    }
  }
}

GEOS_HOST_DEVICE
inline
void
CFLFluxKernel::calculateMeanDensity( integer const checkPhasePresenceInGravity, integer const ip, localIndex const stencilSize,
                                     arraySlice1d< localIndex const > const seri,
                                     arraySlice1d< localIndex const > const sesri,
                                     arraySlice1d< localIndex const > const sei,
                                     ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                                     ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                     real64 & densMean )
{
  integer denom = 0;
  for( localIndex i = 0; i < stencilSize; ++i )
  {
    localIndex const er = seri[i];
    localIndex const esr = sesri[i];
    localIndex const ei = sei[i];

    bool const phaseExists = (phaseVolFrac[er][esr][ei][ip] > 0);
    if( checkPhasePresenceInGravity && !phaseExists )
    {
      continue;
    }

    // average density across the face
    densMean += phaseMassDens[er][esr][ei][0][ip];
    denom++;
  }
  if( denom > 1 )
  {
    densMean /= denom;
  }
}
template< integer NC, typename STENCILWRAPPER_TYPE >
void CFLFluxKernel::
  launch( integer const numPhases,
          integer const checkPhasePresenceInGravity,
          real64 const dt,
          STENCILWRAPPER_TYPE const & stencilWrapper,
          MeshLevel & mesh )
{
  CompFlowAccessors compFlowAccessors( mesh.getElemManager(), "compFlowAccessors" );
  MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), "multiFluidAccessors" );
  PermeabilityAccessors permeabilityAccessors( mesh.getElemManager(), "permeabilityAccessors" );
  RelPermAccessors relPermAccessors( mesh.getElemManager(), "relPermAccessors" );

  // TODO: find a way to compile with this modifiable accessors in CompFlowAccessors, and remove them from here
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64, compflow::USD_PHASE > > const phaseOutfluxAccessor =
    mesh.getElemManager().constructViewAccessor< array2d< real64, compflow::LAYOUT_PHASE >,
                                                 arrayView2d< real64, compflow::USD_PHASE > >( fields::flow::phaseOutflux::key() );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64, compflow::USD_COMP > > const compOutfluxAccessor =
    mesh.getElemManager().constructViewAccessor< array2d< real64, compflow::LAYOUT_COMP >,
                                                 arrayView2d< real64, compflow::USD_COMP > >( fields::flow::componentOutflux::key() );

  ElementViewConst< arrayView1d< real64 const > > const & pres = compFlowAccessors.get( fields::flow::pressure{} );
  ElementViewConst< arrayView1d< real64 const > > const & gravCoef = compFlowAccessors.get( fields::flow::gravityCoefficient{} );
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac = compFlowAccessors.get( fields::flow::phaseVolumeFraction{} );
  ElementViewConst< arrayView3d< real64 const > > const & permeability = permeabilityAccessors.get( fields::permeability::permeability{} );
  ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres = permeabilityAccessors.get( fields::permeability::dPerm_dPressure{} );
  ElementViewConst< arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > > const & phaseRelPerm = relPermAccessors.get( fields::relperm::phaseRelPerm{} );
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseVisc = multiFluidAccessors.get( fields::multifluid::phaseViscosity{} );
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens = multiFluidAccessors.get( fields::multifluid::phaseDensity{} );
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens = multiFluidAccessors.get( fields::multifluid::phaseMassDensity{} );
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac = multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} );
  ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux = phaseOutfluxAccessor.toNestedView();
  ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux = compOutfluxAccessor.toNestedView();

  typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    // compute transmissibility
    real64 transmissibility[STENCILWRAPPER_TYPE::maxNumConnections][2];
    real64 dTrans_dPres[STENCILWRAPPER_TYPE::maxNumConnections][2];

    stencilWrapper.computeWeights( iconn,
                                   permeability,
                                   dPerm_dPres,
                                   transmissibility,
                                   dTrans_dPres );

    CFLFluxKernel::compute< NC >( numPhases,
                                  checkPhasePresenceInGravity,
                                  sei[iconn].size(),
                                  dt,
                                  seri[iconn],
                                  sesri[iconn],
                                  sei[iconn],
                                  transmissibility[0],
                                  pres,
                                  gravCoef,
                                  phaseVolFrac,
                                  phaseRelPerm,
                                  phaseVisc,
                                  phaseDens,
                                  phaseMassDens,
                                  phaseCompFrac,
                                  phaseOutflux,
                                  compOutflux );
  } );
}

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_CFLFLUXKERNEL_HPP
