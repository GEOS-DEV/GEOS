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
 * @file CompositionalMultiphaseFVMUpwindUtilities.hpp
 */
#ifndef GEOSX_COMPOSITIONALMULTIPHASEFVMUPWINDUTILITIES_HPP
#define GEOSX_COMPOSITIONALMULTIPHASEFVMUPWINDUTILITIES_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geosx
{
namespace CompositionalMultiphaseFVMUpwindUtilities
{

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;



/*** struct holding actual static's ***/
struct UpwindHelpers
{

  using Deriv = constitutive::multifluid::DerivativeOffset;

  /// Minimal mobility for fractional flow construction
  static constexpr real64 minTotMob = 1e-12;

  /***
   * @brief section grouping pure utilities functions
   */

  /**
   * @brief Form the PhasePotentialUpwind from pressure gradient and gravitational heads (legacy)
   * @tparam numComp number of components
   * @tparam NUM_ELEMS numberof elements involve in the stencil connexion
   * @tparam MAX_STENCIL maximum number of points in the stencil
   * @param numPhase number of phases
   * @param ip index of the treated phase
   * @param stencilSize number of points in the stencil
   * @param seri arraySlice of the stencil implied element region index
   * @param sesri arraySlice of the stencil implied element subregion index
   * @param sei arraySlice of the stencil implied element index
   * @param stencilWeights weights associated with elements in the stencil
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOSX_HOST_DEVICE
  static void computePPUPhaseFlux( localIndex const numPhase,
                                   localIndex const ip,
                                   const localIndex (& seri)[numFluxSupportPoints],
                                   const localIndex (& sesri)[numFluxSupportPoints],
                                   const localIndex (& sei)[numFluxSupportPoints],
                                   real64 const (&transmissibility)[2],
                                   real64 const (&dTrans_dPres)[2],
                                   ElementViewConst< arrayView1d< real64 const > > const & pres,
                                   ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                   ElementViewConst< arrayView2d< real64, compflow::USD_PHASE > > const & phaseMob,
                                   ElementViewConst< arrayView3d< real64, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                   ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                   ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                   ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                   ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                   ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                   ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                   integer const capPressureFlag,
                                   localIndex( &k_up ),
                                   real64 & potGrad,
                                   real64( &phaseFlux ),
                                   real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
                                   real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp] )
  {

    real64 densMean{};
    real64 dDensMean_dP[numFluxSupportPoints]{};
    real64 dDensMean_dC[numFluxSupportPoints][numComp]{};

    real64 presGrad{};
    real64 dPresGrad_dP[numFluxSupportPoints]{};
    real64 dPresGrad_dC[numFluxSupportPoints][numComp]{};

    real64 gravHead{};
    real64 dGravHead_dP[numFluxSupportPoints]{};
    real64 dGravHead_dC[numFluxSupportPoints][numComp]{};

    real64 dCapPressure_dC[numComp]{};
// Working array
    real64 dProp_dC[numComp]{};

// calculate quantities on primary connected cells
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      // density
      real64 const density = phaseMassDens[er][esr][ei][0][ip];
      real64 const dDens_dP = dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];

      applyChainRule( numComp,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens[er][esr][ei][0][ip],
                      dProp_dC,
                      Deriv::dC
                      );

// average density and derivatives
      densMean += 0.5 * density;
      dDensMean_dP[i] = 0.5 * dDens_dP;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dC[i][jc] = 0.5 * dProp_dC[jc];
      }
    }
// compute potential difference
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      // capillary pressure
      real64 capPressure = 0.0;
      real64 dCapPressure_dP = 0.0;

      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        dCapPressure_dC[ic] = 0.0;
      }

      if( capPressureFlag )
      {
        capPressure = phaseCapPressure[er][esr][ei][0][ip];

        for( localIndex jp = 0; jp < numPhase; ++jp )
        {
          real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dP];

          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dC + jc];
          }
        }
      }

      presGrad += transmissibility[i] * (pres[er][esr][ei] - capPressure);
      dPresGrad_dP[i] += transmissibility[i] * (1 - dCapPressure_dP)
                         + dTrans_dPres[i] * (pres[er][esr][ei] - capPressure);
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dPresGrad_dC[i][jc] += -transmissibility[i] * dCapPressure_dC[jc];
      }

      real64 const gravD = transmissibility[i] * gravCoef[er][esr][ei];
      real64 const dGravD_dP = dTrans_dPres[i] * gravCoef[er][esr][ei];

      // the density used in the potential difference is always a mass density
      // unlike the density used in the phase mobility, which is a mass density
      // if useMass == 1 and a molar density otherwise
      gravHead += densMean * gravD;

      // need to add contributions from both cells the mean density depends on
      for( localIndex j = 0; j < numFluxSupportPoints; ++j )
      {
        dGravHead_dP[j] += dDensMean_dP[j] * gravD + densMean * dGravD_dP;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
        }
      }
    }

    // compute phase potential gradient
    potGrad = presGrad - gravHead;

    // choose upstream cell
    k_up = (potGrad >= 0) ? 0 : 1;

    localIndex er_up = seri[k_up];
    localIndex esr_up = sesri[k_up];
    localIndex ei_up = sei[k_up];

    real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];
// pressure gradient depends on all points in the stencil
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc];
      }

    }

    // gravitational head depends only on the two cells connected (same as mean density)
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] -= dGravHead_dP[ke];
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[ke][jc] -= dGravHead_dC[ke][jc];
      }
    }

    // compute the phase flux and derivatives using upstream cell mobility
    phaseFlux = mobility * potGrad;

    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] *= mobility;

      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[ke][jc] *= mobility;
      }
    }

    real64 const dMob_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
    arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 >
    const dMob_dC = dPhaseMob[er_up][esr_up][ei_up][ip];

    // add contribution from upstream cell mobility derivatives
    dPhaseFlux_dP[k_up] += dMob_dP * potGrad;

    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dPhaseFlux_dC[k_up][jc] += dMob_dC[Deriv::dC + jc] * potGrad;
    }

  }


  /**
   * @brief Distribute phaseFlux onto component fluxes
   * @tparam numComp number of components
   * @tparam MAX_STENCIL maximum number of points in the stencil
   * @param ip concerned phase index
   * @param k_up upwind direction of the phase
   * @param seri arraySlice of the stencil implied element region index
   * @param sesri arraySlice of the stencil implied element subregion index
   * @param sei arraySlice of the stencil implied element index
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOSX_HOST_DEVICE
  static void computePhaseComponentFlux( localIndex const ip,
                                         localIndex const k_up,
                                         const localIndex (& seri)[numFluxSupportPoints],
                                         const localIndex (& sesri)[numFluxSupportPoints],
                                         const localIndex (& sei)[numFluxSupportPoints],
                                         ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                                         ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > >
                                         const & dPhaseCompFrac,
                                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                         real64 const & phaseFlux,
                                         real64 const (&dPhaseFlux_dPres)[numFluxSupportPoints],
                                         real64 const (&dPhaseFlux_dComp)[numFluxSupportPoints][numComp],
                                         real64 (& compFlux)[numComp],
                                         real64 (& dCompFlux_dPres)[numFluxSupportPoints][numComp],
                                         real64 (& dCompFlux_dComp)[numFluxSupportPoints][numComp][numComp] )
  {
    /*update phaseComp from grav part*/
    localIndex const er_up = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up = sei[k_up];

    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 3 > phaseCompFracSub =
      phaseCompFrac[er_up][esr_up][ei_up][0][ip];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 3 > dPhaseCompFracSub =
      dPhaseCompFrac[er_up][esr_up][ei_up][0][ip];

    real64 dProp_dC[numComp]{};

    // compute component fluxes and derivatives using upstream cell composition
    for( localIndex ic = 0; ic < numComp; ++ic )
    {
      real64 const ycp = phaseCompFracSub[ic];
      compFlux[ic] += phaseFlux * ycp;

      // derivatives stemming from phase flux
      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dCompFlux_dPres[ke][ic] += dPhaseFlux_dPres[ke] * ycp;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dCompFlux_dComp[ke][ic][jc] += dPhaseFlux_dComp[ke][jc] * ycp;
        }
      }

      // additional derivatives stemming from upstream cell phase composition
      dCompFlux_dPres[k_up][ic] += phaseFlux * dPhaseCompFracSub[ic][Deriv::dP];

      // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component
      // densities
      applyChainRule( numComp, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseCompFracSub[ic],
                      dProp_dC, Deriv::dC );
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dCompFlux_dComp[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
      }
    }
  }


  /***
   * @brief section grouping acting functions that are used for fractional flow formulation
   */

  /**
   * @brief Function returning upwinded mobility (and derivatives)  of specified phase as well as uppwind direction
   *        according to specified upwind scheme
   * @tparam NC number of components
   * @tparam UpwindScheme Desciption of how to construct potential used to decide upwind direction
   * @param numPhase total number of phases
   * @param ip concerned phase index
   * @param stencilSizei number of points in the stencil
   * @param seri arraySlice of the stencil implied element region index
   * @param sesri arraySlice of the stencil implied element subregion index
   * @param sei arraySlice of the stencil implied element index
   * @param stencilWeights weights associated with elements in the stencil
   * @param totFlux total flux signed value
   */
  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  static void
  upwindMobilityViscous( localIndex const numPhase,
                         localIndex const ip,
                         localIndex const (&seri)[numFluxSupportPoints],
                         localIndex const (&sesri)[numFluxSupportPoints],
                         localIndex const (&sei)[numFluxSupportPoints],
                         real64 const (&transmissibility)[2],
                         real64 const (&dTrans_dPres)[2],
                         real64 const totFlux,  //in fine should be a ElemnetViewConst once seq form are in place
                         ElementViewConst< arrayView1d< real64 const > > const & pres,
                         ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                         ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                         ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                         ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                         ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                         integer const capPressureFlag,
                         localIndex & upwindDir,
                         real64 & mobility,
                         real64( &dMobility_dP),
                         real64 ( & dMobility_dC)[numComp]
                         )
  {

    //reinit
    mobility = 0.0;
    dMobility_dP = 0.0;
    for( localIndex ic = 0; ic < numComp; ++ic )
    {
      dMobility_dC[ic] = 0.0;
    }

    UPWIND scheme;
    scheme.template getUpwindDirectionViscous< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                        ip,
                                                                                        seri,
                                                                                        sesri,
                                                                                        sei,
                                                                                        transmissibility,
                                                                                        dTrans_dPres,
                                                                                        totFlux,
                                                                                        pres,
                                                                                        gravCoef,
                                                                                        phaseMob,
                                                                                        dCompFrac_dCompDens,
                                                                                        phaseMassDens,
                                                                                        dPhaseMassDens,
                                                                                        dPhaseVolFrac,
                                                                                        phaseCapPressure,
                                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                                        capPressureFlag,
                                                                                        upwindDir );

    localIndex const er_up = seri[upwindDir];
    localIndex const esr_up = sesri[upwindDir];
    localIndex const ei_up = sei[upwindDir];

    if( std::fabs( phaseMob[er_up][esr_up][ei_up][ip] ) > 1e-20 )
    {
      mobility = phaseMob[er_up][esr_up][ei_up][ip];
      dMobility_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        dMobility_dC[ic] = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dC + ic];
      }
    }
  }

  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  static void
  upwindMobilityGravity( localIndex const numPhase,
                         localIndex const ip,
                         localIndex const (&seri)[numFluxSupportPoints],
                         localIndex const (&sesri)[numFluxSupportPoints],
                         localIndex const (&sei)[numFluxSupportPoints],
                         real64 const (&transmissibility)[2],
                         real64 const (&dTrans_dPres)[2],
                         real64 const totFlux,    //in fine should be a ElemnetViewConst once seq form are in place
                         ElementViewConst< arrayView1d< real64 const > > const & pres,
                         ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                         ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                         ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                         ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                         ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                         integer const capPressureFlag,
                         localIndex & upwindDir,
                         real64 & mobility,
                         real64( &dMobility_dP),
                         real64 ( & dMobility_dC)[numComp]
                         )
  {

    //reinit
    mobility = 0.0;
    dMobility_dP = 0.0;
    for( localIndex ic = 0; ic < numComp; ++ic )
    {
      dMobility_dC[ic] = 0.0;
    }

    UPWIND scheme;
    scheme.template getUpwindDirectionGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                        ip,
                                                                                        seri,
                                                                                        sesri,
                                                                                        sei,
                                                                                        transmissibility,
                                                                                        dTrans_dPres,
                                                                                        totFlux,
                                                                                        pres,
                                                                                        gravCoef,
                                                                                        phaseMob,
                                                                                        dCompFrac_dCompDens,
                                                                                        phaseMassDens,
                                                                                        dPhaseMassDens,
                                                                                        dPhaseVolFrac,
                                                                                        phaseCapPressure,
                                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                                        capPressureFlag,
                                                                                        upwindDir );

    localIndex const er_up = seri[upwindDir];
    localIndex const esr_up = sesri[upwindDir];
    localIndex const ei_up = sei[upwindDir];

    if( std::fabs( phaseMob[er_up][esr_up][ei_up][ip] ) > 1e-20 )
    {
      mobility = phaseMob[er_up][esr_up][ei_up][ip];
      dMobility_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        dMobility_dC[ic] = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dC + ic];
      }
    }
  }

  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  static void
  upwindMobilityCapillary( localIndex const numPhase,
                           localIndex const ip,
                           localIndex const (&seri)[numFluxSupportPoints],
                           localIndex const (&sesri)[numFluxSupportPoints],
                           localIndex const (&sei)[numFluxSupportPoints],
                           real64 const (&transmissibility)[2],
                           real64 const (&dTrans_dPres)[2],
                           real64 const totFlux,  //in fine should be a ElemnetViewConst once seq form are in place
                           ElementViewConst< arrayView1d< real64 const > > const & pres,
                           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                           ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                           ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                           integer const capPressureFlag,
                           localIndex & upwindDir,
                           real64 & mobility,
                           real64( &dMobility_dP),
                           real64 ( & dMobility_dC)[numComp]
                           )
  {

    //reinit
    mobility = 0.0;
    dMobility_dP = 0.0;
    for( localIndex ic = 0; ic < numComp; ++ic )
    {
      dMobility_dC[ic] = 0.0;
    }

    UPWIND scheme;
    scheme.template getUpwindDirectionCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                          ip,
                                                                                          seri,
                                                                                          sesri,
                                                                                          sei,
                                                                                          transmissibility,
                                                                                          dTrans_dPres,
                                                                                          totFlux,
                                                                                          pres,
                                                                                          gravCoef,
                                                                                          phaseMob,
                                                                                          dCompFrac_dCompDens,
                                                                                          phaseMassDens,
                                                                                          dPhaseMassDens,
                                                                                          dPhaseVolFrac,
                                                                                          phaseCapPressure,
                                                                                          dPhaseCapPressure_dPhaseVolFrac,
                                                                                          capPressureFlag,
                                                                                          upwindDir );

    localIndex const er_up = seri[upwindDir];
    localIndex const esr_up = sesri[upwindDir];
    localIndex const ei_up = sei[upwindDir];

    if( std::fabs( phaseMob[er_up][esr_up][ei_up][ip] ) > 1e-20 )
    {
      mobility = phaseMob[er_up][esr_up][ei_up][ip];
      dMobility_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        dMobility_dC[ic] = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dC + ic];
      }
    }
  }


  /**
   * @brief Function returning upwinded fractional flow (and derivatives) as well as upwind direction of specified phase
   *        according to specified upwind scheme
   * @tparam NC number of components
   * @tparam UpwindScheme Desciption of how to construct potential used to decide upwind direction
   * @param numPhase total number of phases
   * @param ip concerned phase index
   * @param stencilSizei number of points in the stencil
   * @param seri arraySlice of the stencil implied element region index
   * @param sesri arraySlice of the stencil implied element subregion index
   * @param sei arraySlice of the stencil implied element index
   * @param stencilWeights weights associated with elements in the stencil
   * @param totFlux total flux signed value
   */
  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  static void
  computeFractionalFlowViscous( localIndex const numPhase,
                                localIndex const ip,
                                localIndex const (&seri)[numFluxSupportPoints],
                                localIndex const (&sesri)[numFluxSupportPoints],
                                localIndex const (&sei)[numFluxSupportPoints],
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const totFlux,  //in fine should be a ElemnetViewConst once seq form are in place
                                ElementViewConst< arrayView1d< real64 const > > const & pres,
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                integer const capPressureFlag,
                                localIndex & k_up_main,
                                real64 & fractionalFlow,
                                real64 ( & dFractionalFlow_dP)[numFluxSupportPoints],
                                real64 ( & dFractionalFlow_dC)[numFluxSupportPoints][numComp]
                                )
  {
    // get var to memorized the numerator mobility properly upwinded
    real64 mainMob{};
    real64 dMMob_dP{};
    real64 dMMob_dC[numComp]{};

    real64 totMob{};
    real64 dTotMob_dP[numFluxSupportPoints]{};
    real64 dTotMob_dC[numFluxSupportPoints][numComp]{};

    //reinit
    //fractional flow too low to let the upstream phase flow
    k_up_main = -1;             //to throw error if unmodified
    fractionalFlow = 0;
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dFractionalFlow_dP[ke] = 0;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dFractionalFlow_dC[ke][jc] = 0;
      }
    }

    //Form totMob
    for( localIndex jp = 0; jp < numPhase; ++jp )
    {

      localIndex k_up;
      real64 mob{};
      real64 dMob_dP{};
      real64 dMob_dC[numComp]{};

      upwindMobilityViscous< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                      jp,
                                                                      seri,
                                                                      sesri,
                                                                      sei,
                                                                      transmissibility,
                                                                      dTrans_dPres,
                                                                      totFlux,
                                                                      pres,
                                                                      gravCoef,
                                                                      dCompFrac_dCompDens,
                                                                      phaseMassDens,
                                                                      dPhaseMassDens,
                                                                      phaseMob,
                                                                      dPhaseMob,
                                                                      dPhaseVolFrac,
                                                                      phaseCapPressure,
                                                                      dPhaseCapPressure_dPhaseVolFrac,
                                                                      capPressureFlag,
                                                                      k_up,
                                                                      mob,
                                                                      dMob_dP,
                                                                      dMob_dC );


      totMob += mob;
      dTotMob_dP[k_up] += dMob_dP;
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        dTotMob_dC[k_up][ic] += dMob_dC[ic];
      }

      if( jp == ip )
      {
        k_up_main = k_up;
        mainMob = mob;
        dMMob_dP = dMob_dP;
        for( localIndex ic = 0; ic < numComp; ++ic )
        {
          dMMob_dC[ic] = dMob_dC[ic];
        }
      }
    }

    //guard against no flow region
    if( std::fabs( mainMob ) > 1e-20 )
    {
      totMob = LvArray::math::max( totMob, minTotMob );
      fractionalFlow = mainMob / totMob;
      dFractionalFlow_dP[k_up_main] = dMMob_dP / totMob;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dFractionalFlow_dC[k_up_main][jc] = dMMob_dC[jc] / totMob;

      }

      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dFractionalFlow_dP[ke] -= fractionalFlow * dTotMob_dP[ke] / totMob;

        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dFractionalFlow_dC[ke][jc] -= fractionalFlow * dTotMob_dC[ke][jc] / totMob;
        }
      }
    }
  }


  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  static void
  computeFractionalFlowGravity( localIndex const numPhase,
                                localIndex const ip,
                                localIndex const (&seri)[numFluxSupportPoints],
                                localIndex const (&sesri)[numFluxSupportPoints],
                                localIndex const (&sei)[numFluxSupportPoints],
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const totFlux,           //in fine should be a ElemnetViewConst once seq form are in place
                                ElementViewConst< arrayView1d< real64 const > > const & pres,
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                integer const capPressureFlag,
                                localIndex & k_up_main,
                                real64 & fractionalFlow,
                                real64 ( & dFractionalFlow_dP)[numFluxSupportPoints],
                                real64 ( & dFractionalFlow_dC)[numFluxSupportPoints][numComp]
                                )
  {
    // get var to memorized the numerator mobility properly upwinded
    real64 mainMob{};
    real64 dMMob_dP{};
    real64 dMMob_dC[numComp]{};

    real64 totMob{};
    real64 dTotMob_dP[numFluxSupportPoints]{};
    real64 dTotMob_dC[numFluxSupportPoints][numComp]{};

    //reinit
    //fractional flow too low to let the upstream phase flow
    k_up_main = -1;                 //to throw error if unmodified
    fractionalFlow = 0;
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dFractionalFlow_dP[ke] = 0;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dFractionalFlow_dC[ke][jc] = 0;
      }
    }

    //Form totMob
    for( localIndex jp = 0; jp < numPhase; ++jp )
    {

      localIndex k_up;
      real64 mob{};
      real64 dMob_dP{};
      real64 dMob_dC[numComp]{};

      upwindMobilityGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                      jp,
                                                                      seri,
                                                                      sesri,
                                                                      sei,
                                                                      transmissibility,
                                                                      dTrans_dPres,
                                                                      totFlux,
                                                                      pres,
                                                                      gravCoef,
                                                                      dCompFrac_dCompDens,
                                                                      phaseMassDens,
                                                                      dPhaseMassDens,
                                                                      phaseMob,
                                                                      dPhaseMob,
                                                                      dPhaseVolFrac,
                                                                      phaseCapPressure,
                                                                      dPhaseCapPressure_dPhaseVolFrac,
                                                                      capPressureFlag,
                                                                      k_up,
                                                                      mob,
                                                                      dMob_dP,
                                                                      dMob_dC );


      totMob += mob;
      dTotMob_dP[k_up] += dMob_dP;
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        dTotMob_dC[k_up][ic] += dMob_dC[ic];
      }

      if( jp == ip )
      {
        k_up_main = k_up;
        mainMob = mob;
        dMMob_dP = dMob_dP;
        for( localIndex ic = 0; ic < numComp; ++ic )
        {
          dMMob_dC[ic] = dMob_dC[ic];
        }
      }
    }

    //guard against no flow region
    if( std::fabs( mainMob ) > 1e-20 )
    {
      totMob = LvArray::math::max( totMob, minTotMob );
      fractionalFlow = mainMob / totMob;
      dFractionalFlow_dP[k_up_main] = dMMob_dP / totMob;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dFractionalFlow_dC[k_up_main][jc] = dMMob_dC[jc] / totMob;

      }

      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dFractionalFlow_dP[ke] -= fractionalFlow * dTotMob_dP[ke] / totMob;

        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dFractionalFlow_dC[ke][jc] -= fractionalFlow * dTotMob_dC[ke][jc] / totMob;
        }
      }
    }
  }


  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  static void
  computeFractionalFlowCapillary( localIndex const numPhase,
                                  localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  real64 const (&transmissibility)[2],
                                  real64 const (&dTrans_dPres)[2],
                                  real64 const totFlux,         //in fine should be a ElemnetViewConst once seq form are in place
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                  integer const capPressureFlag,
                                  localIndex & k_up_main,
                                  real64 & fractionalFlow,
                                  real64 ( & dFractionalFlow_dP)[numFluxSupportPoints],
                                  real64 ( & dFractionalFlow_dC)[numFluxSupportPoints][numComp]
                                  )
  {
    // get var to memorized the numerator mobility properly upwinded
    real64 mainMob{};
    real64 dMMob_dP{};
    real64 dMMob_dC[numComp]{};

    real64 totMob{};
    real64 dTotMob_dP[numFluxSupportPoints]{};
    real64 dTotMob_dC[numFluxSupportPoints][numComp]{};

    //reinit
    //fractional flow too low to let the upstream phase flow
    k_up_main = -1;                 //to throw error if unmodified
    fractionalFlow = 0;
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dFractionalFlow_dP[ke] = 0;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dFractionalFlow_dC[ke][jc] = 0;
      }
    }

    //Form totMob
    for( localIndex jp = 0; jp < numPhase; ++jp )
    {

      localIndex k_up;
      real64 mob{};
      real64 dMob_dP{};
      real64 dMob_dC[numComp]{};

      upwindMobilityCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                        jp,
                                                                        seri,
                                                                        sesri,
                                                                        sei,
                                                                        transmissibility,
                                                                        dTrans_dPres,
                                                                        totFlux,
                                                                        pres,
                                                                        gravCoef,
                                                                        dCompFrac_dCompDens,
                                                                        phaseMassDens,
                                                                        dPhaseMassDens,
                                                                        phaseMob,
                                                                        dPhaseMob,
                                                                        dPhaseVolFrac,
                                                                        phaseCapPressure,
                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                        capPressureFlag,
                                                                        k_up,
                                                                        mob,
                                                                        dMob_dP,
                                                                        dMob_dC );


      totMob += mob;
      dTotMob_dP[k_up] += dMob_dP;
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        dTotMob_dC[k_up][ic] += dMob_dC[ic];
      }

      if( jp == ip )
      {
        k_up_main = k_up;
        mainMob = mob;
        dMMob_dP = dMob_dP;
        for( localIndex ic = 0; ic < numComp; ++ic )
        {
          dMMob_dC[ic] = dMob_dC[ic];
        }
      }
    }

    //guard against no flow region
    if( std::fabs( mainMob ) > 1e-20 )
    {
      totMob = LvArray::math::max( totMob, minTotMob );
      fractionalFlow = mainMob / totMob;
      dFractionalFlow_dP[k_up_main] = dMMob_dP / totMob;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dFractionalFlow_dC[k_up_main][jc] = dMMob_dC[jc] / totMob;

      }

      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dFractionalFlow_dP[ke] -= fractionalFlow * dTotMob_dP[ke] / totMob;

        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dFractionalFlow_dC[ke][jc] -= fractionalFlow * dTotMob_dC[ke][jc] / totMob;
        }
      }
    }
  }

  /**
   * @brief  Struct defining formation of potential from different Physics (flagged by enum type T) to be used
   *            in Upwind discretization schemes
   * @tparam numComp number of component
   * @tparam T the concerned physics (Viscou,Gravity or Capillary)
   * @tparam numFluxSupportPoints number of point in the stencil
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  struct computePotentialViscous
  {

    /**
     * @brief functor member computing the potential to be specialized (here empty)
     * @param numPhase total number of fluid phases
     * @param ip index of the selected fluid phase
     * @param seri cell region index array for the stencil
     * @param sesri cell subregion index array for the stencil
     * @param sei cell element region index array for the stencil
     * @param transmissibility pair of half transmissibility for the selected interface
     * @param dTrans_dPres pair of half transmissibility for the selected interface
     * @param totFlux total flux needed for fractional flow formulation
     * @param gravCoeff condensed gravity projected contribution
     * @param dCompFrac_dCompDens derivatives of component fraction wrt component densities
     * @param phaseMassDens phase mass densities
     * @param dPhaseMassDens table of stacked derivatives of phase mass densities
     * @param dPhaseVolFrac table of stacked derivatives of phase volume fractions
     * @param phaseCapPressure phase capillary pressure
     * @param dPhaseCapPressure_dPhaseVolFrac derivative of phase capillary pressure wrt phase volume fraction
     * @param[out] potentialHead head difference for the selected driving force potential
     * @param[out] dPotentialHead_dP pressure derivative of the head difference for the selected driving force potential
     * @param[out] dPotentialHead_dC component wise derivative of head difference for the selected driving force potential
     * @param[out] dProp_dComp resuting properties derivative through chain rule
     */
    GEOSX_HOST_DEVICE
    static void compute( localIndex const GEOSX_UNUSED_PARAM( numPhase ),
                         localIndex const GEOSX_UNUSED_PARAM( ip ),
                         localIndex const (&GEOSX_UNUSED_PARAM( seri ))[numFluxSupportPoints],
                         localIndex const (&GEOSX_UNUSED_PARAM( sesri ))[numFluxSupportPoints],
                         localIndex const (&GEOSX_UNUSED_PARAM( sei ))[numFluxSupportPoints],
                         real64 const (&GEOSX_UNUSED_PARAM( transmissibility ))[2],
                         real64 const (&GEOSX_UNUSED_PARAM( dTrans_dPres ))[2],
                         real64 const totFlux,
                         ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravCoef ),
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const &
                         GEOSX_UNUSED_PARAM( dCompFrac_dCompDens ),
                         ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const &
                         GEOSX_UNUSED_PARAM( phaseMassDens ),
                         ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const &
                         GEOSX_UNUSED_PARAM( dPhaseMassDens ),
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const &
                         GEOSX_UNUSED_PARAM( dPhaseVolFrac ),
                         ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const &
                         GEOSX_UNUSED_PARAM( phaseCapPressure ),
                         ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const &
                         GEOSX_UNUSED_PARAM( dPhaseCapPressure_dPhaseVolFrac ),
                         real64 & pot,
                         real64( &GEOSX_UNUSED_PARAM( dPot_dPres ))[numFluxSupportPoints],
                         real64( &GEOSX_UNUSED_PARAM( dPot_dComp ))[numFluxSupportPoints][numComp],
                         real64( &GEOSX_UNUSED_PARAM( dProp_dComp ))[numComp] )
    {
      pot = totFlux;
      //could be relevant for symmetry to include derivative

    }
  };

  /*! @copydoc computePotential
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  struct computePotentialGravity
  {
    /*! @copydoc computePotential::compute
     *
     * @brief specialization for gravitational driving forces which only relies on total flux
     */
    GEOSX_HOST_DEVICE
    static void compute( localIndex const GEOSX_UNUSED_PARAM( numPhase ),
                         localIndex const ip,
                         localIndex const (&seri)[numFluxSupportPoints],
                         localIndex const (&sesri)[numFluxSupportPoints],
                         localIndex const (&sei)[numFluxSupportPoints],
                         real64 const (&transmissibility)[2],
                         real64 const (&dTrans_dPres)[2],
                         real64 const GEOSX_UNUSED_PARAM( totFlux ),
                         ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                         ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                         ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const &
                         GEOSX_UNUSED_PARAM( dPhaseVolFrac ),
                         ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const &
                         GEOSX_UNUSED_PARAM( phaseCapPressure ),
                         ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const &
                         GEOSX_UNUSED_PARAM( dPhaseCapPressure_dPhaseVolFrac ),
                         real64 & pot,
                         real64 ( & dPot_dPres )[numFluxSupportPoints],
                         real64 (& dPot_dComp )[numFluxSupportPoints][numComp],
                         real64 ( & dProp_dComp )[numComp] )
    {
      //working arrays
      real64 densMean{};
      real64 dDensMean_dPres[numFluxSupportPoints]{};
      real64 dDensMean_dComp[numFluxSupportPoints][numComp]{};

      //init
      pot = 0.0;
      for( localIndex i = 0; i < numFluxSupportPoints; ++i )
      {
        dPot_dPres[i] = 0.0;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPot_dComp[i][jc] = 0.0;
          dProp_dComp[jc] = 0.0;
        }
      }

      //inner loop to get average density
      for( localIndex i = 0; i < numFluxSupportPoints; ++i )
      {
        localIndex const er = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei = sei[i];

        // density
        real64 const density = phaseMassDens[er][esr][ei][0][ip];
        real64 const dDens_dPres = dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];

        applyChainRule( numComp,
                        dCompFrac_dCompDens[er][esr][ei],
                        dPhaseMassDens[er][esr][ei][0][ip],
                        dProp_dComp,
                        Deriv::dC );

        // average density and derivatives
        densMean += 0.5 * density;
        dDensMean_dPres[i] = 0.5 * dDens_dPres;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dDensMean_dComp[i][jc] = 0.5 * dProp_dComp[jc];
        }
      }

      // compute potential difference MPFA-style
      for( localIndex i = 0; i < numFluxSupportPoints; ++i )
      {
        localIndex const er = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei = sei[i];

        real64 const gravD = transmissibility[i] * gravCoef[er][esr][ei];
        real64 const dGravD_dP = dTrans_dPres[i] * gravCoef[er][esr][ei];
        pot += densMean * gravD;

        // need to add contributions from both cells the mean density depends on
        for( localIndex j = 0; j < numFluxSupportPoints; ++j )
        {
          dPot_dPres[j] += dDensMean_dPres[j] * gravD + densMean * dGravD_dP;
          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dPot_dComp[j][jc] += dDensMean_dComp[j][jc] * gravD;
          }
        }
      }

    }
  };

  /*! @copydoc computePotential
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  struct computePotentialCapillary
  {
    /*! @copydoc computePotential::compute
     *
     * @brief specialization for capillary driving forces which only relies on total flux
     */
    GEOSX_HOST_DEVICE
    static void compute( localIndex const numPhase,
                         localIndex const ip,
                         localIndex const (&seri)[numFluxSupportPoints],
                         localIndex const (&sesri)[numFluxSupportPoints],
                         localIndex const (&sei)[numFluxSupportPoints],
                         real64 const (&transmissibility)[2],
                         real64 const (&dTrans_dPres)[2],
                         real64 const GEOSX_UNUSED_PARAM( totFlux ),
                         ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravCoef ),
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const &
                         GEOSX_UNUSED_PARAM( dCompFrac_dCompDens ),
                         ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const &
                         GEOSX_UNUSED_PARAM( phaseMassDens ),
                         ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const &
                         GEOSX_UNUSED_PARAM( dPhaseMassDens ),
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                         ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                         ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                         real64 & pot,
                         real64 ( & dPot_dPres)[numFluxSupportPoints],
                         real64 (& dPot_dComp)[numFluxSupportPoints][numComp],
                         real64( &GEOSX_UNUSED_PARAM( dProp_dComp ))[numComp] )
    {


      for( localIndex i = 0; i < numFluxSupportPoints; ++i )
      {
        localIndex const er = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei = sei[i];

        pot += transmissibility[i] * phaseCapPressure[er][esr][ei][0][ip];
        // need to add contributions from both cells
        for( localIndex jp = 0; jp < numPhase; ++jp )
        {

          real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dPot_dPres[i] +=
            transmissibility[i] * dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dP]
            + dTrans_dPres[i] * phaseCapPressure[er][esr][ei][0][jp];

          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dPot_dComp[i][jc] += transmissibility[i] * dCapPressure_dS *
                                 dPhaseVolFrac[er][esr][ei][jp][Deriv::dC + jc];
          }

        }

      }

    }
  };


  /// Form potential-related parts of fluxes
  /***
   *  Compute the flux resulting from the potential given by the driving force
   * @tparam numComp
   * @tparam T
   * @tparam numFluxSupportPoints
   * @tparam UPWIND Upwind Scheme selected
   * @param numPhase
   * @param ip
   * @param seri
   * @param sesri
   * @param sei
   * @param transmissibility
   * @param dTrans_dPres
   * @param totFlux
   * @param pres
   * @param gravCoef
   * @param phaseMob
   * @param dPhaseMob
   * @param dPhaseVolFrac
   * @param dCompFrac_dCompDens
   * @param phaseMassDens
   * @param dPhaseMassDens
   * @param phaseCapPressure
   * @param dPhaseCapPressure_dPhaseVolFrac
   * @param capPressureFlag
   * @param k_up
   * @param k_up_o
   * @param phaseFlux
   * @param dPhaseFlux_dP
   * @param dPhaseFlux_dC
   */
  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  static void computePotentialFluxesGravity( localIndex const numPhase,
                                             localIndex const ip,
                                             localIndex const (&seri)[numFluxSupportPoints],
                                             localIndex const (&sesri)[numFluxSupportPoints],
                                             localIndex const (&sei)[numFluxSupportPoints],
                                             real64 const (&transmissibility)[2],
                                             real64 const (&dTrans_dPres)[2],
                                             real64 const totFlux,
                                             ElementViewConst< arrayView1d< real64 const > > const & pres,
                                             ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                             ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                             ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                             ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                             ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                             ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                             ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                             ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                             ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                             localIndex const capPressureFlag,
                                             localIndex( &k_up),
                                             localIndex (&k_up_o),
                                             real64 & phaseFlux,
                                             real64 (& dPhaseFlux_dP)[numFluxSupportPoints],
                                             real64 ( & dPhaseFlux_dC)[numFluxSupportPoints][numComp] )
  {

    real64 fflow{};
    real64 dFflow_dP[numFluxSupportPoints]{};
    real64 dFflow_dC[numFluxSupportPoints][numComp]{};

    real64 pot{};
    real64 dPot_dP[numFluxSupportPoints]{};
    real64 dPot_dC[numFluxSupportPoints][numComp]{};
    real64 dProp_dC[numComp]{};

    //
    UpwindHelpers::computePotentialGravity< numComp, numFluxSupportPoints >::compute( numPhase,
                                                                                      ip,
                                                                                      seri,
                                                                                      sesri,
                                                                                      sei,
                                                                                      transmissibility,
                                                                                      dTrans_dPres,
                                                                                      totFlux,
                                                                                      gravCoef,
                                                                                      dCompFrac_dCompDens,
                                                                                      phaseMassDens,
                                                                                      dPhaseMassDens,
                                                                                      dPhaseVolFrac,
                                                                                      phaseCapPressure,
                                                                                      dPhaseCapPressure_dPhaseVolFrac,
                                                                                      pot,
                                                                                      dPot_dP,
                                                                                      dPot_dC,
                                                                                      dProp_dC );

    // and the fractional flow for gravitational part as \lambda_i^{up}/\sum_{numPhase}(\lambda_k^{up}) with up decided upon
    // the Upwind strategy
    UpwindHelpers::computeFractionalFlowGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                          ip,
                                                                                          seri,
                                                                                          sesri,
                                                                                          sei,
                                                                                          transmissibility,
                                                                                          dTrans_dPres,
                                                                                          totFlux,
                                                                                          pres,
                                                                                          gravCoef,
                                                                                          dCompFrac_dCompDens,
                                                                                          phaseMassDens,
                                                                                          dPhaseMassDens,
                                                                                          phaseMob,
                                                                                          dPhaseMob,
                                                                                          dPhaseVolFrac,
                                                                                          phaseCapPressure,
                                                                                          dPhaseCapPressure_dPhaseVolFrac,
                                                                                          capPressureFlag,
                                                                                          k_up,
                                                                                          fflow,
                                                                                          dFflow_dP,
                                                                                          dFflow_dC );


    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( ip != jp )
      {

        real64 potOther{};
        real64 dPotOther_dP[numFluxSupportPoints]{};
        real64 dPotOther_dC[numFluxSupportPoints][numComp]{};
        real64 dPropOther_dC[numComp]{};

        //Fetch pot for phase j!=i defined as \rho_j g dz/dx
        UpwindHelpers::computePotentialGravity< numComp, numFluxSupportPoints >::compute( numPhase,
                                                                                          jp,
                                                                                          seri,
                                                                                          sesri,
                                                                                          sei,
                                                                                          transmissibility,
                                                                                          dTrans_dPres,
                                                                                          totFlux,
                                                                                          gravCoef,
                                                                                          dCompFrac_dCompDens,
                                                                                          phaseMassDens,
                                                                                          dPhaseMassDens,
                                                                                          dPhaseVolFrac,
                                                                                          phaseCapPressure,
                                                                                          dPhaseCapPressure_dPhaseVolFrac,
                                                                                          potOther,
                                                                                          dPotOther_dP,
                                                                                          dPotOther_dC,
                                                                                          dPropOther_dC );

        //Eventually get the mobility of the second phase
        real64 mobOther{};
        real64 dMobOther_dP{};
        real64 dMobOther_dC[numComp]{};

        // and the other mobility for gravitational part as \lambda_j^{up} with up decided upon
        // the Upwind strategy - Note that it should be the same as the gravitational fractional flow

        UpwindHelpers::upwindMobilityGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                       jp,
                                                                                       seri,
                                                                                       sesri,
                                                                                       sei,
                                                                                       transmissibility,
                                                                                       dTrans_dPres,
                                                                                       totFlux,
                                                                                       pres,
                                                                                       gravCoef,
                                                                                       dCompFrac_dCompDens,
                                                                                       phaseMassDens,
                                                                                       dPhaseMassDens,
                                                                                       phaseMob,
                                                                                       dPhaseMob,
                                                                                       dPhaseVolFrac,
                                                                                       phaseCapPressure,
                                                                                       dPhaseCapPressure_dPhaseVolFrac,
                                                                                       capPressureFlag,
                                                                                       k_up_o,
                                                                                       mobOther,
                                                                                       dMobOther_dP,
                                                                                       dMobOther_dC );


        // Assembling gravitational flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (G_i - G_k)
        phaseFlux -= fflow * mobOther * (pot - potOther);
        dPhaseFlux_dP[k_up_o] -= fflow * dMobOther_dP * (pot - potOther);
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[k_up_o][jc] -= fflow * dMobOther_dC[jc] * (pot - potOther);
        }

        //mob related part of dFflow_dP is only upstream defined but totMob related is defined everywhere
        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dPhaseFlux_dP[ke] -= dFflow_dP[ke] * mobOther * (pot - potOther);

          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dPhaseFlux_dC[ke][jc] -= dFflow_dC[ke][jc] * mobOther * (pot - potOther);
          }
        }

        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dPhaseFlux_dP[ke] -= fflow * mobOther * (dPot_dP[ke] - dPotOther_dP[ke]);
          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dPhaseFlux_dC[ke][jc] -= fflow * mobOther * (dPot_dC[ke][jc] - dPotOther_dC[ke][jc]);
          }
        }
      }
    }

  }
  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  static void computePotentialFluxesCapillary( localIndex const numPhase,
                                               localIndex const ip,
                                               localIndex const (&seri)[numFluxSupportPoints],
                                               localIndex const (&sesri)[numFluxSupportPoints],
                                               localIndex const (&sei)[numFluxSupportPoints],
                                               real64 const (&transmissibility)[2],
                                               real64 const (&dTrans_dPres)[2],
                                               real64 const totFlux,
                                               ElementViewConst< arrayView1d< real64 const > > const & pres,
                                               ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                               ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                               ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                               ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                               ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                               ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                               ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                               ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                               ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                               localIndex const capPressureFlag,
                                               localIndex( &k_up),
                                               localIndex (&k_up_o),
                                               real64 & phaseFlux,
                                               real64 (& dPhaseFlux_dP)[numFluxSupportPoints],
                                               real64 ( & dPhaseFlux_dC)[numFluxSupportPoints][numComp] )
  {

    real64 fflow{};
    real64 dFflow_dP[numFluxSupportPoints]{};
    real64 dFflow_dC[numFluxSupportPoints][numComp]{};

    real64 pot{};
    real64 dPot_dP[numFluxSupportPoints]{};
    real64 dPot_dC[numFluxSupportPoints][numComp]{};
    real64 dProp_dC[numComp]{};

    //
    UpwindHelpers::computePotentialCapillary< numComp, numFluxSupportPoints >::compute( numPhase,
                                                                                        ip,
                                                                                        seri,
                                                                                        sesri,
                                                                                        sei,
                                                                                        transmissibility,
                                                                                        dTrans_dPres,
                                                                                        totFlux,
                                                                                        gravCoef,
                                                                                        dCompFrac_dCompDens,
                                                                                        phaseMassDens,
                                                                                        dPhaseMassDens,
                                                                                        dPhaseVolFrac,
                                                                                        phaseCapPressure,
                                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                                        pot,
                                                                                        dPot_dP,
                                                                                        dPot_dC,
                                                                                        dProp_dC );

    // and the fractional flow for gravitational part as \lambda_i^{up}/\sum_{numPhase}(\lambda_k^{up}) with up decided upon
    // the Upwind strategy
    UpwindHelpers::computeFractionalFlowCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                            ip,
                                                                                            seri,
                                                                                            sesri,
                                                                                            sei,
                                                                                            transmissibility,
                                                                                            dTrans_dPres,
                                                                                            totFlux,
                                                                                            pres,
                                                                                            gravCoef,
                                                                                            dCompFrac_dCompDens,
                                                                                            phaseMassDens,
                                                                                            dPhaseMassDens,
                                                                                            phaseMob,
                                                                                            dPhaseMob,
                                                                                            dPhaseVolFrac,
                                                                                            phaseCapPressure,
                                                                                            dPhaseCapPressure_dPhaseVolFrac,
                                                                                            capPressureFlag,
                                                                                            k_up,
                                                                                            fflow,
                                                                                            dFflow_dP,
                                                                                            dFflow_dC );


    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( ip != jp )
      {

        real64 potOther{};
        real64 dPotOther_dP[numFluxSupportPoints]{};
        real64 dPotOther_dC[numFluxSupportPoints][numComp]{};
        real64 dPropOther_dC[numComp]{};

        //Fetch pot for phase j!=i defined as \rho_j g dz/dx
        UpwindHelpers::computePotentialCapillary< numComp, numFluxSupportPoints >::compute( numPhase,
                                                                                            jp,
                                                                                            seri,
                                                                                            sesri,
                                                                                            sei,
                                                                                            transmissibility,
                                                                                            dTrans_dPres,
                                                                                            totFlux,
                                                                                            gravCoef,
                                                                                            dCompFrac_dCompDens,
                                                                                            phaseMassDens,
                                                                                            dPhaseMassDens,
                                                                                            dPhaseVolFrac,
                                                                                            phaseCapPressure,
                                                                                            dPhaseCapPressure_dPhaseVolFrac,
                                                                                            potOther,
                                                                                            dPotOther_dP,
                                                                                            dPotOther_dC,
                                                                                            dPropOther_dC );

        //Eventually get the mobility of the second phase
        real64 mobOther{};
        real64 dMobOther_dP{};
        real64 dMobOther_dC[numComp]{};

        // and the other mobility for gravitational part as \lambda_j^{up} with up decided upon
        // the Upwind strategy - Note that it should be the same as the gravitational fractional flow

        UpwindHelpers::upwindMobilityCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                         jp,
                                                                                         seri,
                                                                                         sesri,
                                                                                         sei,
                                                                                         transmissibility,
                                                                                         dTrans_dPres,
                                                                                         totFlux,
                                                                                         pres,
                                                                                         gravCoef,
                                                                                         dCompFrac_dCompDens,
                                                                                         phaseMassDens,
                                                                                         dPhaseMassDens,
                                                                                         phaseMob,
                                                                                         dPhaseMob,
                                                                                         dPhaseVolFrac,
                                                                                         phaseCapPressure,
                                                                                         dPhaseCapPressure_dPhaseVolFrac,
                                                                                         capPressureFlag,
                                                                                         k_up_o,
                                                                                         mobOther,
                                                                                         dMobOther_dP,
                                                                                         dMobOther_dC );


        // Assembling gravitational flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (G_i - G_k)
        phaseFlux -= fflow * mobOther * (pot - potOther);
        dPhaseFlux_dP[k_up_o] -= fflow * dMobOther_dP * (pot - potOther);
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[k_up_o][jc] -= fflow * dMobOther_dC[jc] * (pot - potOther);
        }

        //mob related part of dFflow_dP is only upstream defined but totMob related is defined everywhere
        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dPhaseFlux_dP[ke] -= dFflow_dP[ke] * mobOther * (pot - potOther);

          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dPhaseFlux_dC[ke][jc] -= dFflow_dC[ke][jc] * mobOther * (pot - potOther);
          }
        }

        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dPhaseFlux_dP[ke] -= fflow * mobOther * (dPot_dP[ke] - dPotOther_dP[ke]);
          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dPhaseFlux_dC[ke][jc] -= fflow * mobOther * (dPot_dC[ke][jc] - dPotOther_dC[ke][jc]);
          }
        }
      }
    }

  }


};        //end of struct UpwindHelpers

/************************* UPWIND ******************/

/**
 * @brief Template base class for different upwind Scheme
 * @tparam T physics concerned by the scheme if specialized
 */
class UpwindScheme
{

public:

  //default ctor
  UpwindScheme() = default;

  //usual copy ctor
  UpwindScheme( UpwindScheme const & scheme ) = default;

  //default move ctor
  UpwindScheme( UpwindScheme && ) = default;

  //deleted copy and move assignement
  UpwindScheme & operator=( UpwindScheme const & ) = delete;

  UpwindScheme & operator=( UpwindScheme && ) = delete;

  virtual ~UpwindScheme() = default;

  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  void getUpwindDirectionViscous( localIndex const numPhase,
                                  localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  real64 const (&transmissibility)[2],
                                  real64 const (&dTrans_dPres)[2],
                                  real64 const totFlux,  //in fine should be a ElemnetViewConst once seq form are in place
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                  integer const capPressureFlag,
                                  localIndex & upwindDir
                                  )
  {
    real64 pot{};

    /// each derived concrete class has to define a computePotential method that is calling UpwindScheme::potential method with a specific
    /// lamda defining how to get these potentials
    UPWIND::template computePotentialViscous< numComp, numFluxSupportPoints >( numPhase,
                                                                               ip,
                                                                               seri,
                                                                               sesri,
                                                                               sei,
                                                                               transmissibility,
                                                                               dTrans_dPres,
                                                                               totFlux,
                                                                               pres,
                                                                               gravCoef,
                                                                               phaseMob,
                                                                               dCompFrac_dCompDens,
                                                                               phaseMassDens,
                                                                               dPhaseMassDens,
                                                                               dPhaseVolFrac,
                                                                               phaseCapPressure,
                                                                               dPhaseCapPressure_dPhaseVolFrac,
                                                                               capPressureFlag,
                                                                               pot );

    //all definition has been changed to fit pot>0 => first cell is upstream
    upwindDir = (pot > 0) ? 0 : 1;
  }


  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  void getUpwindDirectionGravity( localIndex const numPhase,
                                  localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  real64 const (&transmissibility)[2],
                                  real64 const (&dTrans_dPres)[2],
                                  real64 const totFlux,           //in fine should be a ElemnetViewConst once seq form are in place
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                  integer const capPressureFlag,
                                  localIndex & upwindDir
                                  )
  {
    real64 pot{};

    /// each derived concrete class has to define a computePotential method that is calling UpwindScheme::potential method with a specific
    /// lamda defining how to get these potentials
    UPWIND::template computePotentialGravity< numComp, numFluxSupportPoints >( numPhase,
                                                                               ip,
                                                                               seri,
                                                                               sesri,
                                                                               sei,
                                                                               transmissibility,
                                                                               dTrans_dPres,
                                                                               totFlux,
                                                                               pres,
                                                                               gravCoef,
                                                                               phaseMob,
                                                                               dCompFrac_dCompDens,
                                                                               phaseMassDens,
                                                                               dPhaseMassDens,
                                                                               dPhaseVolFrac,
                                                                               phaseCapPressure,
                                                                               dPhaseCapPressure_dPhaseVolFrac,
                                                                               capPressureFlag,
                                                                               pot );

    //all definition has been changed to fit pot>0 => first cell is upstream
    upwindDir = (pot > 0) ? 0 : 1;
  }


  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOSX_HOST_DEVICE
  void getUpwindDirectionCapillary( localIndex const numPhase,
                                    localIndex const ip,
                                    localIndex const (&seri)[numFluxSupportPoints],
                                    localIndex const (&sesri)[numFluxSupportPoints],
                                    localIndex const (&sei)[numFluxSupportPoints],
                                    real64 const (&transmissibility)[2],
                                    real64 const (&dTrans_dPres)[2],
                                    real64 const totFlux,         //in fine should be a ElemnetViewConst once seq form are in place
                                    ElementViewConst< arrayView1d< real64 const > > const & pres,
                                    ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                    ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                    ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                    ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                    ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                    ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                    ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                    ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                    integer const capPressureFlag,
                                    localIndex & upwindDir
                                    )
  {
    real64 pot{};

    /// each derived concrete class has to define a computePotential method that is calling UpwindScheme::potential method with a specific
    /// lamda defining how to get these potentials
    UPWIND::template computePotentialCapillary< numComp, numFluxSupportPoints >( numPhase,
                                                                                 ip,
                                                                                 seri,
                                                                                 sesri,
                                                                                 sei,
                                                                                 transmissibility,
                                                                                 dTrans_dPres,
                                                                                 totFlux,
                                                                                 pres,
                                                                                 gravCoef,
                                                                                 phaseMob,
                                                                                 dCompFrac_dCompDens,
                                                                                 phaseMassDens,
                                                                                 dPhaseMassDens,
                                                                                 dPhaseVolFrac,
                                                                                 phaseCapPressure,
                                                                                 dPhaseCapPressure_dPhaseVolFrac,
                                                                                 capPressureFlag,
                                                                                 pot );

    //all definition has been changed to fit pot>0 => first cell is upstream
    upwindDir = (pot > 0) ? 0 : 1;
  }



  /// templated way of evaluating the potential (to the exception of viscous one) which relies on
  /// up-or-downwinded mobility terms pre-multiplying potential differences
  template< localIndex numComp, localIndex numFluxSupportPoints, typename LAMBDA >
  GEOSX_HOST_DEVICE
  static void potential( localIndex numPhase,
                         localIndex ip,
                         localIndex const (&seri)[numFluxSupportPoints],
                         localIndex const (&sesri)[numFluxSupportPoints],
                         localIndex const (&sei)[numFluxSupportPoints],
                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                         real64 & weightedPotential,
                         LAMBDA && fn )
  {
    //getPhase Pot
    real64 pot{};
    real64 pot_dP[numFluxSupportPoints]{};
    real64 pot_dC[numFluxSupportPoints][numComp]{};
    real64 dProp_dC[numComp]{};

    fn( ip, pot, pot_dP, pot_dC, dProp_dC );

    localIndex const k_up = 0;
    localIndex const k_dw = 1;

    //loop other other phases to form
    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( jp != ip )
      {
        localIndex const er_up = seri[k_up];
        localIndex const esr_up = sesri[k_up];
        localIndex const ei_up = sei[k_up];

        localIndex const er_dw = seri[k_dw];
        localIndex const esr_dw = sesri[k_dw];
        localIndex const ei_dw = sei[k_dw];

        real64 potOther{};
        real64 potOther_dP[numFluxSupportPoints]{};
        real64 potOther_dC[numFluxSupportPoints][numComp]{};
        real64 dPropOther_dC[numComp]{};

        fn( jp, potOther, potOther_dP, potOther_dC, dPropOther_dC );

        real64 const mob_up = phaseMob[er_up][esr_up][ei_up][jp];
        real64 const mob_dw = phaseMob[er_dw][esr_dw][ei_dw][jp];

        weightedPotential += (pot - potOther >= 0) ? mob_dw * (potOther - pot) : mob_up * (potOther - pot);

      }
    }
  }

};

/**
 * @brief  Class describing the Hybrid Upwind scheme as defined in "Consistent upwinding for sequential fully implicit
 *         multiscale compositional simulation" (Moncorge,2020)
 */
class HybridUpwind : public UpwindScheme
{

public:
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOSX_HOST_DEVICE
  static
  void computePotentialViscous( localIndex const numPhase,
                                localIndex const ip,
                                localIndex const (&seri)[numFluxSupportPoints],
                                localIndex const (&sesri)[numFluxSupportPoints],
                                localIndex const (&sei)[numFluxSupportPoints],
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const totalFlux,
                                ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const &
                                GEOSX_UNUSED_PARAM( phaseMob ),
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                integer const GEOSX_UNUSED_PARAM( capPressureFlag ),
                                real64 & potential
                                )
  {
    real64 dPot_dP[numFluxSupportPoints]{};
    real64 dPot_dC[numFluxSupportPoints][numComp]{};
    real64 dProp_dC[numComp]{};


    UpwindHelpers::computePotentialViscous< numComp, numFluxSupportPoints >::compute(
      numPhase,
      ip,
      seri,
      sesri,
      sei,
      transmissibility,
      dTrans_dPres,
      totalFlux,
      gravCoef,
      dCompFrac_dCompDens,
      phaseMassDens,
      dPhaseMassDens,
      dPhaseVolFrac,
      phaseCapPressure,
      dPhaseCapPressure_dPhaseVolFrac,
      potential,
      dPot_dP,
      dPot_dC,
      dProp_dC );
  }

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOSX_HOST_DEVICE
  static
  void computePotentialGravity( localIndex const numPhase,
                                localIndex const ip,
                                localIndex const (&seri)[numFluxSupportPoints],
                                localIndex const (&sesri)[numFluxSupportPoints],
                                localIndex const (&sei)[numFluxSupportPoints],
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const totalFlux,
                                ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                integer const GEOSX_UNUSED_PARAM( capPressureFlag ),
                                real64 & potential
                                )
  {

    //Form total velocity
    potential = 0;

    //the arg lambda allows us to access some genericity
    UpwindScheme::template potential< numComp, numFluxSupportPoints >( numPhase, ip, seri, sesri, sei,
                                                                       phaseMob, potential,
                                                                       [&]( localIndex ipp,
                                                                            real64 & potential_,
                                                                            real64 (& dPotential_dP_)[numFluxSupportPoints],
                                                                            real64 (& dPotential_dC_)[numFluxSupportPoints][numComp],
                                                                            real64 (& dProp_dC)[numComp] ) {

      UpwindHelpers::computePotentialGravity< numComp, numFluxSupportPoints >::compute(
        numPhase,
        ipp,
        seri,
        sesri,
        sei,
        transmissibility,
        dTrans_dPres,
        totalFlux,
        gravCoef,
        dCompFrac_dCompDens,
        phaseMassDens,
        dPhaseMassDens,
        dPhaseVolFrac,
        phaseCapPressure,
        dPhaseCapPressure_dPhaseVolFrac,
        potential_,
        dPotential_dP_,
        dPotential_dC_,
        dProp_dC );

    } );
  }


  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOSX_HOST_DEVICE
  static
  void computePotentialCapillary( localIndex const numPhase,
                                  localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  real64 const (&transmissibility)[2],
                                  real64 const (&dTrans_dPres)[2],
                                  real64 const totalFlux,
                                  ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const &
                                  phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                  integer const GEOSX_UNUSED_PARAM( capPressureFlag ),
                                  real64 & potential
                                  )
  {

    //Form total velocity
    potential = 0;

    //the arg lambda allows us to access some genericity
    UpwindScheme::template potential< numComp, numFluxSupportPoints >( numPhase, ip, seri, sesri, sei,
                                                                       phaseMob, potential,
                                                                       [&]( localIndex ipp,
                                                                            real64 & potential_,
                                                                            real64 (& dPotential_dP_)[numFluxSupportPoints],
                                                                            real64 (& dPotential_dC_)[numFluxSupportPoints][numComp],
                                                                            real64 (& dProp_dC)[numComp] ) {

      UpwindHelpers::computePotentialCapillary< numComp, numFluxSupportPoints >::compute(
        numPhase,
        ipp,
        seri,
        sesri,
        sei,
        transmissibility,
        dTrans_dPres,
        totalFlux,
        gravCoef,
        dCompFrac_dCompDens,
        phaseMassDens,
        dPhaseMassDens,
        dPhaseVolFrac,
        phaseCapPressure,
        dPhaseCapPressure_dPhaseVolFrac,
        potential_,
        dPotential_dP_,
        dPotential_dC_,
        dProp_dC );


    } );
  }

};
}    //namespace CompositionalMultiPhaseUpwindUtilities
}//namespace geosx


#endif //GEOSX_COMPOSITIONALMULTIPHASEFVMUPWINDUTILITIES_HPP
