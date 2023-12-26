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
 * @file IsothermalCompositionalMultiphaseFVMKernelUtilities.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEFVMKERNELUTILITIES_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEFVMKERNELUTILITIES_HPP_

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernelUtilities
{

// TODO make input parameter
static constexpr real64 epsC1PPU = 5000;

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

using Deriv = constitutive::multifluid::DerivativeOffset;

struct PotGrad
{
  template< integer numComp, integer numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute ( integer const numPhase,
            integer const ip,
            integer const hasCapPressure,
            localIndex const ( &seri )[numFluxSupportPoints],
            localIndex const ( &sesri )[numFluxSupportPoints],
            localIndex const ( &sei )[numFluxSupportPoints],
            real64 const ( &trans )[numFluxSupportPoints],
            real64 const ( &dTrans_dPres )[numFluxSupportPoints],
            ElementViewConst< arrayView1d< real64 const > > const & pres,
            ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
            ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
            ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
            ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
            ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
            ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
            ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
            real64 & potGrad,
            real64 ( & dPresGrad_dP )[numFluxSupportPoints],
            real64 ( & dPresGrad_dC )[numFluxSupportPoints][numComp],
            real64 ( & dGravHead_dP )[numFluxSupportPoints],
            real64 ( & dGravHead_dC )[numFluxSupportPoints][numComp] )
  {
    // assign derivatives arrays to zero
    for( integer i = 0; i < numFluxSupportPoints; ++i )
    {
      dPresGrad_dP[i] = 0.0;
      dGravHead_dP[i] = 0.0;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPresGrad_dC[i][jc] = 0.0;
        dGravHead_dC[i][jc] = 0.0;
      }
    }

    // create local work arrays
    real64 densMean = 0.0;
    real64 dDensMean_dP[numFluxSupportPoints]{};
    real64 dDensMean_dC[numFluxSupportPoints][numComp]{};

    real64 presGrad = 0.0;
    real64 gravHead = 0.0;
    real64 dCapPressure_dC[numComp]{};

    real64 dProp_dC[numComp]{};

    // calculate quantities on primary connected cells
    for( integer i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      // density
      real64 const density  = phaseMassDens[er][esr][ei][0][ip];
      real64 const dDens_dP = dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];

      applyChainRule( numComp,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens[er][esr][ei][0][ip],
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

      if( hasCapPressure )
      {
        capPressure = phaseCapPressure[er][esr][ei][0][ip];

        for( integer jp = 0; jp < numPhase; ++jp )
        {
          real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dP];

          for( integer jc = 0; jc < numComp; ++jc )
          {
            dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
          }
        }
      }

      presGrad += trans[i] * (pres[er][esr][ei] - capPressure);
      dPresGrad_dP[i] += trans[i] * (1 - dCapPressure_dP)
                         + dTrans_dPres[i] * (pres[er][esr][ei] - capPressure);
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPresGrad_dC[i][jc] += -trans[i] * dCapPressure_dC[jc];
      }

      real64 const gravD     = trans[i] * gravCoef[er][esr][ei];
      real64 const dGravD_dP = dTrans_dPres[i] * gravCoef[er][esr][ei];

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

    // compute phase potential gradient
    potGrad = presGrad - gravHead;

  }

};

struct PPUPhaseFlux
{
  /**
   * @brief Form the PhasePotentialUpwind from pressure gradient and gravitational head
   * @tparam numComp number of components
   * @tparam numFluxSupportPoints number of flux support points
   * @param numPhase number of phases
   * @param ip phase index
   * @param hasCapPressure flag indicating if there is capillary pressure
   * @param seri arraySlice of the stencil-implied element region index
   * @param sesri arraySlice of the stencil-implied element subregion index
   * @param sei arraySlice of the stencil-implied element index
   * @param trans transmissibility at the connection
   * @param dTrans_dPres derivative of transmissibility wrt pressure
   * @param pres pressure
   * @param gravCoef gravitational coefficient
   * @param phaseMob phase mobility
   * @param dPhaseMob derivative of phase mobility wrt pressure, temperature, comp density
   * @param dPhaseVolFrac derivative of phase volume fraction wrt pressure, temperature, comp density
   * @param dCompFrac_dCompDens derivative of component fraction wrt component density
   * @param phaseMassDens phase mass density
   * @param dPhaseMassDens derivative of phase mass density wrt pressure, temperature, comp fraction
   * @param phaseCapPressure phase capillary pressure
   * @param dPhaseCapPressure_dPhaseVolFrac derivative of phase capillary pressure wrt phase volume fraction
   * @param k_up uptream index for this phase
   * @param potGrad potential gradient for this phase
   * @param phaseFlux phase flux
   * @param dPhaseFlux_dP derivative of phase flux wrt pressure
   * @param dPhaseFlux_dC derivative of phase flux wrt comp density
   */
  template< integer numComp, integer numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute( integer const numPhase,
           integer const ip,
           integer const hasCapPressure,
           localIndex const ( &seri )[numFluxSupportPoints],
           localIndex const ( &sesri )[numFluxSupportPoints],
           localIndex const ( &sei )[numFluxSupportPoints],
           real64 const ( &trans )[2],
           real64 const ( &dTrans_dPres )[2],
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
           localIndex & k_up,
           real64 & potGrad,
           real64 ( &phaseFlux ),
           real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
           real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp] )
  {
    real64 dPresGrad_dP[numFluxSupportPoints]{};
    real64 dPresGrad_dC[numFluxSupportPoints][numComp]{};
    real64 dGravHead_dP[numFluxSupportPoints]{};
    real64 dGravHead_dC[numFluxSupportPoints][numComp]{};
    PotGrad::compute< numComp, numFluxSupportPoints >( numPhase, ip, hasCapPressure, seri, sesri, sei, trans, dTrans_dPres, pres,
                                                       gravCoef, dPhaseVolFrac, dCompFrac_dCompDens, phaseMassDens, dPhaseMassDens,
                                                       phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac, potGrad, dPresGrad_dP,
                                                       dPresGrad_dC, dGravHead_dP, dGravHead_dC );

    // *** upwinding ***

    // choose upstream cell
    k_up = (potGrad >= 0) ? 0 : 1;

    localIndex const er_up  = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up  = sei[k_up];

    real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];

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

    real64 const dMob_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
    arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 > dPhaseMobSub =
      dPhaseMob[er_up][esr_up][ei_up][ip];

    // add contribution from upstream cell mobility derivatives
    dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
    for( integer jc = 0; jc < numComp; ++jc )
    {
      dPhaseFlux_dC[k_up][jc] += dPhaseMobSub[Deriv::dC+jc] * potGrad;
    }
  }
};

struct C1PPUPhaseFlux
{
  /**
   * @brief Form the PhasePotentialUpwind from pressure gradient and gravitational head
   * @tparam numComp number of components
   * @tparam numFluxSupportPoints number of flux support points
   * @param numPhase number of phases
   * @param ip phase index
   * @param hasCapPressure flag indicating if there is capillary pressure
   * @param seri arraySlice of the stencil-implied element region index
   * @param sesri arraySlice of the stencil-implied element subregion index
   * @param sei arraySlice of the stencil-implied element index
   * @param trans transmissibility at the connection
   * @param dTrans_dPres derivative of transmissibility wrt pressure
   * @param pres pressure
   * @param gravCoef gravitational coefficient
   * @param phaseMob phase mobility
   * @param dPhaseMob derivative of phase mobility wrt pressure, temperature, comp density
   * @param dPhaseVolFrac derivative of phase volume fraction wrt pressure, temperature, comp density
   * @param dCompFrac_dCompDens derivative of component fraction wrt component density
   * @param phaseMassDens phase mass density
   * @param dPhaseMassDens derivative of phase mass density wrt pressure, temperature, comp fraction
   * @param phaseCapPressure phase capillary pressure
   * @param dPhaseCapPressure_dPhaseVolFrac derivative of phase capillary pressure wrt phase volume fraction
   * @param k_up uptream index for this phase
   * @param potGrad potential gradient for this phase
   * @param phaseFlux phase flux
   * @param dPhaseFlux_dP derivative of phase flux wrt pressure
   * @param dPhaseFlux_dC derivative of phase flux wrt comp density
   */
  template< integer numComp, integer numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute( integer const numPhase,
           integer const ip,
           integer const hasCapPressure,
           //real64 const epsC1PPU,
           localIndex const ( &seri )[numFluxSupportPoints],
           localIndex const ( &sesri )[numFluxSupportPoints],
           localIndex const ( &sei )[numFluxSupportPoints],
           real64 const ( &trans )[2],
           real64 const ( &dTrans_dPres )[2],
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
           localIndex & k_up,
           real64 & potGrad,
           real64 ( &phaseFlux ),
           real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
           real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp] )
  {
    real64 dPresGrad_dP[numFluxSupportPoints]{};
    real64 dPresGrad_dC[numFluxSupportPoints][numComp]{};
    real64 dGravHead_dP[numFluxSupportPoints]{};
    real64 dGravHead_dC[numFluxSupportPoints][numComp]{};
    PotGrad::compute< numComp, numFluxSupportPoints >( numPhase, ip, hasCapPressure, seri, sesri, sei, trans, dTrans_dPres, pres,
                                                       gravCoef, dPhaseVolFrac, dCompFrac_dCompDens, phaseMassDens, dPhaseMassDens,
                                                       phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac, potGrad, dPresGrad_dP,
                                                       dPresGrad_dC, dGravHead_dP, dGravHead_dC );

    // gravity head
    real64 gravHead = 0.0;
    for( integer i = 0; i < numFluxSupportPoints; i++ )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      real64 const gravD = trans[i] * gravCoef[er][esr][ei];

      gravHead += gravD;
    }

    // *** upwinding ***

    // phase flux and derivatives

    // assuming TPFA in the code below

    real64 Ttrans = fabs( trans[0] );
    potGrad = potGrad / Ttrans;

    real64 const mobility_i = phaseMob[seri[0]][sesri[0]][sei[0]][ip];
    real64 const mobility_j = phaseMob[seri[1]][sesri[1]][sei[1]][ip];

    // compute phase flux, see Eqs. (66) and (69) from the reference above
    real64 smoEps = epsC1PPU;
    if( fabs( gravHead ) <= 1e-20 )
      smoEps = 1000;
    real64 const tmpSqrt = sqrt( potGrad * potGrad + smoEps * smoEps );
    real64 const smoMax = 0.5 * (-potGrad + tmpSqrt);

    phaseFlux = Ttrans * ( potGrad * mobility_i - smoMax * (mobility_j - mobility_i) );

    // derivativess

    // first part, mobility derivative

    // dP
    {
      real64 const dMob_dP = dPhaseMob[seri[0]][sesri[0]][sei[0]][ip][Deriv::dP];
      dPhaseFlux_dP[0] += Ttrans * potGrad * dMob_dP;
    }

    // dC
    {
      arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 >
      dPhaseMobSub = dPhaseMob[seri[0]][sesri[0]][sei[0]][ip];
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[0][jc] += Ttrans * potGrad * dPhaseMobSub[Deriv::dC + jc];
      }
    }

    real64 const tmpInv = 1.0 / tmpSqrt;
    real64 const dSmoMax_x = 0.5 * (1.0 - potGrad * tmpInv);

    // pressure gradient and mobility difference depend on all points in the stencil
    real64 const dMobDiff_sign[numFluxSupportPoints] = {-1.0, 1.0};
    for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      // dP

      real64 const dPotGrad_dP = dPresGrad_dP[ke] - dGravHead_dP[ke];

      // first part
      dPhaseFlux_dP[ke] += dPotGrad_dP * mobility_i;

      // second part
      real64 const dSmoMax_dP = -dPotGrad_dP * dSmoMax_x;
      dPhaseFlux_dP[ke] += -dSmoMax_dP * (mobility_j - mobility_i);

      real64 const dMob_dP = dPhaseMob[seri[ke]][sesri[ke]][sei[ke]][ip][Deriv::dP];
      dPhaseFlux_dP[ke] += -Ttrans * smoMax * dMobDiff_sign[ke] * dMob_dP;

      // dC

      arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 >
      dPhaseMobSub = dPhaseMob[seri[ke]][sesri[ke]][sei[ke]][ip];

      for( integer jc = 0; jc < numComp; ++jc )
      {
        real64 const dPotGrad_dC = dPresGrad_dC[ke][jc] - dGravHead_dC[ke][jc];

        // first part
        dPhaseFlux_dC[ke][jc] += dPotGrad_dC * mobility_i;

        // second part
        real64 const dSmoMax_dC = -dPotGrad_dC * dSmoMax_x;
        dPhaseFlux_dC[ke][jc] += -dSmoMax_dC * (mobility_j - mobility_i);
        dPhaseFlux_dC[ke][jc] += -Ttrans * smoMax * dMobDiff_sign[ke] * dPhaseMobSub[Deriv::dC + jc];
      }
    }

    potGrad = potGrad * Ttrans;

    // choose upstream cell for composition upwinding
    k_up = (phaseFlux >= 0) ? 0 : 1;

  }
};

struct PhaseComponentFlux
{
  /**
   * @brief Compute the component flux for a given phase
   * @tparam numComp number of components
   * @tparam numFluxSupportPoints number of flux support points
   * @param ip phase index
   * @param k_up uptream index for this phase
   * @param seri arraySlice of the stencil-implied element region index
   * @param sesri arraySlice of the stencil-implied element subregion index
   * @param sei arraySlice of the stencil-implied element index
   * @param phaseCompFrac phase component fraction
   * @param dPhaseCompFrac derivative of phase component fraction wrt pressure, temperature, component fraction
   * @param dCompFrac_dCompDens derivative of component fraction wrt component density
   * @param phaseFlux phase flux
   * @param dPhaseFlux_dP derivative of phase flux wrt pressure
   * @param dPhaseFlux_dC derivative of phase flux wrt comp density
   * @param compFlux component flux
   * @param dCompFlux_dP derivative of phase flux wrt pressure
   * @param dCompFlux_dC derivative of phase flux wrt comp density
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute( localIndex const ip,
           localIndex const k_up,
           localIndex const ( &seri )[numFluxSupportPoints],
           localIndex const ( &sesri )[numFluxSupportPoints],
           localIndex const ( &sei )[numFluxSupportPoints],
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           real64 const & phaseFlux,
           real64 const ( &dPhaseFlux_dP )[numFluxSupportPoints],
           real64 const ( &dPhaseFlux_dC )[numFluxSupportPoints][numComp],
           real64 ( & compFlux )[numComp],
           real64 ( & dCompFlux_dP )[numFluxSupportPoints][numComp],
           real64 ( & dCompFlux_dC )[numFluxSupportPoints][numComp][numComp] )
  {
    localIndex const er_up = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up = sei[k_up];

    real64 dProp_dC[numComp]{};

    // slice some constitutive arrays to avoid too much indexing in component loop
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE_COMP-3 > phaseCompFracSub =
      phaseCompFrac[er_up][esr_up][ei_up][0][ip];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC-3 > dPhaseCompFracSub =
      dPhaseCompFrac[er_up][esr_up][ei_up][0][ip];

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
                      dCompFrac_dCompDens[er_up][esr_up][ei_up],
                      dPhaseCompFracSub[ic],
                      dProp_dC,
                      Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dCompFlux_dC[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
      }
    }
  }

};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geosx


#endif // GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEFVMKERNELUTILITIES_HPP_
