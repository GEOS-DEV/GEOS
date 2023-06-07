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
 * @file CompositionalMultiphaseHybridFVMKernels.cpp
 */

#include "UpwindingHelper.hpp"

#include "KernelUtilities.hpp"


namespace geos
{

namespace compositionalMultiphaseHybridFVMKernels
{

/******************************** UpwindingHelper ********************************/

template< integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  upwindViscousCoefficient( localIndex const (&localIds)[ 3 ],
                            localIndex const (&neighborIds)[ 3 ],
                            ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
                            ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
                            ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                            ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                            ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                            ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                            ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                            ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                            real64 const & oneSidedVolFlux,
                            real64 ( & upwPhaseViscCoef )[ NP ][ NC ],
                            real64 ( & dUpwPhaseViscCoef_dPres )[ NP ][ NC ],
                            real64 ( & dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ],
                            globalIndex & upwViscDofNumber )
{
  using Deriv = multifluid::DerivativeOffset;

  real64 dUpwMobRatio_dCompDens[ NC ]{};
  real64 dUpwDensMobRatio_dCompDens[ NC ]{};
  real64 dPhaseDens_dC[ NC ]{};
  real64 dPhaseCompFrac_dC[ NC ]{};

  // 1) Upwind
  localIndex const er  = ( oneSidedVolFlux > 0 ) ? localIds[0] : neighborIds[0];
  localIndex const esr = ( oneSidedVolFlux > 0 ) ? localIds[1] : neighborIds[1];
  localIndex const ei  = ( oneSidedVolFlux > 0 ) ? localIds[2] : neighborIds[2];

  // 2) Compute total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
  real64 totalMob = 0;
  real64 dTotalMob_dPres = 0;
  real64 dTotalMob_dCompDens[ NC ]{};
  for( integer ip = 0; ip < NP; ++ip )
  {
    totalMob = totalMob + phaseMob[er][esr][ei][ip];
    dTotalMob_dPres = dTotalMob_dPres + dPhaseMob[er][esr][ei][ip][Deriv::dP];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dTotalMob_dCompDens[ic] = dTotalMob_dCompDens[ic] + dPhaseMob[er][esr][ei][ip][Deriv::dC+ic];
    }
  }
  real64 const totalMobInv = 1.0 / totalMob;

  for( integer ip = 0; ip < NP; ++ip )
  {
    // 3) Compute viscous mobility ratio: \frac{\lambda_{\ell}}{\lambda_T}
    real64 const upwMobRatio = phaseMob[er][esr][ei][ip] * totalMobInv;
    real64 const dUpwMobRatio_dPres = ( dPhaseMob[er][esr][ei][ip][Deriv::dP] - upwMobRatio * dTotalMob_dPres )
                                      * totalMobInv;
    for( integer ic = 0; ic < NC; ++ic )
    {
      dUpwMobRatio_dCompDens[ic] = ( dPhaseMob[er][esr][ei][ip][Deriv::dC+ic] - upwMobRatio * dTotalMob_dCompDens[ic] )
                                   * totalMobInv;
    }

    // 4) Multiply mobility ratio by phase density: \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
    applyChainRule( NC,
                    dCompFrac_dCompDens[er][esr][ei],
                    dPhaseDens[er][esr][ei][0][ip],
                    dPhaseDens_dC,
                    Deriv::dC );
    real64 const upwDensMobRatio = phaseDens[er][esr][ei][0][ip] * upwMobRatio;
    real64 const dUpwDensMobRatio_dPres = dPhaseDens[er][esr][ei][0][ip][Deriv::dP] * upwMobRatio
                                          + phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dPres;
    for( integer ic = 0; ic < NC; ++ic )
    {
      dUpwDensMobRatio_dCompDens[ic] = dPhaseDens_dC[ic] * upwMobRatio
                                       + phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dCompDens[ic];
    }

    // 5) Multiply density mobility ratio by phase comp fraction: x_{c,\ell} \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
    for( integer ic = 0; ic < NC; ++ic )
    {
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseCompFrac[er][esr][ei][0][ip][ic],
                      dPhaseCompFrac_dC,
                      Deriv::dC );
      upwPhaseViscCoef[ip][ic] = phaseCompFrac[er][esr][ei][0][ip][ic] * upwDensMobRatio;
      dUpwPhaseViscCoef_dPres[ip][ic] = dPhaseCompFrac[er][esr][ei][0][ip][ic][Deriv::dP] * upwDensMobRatio
                                        + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dPres;
      for( integer jc = 0; jc < NC; ++jc )
      {
        dUpwPhaseViscCoef_dCompDens[ip][ic][jc] = dPhaseCompFrac_dC[jc] * upwDensMobRatio
                                                  + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dCompDens[jc];
      }
    }
  }
  // 6) Save the dof number of the upwind cell
  upwViscDofNumber = elemDofNumber[er][esr][ei];

}

template< integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  upwindBuoyancyCoefficient( localIndex const (&localIds)[ 3 ],
                             localIndex const (&neighborIds)[ 3 ],
                             real64 const & transGravCoef,
                             ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
                             ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                             ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                             ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                             ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                             ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                             real64 ( & phaseGravTerm )[ NP ][ NP-1 ],
                             real64 ( & dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                             real64 ( & dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ],
                             real64 ( & upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ],
                             real64 ( & dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ],
                             real64 ( & dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] )
{
  using Deriv = multifluid::DerivativeOffset;

  // 1) Compute the driving force: T ( \rho^{avg}_{\ell} - \rho^{avg}_m ) g \Delta z
  computePhaseGravTerm( localIds,
                        neighborIds,
                        transGravCoef,
                        phaseMassDens,
                        dPhaseMassDens,
                        dCompFrac_dCompDens,
                        phaseGravTerm,
                        dPhaseGravTerm_dPres,
                        dPhaseGravTerm_dCompDens );

  // 2) Compute the total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
  real64 totalMob = 0.0;
  real64 dTotalMob_dPres[ 2 ]{};
  real64 dTotalMob_dCompDens[ 2 ][ NC ]{};

  computeUpwindedTotalMobility( localIds,
                                neighborIds,
                                phaseMob,
                                dPhaseMob,
                                phaseGravTerm,
                                totalMob,
                                dTotalMob_dPres,
                                dTotalMob_dCompDens );
  real64 const totalMobInv = 1.0 / totalMob;

  // 3) Compute the quantities \x_{up}_{c,p} \rho_p \frac{\lambda_p \lambda_m}{\lambda_T}
  real64 dMobRatio_dPres[ 2 ]{};
  real64 dMobRatio_dCompDens[ 2 ][ NC ]{};
  real64 dDensMobRatio_dPres[ 2 ]{};
  real64 dDensMobRatio_dCompDens[ 2 ][ NC ]{};

  real64 dPhaseDens_dC[ NC ]{};
  real64 dPhaseCompFrac_dC[ NC ]{};

  for( integer ip = 0; ip < NP; ++ip )
  {
    localIndex k = 0;
    for( integer jp = 0; jp < NP; ++jp )
    {
      if( ip == jp )
      {
        continue;
      }

      // 3.a) Upwinding using the gravity term
      localIndex eru, esru, eiu, posu; // upwind
      localIndex erd, esrd, eid, posd; // downwind
      setIndicesForMobilityRatioUpwinding( localIds, neighborIds,
                                           phaseGravTerm[ip][k],
                                           eru, esru, eiu, posu,
                                           erd, esrd, eid, posd );

      // 3.b) Compute mobility ratio \frac{\lambda_l \lambda_m}{\lambda_T}
      real64 const mobRatio = phaseMob[eru][esru][eiu][ip] * phaseMob[erd][esrd][eid][jp]
                              * totalMobInv;
      dMobRatio_dPres[posu] = ( dPhaseMob[eru][esru][eiu][ip][Deriv::dP] * phaseMob[erd][esrd][eid][jp]
                                - mobRatio * dTotalMob_dPres[posu] ) * totalMobInv;
      dMobRatio_dPres[posd] = ( dPhaseMob[erd][esrd][eid][jp][Deriv::dP] * phaseMob[eru][esru][eiu][ip]
                                - mobRatio * dTotalMob_dPres[posd] ) * totalMobInv;

      for( integer ic = 0; ic < NC; ++ic )
      {
        dMobRatio_dCompDens[posu][ic] = ( dPhaseMob[eru][esru][eiu][ip][Deriv::dC+ic] * phaseMob[erd][esrd][eid][jp]
                                          - mobRatio * dTotalMob_dCompDens[posu][ic] ) * totalMobInv;
        dMobRatio_dCompDens[posd][ic] = ( dPhaseMob[erd][esrd][eid][jp][Deriv::dC+ic] * phaseMob[eru][esru][eiu][ip]
                                          - mobRatio * dTotalMob_dCompDens[posd][ic] ) * totalMobInv;
      }

      // 3.c) Compute mobility ratio multiplied by upwinded phase density \rho_l \frac{\lambda_l \lambda_m}{\lambda_T}
      applyChainRule( NC,
                      dCompFrac_dCompDens[eru][esru][eiu],
                      dPhaseDens[eru][esru][eiu][0][ip],
                      dPhaseDens_dC,
                      Deriv::dC );
      real64 const densMobRatio = phaseDens[eru][esru][eiu][0][ip] * mobRatio;
      dDensMobRatio_dPres[posu] = dPhaseDens[eru][esru][eiu][0][ip][Deriv::dP] * mobRatio
                                  + phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posu];
      dDensMobRatio_dPres[posd] = phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posd];
      for( integer ic = 0; ic < NC; ++ic )
      {
        dDensMobRatio_dCompDens[posu][ic] = dPhaseDens_dC[ic] * mobRatio
                                            + phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posu][ic];
        dDensMobRatio_dCompDens[posd][ic] = phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posd][ic];
      }

      // 3.d) Compute the final gravity coefficient \x_{up}_{c,p} \rho_p \frac{\lambda_l \lambda_m}{\lambda_T}
      for( integer ic = 0; ic < NC; ++ic )
      {
        applyChainRule( NC,
                        dCompFrac_dCompDens[eru][esru][eiu],
                        dPhaseCompFrac[eru][esru][eiu][0][ip][ic],
                        dPhaseCompFrac_dC,
                        Deriv::dC );
        upwPhaseGravCoef[ip][k][ic] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * densMobRatio;
        dUpwPhaseGravCoef_dPres[ip][k][ic][posu] = dPhaseCompFrac[eru][esru][eiu][0][ip][ic][Deriv::dP] * densMobRatio
                                                   + phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posu];
        dUpwPhaseGravCoef_dPres[ip][k][ic][posd] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posd];

        for( integer jc = 0; jc < NC; ++jc )
        {
          dUpwPhaseGravCoef_dCompDens[ip][k][ic][posu][jc] = dPhaseCompFrac_dC[jc] * densMobRatio
                                                             + phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dCompDens[posu][jc];
          dUpwPhaseGravCoef_dCompDens[ip][k][ic][posd][jc] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dCompDens[posd][jc];
        }
      }
      ++k;
    }
  }

}

template< integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  computePhaseGravTerm( localIndex const (&localIds)[ 3 ],
                        localIndex const (&neighborIds)[ 3 ],
                        real64 const & transGravCoef,
                        ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                        ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                        real64 ( & phaseGravTerm )[ NP ][ NP-1 ],
                        real64 ( & dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                        real64 ( & dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] )
{
  using Deriv = multifluid::DerivativeOffset;

  localIndex const er   = localIds[0];
  localIndex const esr  = localIds[1];
  localIndex const ei   = localIds[2];
  localIndex const ern  = neighborIds[0];
  localIndex const esrn = neighborIds[1];
  localIndex const ein  = neighborIds[2];

  real64 dPhaseMassDens_dCLoc[ NC ]{};
  real64 dPhaseMassDens_dCNeighbor[ NC ]{};
  real64 dPhaseMassDens_dC[ NC ]{};

  for( integer ip = 0; ip < NP; ++ip )
  {
    applyChainRule( NC,
                    dCompFrac_dCompDens[er][esr][ei],
                    dPhaseMassDens[er][esr][ei][0][ip],
                    dPhaseMassDens_dCLoc,
                    Deriv::dC );
    applyChainRule( NC,
                    dCompFrac_dCompDens[ern][esrn][ein],
                    dPhaseMassDens[ern][esrn][ein][0][ip],
                    dPhaseMassDens_dCNeighbor,
                    Deriv::dC );

    localIndex k = 0;
    for( integer jp = 0; jp < NP; ++jp )
    {
      if( ip == jp )
      {
        continue;
      }

      phaseGravTerm[ip][k] = -( phaseMassDens[er][esr][ei][0][ip] + phaseMassDens[ern][esrn][ein][0][ip] );
      phaseGravTerm[ip][k] += ( phaseMassDens[er][esr][ei][0][jp] + phaseMassDens[ern][esrn][ein][0][jp] );
      phaseGravTerm[ip][k] *= 0.5 * transGravCoef;

      dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] = ( -dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP] + dPhaseMassDens[er][esr][ei][0][jp][Deriv::dP] );
      dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] *= 0.5 * transGravCoef;

      dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] = ( -dPhaseMassDens[ern][esrn][ein][0][ip][Deriv::dP] + dPhaseMassDens[ern][esrn][ein][0][jp][Deriv::dP] );
      dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] *= 0.5 * transGravCoef;

      for( integer ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCLoc[ic];
        dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCNeighbor[ic];
      }
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens[er][esr][ei][0][jp],
                      dPhaseMassDens_dC,
                      Deriv::dC );
      for( integer ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
      }
      applyChainRule( NC,
                      dCompFrac_dCompDens[ern][esrn][ein],
                      dPhaseMassDens[ern][esrn][ein][0][jp],
                      dPhaseMassDens_dC,
                      Deriv::dC );
      for( integer ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
      }
      ++k;
    }
  }
}

template< integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  computeUpwindedTotalMobility( localIndex const (&localIds)[ 3 ],
                                localIndex const (&neighborIds)[ 3 ],
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                                real64 & totalMob,
                                real64 ( & dTotalMob_dPres )[ 2 ],
                                real64 ( & dTotalMob_dCompDens )[ 2 ][ NC ] )
{
  using Deriv = multifluid::DerivativeOffset;

  localIndex totalMobIds[ NP ][ 3 ]{};
  localIndex totalMobPos[ NP ]{};
  setIndicesForTotalMobilityUpwinding< NP >( localIds,
                                             neighborIds,
                                             phaseGravTerm,
                                             totalMobIds,
                                             totalMobPos );
  for( integer ip = 0; ip < NP; ++ip )
  {
    localIndex const er  = totalMobIds[ip][0];
    localIndex const esr = totalMobIds[ip][1];
    localIndex const ei  = totalMobIds[ip][2];
    localIndex const pos = totalMobPos[ip];
    totalMob = totalMob + phaseMob[er][esr][ei][ip];
    dTotalMob_dPres[pos] = dTotalMob_dPres[pos] + dPhaseMob[er][esr][ei][pos][Deriv::dP];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dTotalMob_dCompDens[pos][ic] = dTotalMob_dCompDens[pos][ic] + dPhaseMob[er][esr][ei][ip][Deriv::dC+ic];
    }
  }
  if( totalMob < 1e-12 )
  {
    totalMob = 1e-12;
  }
}

GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  setIndicesForMobilityRatioUpwinding( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const & gravTerm,
                                       localIndex & eru, localIndex & esru, localIndex & eiu, localIndex & posu,
                                       localIndex & erd, localIndex & esrd, localIndex & eid, localIndex & posd )
{
  if( gravTerm > 0 )
  {
    eru  = localIds[0];
    esru = localIds[1];
    eiu  = localIds[2];
    posu = Pos::LOCAL;
    erd  = neighborIds[0];
    esrd = neighborIds[1];
    eid  = neighborIds[2];
    posd = Pos::NEIGHBOR;
  }
  else
  {
    eru  = neighborIds[0];
    esru = neighborIds[1];
    eiu  = neighborIds[2];
    posu = Pos::NEIGHBOR;
    erd  = localIds[0];
    esrd = localIds[1];
    eid  = localIds[2];
    posd = Pos::LOCAL;
  }
}

template< integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  setIndicesForTotalMobilityUpwinding( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const (&gravTerm)[ NP ][ NP-1 ],
                                       localIndex ( & totalMobIds )[ NP ][ 3 ],
                                       localIndex ( & totalMobPos )[ NP ] )
{
  if( NP == 2 )
  {
    if( gravTerm[0][0] > 0 )
    {
      totalMobIds[0][0] = localIds[0];
      totalMobIds[0][1] = localIds[1];
      totalMobIds[0][2] = localIds[2];
      totalMobPos[0] = Pos::LOCAL;
      totalMobIds[1][0] = neighborIds[0];
      totalMobIds[1][1] = neighborIds[1];
      totalMobIds[1][2] = neighborIds[2];
      totalMobPos[1] = Pos::NEIGHBOR;
    }
    else
    {
      totalMobIds[0][0] = neighborIds[0];
      totalMobIds[0][1] = neighborIds[1];
      totalMobIds[0][2] = neighborIds[2];
      totalMobPos[0] = Pos::NEIGHBOR;
      totalMobIds[1][0] = localIds[0];
      totalMobIds[1][1] = localIds[1];
      totalMobIds[1][2] = localIds[2];
      totalMobPos[1] = Pos::LOCAL;
    }
  }
  else if( NP == 3 )
  {
    // TODO Francois: this should be improved
    // currently this implements the algorithm proposed by SH Lee
    for( integer ip = 0; ip < NP; ++ip )
    {
      if( ( gravTerm[ip][0] >= 0 && gravTerm[ip][1] >= 0 ) || // includes the no-buoyancy case
          ( ( LvArray::math::abs( gravTerm[ip][0] ) >= LvArray::math::abs( gravTerm[ip][1] ) ) && gravTerm[ip][1] >= 0 ) ||
          ( ( LvArray::math::abs( gravTerm[ip][1] ) >= LvArray::math::abs( gravTerm[ip][0] ) ) && gravTerm[ip][0] >= 0 ) )
      {
        totalMobIds[ip][0] = localIds[0];
        totalMobIds[ip][1] = localIds[1];
        totalMobIds[ip][2] = localIds[2];
        totalMobPos[ip] = Pos::LOCAL;
      }
      else
      {
        totalMobIds[ip][0] = neighborIds[0];
        totalMobIds[ip][1] = neighborIds[1];
        totalMobIds[ip][2] = neighborIds[2];
        totalMobPos[ip] = Pos::NEIGHBOR;
      }
    }
  }
}


#define INST_UpwindingHelperNCNP( NC, NP ) \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    upwindViscousCoefficient< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                        localIndex const (&neighborIds)[ 3 ], \
                                        ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                                        ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                        ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                        ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                        ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                                        ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                        real64 const & oneSidedVolFlux, \
                                        real64 ( &upwPhaseViscCoef )[ NP ][ NC ], \
                                        real64 ( &dUpwPhaseViscCoef_dPres )[ NP ][ NC ], \
                                        real64 ( &dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ], \
                                        globalIndex & upwViscDofNumber ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    upwindBuoyancyCoefficient< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                         localIndex const (&neighborIds)[ 3 ], \
                                         real64 const & transGravCoef,  \
                                         ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                                         ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                         ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                                         real64 ( &phaseGravTerm )[ NP ][ NP-1 ], \
                                         real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ], \
                                         real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ], \
                                         real64 ( &upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ], \
                                         real64 ( &dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ], \
                                         real64 ( &dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    computePhaseGravTerm< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                    localIndex const (&neighborIds)[ 3 ], \
                                    real64 const & transGravCoef, \
                                    ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                    ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                                    ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                    real64 ( &phaseGravTerm )[ NP ][ NP-1 ], \
                                    real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ], \
                                    real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    computeUpwindedTotalMobility< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                            localIndex const (&neighborIds)[ 3 ], \
                                            ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                            ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                            real64 const (&phaseGravTerm)[ NP ][ NP-1 ], \
                                            real64 & totalMob,    \
                                            real64 ( &dTotalMob_dPres )[ 2 ], \
                                            real64 ( &dTotalMob_dCompDens )[ 2 ][ NC ] )

INST_UpwindingHelperNCNP( 1, 2 );
INST_UpwindingHelperNCNP( 2, 2 );
INST_UpwindingHelperNCNP( 3, 2 );
INST_UpwindingHelperNCNP( 4, 2 );
INST_UpwindingHelperNCNP( 5, 2 );

INST_UpwindingHelperNCNP( 1, 3 );
INST_UpwindingHelperNCNP( 2, 3 );
INST_UpwindingHelperNCNP( 3, 3 );
INST_UpwindingHelperNCNP( 4, 3 );
INST_UpwindingHelperNCNP( 5, 3 );

#undef INST_UpwindingHelperNCNP

#define INST_UpwindingHelperNP( NP ) \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    setIndicesForTotalMobilityUpwinding< NP >( localIndex const (&localIds)[ 3 ], \
                                               localIndex const (&neighborIds)[ 3 ], \
                                               real64 const (&gravTerm)[ NP ][ NP-1 ], \
                                               localIndex ( &totalMobIds )[ NP ][ 3 ], \
                                               localIndex ( &totalMobPos )[ NP ] )

INST_UpwindingHelperNP( 2 );
INST_UpwindingHelperNP( 3 );

#undef INST_UpwindingHelperNP

}

}
