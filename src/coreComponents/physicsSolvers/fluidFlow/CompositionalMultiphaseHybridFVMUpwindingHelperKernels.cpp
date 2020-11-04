/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseHybridFVMKernels.cpp
 */

#include "CompositionalMultiphaseHybridFVMKernels.hpp"

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"


namespace geosx
{

namespace CompositionalMultiphaseHybridFVMKernels
{

/******************************** UpwindingHelper ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
UpwindingHelper::UpwindViscousCoefficient( localIndex const (&localIds)[ 3 ],
                                           localIndex const (&neighborIds)[ 3 ],
                                           ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                                           ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                                           ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                                           ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                                           ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                                           ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
                                           ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                                           ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
                                           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                                           real64 const & oneSidedVolFlux,
                                           real64 (& upwPhaseViscCoef)[ NP ][ NC ],
                                           real64 (& dUpwPhaseViscCoef_dPres)[ NP ][ NC ],
                                           real64 (& dUpwPhaseViscCoef_dCompDens)[ NP ][ NC ][ NC ],
                                           globalIndex & upwViscDofNumber )
{
  real64 dUpwMobRatio_dCompDens[ NC ] = { 0.0 };
  real64 dUpwDensMobRatio_dCompDens[ NC ] = { 0.0 };
  real64 dPhaseDens_dC[ NC ] = { 0.0 };
  real64 dPhaseCompFrac_dC[ NC ] = { 0.0 };

  // 1) Upwind
  localIndex const er  = ( oneSidedVolFlux > 0 ) ? localIds[0] : neighborIds[0];
  localIndex const esr = ( oneSidedVolFlux > 0 ) ? localIds[1] : neighborIds[1];
  localIndex const ei  = ( oneSidedVolFlux > 0 ) ? localIds[2] : neighborIds[2];

  // 2) Compute total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
  real64 totalMob = 0;
  real64 dTotalMob_dPres = 0;
  real64 dTotalMob_dCompDens[ NC ] = { 0.0 };
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    totalMob = totalMob + phaseMob[er][esr][ei][ip];
    dTotalMob_dPres = dTotalMob_dPres + dPhaseMob_dPres[er][esr][ei][ip];
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dTotalMob_dCompDens[ic] = dTotalMob_dCompDens[ic] + dPhaseMob_dCompDens[er][esr][ei][ip][ic];
    }
  }
  real64 const totalMobInv = 1.0 / totalMob;

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    // 3) Compute viscous mobility ratio: \frac{\lambda_{\ell}}{\lambda_T}
    real64 const upwMobRatio = phaseMob[er][esr][ei][ip] * totalMobInv;
    real64 const dUpwMobRatio_dPres = ( dPhaseMob_dPres[er][esr][ei][ip] - upwMobRatio * dTotalMob_dPres )
                                      * totalMobInv;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dUpwMobRatio_dCompDens[ic] = ( dPhaseMob_dCompDens[er][esr][ei][ip][ic] - upwMobRatio * dTotalMob_dCompDens[ic] )
                                   * totalMobInv;
    }

    // 4) Multiply mobility ratio by phase density: \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
    applyChainRule( NC,
                    dCompFrac_dCompDens[er][esr][ei],
                    dPhaseDens_dCompFrac[er][esr][ei][0][ip],
                    dPhaseDens_dC );
    real64 const upwDensMobRatio = phaseDens[er][esr][ei][0][ip] * upwMobRatio;
    real64 const dUpwDensMobRatio_dPres = dPhaseDens_dPres[er][esr][ei][0][ip] * upwMobRatio
                                          + phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dPres;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dUpwDensMobRatio_dCompDens[ic] = dPhaseDens_dC[ic] * upwMobRatio
                                       + phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dCompDens[ic];
    }

    // 5) Multiply density mobility ratio by phase comp fraction: x_{c,\ell} \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseCompFrac_dCompFrac[er][esr][ei][0][ip][ic],
                      dPhaseCompFrac_dC );
      upwPhaseViscCoef[ip][ic] = phaseCompFrac[er][esr][ei][0][ip][ic] * upwDensMobRatio;
      dUpwPhaseViscCoef_dPres[ip][ic] = dPhaseCompFrac_dPres[er][esr][ei][0][ip][ic] * upwDensMobRatio
                                        + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dPres;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dUpwPhaseViscCoef_dCompDens[ip][ic][jc] = dPhaseCompFrac_dC[jc] * upwDensMobRatio
                                                  + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dCompDens[jc];
      }
    }
  }
  // 6) Save the dof number of the upwind cell
  upwViscDofNumber = elemDofNumber[er][esr][ei];
}


template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
UpwindingHelper::UpwindBuoyancyCoefficient( localIndex const (&localIds)[ 3 ],
                                            localIndex const (&neighborIds)[ 3 ],
                                            real64 const & transGravCoef,
                                            ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                                            ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                                            ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                                            ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                                            ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                                            ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                                            ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                                            ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
                                            ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                                            ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
                                            real64 (& phaseGravTerm)[ NP ][ NP-1 ],
                                            real64 (& dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ],
                                            real64 (& dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ],
                                            real64 (& upwPhaseGravCoef)[ NP ][ NP-1 ][ NC ],
                                            real64 (& dUpwPhaseGravCoef_dPres)[ NP ][ NP-1 ][ NC ][ 2 ],
                                            real64 (& dUpwPhaseGravCoef_dCompDens)[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] )
{
  // 1) Compute the driving force: T ( \rho^{avg}_{\ell} - \rho^{avg}_m ) g \Delta z
  ComputePhaseGravTerm( localIds,
                        neighborIds,
                        transGravCoef,
                        phaseDens,
                        dPhaseDens_dPres,
                        dPhaseDens_dCompFrac,
                        dCompFrac_dCompDens,
                        phaseGravTerm,
                        dPhaseGravTerm_dPres,
                        dPhaseGravTerm_dCompDens );

  // 2) Compute the total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
  real64 totalMob = 0.0;
  real64 dTotalMob_dPres[ 2 ] = { 0.0 };
  real64 dTotalMob_dCompDens[ 2 ][ NC ]  = {{ 0.0 }};

  ComputeUpwindedTotalMobility( localIds,
                                neighborIds,
                                phaseMob,
                                dPhaseMob_dPres,
                                dPhaseMob_dCompDens,
                                phaseGravTerm,
                                totalMob,
                                dTotalMob_dPres,
                                dTotalMob_dCompDens );
  real64 const totalMobInv = 1.0 / totalMob;

  // 3) Compute the quantities \x_{up}_{c,p} \rho_p \frac{\lambda_p \lambda_m}{\lambda_T}
  real64 dMobRatio_dPres[ 2 ] = { 0.0 };
  real64 dMobRatio_dCompDens[ 2 ][ NC ] = {{ 0.0 }};
  real64 dDensMobRatio_dPres[ 2 ] = { 0.0 };
  real64 dDensMobRatio_dCompDens[ 2 ][ NC ] = {{ 0.0 }};

  real64 dPhaseDens_dC[ NC ] = { 0.0 };
  real64 dPhaseCompFrac_dC[ NC ] = { 0.0 };

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    localIndex k = 0;
    for( localIndex jp = 0; jp < NP; ++jp )
    {
      if( ip == jp )
      {
        continue;
      }

      // 3.a) Upwinding using the gravity term
      localIndex eru, esru, eiu, posu; // upwind
      localIndex erd, esrd, eid, posd; // downwind
      SetIndicesForMobilityRatioUpwinding( localIds, neighborIds,
                                           phaseGravTerm[ip][k],
                                           eru, esru, eiu, posu,
                                           erd, esrd, eid, posd );

      // 3.b) Compute mobility ratio \frac{\lambda_l \lambda_m}{\lambda_T}
      real64 const mobRatio = phaseMob[eru][esru][eiu][ip] * phaseMob[erd][esrd][eid][jp]
                              * totalMobInv;
      dMobRatio_dPres[posu] = ( dPhaseMob_dPres[eru][esru][eiu][ip] * phaseMob[erd][esrd][eid][jp]
                                - mobRatio * dTotalMob_dPres[posu] ) * totalMobInv;
      dMobRatio_dPres[posd] = ( dPhaseMob_dPres[erd][esrd][eid][jp] * phaseMob[eru][esru][eiu][ip]
                                - mobRatio * dTotalMob_dPres[posd] ) * totalMobInv;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dMobRatio_dCompDens[posu][ic] = ( dPhaseMob_dCompDens[eru][esru][eiu][ip][ic] * phaseMob[erd][esrd][eid][jp]
                                          - mobRatio * dTotalMob_dCompDens[posu][ic] ) * totalMobInv;
        dMobRatio_dCompDens[posd][ic] = ( dPhaseMob_dCompDens[erd][esrd][eid][jp][ic] * phaseMob[eru][esru][eiu][ip]
                                          - mobRatio * dTotalMob_dCompDens[posd][ic] ) * totalMobInv;
      }

      // 3.c) Compute mobility ratio multiplied by upwinded phase density \rho_l \frac{\lambda_l \lambda_m}{\lambda_T}
      applyChainRule( NC,
                      dCompFrac_dCompDens[eru][esru][eiu],
                      dPhaseDens_dCompFrac[eru][esru][eiu][0][ip],
                      dPhaseDens_dC );
      real64 const densMobRatio = phaseDens[eru][esru][eiu][0][ip] * mobRatio;
      dDensMobRatio_dPres[posu] = dPhaseDens_dPres[eru][esru][eiu][0][ip] * mobRatio
                                  + phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posu];
      dDensMobRatio_dPres[posd] = phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posd];
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dDensMobRatio_dCompDens[posu][ic] = dPhaseDens_dC[ic] * mobRatio
                                            + phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posu][ic];
        dDensMobRatio_dCompDens[posd][ic] = phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posd][ic];
      }

      // 3.d) Compute the final gravity coefficient \x_{up}_{c,p} \rho_p \frac{\lambda_l \lambda_m}{\lambda_T}
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        applyChainRule( NC,
                        dCompFrac_dCompDens[eru][esru][eiu],
                        dPhaseCompFrac_dCompFrac[eru][esru][eiu][0][ip][ic],
                        dPhaseCompFrac_dC );
        upwPhaseGravCoef[ip][k][ic] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * densMobRatio;
        dUpwPhaseGravCoef_dPres[ip][k][ic][posu] = dPhaseCompFrac_dPres[eru][esru][eiu][0][ip][ic] * densMobRatio
                                                   + phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posu];
        dUpwPhaseGravCoef_dPres[ip][k][ic][posd] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posd];

        for( localIndex jc = 0; jc < NC; ++jc )
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

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
UpwindingHelper::ComputePhaseGravTerm( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const & transGravCoef,
                                       ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                                       ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                                       ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                                       ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                                       real64 (& phaseGravTerm)[ NP ][ NP-1 ],
                                       real64 (& dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ],
                                       real64 (& dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ] )
{
  localIndex const er   = localIds[0];
  localIndex const esr  = localIds[1];
  localIndex const ei   = localIds[2];
  localIndex const ern  = neighborIds[0];
  localIndex const esrn = neighborIds[1];
  localIndex const ein  = neighborIds[2];

  real64 dPhaseDens_dCLoc[ NC ] = { 0.0 };
  real64 dPhaseDens_dCNeighbor[ NC ] = { 0.0 };
  real64 dPhaseDens_dC[ NC ] = { 0.0 };

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    applyChainRule( NC,
                    dCompFrac_dCompDens[er][esr][ei],
                    dPhaseDens_dCompFrac[er][esr][ei][0][ip],
                    dPhaseDens_dCLoc );
    applyChainRule( NC,
                    dCompFrac_dCompDens[ern][esrn][ein],
                    dPhaseDens_dCompFrac[ern][esrn][ein][0][ip],
                    dPhaseDens_dCNeighbor );

    localIndex k = 0;
    for( localIndex jp = 0; jp < NP; ++jp )
    {
      if( ip == jp )
      {
        continue;
      }

      phaseGravTerm[ip][k] = -( phaseDens[er][esr][ei][0][ip] + phaseDens[ern][esrn][ein][0][ip] );
      phaseGravTerm[ip][k] += ( phaseDens[er][esr][ei][0][jp] + phaseDens[ern][esrn][ein][0][jp] );
      phaseGravTerm[ip][k] *= 0.5 * transGravCoef;

      dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] = ( -dPhaseDens_dPres[er][esr][ei][0][ip] + dPhaseDens_dPres[er][esr][ei][0][jp] );
      dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] *= 0.5 * transGravCoef;

      dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] = ( -dPhaseDens_dPres[ern][esrn][ein][0][ip] + dPhaseDens_dPres[ern][esrn][ein][0][jp] );
      dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] *= 0.5 * transGravCoef;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] = -0.5 * transGravCoef * dPhaseDens_dCLoc[ic];
        dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] = -0.5 * transGravCoef * dPhaseDens_dCNeighbor[ic];
      }
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseDens_dCompFrac[er][esr][ei][0][jp],
                      dPhaseDens_dC );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] += 0.5 * transGravCoef * dPhaseDens_dC[ic];
      }
      applyChainRule( NC,
                      dCompFrac_dCompDens[ern][esrn][ein],
                      dPhaseDens_dCompFrac[ern][esrn][ein][0][jp],
                      dPhaseDens_dC );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] += 0.5 * transGravCoef * dPhaseDens_dC[ic];
      }
      ++k;
    }
  }
}


template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
UpwindingHelper::ComputeUpwindedTotalMobility( localIndex const (&localIds)[ 3 ],
                                               localIndex const (&neighborIds)[ 3 ],
                                               ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                                               ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                                               ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                                               real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                                               real64 & totalMob,
                                               real64 (& dTotalMob_dPres)[ 2 ],
                                               real64 (& dTotalMob_dCompDens)[ 2 ][ NC ] )
{

  localIndex totalMobIds[ NP ][ 3 ] = {{ 0 }};
  localIndex totalMobPos[ NP ] = { 0 };
  SetIndicesForTotalMobilityUpwinding< NP >( localIds,
                                             neighborIds,
                                             phaseGravTerm,
                                             totalMobIds,
                                             totalMobPos );
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    localIndex const er  = totalMobIds[ip][0];
    localIndex const esr = totalMobIds[ip][1];
    localIndex const ei  = totalMobIds[ip][2];
    localIndex const pos = totalMobPos[ip];
    totalMob = totalMob + phaseMob[er][esr][ei][ip];
    dTotalMob_dPres[pos] = dTotalMob_dPres[pos] + dPhaseMob_dPres[er][esr][ei][pos];
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dTotalMob_dCompDens[pos][ic] = dTotalMob_dCompDens[pos][ic] + dPhaseMob_dCompDens[er][esr][ei][ip][ic];
    }
  }
  if( totalMob < 1e-12 )
  {
    totalMob = 1e-12;
  }
}


GEOSX_HOST_DEVICE
void
UpwindingHelper::SetIndicesForMobilityRatioUpwinding( localIndex const (&localIds)[ 3 ],
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

template<>
GEOSX_HOST_DEVICE
void
UpwindingHelper::SetIndicesForTotalMobilityUpwinding< 2 >( localIndex const (&localIds)[ 3 ],
                                                           localIndex const (&neighborIds)[ 3 ],
                                                           real64 const (&gravTerm)[ 2 ][ 1 ],
                                                           localIndex (& totalMobIds)[ 2 ][ 3 ],
                                                           localIndex (& totalMobPos)[ 2 ] )
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

template<>
GEOSX_HOST_DEVICE
void
UpwindingHelper::SetIndicesForTotalMobilityUpwinding< 3 >( localIndex const (&localIds)[ 3 ],
                                                           localIndex const (&neighborIds)[ 3 ],
                                                           real64 const (&gravTerm)[ 3 ][ 2 ],
                                                           localIndex (& totalMobIds)[ 3 ][ 3 ],
                                                           localIndex (& totalMobPos)[ 3 ] )
{
  GEOSX_UNUSED_VAR( localIds );
  GEOSX_UNUSED_VAR( neighborIds );
  GEOSX_UNUSED_VAR( gravTerm );
  GEOSX_UNUSED_VAR( totalMobIds );
  GEOSX_UNUSED_VAR( totalMobPos );
}

#define INST_UpwindingHelper( NC, NP ) \
  template \
  void \
  UpwindingHelper::UpwindViscousCoefficient< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                                       localIndex const (&neighborIds)[ 3 ], \
                                                       ElementViewConst< arrayView3d< real64 const > > const & phaseDens, \
                                                       ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres, \
                                                       ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac, \
                                                       ElementViewConst< arrayView2d< real64 const > > const & phaseMob, \
                                                       ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres, \
                                                       ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens, \
                                                       ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens, \
                                                       ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac, \
                                                       ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres, \
                                                       ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac, \
                                                       ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                                       real64 const & oneSidedVolFlux, \
                                                       real64 ( &upwPhaseViscCoef )[ NP ][ NC ], \
                                                       real64 ( &dUpwPhaseViscCoef_dPres )[ NP ][ NC ], \
                                                       real64 ( &dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ], \
                                                       globalIndex & upwViscDofNumber ); \
  template \
  void \
  UpwindingHelper::UpwindBuoyancyCoefficient< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                                        localIndex const (&neighborIds)[ 3 ], \
                                                        real64 const & transGravCoef, \
                                                        ElementViewConst< arrayView3d< real64 const > > const & phaseDens, \
                                                        ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres, \
                                                        ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac, \
                                                        ElementViewConst< arrayView2d< real64 const > > const & phaseMob, \
                                                        ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres, \
                                                        ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens, \
                                                        ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens, \
                                                        ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac, \
                                                        ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres, \
                                                        ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac, \
                                                        real64 ( &phaseGravTerm )[ NP ][ NP-1 ], \
                                                        real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ], \
                                                        real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ], \
                                                        real64 ( &upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ], \
                                                        real64 ( &dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ], \
                                                        real64 ( &dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] ); \
  template \
  void \
  UpwindingHelper::ComputePhaseGravTerm< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                                   localIndex const (&neighborIds)[ 3 ], \
                                                   real64 const & transGravCoef, \
                                                   ElementViewConst< arrayView3d< real64 const > > const & phaseDens, \
                                                   ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres, \
                                                   ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac, \
                                                   ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens, \
                                                   real64 ( &phaseGravTerm )[ NP ][ NP-1 ], \
                                                   real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ], \
                                                   real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] ); \
  template \
  void \
  UpwindingHelper::ComputeUpwindedTotalMobility< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                                           localIndex const (&neighborIds)[ 3 ], \
                                                           ElementViewConst< arrayView2d< real64 const > > const & phaseMob, \
                                                           ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres, \
                                                           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens, \
                                                           real64 const (&phaseGravTerm)[ NP ][ NP-1 ], \
                                                           real64 & totalMob, \
                                                           real64 ( &dTotalMob_dPres )[ 2 ], \
                                                           real64 ( &dTotalMob_dCompDens )[ 2 ][ NC ] )


INST_UpwindingHelper( 2, 2 );
INST_UpwindingHelper( 3, 2 );
INST_UpwindingHelper( 4, 2 );
INST_UpwindingHelper( 5, 2 );
INST_UpwindingHelper( 2, 3 );
INST_UpwindingHelper( 3, 3 );
INST_UpwindingHelper( 4, 3 );
INST_UpwindingHelper( 5, 3 );

#undef INST_UpwindingHelper

} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx
