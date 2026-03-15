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
 * @file CompositionalMultiphaseHybridFVMKernels_impl.hpp
 */

#include "CompositionalMultiphaseHybridFVMKernels.hpp"

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"

#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/kernels/HybridFVMHelperKernels.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"

namespace geos
{
using namespace constitutive;

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
  if constexpr ( NP == 2 )
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
  else if constexpr ( NP == 3 )
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

/******************************** AssemblerKernelHelper ********************************/

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  applyGradient( arrayView1d< real64 const > const & facePres,
                 arrayView1d< real64 const > const & faceGravCoef,
                 arraySlice1d< localIndex const > const & elemToFaces,
                 real64 const & elemPres,
                 real64 const & elemGravCoef,
                 arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & elemPhaseMassDens,
                 arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dElemPhaseMassDens,
                 arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & elemPhaseMob,
                 arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dElemPhaseMob,
                 arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dElemCompFrac_dCompDens,
                 arraySlice2d< real64 const > const & transMatrix,
                 real64 ( & oneSidedVolFlux )[ NF ],
                 real64 ( & dOneSidedVolFlux_dPres )[ NF ],
                 real64 ( & dOneSidedVolFlux_dFacePres )[ NF ][ NF ],
                 real64 ( & dOneSidedVolFlux_dCompDens )[ NF ][ NC ] )
{
  using Deriv = multifluid::DerivativeOffset;

  real64 dPhaseMassDens_dC[ NP ][ NC ]{};
  real64 dPresDif_dCompDens[ NC ]{};
  real64 dPhaseGravDif_dCompDens[ NC ]{};
  real64 dPhaseMobPotDif_dCompDens[ NC ]{};

  // 0) precompute dPhaseDens_dC since it is always computed at the element center
  for( integer ip = 0; ip < NP; ++ip )
  {
    applyChainRule( NC,
                    dElemCompFrac_dCompDens,
                    dElemPhaseMassDens[ip],
                    dPhaseMassDens_dC[ip],
                    Deriv::dC );
  }

  for( integer ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
    for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {

      // depth difference between element center and face center
      real64 const ccGravCoef = elemGravCoef;
      real64 const fGravCoef = faceGravCoef[elemToFaces[jfaceLoc]];
      real64 const gravCoefDif = ccGravCoef - fGravCoef;

      for( integer ip = 0; ip < NP; ++ip )
      {

        // 1) compute the potential diff between the cell center and the face center
        real64 const ccPres = elemPres;
        real64 const fPres  = facePres[elemToFaces[jfaceLoc]];

        // pressure difference
        real64 const presDif = ccPres - fPres;
        real64 const dPresDif_dPres = 1;
        real64 const dPresDif_dFacePres = -1;
        for( integer ic = 0; ic < NC; ++ic )
        {
          dPresDif_dCompDens[ic] = 0.0; // no capillary pressure
        }

        // gravity term
        real64 const phaseGravDif = elemPhaseMassDens[ip] * gravCoefDif;
        real64 const dPhaseGravDif_dPres = dElemPhaseMassDens[ip][Deriv::dP] * gravCoefDif;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPhaseGravDif_dCompDens[ic] = dPhaseMassDens_dC[ip][ic] * gravCoefDif;
        }
        // no density evaluated at the face center

        // potential difference
        real64 const phasePotDif = presDif - phaseGravDif;
        real64 const phaseMobPotDif = elemPhaseMob[ip] * phasePotDif;
        real64 const dPhaseMobPotDif_dPres = dElemPhaseMob[ip][Deriv::dP] * phasePotDif
                                             + elemPhaseMob[ip] * (dPresDif_dPres - dPhaseGravDif_dPres);
        real64 const dPhaseMobPotDif_dFacePres = elemPhaseMob[ip] * dPresDif_dFacePres;
        for( integer ic = 0; ic < NC; ++ic )
        {
          dPhaseMobPotDif_dCompDens[ic] = dElemPhaseMob[ip][Deriv::dC+ic] * phasePotDif
                                          + elemPhaseMob[ip] * (dPresDif_dCompDens[ic] - dPhaseGravDif_dCompDens[ic]);
        }

        // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
        oneSidedVolFlux[ifaceLoc] = oneSidedVolFlux[ifaceLoc]
                                    + transMatrix[ifaceLoc][jfaceLoc] * phaseMobPotDif;
        dOneSidedVolFlux_dPres[ifaceLoc] = dOneSidedVolFlux_dPres[ifaceLoc]
                                           + transMatrix[ifaceLoc][jfaceLoc] * dPhaseMobPotDif_dPres;
        dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc] = dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc]
                                                         + transMatrix[ifaceLoc][jfaceLoc] * dPhaseMobPotDif_dFacePres;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dOneSidedVolFlux_dCompDens[ifaceLoc][ic] = dOneSidedVolFlux_dCompDens[ifaceLoc][ic]
                                                     + transMatrix[ifaceLoc][jfaceLoc] * dPhaseMobPotDif_dCompDens[ic];
        }
      }
    }
  }
}

template< integer NF, integer NC, integer NP, typename MATRIX_VIEW >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  assembleFluxDivergence( localIndex const (&localIds)[ 3 ],
                          globalIndex const rankOffset,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          SortedArrayView< localIndex const > const & regionFilter,
                          arrayView1d< globalIndex const > const & faceDofNumber,
                          arrayView1d< integer const > const & isBoundaryFace,
                          arrayView1d< real64 const > const & mimFaceGravCoef,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          real64 const & elemGravCoef,
                          integer const useTotalMassEquation,
                          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseDens,
                          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseDens,
                          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                          ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                          arraySlice2d< real64 const > const & transMatrixGrav,
                          real64 const (&oneSidedVolFlux)[ NF ],
                          real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                          real64 const & dt,
                          MATRIX_VIEW const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
{
  using namespace compositionalMultiphaseUtilities;
  integer constexpr NDOF = NC+1;

  // dof numbers
  globalIndex dofColIndicesElemVars[ NDOF*(NF+1) ]{};
  globalIndex dofColIndicesFaceVars[ NF ]{};
  integer faceIndexMap[ NF ]{}; // Maps local face index to position in compact DOF array
  integer numNonBoundaryFaces = 0;
  for( integer idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[idof] = elemDofNumber[localIds[0]][localIds[1]][localIds[2]] + idof;
  }

  // divergence of fluxes
  real64 divMassFluxes[ NC ]{};
  real64 dDivMassFluxes_dElemVars[ NC ][ NDOF*(NF+1) ]{};
  real64 dDivMassFluxes_dFaceVars[ NC ][ NF ]{};

  // auxiliary variables for upwinding

  // upwinding phase buoyancy transport coefficients
  real64 upwPhaseViscCoef[ NP ][ NC ]{};
  real64 dUpwPhaseViscCoef_dPres[ NP ][ NC ]{};
  real64 dUpwPhaseViscCoef_dCompDens[ NP ][ NC ][ NC ]{};
  globalIndex upwViscDofNumber = 0;

  // gravity term: ( \rho_l - \rho_m ) g \Delta z
  real64 phaseGravTerm[ NP ][ NP-1 ]{};
  real64 dPhaseGravTerm_dPres[ NP ][ NP-1 ][ 2 ]{};
  real64 dPhaseGravTerm_dCompDens[ NP ][ NP-1 ][ 2 ][ NC ]{};

  // upwinding phase buoyancy transport coefficients
  real64 upwPhaseGravCoef[ NP ][ NP-1 ][ NC ]{};
  real64 dUpwPhaseGravCoef_dPres[ NP ][ NP-1 ][ NC ][ 2 ]{};
  real64 dUpwPhaseGravCoef_dCompDens[ NP ][ NP-1 ][ NC ][ 2 ][ NC ]{};

  // for each element, loop over the one-sided faces
  for( integer ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {

    // Contribute and Collect face DOF number only if this is not a boundary face
    if( isBoundaryFace[elemToFaces[ifaceLoc]] == 0 )
    {

      faceIndexMap[ifaceLoc] = numNonBoundaryFaces;
      dofColIndicesFaceVars[numNonBoundaryFaces] = faceDofNumber[elemToFaces[ifaceLoc]];
      numNonBoundaryFaces++;


      // 1) Find if there is a neighbor, and if there is, grab the indices of the neighbor element
      localIndex neighborIds[ 3 ] = { localIds[0], localIds[1], localIds[2] };
      hybridFVMKernels::CellConnectivity::isNeighborFound( localIds,
                                                           ifaceLoc,
                                                           elemRegionList,
                                                           elemSubRegionList,
                                                           elemList,
                                                           regionFilter,
                                                           elemToFaces,
                                                           neighborIds );
      localIndex const neighborDofNumber = elemDofNumber[neighborIds[0]][neighborIds[1]][neighborIds[2]];

      // 2) *************** Assemble viscous terms ******************

      // 2.a) Compute the upwinded x_{c, \ell} \rho_{\ell} \frac{\lambda_{\ell}}{\lambda_T} for each phase at this face
      UpwindingHelper::upwindViscousCoefficient< NC, NP >( localIds,
                                                           neighborIds,
                                                           phaseDens,
                                                           dPhaseDens,
                                                           phaseMob,
                                                           dPhaseMob,
                                                           dCompFrac_dCompDens,
                                                           phaseCompFrac,
                                                           dPhaseCompFrac,
                                                           elemDofNumber,
                                                           oneSidedVolFlux[ifaceLoc],
                                                           upwPhaseViscCoef,
                                                           dUpwPhaseViscCoef_dPres,
                                                           dUpwPhaseViscCoef_dCompDens,
                                                           upwViscDofNumber );

      // 2.b) Add the \x_{c,\ell} \rho_{\ell} \frac{\lambda_{\ell}}{\lambda_T} q_T of this face to the divergence of the flux in this cell
      assembleViscousFlux< NF, NC, NP >( ifaceLoc,
                                         oneSidedVolFlux,
                                         dOneSidedVolFlux_dPres,
                                         dOneSidedVolFlux_dFacePres,
                                         dOneSidedVolFlux_dCompDens,
                                         upwPhaseViscCoef,
                                         dUpwPhaseViscCoef_dPres,
                                         dUpwPhaseViscCoef_dCompDens,
                                         elemDofNumber[localIds[0]][localIds[1]][localIds[2]],
                                         neighborDofNumber,
                                         upwViscDofNumber,
                                         dt,
                                         divMassFluxes,
                                         dDivMassFluxes_dElemVars,
                                         dDivMassFluxes_dFaceVars,
                                         dofColIndicesElemVars );


      // 3) *************** Assemble buoyancy terms ******************

      real64 const transGravCoef = (localIds[0] != neighborIds[0] || localIds[1] != neighborIds[1] || localIds[2] != neighborIds[2])
                                   * transMatrixGrav[ifaceLoc][ifaceLoc] * (elemGravCoef - mimFaceGravCoef[elemToFaces[ifaceLoc]]);

      // 3.a) Compute the upwinded x_{c, \ell} \rho_{\ell} \frac{\lambda_{\ell}\lambda_m}{\lambda_T}
      //      and (\rho_{\ell} - \rho_m) g \Delta z for each phase at this face
      UpwindingHelper::upwindBuoyancyCoefficient< NC, NP >( localIds,
                                                            neighborIds,
                                                            transGravCoef,
                                                            phaseDens,
                                                            dPhaseDens,
                                                            phaseMassDens,
                                                            dPhaseMassDens,
                                                            phaseMob,
                                                            dPhaseMob,
                                                            dCompFrac_dCompDens,
                                                            phaseCompFrac,
                                                            dPhaseCompFrac,
                                                            phaseGravTerm,
                                                            dPhaseGravTerm_dPres,
                                                            dPhaseGravTerm_dCompDens,
                                                            upwPhaseGravCoef,
                                                            dUpwPhaseGravCoef_dPres,
                                                            dUpwPhaseGravCoef_dCompDens );

      // 3.b) Add the buoyancy term of this face to the divergence of the flux in this cell
      assembleBuoyancyFlux< NF, NC, NP >( ifaceLoc,
                                          phaseGravTerm,
                                          dPhaseGravTerm_dPres,
                                          dPhaseGravTerm_dCompDens,
                                          upwPhaseGravCoef,
                                          dUpwPhaseGravCoef_dPres,
                                          dUpwPhaseGravCoef_dCompDens,
                                          dt,
                                          divMassFluxes,
                                          dDivMassFluxes_dElemVars );

    }
    else
    {
      faceIndexMap[ifaceLoc] = -1; // Mark boundary faces
    }

  }

  if( useTotalMassEquation > 0 )
  {
    // Apply equation/variable change transformation(s)
    real64 work[NDOF * ( NF + 1 )];
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NDOF * ( NF + 1 ), dDivMassFluxes_dElemVars, work );
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NF, dDivMassFluxes_dFaceVars, work );
    shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, divMassFluxes );
  }

  // we are ready to assemble the local flux and its derivatives
  // no need for atomic adds - each row is assembled by a single thread

  for( integer ic = 0; ic < NC; ++ic )
  {
    localIndex const eqnRowLocalIndex =
      LvArray::integerConversion< localIndex >( elemDofNumber[localIds[0]][localIds[1]][localIds[2]] + ic - rankOffset );

    GEOS_ASSERT_GE( eqnRowLocalIndex, 0 );
    GEOS_ASSERT_GT( localMatrix.numRows(), eqnRowLocalIndex );

    // residual
    localRhs[eqnRowLocalIndex] = localRhs[eqnRowLocalIndex] + divMassFluxes[ic];

    // jacobian -- derivative wrt elem centered vars
    // Need to compact both DOF indices and derivatives to only include non-boundary neighbors
    globalIndex compactElemDofs[ NDOF*(NF+1) ];
    real64 compactElemDerivs[ NDOF*(NF+1) ];

    // Copy current element DOFs and derivatives
    for( integer idof = 0; idof < NDOF; ++idof )
    {
      compactElemDofs[idof] = dofColIndicesElemVars[idof];
      compactElemDerivs[idof] = dDivMassFluxes_dElemVars[ic][idof];
    }

    // Copy neighbor DOFs and derivatives only for non-boundary faces
    integer compactIdx = NDOF;
    for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      if( faceIndexMap[jfaceLoc] >= 0 )  // non-boundary face
      {
        integer const neighborOffset = NDOF * (jfaceLoc + 1);
        for( integer idof = 0; idof < NDOF; ++idof )
        {
          compactElemDofs[compactIdx] = dofColIndicesElemVars[neighborOffset + idof];
          compactElemDerivs[compactIdx] = dDivMassFluxes_dElemVars[ic][neighborOffset + idof];
          compactIdx++;
        }
      }
    }

    localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                              compactElemDofs,
                                                              compactElemDerivs,
                                                              compactIdx );

    // jacobian -- derivatives wrt face centered vars (only non-boundary faces)
    // Need to compact the derivatives array to only include non-boundary faces
    if( numNonBoundaryFaces > 0 )
    {
      real64 compactFaceDerivs[ NF ]{};
      for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        if( faceIndexMap[jfaceLoc] >= 0 )
        {
          compactFaceDerivs[faceIndexMap[jfaceLoc]] = dDivMassFluxes_dFaceVars[ic][jfaceLoc];
        }
      }

      localMatrix.template addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                                &dofColIndicesFaceVars[0],
                                                                compactFaceDerivs,
                                                                numNonBoundaryFaces );
    }
  }
}

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  assembleViscousFlux( localIndex const ifaceLoc,
                       real64 const (&oneSidedVolFlux)[ NF ],
                       real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                       real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                       real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                       real64 const (&upwPhaseViscCoef)[ NP ][ NC ],
                       real64 const (&dUpwPhaseViscCoef_dPres)[ NP ][ NC ],
                       real64 const (&dUpwPhaseViscCoef_dCompDens)[ NP ][ NC ][ NC ],
                       globalIndex const elemDofNumber,
                       globalIndex const neighborDofNumber,
                       globalIndex const upwViscDofNumber,
                       real64 const & dt,
                       real64 ( & divMassFluxes )[ NC ],
                       real64 ( & dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ],
                       real64 ( & dDivMassFluxes_dFaceVars )[ NC ][ NF ],
                       globalIndex ( & dofColIndicesElemVars )[ (NC+1)*(NF+1) ] )
{
  integer constexpr NDOF = NC+1;
  localIndex const elemVarsOffset = NDOF*(ifaceLoc+1);

  for( integer ip = 0; ip < NP; ++ip )
  {
    for( integer ic = 0; ic < NC; ++ic )
    {
      // compute the mass flux at the one-sided face plus its derivatives
      // add the newly computed flux to the sum

      real64 const dt_upwPhaseViscCoef = dt * upwPhaseViscCoef[ip][ic];

      // residual
      divMassFluxes[ic] = divMassFluxes[ic] + dt_upwPhaseViscCoef * oneSidedVolFlux[ifaceLoc];

      // local derivatives
      dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                        + dt_upwPhaseViscCoef * dOneSidedVolFlux_dPres[ifaceLoc];
      dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                        + ( elemDofNumber == upwViscDofNumber )
                                        * dt * dUpwPhaseViscCoef_dPres[ip][ic] * oneSidedVolFlux[ifaceLoc];
      for( integer jc = 0; jc < NC; ++jc )
      {
        dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                             + dt_upwPhaseViscCoef * dOneSidedVolFlux_dCompDens[ifaceLoc][jc];
        dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                             + ( elemDofNumber == upwViscDofNumber )
                                             * dt * dUpwPhaseViscCoef_dCompDens[ip][ic][jc] * oneSidedVolFlux[ifaceLoc];
      }

      // neighbor derivatives
      dDivMassFluxes_dElemVars[ic][elemVarsOffset] = dDivMassFluxes_dElemVars[ic][elemVarsOffset]
                                                     + ( elemDofNumber != upwViscDofNumber )
                                                     * dt * dUpwPhaseViscCoef_dPres[ip][ic] * oneSidedVolFlux[ifaceLoc];

      for( integer jc = 0; jc < NC; ++jc )
      {
        dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                            + ( elemDofNumber != upwViscDofNumber )
                                                            * dt * dUpwPhaseViscCoef_dCompDens[ip][ic][jc] * oneSidedVolFlux[ifaceLoc];
      }

      for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        dDivMassFluxes_dFaceVars[ic][jfaceLoc] = dDivMassFluxes_dFaceVars[ic][jfaceLoc]
                                                 + dt_upwPhaseViscCoef * dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc];
      }
    }
  }

  // collect the relevant dof numbers
  for( integer idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[elemVarsOffset+idof] = neighborDofNumber + idof;
  }
}

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  assembleBuoyancyFlux( localIndex const ifaceLoc,
                        real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                        real64 const (&dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ],
                        real64 const (&dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ],
                        real64 const (&upwPhaseGravCoef)[ NP ][ NP-1 ][ NC ],
                        real64 const (&dUpwPhaseGravCoef_dPres)[ NP ][ NP-1 ][ NC ][ 2 ],
                        real64 const (&dUpwPhaseGravCoef_dCompDens)[ NP ][ NP-1 ][ NC ][ 2 ][ NC ],
                        real64 const & dt,
                        real64 ( & divMassFluxes )[ NC ],
                        real64 ( & dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ] )
{
  integer constexpr NDOF = NC+1;
  localIndex const elemVarsOffset = NDOF*(ifaceLoc+1);

  for( integer ip = 0; ip < NP; ++ip )
  {
    for( integer jp = 0; jp < NP - 1; ++jp )
    {
      for( integer ic = 0; ic < NC; ++ic )
      {
        real64 const dt_upwPhaseGravCoef = dt * upwPhaseGravCoef[ip][jp][ic];

        // residual
        divMassFluxes[ic] = divMassFluxes[ic] + dt_upwPhaseGravCoef * phaseGravTerm[ip][jp];

        // local derivatives
        dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                          + dt * dUpwPhaseGravCoef_dPres[ip][jp][ic][Pos::LOCAL] * phaseGravTerm[ip][jp];
        dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                          + dt_upwPhaseGravCoef * dPhaseGravTerm_dPres[ip][jp][Pos::LOCAL];

        for( integer jc = 0; jc < NC; ++jc )
        {
          dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                               + dt * dUpwPhaseGravCoef_dCompDens[ip][jp][ic][Pos::LOCAL][jc] * phaseGravTerm[ip][jp];
          dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                               + dt_upwPhaseGravCoef * dPhaseGravTerm_dCompDens[ip][jp][Pos::LOCAL][jc];
        }

        // neighbor derivatives
        dDivMassFluxes_dElemVars[ic][elemVarsOffset] = dDivMassFluxes_dElemVars[ic][elemVarsOffset]
                                                       + dt * dUpwPhaseGravCoef_dPres[ip][jp][ic][Pos::NEIGHBOR] * phaseGravTerm[ip][jp];
        dDivMassFluxes_dElemVars[ic][elemVarsOffset] = dDivMassFluxes_dElemVars[ic][elemVarsOffset]
                                                       + dt_upwPhaseGravCoef * dPhaseGravTerm_dPres[ip][jp][Pos::NEIGHBOR];

        for( integer jc = 0; jc < NC; ++jc )
        {
          dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                              + dt * dUpwPhaseGravCoef_dCompDens[ip][jp][ic][Pos::NEIGHBOR][jc] * phaseGravTerm[ip][jp];
          dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                              + dt_upwPhaseGravCoef * dPhaseGravTerm_dCompDens[ip][jp][Pos::NEIGHBOR][jc];
        }
      }
    }
  }
}

template< integer NF, integer NC, typename MATRIX_VIEW >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  assembleFaceConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                           arrayView1d< integer const > const & faceGhostRank,
                           arrayView1d< integer const > const & isBoundaryFace,
                           arraySlice1d< localIndex const > const & elemToFaces,
                           globalIndex const elemDofNumber,
                           globalIndex const rankOffset,
                           real64 const (&oneSidedVolFlux)[ NF ],
                           real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                           real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                           real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                           MATRIX_VIEW const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
{
  integer constexpr NDOF = NC+1;

  // fluxes
  real64 dFlux_dElemVars[ NDOF ]{};
  real64 dFlux_dFaceVars[ NF ]{};

  // dof numbers
  globalIndex dofColIndicesElemVars[ NDOF ]{};
  globalIndex dofColIndicesFaceVars[ NF ]{};
  for( integer idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[idof] = elemDofNumber + idof;
  }

  // for each element, loop over the local (one-sided) faces
  for( integer ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    localIndex const kf = elemToFaces[ifaceLoc];

    // Skip ghost faces
    if( faceGhostRank[kf] >= 0 )
    {
      continue;
    }

    // Skip boundary faces - they have Dirichlet constraints, not flux continuity
    if( isBoundaryFace[kf] > 0 )
    {
      continue;
    }

    // flux at this face
    real64 const flux = oneSidedVolFlux[ifaceLoc];
    dFlux_dElemVars[0] = dOneSidedVolFlux_dPres[ifaceLoc];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dFlux_dElemVars[ic+1] = dOneSidedVolFlux_dCompDens[ifaceLoc][ic];
    }

    for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      dFlux_dFaceVars[jfaceLoc] = dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc];
      dofColIndicesFaceVars[jfaceLoc] = faceDofNumber[elemToFaces[jfaceLoc]];
    }

    // dof number of this face constraint
    localIndex const eqnLocalRowIndex = LvArray::integerConversion< localIndex >( faceDofNumber[elemToFaces[ifaceLoc]] - rankOffset );

    GEOS_ASSERT_GE( eqnLocalRowIndex, 0 );
    GEOS_ASSERT_GT( localMatrix.numRows(), eqnLocalRowIndex );

    // residual
    RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnLocalRowIndex], flux );

    // jacobian -- derivatives wrt elem-centered terms
    localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                      &dofColIndicesElemVars[0],
                                                                      &dFlux_dElemVars[0],
                                                                      NDOF );

    // jacobian -- derivatives wrt face pressure terms
    localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                      &dofColIndicesFaceVars[0],
                                                                      &dFlux_dFaceVars[0],
                                                                      NF );
  }
}

/******************************** AssemblerKernel ********************************/

template< integer NF, integer NC, integer NP, typename MATRIX_VIEW >
GEOS_HOST_DEVICE
inline
void
AssemblerKernel::
  compute( localIndex const er, localIndex const esr, localIndex const ei,
           SortedArrayView< localIndex const > const & regionFilter,
           arrayView2d< localIndex const > const & elemRegionList,
           arrayView2d< localIndex const > const & elemSubRegionList,
           arrayView2d< localIndex const > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< integer const > const & faceGhostRank,
           arrayView1d< integer const > const & isBoundaryFace,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arrayView1d< real64 const > const & mimFaceGravCoef,
           arraySlice1d< localIndex const > const & elemToFaces,
           real64 const & elemPres,
           real64 const & elemGravCoef,
           integer const useTotalMassEquation,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           arraySlice2d< real64 const > const & transMatrixGrav,
           MATRIX_VIEW const & localMatrix,
           arrayView1d< real64 > const & localRhs )
{
  // one sided flux
  real64 oneSidedVolFlux[ NF ]{};
  real64 dOneSidedVolFlux_dPres[ NF ]{};
  real64 dOneSidedVolFlux_dFacePres[ NF ][ NF ]{};
  real64 dOneSidedVolFlux_dCompDens[ NF ][ NC ]{};

  localIndex const localIds[3] = { er, esr, ei };

  /*
   * compute auxiliary quantities at the one sided faces of this element:
   * 1) One-sided volumetric fluxes
   * 2) Upwinded mobilities
   */

  // for each one-sided face of the elem,
  // compute the volumetric flux using transMatrix
  AssemblerKernelHelper::applyGradient< NF, NC, NP >( facePres,
                                                      faceGravCoef,
                                                      elemToFaces,
                                                      elemPres,
                                                      elemGravCoef,
                                                      phaseMassDens[er][esr][ei][0],
                                                      dPhaseMassDens[er][esr][ei][0],
                                                      phaseMob[er][esr][ei],
                                                      dPhaseMob[er][esr][ei],
                                                      dCompFrac_dCompDens[er][esr][ei],
                                                      transMatrix,
                                                      oneSidedVolFlux,
                                                      dOneSidedVolFlux_dPres,
                                                      dOneSidedVolFlux_dFacePres,
                                                      dOneSidedVolFlux_dCompDens );

  // at this point, we know the local flow direction in the element
  // so we can upwind the transport coefficients (mobilities) at the one sided faces
  // ** this function needs non-local information **
  if( elemGhostRank < 0 )
  {
    /*
     * perform assembly in this element in two steps:
     * 1) mass conservation equations
     * 2) face constraints
     */

    // use the computed one sided vol fluxes and the upwinded mobilities
    // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
    AssemblerKernelHelper::assembleFluxDivergence< NF, NC, NP >( localIds,
                                                                 rankOffset,
                                                                 elemRegionList,
                                                                 elemSubRegionList,
                                                                 elemList,
                                                                 regionFilter,
                                                                 faceDofNumber,
                                                                 isBoundaryFace,
                                                                 mimFaceGravCoef,
                                                                 elemToFaces,
                                                                 elemGravCoef,
                                                                 useTotalMassEquation,
                                                                 phaseDens,
                                                                 dPhaseDens,
                                                                 phaseMassDens,
                                                                 dPhaseMassDens,
                                                                 phaseMob,
                                                                 dPhaseMob,
                                                                 dCompFrac_dCompDens,
                                                                 phaseCompFrac,
                                                                 dPhaseCompFrac,
                                                                 elemDofNumber,
                                                                 transMatrixGrav,
                                                                 oneSidedVolFlux,
                                                                 dOneSidedVolFlux_dPres,
                                                                 dOneSidedVolFlux_dFacePres,
                                                                 dOneSidedVolFlux_dCompDens,
                                                                 dt,
                                                                 localMatrix,
                                                                 localRhs );
  }

  // use the computed one sided vol fluxes to assemble the constraints
  // enforcing flux continuity at this element's faces
  AssemblerKernelHelper::assembleFaceConstraints< NF, NC >( faceDofNumber,
                                                            faceGhostRank,
                                                            isBoundaryFace,
                                                            elemToFaces,
                                                            elemDofNumber[er][esr][ei],
                                                            rankOffset,
                                                            oneSidedVolFlux,
                                                            dOneSidedVolFlux_dPres,
                                                            dOneSidedVolFlux_dFacePres,
                                                            dOneSidedVolFlux_dCompDens,
                                                            localMatrix,
                                                            localRhs );
}

/******************************** FluxKernel ********************************/

template< integer NF, integer NC, integer NP, typename IP_TYPE, typename MATRIX_VIEW >
void
FluxKernel::
  launch( localIndex er, localIndex esr,
          CellElementSubRegion const & subRegion,
          constitutive::PermeabilityBase const & permeabilityModel,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView1d< globalIndex const > const & faceDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView1d< integer const > const & isBoundaryFace,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & faceGravCoef,
          arrayView1d< real64 const > const & mimFaceGravCoef,
          arrayView1d< real64 const > const & transMultiplier,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
          globalIndex const rankOffset,
          real64 const lengthTolerance,
          real64 const dt,
          integer const useTotalMassEquation,
          MATRIX_VIEW const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  // get the cell-centered DOF numbers and ghost rank for the assembly
  arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

  // get the map from elem to faces
  arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

  // get the cell-centered pressures
  arrayView1d< real64 const > const & elemPres  =
    subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );

  // get the element data needed for transmissibility computation
  arrayView2d< real64 const > const & elemCenter =
    subRegion.getReference< array2d< real64 > >( CellElementSubRegion::viewKeyStruct::elementCenterString() );
  arrayView1d< real64 const > const & elemVolume =
    subRegion.getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString() );

  // TODO add this dependency to the compute function
  //arrayView3d< real64 const > const elemdPermdPres = permeabilityModel.dPerm_dPressure();

  arrayView3d< real64 const > const & elemPerm = permeabilityModel.permeability();

  // get the cell-centered depth
  arrayView1d< real64 const > const & elemGravCoef =
    subRegion.getReference< array1d< real64 > >( fields::flow::gravityCoefficient::key() );

  // assemble the residual and Jacobian element by element
  // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const ei )
  {
    // transmissibility matrix
    stackArray2d< real64, NF *NF > transMatrix( NF, NF );
    stackArray2d< real64, NF *NF > transMatrixGrav( NF, NF );

    real64 const perm[ 3 ] = { elemPerm[ei][0][0], elemPerm[ei][0][1], elemPerm[ei][0][2] };

    // recompute the local transmissibility matrix at each iteration
    // we can decide later to precompute transMatrix if needed
    IP_TYPE::template compute< NF >( nodePosition,
                                     transMultiplier,
                                     faceToNodes,
                                     elemToFaces[ei],
                                     elemCenter[ei],
                                     elemVolume[ei],
                                     perm,
                                     lengthTolerance,
                                     transMatrix );

    // currently the gravity term in the transport scheme is treated as in MRST, that is, always with TPFA
    // this is why below we have to recompute the TPFA transmissibility in addition to the transmissibility matrix above
    // TODO: treat the gravity term with a consistent inner product
    mimeticInnerProduct::TPFAInnerProduct::compute< NF >( nodePosition,
                                                          transMultiplier,
                                                          faceToNodes,
                                                          elemToFaces[ei],
                                                          elemCenter[ei],
                                                          elemVolume[ei],
                                                          perm,
                                                          lengthTolerance,
                                                          transMatrixGrav );

    // perform flux assembly in this element
    compositionalMultiphaseHybridFVMKernels::AssemblerKernel::compute< NF, NC, NP >( er, esr, ei,
                                                                                     regionFilter,
                                                                                     elemRegionList,
                                                                                     elemSubRegionList,
                                                                                     elemList,
                                                                                     faceDofNumber,
                                                                                     faceGhostRank,
                                                                                     isBoundaryFace,
                                                                                     facePres,
                                                                                     faceGravCoef,
                                                                                     mimFaceGravCoef,
                                                                                     elemToFaces[ei],
                                                                                     elemPres[ei],
                                                                                     elemGravCoef[ei],
                                                                                     useTotalMassEquation,
                                                                                     phaseDens,
                                                                                     dPhaseDens,
                                                                                     phaseMassDens,
                                                                                     dPhaseMassDens,
                                                                                     phaseMob,
                                                                                     dPhaseMob,
                                                                                     dCompFrac_dCompDens,
                                                                                     phaseCompFrac,
                                                                                     dPhaseCompFrac,
                                                                                     elemDofNumber,
                                                                                     elemGhostRank[ei],
                                                                                     rankOffset,
                                                                                     dt,
                                                                                     transMatrix,
                                                                                     transMatrixGrav,
                                                                                     localMatrix,
                                                                                     localRhs );
  } );
}

/******************************** DirichletFluxKernel ********************************/

template< integer NF, integer NC, integer NP, typename IP_TYPE, typename MATRIX_VIEW >
void
DirichletFluxKernel::
  launch( integer const numPhases,
          localIndex const er,
          localIndex const esr,
          CellElementSubRegion const & subRegion,
          constitutive::PermeabilityBase const & permeabilityModel,
          SortedArrayView< localIndex const > const & boundaryFaceSet,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & faceTemp,
          arrayView2d< real64 const, compflow::USD_PHASE > const & facePhaseMob,
          arrayView2d< real64 const, compflow::USD_PHASE > const & facePhaseMassDens,
          arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & facePhaseCompFrac,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          arrayView1d< globalIndex const > const & faceDofNumber,
          arrayView1d< real64 const > const & faceGravCoef,
          arrayView1d< real64 const > const & transMultiplier,
          real64 const lengthTolerance,
          real64 const dt,
          globalIndex const rankOffset,
          integer const useTotalMassEquation,
          CompFlowAccessors const & compFlowAccessors,
          MultiFluidAccessors const & multiFluidAccessors,
          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
          MATRIX_VIEW const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( faceTemp );
  using Deriv = multifluid::DerivativeOffset;

  // Get element data
  arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();
  arrayView1d< real64 const > const & elemPres = subRegion.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const & elemGravCoef = subRegion.getField< fields::flow::gravityCoefficient >();
  arrayView2d< real64 const > const & elemCenter = subRegion.getReference< array2d< real64 > >( CellElementSubRegion::viewKeyStruct::elementCenterString() );
  arrayView1d< real64 const > const & elemVolume = subRegion.getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString() );
  arrayView3d< real64 const > const & elemPerm = permeabilityModel.permeability();
  arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

  // Get compositional flow fields
  auto phaseMob = compFlowAccessors.get( fields::flow::phaseMobility{} );
  auto dPhaseMob = compFlowAccessors.get( fields::flow::dPhaseMobility{} );
  auto dCompFrac_dCompDens = compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity{} );

  // Get multifluid fields
  auto phaseMassDens = multiFluidAccessors.get( fields::multifluid::phaseMassDensity{} );
  auto dPhaseMassDens = multiFluidAccessors.get( fields::multifluid::dPhaseMassDensity{} );
  auto phaseCompFrac = multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} );
  auto dPhaseCompFrac = multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction{} );

  // Loop over boundary faces
  forAll< parallelDevicePolicy<> >( boundaryFaceSet.size(), [=] GEOS_DEVICE ( localIndex const iset )
  {
    localIndex const kf = boundaryFaceSet[iset];

    // Find the element adjacent to this boundary face
    localIndex erAdj = -1, esrAdj = -1, eiAdj = -1;
    for( integer ke = 0; ke < elemRegionList.size( 1 ); ++ke )
    {
      if( elemRegionList[kf][ke] == er && elemSubRegionList[kf][ke] == esr )
      {
        erAdj = elemRegionList[kf][ke];
        esrAdj = elemSubRegionList[kf][ke];
        eiAdj = elemList[kf][ke];
        break;
      }
    }

    // Skip if no adjacent element in target region
    if( eiAdj < 0 || elemGhostRank[eiAdj] >= 0 )
    {
      return;
    }

    // Compute one-sided transmissibility
    stackArray2d< real64, NF * NF > transMatrix( NF, NF );
    real64 const perm[3] = { elemPerm[eiAdj][0][0], elemPerm[eiAdj][0][1], elemPerm[eiAdj][0][2] };

    IP_TYPE::template compute< NF >( nodePosition,
                                     transMultiplier,
                                     faceToNodes,
                                     elemToFaces[eiAdj],
                                     elemCenter[eiAdj],
                                     elemVolume[eiAdj],
                                     perm,
                                     lengthTolerance,
                                     transMatrix );

    // Find local face index
    integer ifaceLoc = -1;
    for( integer j = 0; j < NF; ++j )
    {
      if( elemToFaces[eiAdj][j] == kf )
      {
        ifaceLoc = j;
        break;
      }
    }
    if( ifaceLoc < 0 )
    {
      return;
    }

    // Get all global face indices of the cell touching the boundary
    stackArray1d< localIndex, NF > cellFaces( NF );
    for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      cellFaces[jfaceLoc] = elemToFaces[eiAdj][jfaceLoc];
    }

    // Component fluxes
    real64 compFlux[NC]{};
    real64 dCompFlux_dP[NC]{};
    real64 dCompFlux_dC[NC][NC]{};
    real64 dCompFlux_dFaceP[NC][NF]{}; // Derivatives w.r.t. face pressures

    // Loop over phases to compute component fluxes
    // Use face-based properties for boundary conditions (upwinding)
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      real64 dDensMean_dC[NC]{};
      real64 dF_dC[NC]{};
      real64 dProp_dC[NC]{};

      // Use average of element and face phase mass density for gravity term
      real64 const elemDens = phaseMassDens[erAdj][esrAdj][eiAdj][0][ip];
      real64 const faceDens = facePhaseMassDens[kf][ip];
      real64 const densMean = 0.5 * (elemDens + faceDens);

      applyChainRule( NC,
                      dCompFrac_dCompDens[erAdj][esrAdj][eiAdj],
                      dPhaseMassDens[erAdj][esrAdj][eiAdj][0][ip],
                      dProp_dC,
                      Deriv::dC );

      real64 const dDensMean_dP = 0.5 * dPhaseMassDens[erAdj][esrAdj][eiAdj][0][ip][Deriv::dP];
      for( integer jc = 0; jc < NC; ++jc )
      {
        dDensMean_dC[jc] = 0.5 * dProp_dC[jc];
      }

      // Compute flux by looping over all connected faces with consistent transmissibility
      real64 f = 0.0;
      real64 dF_dP = 0.0;
      real64 dF_dFaceP[NF]{}; // Derivatives of flux w.r.t. face pressures
      for( integer jc = 0; jc < NC; ++jc )
      {
        dF_dC[jc] = 0.0;
      }

      // Sum contributions from all connected faces
      for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        // Each stencil face j uses its own face pressure
        localIndex const kfj = cellFaces[jfaceLoc];
        real64 const gravTimesDz_j = elemGravCoef[eiAdj] - faceGravCoef[kfj];
        real64 const Tij = transMatrix[ifaceLoc][jfaceLoc];
        real64 const potDif_j = elemPres[eiAdj] - facePres[kfj] - densMean * gravTimesDz_j;

        // Accumulate flux and derivatives
        f     += Tij * potDif_j;
        dF_dP += Tij * ( 1.0 - dDensMean_dP * gravTimesDz_j );

        // Correct Jacobian: d(potDif_j)/d(facePres[kfj]) = -1, so dF/d(facePres_j) = -T_ij
        dF_dFaceP[jfaceLoc] = -Tij;

        for( integer jc = 0; jc < NC; ++jc )
        {
          dF_dC[jc] += -Tij * dDensMean_dC[jc] * gravTimesDz_j;
        }
      }

      // Use element mobility (simplified upwinding for Dirichlet BC)
      // Upwind phase component fraction based on flow direction
      // Use the total flux f for upwinding decision
      real64 const beta = (f > 0.0) ? 1.0  : 0.0;
      real64 const facePhaseMobility = facePhaseMob[kf][ip];
      real64 const phaseMobility = beta * phaseMob[erAdj][esrAdj][eiAdj][ip] + (1.0 - beta) * facePhaseMobility;
      real64 const phaseFlux = phaseMobility * f;
      real64 const dPhaseFlux_dP = phaseMobility * dF_dP + beta * dPhaseMob[erAdj][esrAdj][eiAdj][ip][Deriv::dP] * f;
      real64 dPhaseFlux_dC[NC];
      real64 dPhaseFlux_dFaceP[NF];
      for( integer jc = 0; jc < NC; ++jc )
      {
        dPhaseFlux_dC[jc] = phaseMobility * dF_dC[jc] + beta * dPhaseMob[erAdj][esrAdj][eiAdj][ip][Deriv::dC+jc] * f;
      }
      for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        dPhaseFlux_dFaceP[jfaceLoc] = phaseMobility * dF_dFaceP[jfaceLoc];
      }

      // Use face-based phase composition for boundary conditions when flow is INTO the domain (potDif < 0)
      // Use element phase composition when flow is OUT of the domain (potDif > 0)
      for( integer ic = 0; ic < NC; ++ic )
      {
        real64 const ycpElem = phaseCompFrac[erAdj][esrAdj][eiAdj][0][ip][ic];
        real64 const ycpFace = facePhaseCompFrac[kf][ip][ic];

        // Upwind phase component fraction based on flow direction
        real64 const ycp = beta * ycpElem + (1.0 - beta) * ycpFace;

        compFlux[ic] += ycp * phaseFlux;

        // For derivatives, use element properties (simplified for boundary conditions)
        dCompFlux_dP[ic] += dPhaseFlux_dP * ycp + beta * phaseFlux * dPhaseCompFrac[erAdj][esrAdj][eiAdj][0][ip][ic][Deriv::dP];

        // Compute dycpElem/dC for this component ic
        real64 dycpElem_dC[NC];
        applyChainRule( NC,
                        dCompFrac_dCompDens[erAdj][esrAdj][eiAdj],
                        dPhaseCompFrac[erAdj][esrAdj][eiAdj][0][ip][ic],
                        dycpElem_dC,
                        Deriv::dC );
        for( integer jc = 0; jc < NC; ++jc )
        {
          // d(ycp * phaseFlux)/dC = dycp/dC * phaseFlux + ycp * dPhaseFlux/dC
          // dycp/dC = beta * dycpElem/dC (ycpFace is BC, independent of C)
          dCompFlux_dC[ic][jc] += ycp * dPhaseFlux_dC[jc] + beta * phaseFlux * dycpElem_dC[jc];
        }

        // Add derivatives w.r.t. face pressures
        for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
        {
          dCompFlux_dFaceP[ic][jfaceLoc] += dPhaseFlux_dFaceP[jfaceLoc] * ycp;
        }
      }
    }

    // Assemble into residual and Jacobian
    real64 localFlux[NC];
    real64 localFluxJacobian[NC][NC+1];

    for( integer ic = 0; ic < NC; ++ic )
    {
      localFlux[ic] = dt * compFlux[ic];
      localFluxJacobian[ic][0] = dt * dCompFlux_dP[ic];
      for( integer jc = 0; jc < NC; ++jc )
      {
        localFluxJacobian[ic][jc+1] = dt * dCompFlux_dC[ic][jc];
      }
    }

    // Apply total mass equation transformation if needed
    if( useTotalMassEquation )
    {
      real64 work[NC+1]{};
      compositionalMultiphaseUtilities::shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC+1, localFluxJacobian, work );
      compositionalMultiphaseUtilities::shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, localFlux );
    }

    // Add to global system
    globalIndex const elemDof = elemDofNumber[erAdj][esrAdj][eiAdj];
    localIndex const localRow = LvArray::integerConversion< localIndex >( elemDof - rankOffset );

    for( integer ic = 0; ic < NC; ++ic )
    {
      RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow + ic], localFlux[ic] );

      globalIndex dofColIndices[NC+1];
      dofColIndices[0] = elemDof;
      for( integer jc = 0; jc < NC; ++jc )
      {
        dofColIndices[jc+1] = elemDof + jc + 1;
      }

      localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + ic,
                                                                        dofColIndices,
                                                                        localFluxJacobian[ic],
                                                                        NC+1 );

      // Add contributions from face pressure derivatives
      for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        localIndex const kfj = cellFaces[jfaceLoc];
        globalIndex const faceDof = faceDofNumber[kfj];

        real64 facePressureJacobian = dt * dCompFlux_dFaceP[ic][jfaceLoc];

        // Apply total mass equation transformation if needed
        if( useTotalMassEquation && ic == 0 )
        {
          // For total mass equation, sum all component derivatives
          facePressureJacobian = 0.0;
          for( integer icc = 0; icc < NC; ++icc )
          {
            facePressureJacobian += dt * dCompFlux_dFaceP[icc][jfaceLoc];
          }
        }

        localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + ic,
                                                                          &faceDof,
                                                                          &facePressureJacobian,
                                                                          1 );
      }

    }

  } );
}

/******************************** EvaluateBCFacePropertiesKernel ********************************/

/**
 * @brief Evaluate constitutive properties at BC face conditions using flash calculations
 * @tparam NC number of components
 * @tparam NP number of phases
 * @param[in] numPhases number of phases in the simulation
 * @param[in] boundaryFaceSet sorted array view of boundary face indices
 * @param[in] facePres pressure values at boundary faces
 * @param[in] faceTemp temperature values at boundary faces
 * @param[in] faceCompFrac component fraction values at boundary faces
 * @param[in] elemRegionList face-to-element region mapping
 * @param[in] elemSubRegionList face-to-element subregion mapping
 * @param[in] elemList face-to-element index mapping
 * @param[in] er target element region index
 * @param[in] esr target element subregion index
 * @param[in] fluid reference to multifluid constitutive model
 * @param[in] relperm reference to relative permeability constitutive model
 * @param[out] facePhaseMob computed phase mobility at BC faces
 * @param[out] facePhaseMassDens computed phase mass density at BC faces
 * @param[out] facePhaseCompFrac computed phase component fractions at BC faces
 *
 * @note This function uses serialPolicy (host execution) because CUDA does not allow
 *       extended device lambdas inside generic lambda expressions. The output arrays
 *       must be explicitly moved to device memory after this function completes.
 */
template< integer NC, integer NP >
void
evaluateBCFaceProperties( integer const numPhases,
                          SortedArrayView< localIndex const > const & boundaryFaceSet,
                          arrayView1d< real64 const > const & facePres,
                          arrayView1d< real64 const > const & faceTemp,
                          arrayView2d< real64 const, compflow::USD_COMP > const & faceCompFrac,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          localIndex const er,
                          localIndex const esr,
                          constitutive::MultiFluidBase & fluid,
                          constitutive::RelativePermeabilityBase & relperm,
                          arrayView2d< real64, compflow::USD_PHASE > const & facePhaseMob,
                          arrayView2d< real64, compflow::USD_PHASE > const & facePhaseMassDens,
                          arrayView3d< real64, compflow::USD_PHASE_COMP > const & facePhaseCompFrac )
{
  // Use constitutiveUpdatePassThru to dispatch to concrete fluid and relperm types
  constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    constitutive::constitutiveUpdatePassThru( relperm, [&] ( auto & castedRelperm )
    {
      using RelPermType = TYPEOFREF( castedRelperm );
      typename RelPermType::KernelWrapper relPermWrapper = castedRelperm.createKernelWrapper();
      GEOS_UNUSED_VAR( relPermWrapper );

      // Loop over BC faces and evaluate properties at BC conditions
      // Note: Using serialPolicy (host) because extended device lambdas cannot be defined inside generic lambdas
      // The output arrays will be moved to device after this completes
      forAll< serialPolicy >( boundaryFaceSet.size(), [=, &facePhaseMob, &facePhaseMassDens, &facePhaseCompFrac] ( localIndex const iset )
      {
        localIndex const kf = boundaryFaceSet[iset];

        // Find adjacent element in target region
        localIndex eiAdj = -1;
        for( integer ke = 0; ke < elemRegionList.size( 1 ); ++ke )
        {
          if( elemRegionList[kf][ke] == er && elemSubRegionList[kf][ke] == esr )
          {
            eiAdj = elemList[kf][ke];
            break;
          }
        }
        if( eiAdj < 0 )
          return;

        // Allocate temporary storage for face constitutive properties
        StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseFrac( 1, 1, numPhases );
        StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseDens( 1, 1, numPhases );
        StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseMassDensLocal( 1, 1, numPhases );
        StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseVisc( 1, 1, numPhases );
        StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseEnthalpy( 1, 1, numPhases );
        StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseInternalEnergy( 1, 1, numPhases );
        StackArray< real64, 4, constitutive::MultiFluidBase::MAX_NUM_PHASES * NC,
                    constitutive::multifluid::LAYOUT_PHASE_COMP > facePhaseCompFracLocal( 1, 1, numPhases, NC );

        // Initialize all temporary arrays to zero to ensure deterministic behavior on GPU
        for( integer ip = 0; ip < numPhases; ++ip )
        {
          facePhaseFrac[0][0][ip] = 0.0;
          facePhaseDens[0][0][ip] = 0.0;
          facePhaseMassDensLocal[0][0][ip] = 0.0;
          facePhaseVisc[0][0][ip] = 0.0;
          facePhaseEnthalpy[0][0][ip] = 0.0;
          facePhaseInternalEnergy[0][0][ip] = 0.0;
          for( integer ic = 0; ic < NC; ++ic )
          {
            facePhaseCompFracLocal[0][0][ip][ic] = 0.0;
          }
        }

        real64 faceTotalDens = 0.0;

        // Evaluate fluid properties at BC face conditions using flash calculation
        constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                    facePres[kf],
                                                                    faceTemp[kf],
                                                                    faceCompFrac[kf],
                                                                    facePhaseFrac[0][0],
                                                                    facePhaseDens[0][0],
                                                                    facePhaseMassDensLocal[0][0],
                                                                    facePhaseVisc[0][0],
                                                                    facePhaseEnthalpy[0][0],
                                                                    facePhaseInternalEnergy[0][0],
                                                                    facePhaseCompFracLocal[0][0],
                                                                    faceTotalDens );

        // Evaluate relative permeability at face saturation from flash calculation
//      relPermWrapper.compute( facePhaseFrac[0][0], 0, 0 );

        // Store computed properties in output arrays
        for( integer ip = 0; ip < numPhases; ++ip )
        {
          // Store phase mass density from flash calculation
          facePhaseMassDens[kf][ip] = facePhaseMassDensLocal[0][0][ip];

          // Compute mobility from relative permeability evaluated at face conditions
          real64 const faceKr = facePhaseFrac[0][0][ip]; // phaseRelPerm[eiAdj][0][ip];
          real64 const mu = facePhaseVisc[0][0][ip];
          // Safety check: ensure faceTotalDens and mu are valid before computing mobility
          facePhaseMob[kf][ip] = (mu > 0 && faceTotalDens > 0) ? faceTotalDens * faceKr / mu : 0.0;

          // Store phase composition from flash calculation
          for( integer ic = 0; ic < NC; ++ic )
          {
            facePhaseCompFrac[kf][ip][ic] = facePhaseCompFracLocal[0][0][ip][ic];
          }
        }
      } );
    } ); // end nested constitutiveUpdatePassThru for relperm
  } ); // end constitutiveUpdatePassThru for fluid
}

} // namespace compositionalMultiphaseHybridFVMKernels

} // namespace geos
