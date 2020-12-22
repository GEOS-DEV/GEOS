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
 * @file CompositionalMultiphaseHybridFVMHelperKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMUPWINDINGHELPERKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMUPWINDINGHELPERKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"

namespace geosx
{

namespace CompositionalMultiphaseHybridFVMKernels
{

// struct to specify local and neighbor derivatives
struct Pos
{
  static constexpr integer LOCAL = 0;
  static constexpr integer NEIGHBOR = 1;
};

/******************************** UpwindingHelper ********************************/

struct UpwindingHelper
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief At a given one-sided face, compute the upwind viscous transport coefficient
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions in the domain wrt pressure (non-local)
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions in the domain wrt component fraction (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] oneSidedVolFlux the total volumetric flux at this face
   * @param[out] upwPhaseViscCoef the upwind viscous transport coef at this face
   * @param[out] dUpwPhaseViscCoef_dPres the derivative of the upwind viscous transport coef wrt pressure at this face
   * @param[out] dUpwPhaseViscCoef_dCompDens the derivatives of the upwind viscous transport coef wrt component density at this face
   * @param[out] upwViscDofNumber the dof number of the upwind cell at this face
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  UpwindViscousCoefficient( localIndex const (&localIds)[ 3 ],
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
                            real64 ( & upwPhaseViscCoef )[ NP ][ NC ],
                            real64 ( & dUpwPhaseViscCoef_dPres )[ NP ][ NC ],
                            real64 ( & dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ],
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

  /**
   * @brief At a given one-sided face, compute the upwind viscous transport coefficient
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] transGravCoef
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions in the domain wrt pressure (non-local)
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions in the domain wrt component fraction (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[inout] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[inout] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   * @param[inout] upwPhaseGravCoef the upwinded buoyancy transport coefficient at this face (ifaceLoc)
   * @param[inout] dUpwPhaseGravCoef_dPres the derivative of the upwinded buoyancy transport coefficient wrt pressure
   * @param[inout] dUpwPhaseGravCoef_dCompDens the derivative of the upwinded buoyancy transport coefficient wrt component density
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  UpwindBuoyancyCoefficient( localIndex const (&localIds)[ 3 ],
                             localIndex const (&neighborIds)[ 3 ],
                             real64 const & transGravCoef,
                             ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                             ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                             ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                             ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                             ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                             ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dCompFrac,
                             ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                             ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                             ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                             ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                             ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
                             ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                             ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
                             real64 ( & phaseGravTerm )[ NP ][ NP-1 ],
                             real64 ( & dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                             real64 ( & dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ],
                             real64 ( & upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ],
                             real64 ( & dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ],
                             real64 ( & dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] )
  {
    // 1) Compute the driving force: T ( \rho^{avg}_{\ell} - \rho^{avg}_m ) g \Delta z
    ComputePhaseGravTerm( localIds,
                          neighborIds,
                          transGravCoef,
                          phaseMassDens,
                          dPhaseMassDens_dPres,
                          dPhaseMassDens_dCompFrac,
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


  /**
   * @brief At a given one-sided face, compute the gravCoef multiplied by the difference in phase densities
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] transGravCoef
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[inout] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[inout] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  ComputePhaseGravTerm( localIndex const (&localIds)[ 3 ],
                        localIndex const (&neighborIds)[ 3 ],
                        real64 const & transGravCoef,
                        ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                        ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                        ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dCompFrac,
                        ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                        real64 ( & phaseGravTerm )[ NP ][ NP-1 ],
                        real64 ( & dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                        real64 ( & dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] )
  {
    localIndex const er   = localIds[0];
    localIndex const esr  = localIds[1];
    localIndex const ei   = localIds[2];
    localIndex const ern  = neighborIds[0];
    localIndex const esrn = neighborIds[1];
    localIndex const ein  = neighborIds[2];

    real64 dPhaseMassDens_dCLoc[ NC ] = { 0.0 };
    real64 dPhaseMassDens_dCNeighbor[ NC ] = { 0.0 };
    real64 dPhaseMassDens_dC[ NC ] = { 0.0 };

    for( localIndex ip = 0; ip < NP; ++ip )
    {
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens_dCompFrac[er][esr][ei][0][ip],
                      dPhaseMassDens_dCLoc );
      applyChainRule( NC,
                      dCompFrac_dCompDens[ern][esrn][ein],
                      dPhaseMassDens_dCompFrac[ern][esrn][ein][0][ip],
                      dPhaseMassDens_dCNeighbor );

      localIndex k = 0;
      for( localIndex jp = 0; jp < NP; ++jp )
      {
        if( ip == jp )
        {
          continue;
        }

        phaseGravTerm[ip][k] = -( phaseMassDens[er][esr][ei][0][ip] + phaseMassDens[ern][esrn][ein][0][ip] );
        phaseGravTerm[ip][k] += ( phaseMassDens[er][esr][ei][0][jp] + phaseMassDens[ern][esrn][ein][0][jp] );
        phaseGravTerm[ip][k] *= 0.5 * transGravCoef;

        dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] = ( -dPhaseMassDens_dPres[er][esr][ei][0][ip] + dPhaseMassDens_dPres[er][esr][ei][0][jp] );
        dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] *= 0.5 * transGravCoef;

        dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] = ( -dPhaseMassDens_dPres[ern][esrn][ein][0][ip] + dPhaseMassDens_dPres[ern][esrn][ein][0][jp] );
        dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] *= 0.5 * transGravCoef;

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCLoc[ic];
          dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCNeighbor[ic];
        }
        applyChainRule( NC,
                        dCompFrac_dCompDens[er][esr][ei],
                        dPhaseMassDens_dCompFrac[er][esr][ei][0][jp],
                        dPhaseMassDens_dC );
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
        }
        applyChainRule( NC,
                        dCompFrac_dCompDens[ern][esrn][ein],
                        dPhaseMassDens_dCompFrac[ern][esrn][ein][0][jp],
                        dPhaseMassDens_dC );
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
        }
        ++k;
      }
    }
  }

  /**
   * @brief At a given one-sided face, compute the upwinded total mobility
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] totalMob the upwinded total mobility
   * @param[inout] dTotalMob_dPres the derivative of the upwinded total mobility wrt pressure
   * @param[inout] dTotalMob_dCompDens the derivative of the upwinded total mobility wrt component density
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  ComputeUpwindedTotalMobility( localIndex const (&localIds)[ 3 ],
                                localIndex const (&neighborIds)[ 3 ],
                                ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                                ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                                ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                                real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                                real64 & totalMob,
                                real64 ( & dTotalMob_dPres )[ 2 ],
                                real64 ( & dTotalMob_dCompDens )[ 2 ][ NC ] )
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

  /**
   * @brief Set the element indices used to evaluate the mobility ratios of the buoyancy term in hybrid upwinding
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] gravTerm the gravCoef multiplied by the difference in phase densities
   * @param[in] totalMob the upwinded total mobility
   * @param[in] eru region index of the upwind element
   * @param[in] esru subRegion index of the upwind element
   * @param[in] eiu element index of the upwind element
   * @param[in] posu position (local or neighbor) of the upwind element
   * @param[in] erd region index of the downwind element
   * @param[in] esrd subRegion index of the downwind element
   * @param[in] eid element index of the downwind element
   * @param[in] posd position (local or neighbor) of the downwind element
   */
  GEOSX_HOST_DEVICE
  static void
  SetIndicesForMobilityRatioUpwinding( localIndex const (&localIds)[ 3 ],
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

  /**
   * @brief Set the element indices used to evaluate the total mobility of the buoyancy term in hybrid upwinding
   * @param[in] localIds triplet of indices for the local element
   * @param[in] neighborIds triplet of indices for the neighbor element
   * @param[in] gravTerm gravity term used to upwind
   * @param[out] totalMobIds for each phase, triplet of indices of the upwind element
   * @param[out] totalMobPos for each phase, flag specifying with the upwind element is local or neighbor
   */
  template< localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  SetIndicesForTotalMobilityUpwinding( localIndex const (&localIds)[ 3 ],
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
      // TODO: implement the upwinding here
    }
  }

};



} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMUPWINDINGHELPERKERNELS_HPP
