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

/******************************** PhaseMobilityKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
PhaseMobilityKernel::
  Compute( arraySlice2d< real64 const > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const > const & phaseVisc,
           arraySlice1d< real64 const > const & dPhaseVisc_dPres,
           arraySlice2d< real64 const > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const > const & phaseRelPerm,
           arraySlice2d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64 > const & phaseMob,
           arraySlice1d< real64 > const & dPhaseMob_dPres,
           arraySlice2d< real64 > const & dPhaseMob_dComp )
{
  real64 dRelPerm_dC[NC];
  real64 dVisc_dC[NC];

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    real64 const viscosity = phaseVisc[ip];
    real64 const dVisc_dP = dPhaseVisc_dPres[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseVisc_dComp[ip], dVisc_dC );

    real64 const relPerm = phaseRelPerm[ip];
    real64 dRelPerm_dP = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dRelPerm_dC[ic] = 0.0;
    }

    for( localIndex jp = 0; jp < NP; ++jp )
    {
      real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
      dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_dPres[jp];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_dComp[jp][jc];
      }
    }

    real64 const mobility = relPerm / viscosity;

    phaseMob[ip] = mobility;
    dPhaseMob_dPres[ip] = dRelPerm_dP / viscosity
                          - mobility * dVisc_dP / viscosity;

    // compositional derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] / viscosity
                                - mobility * dVisc_dC[jc] / viscosity;
    }
  }
}

template< localIndex NC, localIndex NP >
void PhaseMobilityKernel::
  Launch( localIndex const size,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    Compute< NC, NP >( dCompFrac_dCompDens[a],
                       phaseVisc[a][0],
                       dPhaseVisc_dPres[a][0],
                       dPhaseVisc_dComp[a][0],
                       phaseRelPerm[a][0],
                       dPhaseRelPerm_dPhaseVolFrac[a][0],
                       dPhaseVolFrac_dPres[a],
                       dPhaseVolFrac_dComp[a],
                       phaseMob[a],
                       dPhaseMob_dPres[a],
                       dPhaseMob_dComp[a] );
  } );
}

template< localIndex NC, localIndex NP >
void PhaseMobilityKernel::
  Launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute< NC, NP >( dCompFrac_dCompDens[a],
                       phaseVisc[a][0],
                       dPhaseVisc_dPres[a][0],
                       dPhaseVisc_dComp[a][0],
                       phaseRelPerm[a][0],
                       dPhaseRelPerm_dPhaseVolFrac[a][0],
                       dPhaseVolFrac_dPres[a],
                       dPhaseVolFrac_dComp[a],
                       phaseMob[a],
                       dPhaseMob_dPres[a],
                       dPhaseMob_dComp[a] );
  } );
}

#define INST_PhaseMobilityKernel( NC, NP ) \
  template \
  void \
  PhaseMobilityKernel:: \
    Launch< NC, NP >( localIndex const size, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseVisc, \
                      arrayView3d< real64 const > const & dPhaseVisc_dPres, \
                      arrayView4d< real64 const > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const > const & phaseRelPerm, \
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64 > const & phaseMob, \
                      arrayView2d< real64 > const & dPhaseMob_dPres, \
                      arrayView3d< real64 > const & dPhaseMob_dComp ); \
  template \
  void \
  PhaseMobilityKernel:: \
    Launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseVisc, \
                      arrayView3d< real64 const > const & dPhaseVisc_dPres, \
                      arrayView4d< real64 const > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const > const & phaseRelPerm, \
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64 > const & phaseMob, \
                      arrayView2d< real64 > const & dPhaseMob_dPres, \
                      arrayView3d< real64 > const & dPhaseMob_dComp )

INST_PhaseMobilityKernel( 1, 1 );
INST_PhaseMobilityKernel( 2, 1 );
INST_PhaseMobilityKernel( 3, 1 );
INST_PhaseMobilityKernel( 4, 1 );
INST_PhaseMobilityKernel( 5, 1 );

INST_PhaseMobilityKernel( 1, 2 );
INST_PhaseMobilityKernel( 2, 2 );
INST_PhaseMobilityKernel( 3, 2 );
INST_PhaseMobilityKernel( 4, 2 );
INST_PhaseMobilityKernel( 5, 2 );

INST_PhaseMobilityKernel( 1, 3 );
INST_PhaseMobilityKernel( 2, 3 );
INST_PhaseMobilityKernel( 3, 3 );
INST_PhaseMobilityKernel( 4, 3 );
INST_PhaseMobilityKernel( 5, 3 );

#undef INST_PhaseMobilityKernel


/******************************** UpwindingHelper ********************************/

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
UpwindingHelper::UpwindViscousTerm( localIndex const er, localIndex const esr, localIndex const ei,
                                    localIndex const ifaceLoc,
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
                                    real64 (& upwPhaseViscCoef)[ NF ][ NP ][ NC ],
                                    real64 (& dUpwPhaseViscCoef_dPres)[ NF ][ NP ][ NC ],
                                    real64 (& dUpwPhaseViscCoef_dCompDens)[ NF ][ NP ][ NC ][ NC ],
                                    globalIndex (& upwViscDofNumber)[ NF ] )
{
  real64 dUpwMobRatio_dCompDens[ NC ] = { 0.0 };
  real64 dUpwDensMobRatio_dCompDens[ NC ] = { 0.0 };
  real64 dPhaseDens_dC[ NC ] = { 0.0 };
  real64 dPhaseCompFrac_dC[ NC ] = { 0.0 };

  // 1) Compute total mobility
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

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    // 2) Compute viscous mobility ratio
    real64 const upwMobRatio = phaseMob[er][esr][ei][ip] / totalMob;
    real64 const dUpwMobRatio_dPres = ( dPhaseMob_dPres[er][esr][ei][ip] - upwMobRatio * dTotalMob_dPres )
                                      / totalMob;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dUpwMobRatio_dCompDens[ic] = ( dPhaseMob_dCompDens[er][esr][ei][ip][ic] - upwMobRatio * dTotalMob_dCompDens[ic] )
                                   / totalMob;
    }

    // 3) Multiply mobility ratio by phase density
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

    // 4) Multiply density mobility ratio by phase comp fraction
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseCompFrac_dCompFrac[er][esr][ei][0][ip][ic],
                      dPhaseCompFrac_dC );
      upwPhaseViscCoef[ifaceLoc][ip][ic] = phaseCompFrac[er][esr][ei][0][ip][ic] * upwDensMobRatio;
      dUpwPhaseViscCoef_dPres[ifaceLoc][ip][ic] = dPhaseCompFrac_dPres[er][esr][ei][0][ip][ic] * upwDensMobRatio
                                                  + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dPres;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dUpwPhaseViscCoef_dCompDens[ifaceLoc][ip][ic][jc] = dPhaseCompFrac_dC[jc] * upwDensMobRatio
                                                            + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dCompDens[jc];
      }
    }
  }
  // 5) Save the dof number of the upwind cell
  upwViscDofNumber[ifaceLoc] = elemDofNumber[er][esr][ei];
}


/******************************** AssemblerKernelHelper ********************************/

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernelHelper::ComputeOneSidedVolFluxes( arrayView1d< real64 const > const & facePres,
                                                 arrayView1d< real64 const > const & dFacePres,
                                                 arrayView1d< real64 const > const & faceGravCoef,
                                                 arraySlice1d< localIndex const > const & elemToFaces,
                                                 real64 const & elemPres,
                                                 real64 const & dElemPres,
                                                 real64 const & elemGravCoef,
                                                 arraySlice1d< real64 const > const & elemPhaseDens,
                                                 arraySlice1d< real64 const > const & dElemPhaseDens_dPres,
                                                 arraySlice2d< real64 const > const & dElemPhaseDens_dCompFrac,
                                                 arraySlice1d< real64 const > const & elemPhaseMob,
                                                 arraySlice1d< real64 const > const & dElemPhaseMob_dPres,
                                                 arraySlice2d< real64 const > const & dElemPhaseMob_dCompDens,
                                                 arraySlice2d< real64 const > const & dElemCompFrac_dCompDens,
                                                 arraySlice2d< real64 const > const & transMatrix,
                                                 real64 (& oneSidedVolFlux)[ NF ],
                                                 real64 (& dOneSidedVolFlux_dPres)[ NF ],
                                                 real64 (& dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                                                 real64 (& dOneSidedVolFlux_dCompDens)[ NF ][ NC ] )
{
  real64 dPhaseDens_dC[ NP ][ NC ] = {{ 0.0 }};
  real64 dPresDif_dCompDens[ NC ] = { 0.0 };
  real64 dPhaseGravDif_dCompDens[ NC ] = { 0.0 };
  real64 dPhaseMobPotDif_dCompDens[ NC ] = { 0.0 };

  // 0) precompute dPhaseDens_dC since it is always computed at the element center
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    applyChainRule( NC,
                    dElemCompFrac_dCompDens,
                    dElemPhaseDens_dCompFrac[ip],
                    dPhaseDens_dC[ip] );
  }

  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {

      // depth difference between element center and face center
      real64 const ccGravCoef = elemGravCoef;
      real64 const fGravCoef = faceGravCoef[elemToFaces[ifaceLoc]];
      real64 const gravCoefDif = ccGravCoef - fGravCoef;

      for( localIndex ip = 0; ip < NP; ++ip )
      {

        // 1) compute the potential diff between the cell center and the face center
        real64 const ccPres = elemPres + dElemPres;
        real64 const fPres  = facePres[elemToFaces[jfaceLoc]] + dFacePres[elemToFaces[jfaceLoc]];

        // pressure difference
        real64 const presDif = ccPres - fPres;
        real64 const dPresDif_dPres = 1;
        real64 const dPresDif_dFacePres = -1;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPresDif_dCompDens[ic] = 0.0; // no capillary pressure
        }

        // gravity term
        real64 const phaseGravDif = elemPhaseDens[ip] * gravCoefDif;
        real64 const dPhaseGravDif_dPres = dElemPhaseDens_dPres[ip] * gravCoefDif;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPhaseGravDif_dCompDens[ic] = dPhaseDens_dC[ip][ic] * gravCoefDif;
        }
        // no density evaluated at the face center

        // potential difference
        real64 const phasePotDif = presDif - phaseGravDif;
        real64 const phaseMobPotDif = elemPhaseMob[ip] * phasePotDif;
        real64 const dPhaseMobPotDif_dPres = dElemPhaseMob_dPres[ip] * phasePotDif
                                             + elemPhaseMob[ip] * (dPresDif_dPres - dPhaseGravDif_dPres);
        real64 const dPhaseMobPotDif_dFacePres = elemPhaseMob[ip] * dPresDif_dFacePres;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPhaseMobPotDif_dCompDens[ic] = dElemPhaseMob_dCompDens[ip][ic] * phasePotDif
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

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernelHelper::UpdateUpwindedCoefficients( localIndex const er, localIndex const esr, localIndex const ei,
                                                   arrayView2d< localIndex const > const & elemRegionList,
                                                   arrayView2d< localIndex const > const & elemSubRegionList,
                                                   arrayView2d< localIndex const > const & elemList,
                                                   SortedArrayView< localIndex const > const & regionFilter,
                                                   arraySlice1d< localIndex const > const & elemToFaces,
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
                                                   real64 const (&oneSidedVolFlux)[ NF ],
                                                   real64 (& upwPhaseViscCoef)[ NF ][ NP ][ NC ],
                                                   real64 (& dUpwPhaseViscCoef_dPres)[ NF ][ NP ][ NC ],
                                                   real64 (& dUpwPhaseViscCoef_dCompDens)[ NF ][ NP ][ NC ][ NC ],
                                                   globalIndex (& upwViscDofNumber)[ NF ] )
{
  // for this element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {

    // we initialize these upw quantities with the values of the local elem
    UpwindingHelper::UpwindViscousTerm< NF, NC, NP >( er, esr, ei, ifaceLoc,
                                                      phaseDens,
                                                      dPhaseDens_dPres,
                                                      dPhaseDens_dCompFrac,
                                                      phaseMob,
                                                      dPhaseMob_dPres,
                                                      dPhaseMob_dCompDens,
                                                      dCompFrac_dCompDens,
                                                      phaseCompFrac,
                                                      dPhaseCompFrac_dPres,
                                                      dPhaseCompFrac_dCompFrac,
                                                      elemDofNumber,
                                                      upwPhaseViscCoef,
                                                      dUpwPhaseViscCoef_dPres,
                                                      dUpwPhaseViscCoef_dCompDens,
                                                      upwViscDofNumber );

    // if the local elem if upstream, we are done, we can proceed to the next one-sided face
    // otherwise, we have to access the properties of the neighbor element
    // this is done on the fly below
    if( oneSidedVolFlux[ifaceLoc] < 0 )
    {

      // the face has at most two adjacent elements
      // one of these two elements is the current element indexed by er, esr, ei
      // but here we are interested in the indices of the other element
      // this other element is "the neighbor" for this one-sided face
      for( localIndex k=0; k<elemRegionList.size( 1 ); ++k )
      {

        localIndex const erNeighbor  = elemRegionList[elemToFaces[ifaceLoc]][k];
        localIndex const esrNeighbor = elemSubRegionList[elemToFaces[ifaceLoc]][k];
        localIndex const eiNeighbor  = elemList[elemToFaces[ifaceLoc]][k];

        // this element is not the current element
        // we have found the neighbor or we are at the boundary
        if( erNeighbor != er || esrNeighbor != esr || eiNeighbor != ei )
        {
          bool const onBoundary       = (erNeighbor == -1 || esrNeighbor == -1 || eiNeighbor == -1);
          bool const neighborInTarget = regionFilter.contains( erNeighbor );

          // if not on boundary, save the mobility and the upwViscDofNumber
          if( !onBoundary && neighborInTarget )
          {
            UpwindingHelper::UpwindViscousTerm< NF, NC, NP >( erNeighbor, esrNeighbor, eiNeighbor, ifaceLoc,
                                                              phaseDens,
                                                              dPhaseDens_dPres,
                                                              dPhaseDens_dCompFrac,
                                                              phaseMob,
                                                              dPhaseMob_dPres,
                                                              dPhaseMob_dCompDens,
                                                              dCompFrac_dCompDens,
                                                              phaseCompFrac,
                                                              dPhaseCompFrac_dPres,
                                                              dPhaseCompFrac_dCompFrac,
                                                              elemDofNumber,
                                                              upwPhaseViscCoef,
                                                              dUpwPhaseViscCoef_dPres,
                                                              dUpwPhaseViscCoef_dCompDens,
                                                              upwViscDofNumber );
          }
          // if the face is on the boundary, use the properties of the local elem
        }
      }
    }
  }
}

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernelHelper::AssembleOneSidedMassFluxes( arrayView1d< globalIndex const > const & faceDofNumber,
                                                   arraySlice1d< localIndex const > const & elemToFaces,
                                                   globalIndex const elemDofNumber,
                                                   globalIndex const rankOffset,
                                                   real64 const (&oneSidedVolFlux)[ NF ],
                                                   real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                                                   real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                                                   real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                                                   real64 const (&upwPhaseViscCoef)[ NF ][ NP ][ NC ],
                                                   real64 const (&dUpwPhaseViscCoef_dPres)[ NF ][ NP ][ NC ],
                                                   real64 const (&dUpwPhaseViscCoef_dCompDens)[ NF ][ NP ][ NC ][ NC ],
                                                   globalIndex const (&upwViscDofNumber)[ NF ],
                                                   real64 const & dt,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  localIndex constexpr NDOF = NC+1;

  // dof numbers
  globalIndex dofColIndicesElemVars[ NDOF * (NF + 1) ] = { 0 };
  globalIndex dofColIndicesFaceVars[ NF ] = { 0 };
  for( localIndex idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[idof] = elemDofNumber + idof;
  }

  // fluxes
  real64 sumLocalMassFluxes[ NC ] = { 0.0 };
  real64 dSumLocalMassFluxes_dElemVars[ NC ][ NDOF * (NF + 1) ] = {{ 0.0 }};
  real64 dSumLocalMassFluxes_dFaceVars[ NC ][ NF ] = {{ 0.0 }};

  // for each element, loop over the one-sided faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {

    localIndex const elemVarsOffset = NDOF*(ifaceLoc+1);

    for( localIndex ip = 0; ip < NP; ++ip )
    {

      for( localIndex ic = 0; ic < NC; ++ic )
      {

        // compute the mass flux at the one-sided face plus its derivatives
        // add the newly computed flux to the sum

        real64 const dt_upwPhaseViscCoef = dt * upwPhaseViscCoef[ifaceLoc][ip][ic];

        sumLocalMassFluxes[ic] = sumLocalMassFluxes[ic] + dt_upwPhaseViscCoef * oneSidedVolFlux[ifaceLoc];
        dSumLocalMassFluxes_dElemVars[ic][0] = dSumLocalMassFluxes_dElemVars[ic][0]
                                               + dt_upwPhaseViscCoef * dOneSidedVolFlux_dPres[ifaceLoc];
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dSumLocalMassFluxes_dElemVars[ic][jc+1] = dSumLocalMassFluxes_dElemVars[ic][jc+1]
                                                    + dt_upwPhaseViscCoef * dOneSidedVolFlux_dCompDens[ifaceLoc][jc];
        }

        dSumLocalMassFluxes_dElemVars[ic][elemVarsOffset] = dSumLocalMassFluxes_dElemVars[ic][elemVarsOffset]
                                                            + dt * dUpwPhaseViscCoef_dPres[ifaceLoc][ip][ic] * oneSidedVolFlux[ifaceLoc];
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dSumLocalMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dSumLocalMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                                   + dt * dUpwPhaseViscCoef_dCompDens[ifaceLoc][ip][ic][jc] * oneSidedVolFlux[ifaceLoc];
        }

        for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
        {
          dSumLocalMassFluxes_dFaceVars[ic][jfaceLoc] = dSumLocalMassFluxes_dFaceVars[ic][jfaceLoc]
                                                        + dt_upwPhaseViscCoef * dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc];
        }
      }
    }

    // collect the relevant dof numbers
    for( localIndex idof = 0; idof < NDOF; ++idof )
    {
      dofColIndicesElemVars[elemVarsOffset+idof] = upwViscDofNumber[ifaceLoc] + idof;
    }
    dofColIndicesFaceVars[ifaceLoc] = faceDofNumber[elemToFaces[ifaceLoc]];
  }

  // we are ready to assemble the local flux and its derivatives
  // no need for atomic adds - each row is assembled by a single thread

  for( localIndex ic = 0; ic < NC; ++ic )
  {
    localIndex const eqnRowLocalIndex = LvArray::integerConversion< localIndex >( elemDofNumber + ic - rankOffset );

    GEOSX_ASSERT_GE( eqnRowLocalIndex, 0 );
    GEOSX_ASSERT_GT( localMatrix.numRows(), eqnRowLocalIndex );

    // residual
    localRhs[eqnRowLocalIndex] = localRhs[eqnRowLocalIndex] + sumLocalMassFluxes[ic];

    // jacobian -- derivative wrt elem centered vars
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                              &dofColIndicesElemVars[0],
                                                              &dSumLocalMassFluxes_dElemVars[0][0] + ic * NDOF * (NF + 1),
                                                              NDOF * (NF + 1) );

    // jacobian -- derivatives wrt face centered vars
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                              &dofColIndicesFaceVars[0],
                                                              &dSumLocalMassFluxes_dFaceVars[0][0] + ic * NF,
                                                              NF );
  }
}

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernelHelper::AssembleConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                                            arrayView1d< integer const > const & faceGhostRank,
                                            arraySlice1d< localIndex const > const & elemToFaces,
                                            globalIndex const elemDofNumber,
                                            globalIndex const rankOffset,
                                            real64 const (&oneSidedVolFlux)[ NF ],
                                            real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                                            real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                                            real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
{
  localIndex constexpr NDOF = NC+1;

  // fluxes
  real64 dFlux_dElemVars[ NDOF ] = { 0.0 };
  real64 dFlux_dFaceVars[ NF ] = { 0.0 };

  // dof numbers
  globalIndex dofColIndicesElemVars[ NDOF ] = { 0 };
  globalIndex dofColIndicesFaceVars[ NF ] = { 0 };
  for( localIndex idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[idof] = elemDofNumber + idof;
  }

  // for each element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    if( faceGhostRank[elemToFaces[ifaceLoc]] >= 0 )
    {
      continue;
    }

    // flux at this face
    real64 const flux = oneSidedVolFlux[ifaceLoc];
    dFlux_dElemVars[0] = dOneSidedVolFlux_dPres[ifaceLoc];
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dFlux_dElemVars[ic+1] = dOneSidedVolFlux_dCompDens[ifaceLoc][ic];
    }

    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      dFlux_dFaceVars[jfaceLoc] = dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc];
      dofColIndicesFaceVars[jfaceLoc] = faceDofNumber[elemToFaces[jfaceLoc]];
    }

    // dof number of this face constraint
    localIndex const eqnLocalRowIndex = LvArray::integerConversion< localIndex >( faceDofNumber[elemToFaces[ifaceLoc]] - rankOffset );

    GEOSX_ASSERT_GE( eqnLocalRowIndex, 0 );
    GEOSX_ASSERT_GT( localMatrix.numRows(), eqnLocalRowIndex );

    // residual
    atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnLocalRowIndex], flux );

    // jacobian -- derivatives wrt elem-centered terms
    localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                      &dofColIndicesElemVars[0],
                                                                      &dFlux_dElemVars[0],
                                                                      NDOF );

    // jacobian -- derivatives wrt face pressure terms
    localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                      &dofColIndicesFaceVars[0],
                                                                      &dFlux_dFaceVars[0],
                                                                      NF );
  }
}


/******************************** AssemblerKernel ********************************/

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernel::Compute( localIndex const er, localIndex const esr, localIndex const ei,
                          SortedArrayView< localIndex const > const & regionFilter,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          arrayView1d< globalIndex const > const & faceDofNumber,
                          arrayView1d< integer const > const & faceGhostRank,
                          arrayView1d< real64 const > const & facePres,
                          arrayView1d< real64 const > const & dFacePres,
                          arrayView1d< real64 const > const & faceGravCoef,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          real64 const & elemPres,
                          real64 const & dElemPres,
                          real64 const & elemGravCoef,
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
                          integer const elemGhostRank,
                          globalIndex const rankOffset,
                          real64 const & dt,
                          arraySlice2d< real64 const > const & transMatrix,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
{
  // one sided flux
  real64 oneSidedVolFlux[ NF ] = { 0.0 };
  real64 dOneSidedVolFlux_dPres[ NF ] = { 0.0 };
  real64 dOneSidedVolFlux_dFacePres[ NF ][ NF ] = {{ 0.0 }};
  real64 dOneSidedVolFlux_dCompDens[ NF ][ NC ] = {{ 0.0 }};

  // upwinded phase viscous transport coefficient
  real64 upwPhaseViscCoef[ NF ][ NP ][ NC ] = {{{ 0.0 }}};
  real64 dUpwPhaseViscCoef_dPres[ NF ][ NP ][ NC ] = {{{ 0.0 }}};
  real64 dUpwPhaseViscCoef_dCompDens[ NF ][ NP ][ NC ][ NC ] = {{{{ 0.0 }}}};
  globalIndex upwViscDofNumber[ NF ] = { 0 };

  /*
   * compute auxiliary quantities at the one sided faces of this element:
   * 1) One-sided volumetric fluxes
   * 2) Upwinded mobilities
   */

  // for each one-sided face of the elem,
  // compute the volumetric flux using transMatrix
  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF, NC, NP >( facePres,
                                                                 dFacePres,
                                                                 faceGravCoef,
                                                                 elemToFaces,
                                                                 elemPres,
                                                                 dElemPres,
                                                                 elemGravCoef,
                                                                 phaseDens[er][esr][ei][0],
                                                                 dPhaseDens_dPres[er][esr][ei][0],
                                                                 dPhaseDens_dCompFrac[er][esr][ei][0],
                                                                 phaseMob[er][esr][ei],
                                                                 dPhaseMob_dPres[er][esr][ei],
                                                                 dPhaseMob_dCompDens[er][esr][ei],
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
    AssemblerKernelHelper::UpdateUpwindedCoefficients< NF, NC, NP >( er, esr, ei,
                                                                     elemRegionList,
                                                                     elemSubRegionList,
                                                                     elemList,
                                                                     regionFilter,
                                                                     elemToFaces,
                                                                     phaseDens,
                                                                     dPhaseDens_dPres,
                                                                     dPhaseDens_dCompFrac,
                                                                     phaseMob,
                                                                     dPhaseMob_dPres,
                                                                     dPhaseMob_dCompDens,
                                                                     dCompFrac_dCompDens,
                                                                     phaseCompFrac,
                                                                     dPhaseCompFrac_dPres,
                                                                     dPhaseCompFrac_dCompFrac,
                                                                     elemDofNumber,
                                                                     oneSidedVolFlux,
                                                                     upwPhaseViscCoef,
                                                                     dUpwPhaseViscCoef_dPres,
                                                                     dUpwPhaseViscCoef_dCompDens,
                                                                     upwViscDofNumber );

    /*
     * perform assembly in this element in two steps:
     * 1) mass conservation equations
     * 2) face constraints
     */

    // use the computed one sided vol fluxes and the upwinded mobilities
    // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
    AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF, NC, NP >( faceDofNumber,
                                                                     elemToFaces,
                                                                     elemDofNumber[er][esr][ei],
                                                                     rankOffset,
                                                                     oneSidedVolFlux,
                                                                     dOneSidedVolFlux_dPres,
                                                                     dOneSidedVolFlux_dFacePres,
                                                                     dOneSidedVolFlux_dCompDens,
                                                                     upwPhaseViscCoef,
                                                                     dUpwPhaseViscCoef_dPres,
                                                                     dUpwPhaseViscCoef_dCompDens,
                                                                     upwViscDofNumber,
                                                                     dt,
                                                                     localMatrix,
                                                                     localRhs );
  }

  // use the computed one sided vol fluxes to assemble the constraints
  // enforcing flux continuity at this element's faces
  AssemblerKernelHelper::AssembleConstraints< NF, NC, NP >( faceDofNumber,
                                                            faceGhostRank,
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

template< localIndex NF, localIndex NC, localIndex NP >
void
FluxKernel::Launch( localIndex er,
                    localIndex esr,
                    CellElementSubRegion const & subRegion,
                    SortedArrayView< localIndex const > const & regionFilter,
                    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                    arrayView2d< localIndex const > const & elemRegionList,
                    arrayView2d< localIndex const > const & elemSubRegionList,
                    arrayView2d< localIndex const > const & elemList,
                    ArrayOfArraysView< localIndex const > const & faceToNodes,
                    arrayView1d< globalIndex const > const & faceDofNumber,
                    arrayView1d< integer const > const & faceGhostRank,
                    arrayView1d< real64 const > const & facePres,
                    arrayView1d< real64 const > const & dFacePres,
                    arrayView1d< real64 const > const & faceGravCoef,
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
                    globalIndex const rankOffset,
                    real64 const lengthTolerance,
                    real64 const dt,
                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                    arrayView1d< real64 > const & localRhs )
{
  // get the cell-centered DOF numbers and ghost rank for the assembly
  arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

  // get the map from elem to faces
  arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

  // get the cell-centered pressures
  arrayView1d< real64 const > const & elemPres  =
    subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dElemPres =
    subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::deltaPressureString );

  // get the element data needed for transmissibility computation
  arrayView2d< real64 const > const & elemCenter =
    subRegion.getReference< array2d< real64 > >( CellBlock::viewKeyStruct::elementCenterString );
  arrayView1d< real64 const > const & elemVolume =
    subRegion.getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );
  arrayView1d< R1Tensor const > const & elemPerm =
    subRegion.getReference< array1d< R1Tensor > >( CompositionalMultiphaseBase::viewKeyStruct::permeabilityString );

  // get the cell-centered depth
  arrayView1d< real64 const > const & elemGravCoef =
    subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::gravityCoefString );

  // assemble the residual and Jacobian element by element
  // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const ei )
  {

    // transmissibility matrix
    stackArray2d< real64, NF *NF > transMatrix( NF, NF );

    real64 const perm[ 3 ] = { elemPerm[ei][0], elemPerm[ei][1], elemPerm[ei][2] };

    // recompute the local transmissibility matrix at each iteration
    // we can decide later to precompute transMatrix if needed
    HybridFVMInnerProduct::QTPFACellInnerProductKernel::Compute< NF >( nodePosition,
                                                                       faceToNodes,
                                                                       elemToFaces[ei],
                                                                       elemCenter[ei],
                                                                       elemVolume[ei],
                                                                       perm,
                                                                       2.0,
                                                                       lengthTolerance,
                                                                       transMatrix );

    // HybridFVMInnerProduct::TPFACellInnerProductKernel::Compute< NF >( nodePosition,
    //                                                                   faceToNodes,
    //                                                                   elemToFaces[ei],
    //                                                                   elemCenter[ei],
    //                                                                   perm,
    //                                                                   lengthTolerance,
    //                                                                   transMatrix );


    // perform flux assembly in this element
    CompositionalMultiphaseHybridFVMKernels::AssemblerKernel::Compute< NF, NC, NP >( er, esr, ei,
                                                                                     regionFilter,
                                                                                     elemRegionList,
                                                                                     elemSubRegionList,
                                                                                     elemList,
                                                                                     faceDofNumber,
                                                                                     faceGhostRank,
                                                                                     facePres,
                                                                                     dFacePres,
                                                                                     faceGravCoef,
                                                                                     elemToFaces[ei],
                                                                                     elemPres[ei],
                                                                                     dElemPres[ei],
                                                                                     elemGravCoef[ei],
                                                                                     phaseDens,
                                                                                     dPhaseDens_dPres,
                                                                                     dPhaseDens_dCompFrac,
                                                                                     phaseMob,
                                                                                     dPhaseMob_dPres,
                                                                                     dPhaseMob_dCompDens,
                                                                                     dCompFrac_dCompDens,
                                                                                     phaseCompFrac,
                                                                                     dPhaseCompFrac_dPres,
                                                                                     dPhaseCompFrac_dCompFrac,
                                                                                     elemDofNumber,
                                                                                     elemGhostRank[ei],
                                                                                     rankOffset,
                                                                                     dt,
                                                                                     transMatrix,
                                                                                     localMatrix,
                                                                                     localRhs );
  } );
}

#define INST_UpwindingHelper( NF, NC, NP ) \
  template \
  void \
  UpwindingHelper::UpwindViscousTerm< NF, NC, NP >( localIndex const er, localIndex const esr, localIndex const ei, \
                                                    localIndex const ifaceLoc, \
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
                                                    real64 ( &upwPhaseViscCoef )[ NF ][ NP ][ NC ], \
                                                    real64 ( &dUpwPhaseViscCoef_dPres )[ NF ][ NP ][ NC ], \
                                                    real64 ( &dUpwPhaseViscCoef_dCompDens )[ NF ][ NP ][ NC ][ NC ], \
                                                    globalIndex ( &upwViscDofNumber )[ NF ] )

INST_UpwindingHelper( 4, 2, 2 );
INST_UpwindingHelper( 4, 3, 2 );
INST_UpwindingHelper( 4, 4, 2 );
INST_UpwindingHelper( 4, 5, 2 );
INST_UpwindingHelper( 4, 2, 3 );
INST_UpwindingHelper( 4, 3, 3 );
INST_UpwindingHelper( 4, 4, 3 );
INST_UpwindingHelper( 4, 5, 3 );
INST_UpwindingHelper( 5, 2, 2 );
INST_UpwindingHelper( 5, 3, 2 );
INST_UpwindingHelper( 5, 4, 2 );
INST_UpwindingHelper( 5, 5, 2 );
INST_UpwindingHelper( 5, 2, 3 );
INST_UpwindingHelper( 5, 3, 3 );
INST_UpwindingHelper( 5, 4, 3 );
INST_UpwindingHelper( 5, 5, 3 );
INST_UpwindingHelper( 6, 2, 2 );
INST_UpwindingHelper( 6, 3, 2 );
INST_UpwindingHelper( 6, 4, 2 );
INST_UpwindingHelper( 6, 5, 2 );
INST_UpwindingHelper( 6, 2, 3 );
INST_UpwindingHelper( 6, 3, 3 );
INST_UpwindingHelper( 6, 4, 3 );
INST_UpwindingHelper( 6, 5, 3 );

#undef INST_UpwindingHelper

#define INST_AssemblerKernelHelper( NF, NC, NP ) \
  template \
  void \
  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF, NC, NP >( arrayView1d< real64 const > const & facePres, \
                                                                 arrayView1d< real64 const > const & dFacePres, \
                                                                 arrayView1d< real64 const > const & faceGravCoef, \
                                                                 arraySlice1d< localIndex const > const & elemToFaces, \
                                                                 real64 const & elemPres, \
                                                                 real64 const & dElemPres, \
                                                                 real64 const & elemGravCoef, \
                                                                 arraySlice1d< real64 const > const & elemPhaseDens, \
                                                                 arraySlice1d< real64 const > const & dElemPhaseDens_dPres, \
                                                                 arraySlice2d< real64 const > const & dElemPhaseDens_dCompFrac, \
                                                                 arraySlice1d< real64 const > const & elemPhaseMob, \
                                                                 arraySlice1d< real64 const > const & dElemPhaseMob_dPres, \
                                                                 arraySlice2d< real64 const > const & dElemPhaseMob_dCompDens, \
                                                                 arraySlice2d< real64 const > const & dElemCompFrac_dCompDens, \
                                                                 arraySlice2d< real64 const > const & transMatrix, \
                                                                 real64 ( &oneSidedVolFlux )[ NF ], \
                                                                 real64 ( &dOneSidedVolFlux_dPres )[ NF ], \
                                                                 real64 ( &dOneSidedVolFlux_dFacePres )[ NF ][ NF ], \
                                                                 real64 ( &dOneSidedVolFlux_dCompDens )[ NF ][ NC ] ); \
  template \
  void \
  AssemblerKernelHelper::UpdateUpwindedCoefficients< NF, NC, NP >( localIndex const er, localIndex const esr, localIndex const ei, \
                                                                   arrayView2d< localIndex const > const & elemRegionList, \
                                                                   arrayView2d< localIndex const > const & elemSubRegionList, \
                                                                   arrayView2d< localIndex const > const & elemList, \
                                                                   SortedArrayView< localIndex const > const & regionFilter, \
                                                                   arraySlice1d< localIndex const > const & elemToFaces, \
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
                                                                   real64 const (&oneSidedVolFlux)[ NF ], \
                                                                   real64 ( &upwPhaseViscCoef )[ NF ][ NP ][ NC ], \
                                                                   real64 ( &dUpwPhaseViscCoef_dPres )[ NF ][ NP ][ NC ], \
                                                                   real64 ( &dUpwPhaseViscCoef_dCompDens )[ NF ][ NP ][ NC ][ NC ], \
                                                                   globalIndex ( &upwViscDofNumber )[ NF ] ); \
  template \
  void \
  AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF, NC, NP >( arrayView1d< globalIndex const > const & faceDofNumber, \
                                                                   arraySlice1d< localIndex const > const & elemToFaces, \
                                                                   globalIndex const elemDofNumber, \
                                                                   globalIndex const rankOffset, \
                                                                   real64 const (&oneSidedVolFlux)[ NF ], \
                                                                   real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                                                   real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                                                   real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                                                   real64 const (&upwPhaseViscCoef)[ NF ][ NP ][ NC ], \
                                                                   real64 const (&dUpwPhaseViscCoef_dPres)[ NF ][ NP ][ NC ], \
                                                                   real64 const (&dUpwPhaseViscCoef_dCompDens)[ NF ][ NP ][ NC ][ NC ], \
                                                                   globalIndex const (&upwViscDofNumber)[ NF ], \
                                                                   real64 const & dt, \
                                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                                                   arrayView1d< real64 > const & localRhs ); \
  template \
  void \
  AssemblerKernelHelper::AssembleConstraints< NF, NC, NP >( arrayView1d< globalIndex const > const & faceDofNumber, \
                                                            arrayView1d< integer const > const & faceGhostRank, \
                                                            arraySlice1d< localIndex const > const & elemToFaces, \
                                                            globalIndex const elemDofNumber, \
                                                            globalIndex const rankOffset, \
                                                            real64 const (&oneSidedVolFlux)[ NF ], \
                                                            real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                                            real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                                            real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                                            CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                                            arrayView1d< real64 > const & localRhs )

INST_AssemblerKernelHelper( 4, 2, 2 );
INST_AssemblerKernelHelper( 4, 3, 2 );
INST_AssemblerKernelHelper( 4, 4, 2 );
INST_AssemblerKernelHelper( 4, 5, 2 );
INST_AssemblerKernelHelper( 4, 2, 3 );
INST_AssemblerKernelHelper( 4, 3, 3 );
INST_AssemblerKernelHelper( 4, 4, 3 );
INST_AssemblerKernelHelper( 4, 5, 3 );
INST_AssemblerKernelHelper( 5, 2, 2 );
INST_AssemblerKernelHelper( 5, 3, 2 );
INST_AssemblerKernelHelper( 5, 4, 2 );
INST_AssemblerKernelHelper( 5, 5, 2 );
INST_AssemblerKernelHelper( 5, 2, 3 );
INST_AssemblerKernelHelper( 5, 3, 3 );
INST_AssemblerKernelHelper( 5, 4, 3 );
INST_AssemblerKernelHelper( 5, 5, 3 );
INST_AssemblerKernelHelper( 6, 2, 2 );
INST_AssemblerKernelHelper( 6, 3, 2 );
INST_AssemblerKernelHelper( 6, 4, 2 );
INST_AssemblerKernelHelper( 6, 5, 2 );
INST_AssemblerKernelHelper( 6, 2, 3 );
INST_AssemblerKernelHelper( 6, 3, 3 );
INST_AssemblerKernelHelper( 6, 4, 3 );
INST_AssemblerKernelHelper( 6, 5, 3 );

#undef INST_AssemblerKernelHelper

#define INST_FluxKernel( NF, NC, NP ) \
  template \
  void \
  FluxKernel::Launch< NF, NC, NP >( localIndex er, \
                                    localIndex esr, \
                                    CellElementSubRegion const & subRegion, \
                                    SortedArrayView< localIndex const > const & regionFilter, \
                                    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
                                    arrayView2d< localIndex const > const & elemRegionList, \
                                    arrayView2d< localIndex const > const & elemSubRegionList, \
                                    arrayView2d< localIndex const > const & elemList, \
                                    ArrayOfArraysView< localIndex const > const & faceToNodes, \
                                    arrayView1d< globalIndex const > const & faceDofNumber, \
                                    arrayView1d< integer const > const & faceGhostRank, \
                                    arrayView1d< real64 const > const & facePres, \
                                    arrayView1d< real64 const > const & dFacePres, \
                                    arrayView1d< real64 const > const & faceGravCoef, \
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
                                    globalIndex const rankOffset, \
                                    real64 const lengthTolerance, \
                                    real64 const dt, \
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                    arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 4, 2, 2 );
INST_FluxKernel( 4, 3, 2 );
INST_FluxKernel( 4, 4, 2 );
INST_FluxKernel( 4, 5, 2 );
INST_FluxKernel( 4, 2, 3 );
INST_FluxKernel( 4, 3, 3 );
INST_FluxKernel( 4, 4, 3 );
INST_FluxKernel( 4, 5, 3 );
INST_FluxKernel( 5, 2, 2 );
INST_FluxKernel( 5, 3, 2 );
INST_FluxKernel( 5, 4, 2 );
INST_FluxKernel( 5, 5, 2 );
INST_FluxKernel( 5, 2, 3 );
INST_FluxKernel( 5, 3, 3 );
INST_FluxKernel( 5, 4, 3 );
INST_FluxKernel( 5, 5, 3 );
INST_FluxKernel( 6, 2, 2 );
INST_FluxKernel( 6, 3, 2 );
INST_FluxKernel( 6, 4, 2 );
INST_FluxKernel( 6, 5, 2 );
INST_FluxKernel( 6, 2, 3 );
INST_FluxKernel( 6, 3, 3 );
INST_FluxKernel( 6, 4, 3 );
INST_FluxKernel( 6, 5, 3 );

#undef INST_FluxKernel



} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx
