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
                          arrayView1d< real64 const > const & mimFaceGravCoef,
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
                          arraySlice2d< real64 const > const & transMatrixGrav,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
{
  // one sided flux
  real64 oneSidedVolFlux[ NF ] = { 0.0 };
  real64 dOneSidedVolFlux_dPres[ NF ] = { 0.0 };
  real64 dOneSidedVolFlux_dFacePres[ NF ][ NF ] = {{ 0.0 }};
  real64 dOneSidedVolFlux_dCompDens[ NF ][ NC ] = {{ 0.0 }};

  localIndex const localIds[3] = { er, esr, ei };

  /*
   * compute auxiliary quantities at the one sided faces of this element:
   * 1) One-sided volumetric fluxes
   * 2) Upwinded mobilities
   */

  // for each one-sided face of the elem,
  // compute the volumetric flux using transMatrix
  AssemblerKernelHelper::ApplyGradient< NF, NC, NP >( facePres,
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
    /*
     * perform assembly in this element in two steps:
     * 1) mass conservation equations
     * 2) face constraints
     */

    // use the computed one sided vol fluxes and the upwinded mobilities
    // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
    AssemblerKernelHelper::AssembleFluxDivergence< NF, NC, NP >( localIds,
                                                                 rankOffset,
                                                                 elemRegionList,
                                                                 elemSubRegionList,
                                                                 elemList,
                                                                 regionFilter,
                                                                 faceDofNumber,
                                                                 mimFaceGravCoef,
                                                                 elemToFaces,
                                                                 elemGravCoef,
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
  AssemblerKernelHelper::AssembleFaceConstraints< NF, NC, NP >( faceDofNumber,
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
                    arrayView1d< real64 const > const & mimFaceGravCoef,
                    arrayView1d< real64 const > const & transMultiplier,
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
    stackArray2d< real64, NF *NF > transMatrixGrav( NF, NF );

    real64 const perm[ 3 ] = { elemPerm[ei][0], elemPerm[ei][1], elemPerm[ei][2] };

    // recompute the local transmissibility matrix at each iteration
    // we can decide later to precompute transMatrix if needed
    HybridFVMInnerProduct::QTPFACellInnerProductKernel::Compute< NF >( nodePosition,
                                                                       transMultiplier,
                                                                       faceToNodes,
                                                                       elemToFaces[ei],
                                                                       elemCenter[ei],
                                                                       elemVolume[ei],
                                                                       perm,
                                                                       1.5,
                                                                       lengthTolerance,
                                                                       transMatrix );

    HybridFVMInnerProduct::TPFACellInnerProductKernel::Compute< NF >( nodePosition,
                                                                      transMultiplier,
                                                                      faceToNodes,
                                                                      elemToFaces[ei],
                                                                      elemCenter[ei],
                                                                      perm,
                                                                      lengthTolerance,
                                                                      transMatrixGrav );


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
                                                                                     mimFaceGravCoef,
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
                                                                                     transMatrixGrav,
                                                                                     localMatrix,
                                                                                     localRhs );
  } );
}

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
                                    arrayView1d< real64 const > const & mimFaceGravCoef, \
                                    arrayView1d< real64 const > const & transMultiplier, \
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
