/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TwoPhaseHybridFVMKernels.cpp
 */

#include "TwoPhaseHybridFVMKernels.hpp"
#include "TwoPhaseBase.hpp"

namespace geosx
{

namespace TwoPhaseHybridFVMKernels
{

/******************************** FluxKernelHelper ********************************/


void FluxKernelHelper::ComputeOneSidedVolFluxes( arrayView2d< real64 const > const & facePotential,
                                                 arrayView2d< real64 const > const & dFacePotential,
                                                 arrayView1d< real64 const > const & faceGravCoef,
                                                 arraySlice1d< localIndex const > const elemToFaces,
                                                 real64 const & elemPres,
                                                 real64 const & dElemPres,
                                                 real64 const & elemGravCoef,
                                                 arraySlice1d< real64 const > const elemDens,
                                                 arraySlice1d< real64 const > const dElemDens_dPres,
                                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                               *HybridFVMInnerProduct::MAX_NUM_FACES > const & transMatrix,
                                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & volFlux,
                                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & dVolFlux_dPres,
                                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & dVolFlux_dSat,
                                                 stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *HybridFVMInnerProduct::MAX_NUM_FACES
                                                               *TwoPhaseBase::NUM_PHASES > & dVolFlux_dFacePres )
{
  localIndex constexpr numPhases   = TwoPhaseBase::NUM_PHASES;
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;

  localIndex const numFacesInElem = elemToFaces.size();

  volFlux = 0;
  dVolFlux_dPres = 0;
  dVolFlux_dSat = 0;
  dVolFlux_dFacePres = 0;

  // this is going to store \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
  stackArray2d< real64, maxNumFaces *numPhases > potDif( numFacesInElem, numPhases );
  stackArray2d< real64, maxNumFaces *numPhases > dPotDif_dPres( numFacesInElem, numPhases );
  stackArray2d< real64, maxNumFaces *numPhases > dPotDif_dSat( numFacesInElem, numPhases );
  stackArray2d< real64, maxNumFaces *numPhases > dPotDif_dFacePres( numFacesInElem, numPhases );
  potDif = 0;
  dPotDif_dPres = 0;
  dPotDif_dSat = 0;
  dPotDif_dFacePres = 0;

  // 1) precompute the potential difference at each one-sided face
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {
    // depth difference between element center and face center
    real64 const ccGravCoef = elemGravCoef;
    real64 const fGravCoef  = faceGravCoef[elemToFaces[ifaceLoc]];
    real64 const depthDif   = ccGravCoef - fGravCoef;

    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      // cell centered phase pressure and gravity terms
      real64 const phasePres = elemPres + dElemPres;
      real64 const dPhasePres_dPres = 1;
      real64 const dPhasePres_dSat  = 0; // for capillary pressure
      real64 const phaseGravTerm    = elemDens[ip] * depthDif;
      real64 const dPhaseGravTerm_dPres = dElemDens_dPres[ip] * depthDif;

      // cell centered phase potential
      real64 const phasePotential = phasePres - phaseGravTerm;
      real64 const dPhasePotential_dPres = dPhasePres_dPres - dPhaseGravTerm_dPres;
      real64 const dPhasePotential_dSat  = dPhasePres_dSat;

      // face phase potential
      real64 const facePhasePotential = facePotential[elemToFaces[ifaceLoc]][ip]
                                        + dFacePotential[elemToFaces[ifaceLoc]][ip];

      // phase potential difference
      potDif[ifaceLoc][ip] = phasePotential - facePhasePotential;
      dPotDif_dPres[ifaceLoc][ip] = dPhasePotential_dPres;
      dPotDif_dSat[ifaceLoc][ip]  = dPhasePotential_dSat;
      dPotDif_dFacePres[ifaceLoc][ip] = -1;
    }
  }

  // 2) multiply the potential difference by the transmissibility
  //    to obtain the total volumetrix flux
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {

    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
    for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
    {

      for( localIndex ip = 0; ip < numPhases; ++ip )
      {
        // this is going to store T (\nabla p_{\ell} - \rho_p_{\ell} g \nabla d)
        volFlux[ifaceLoc][ip]        += transMatrix[ifaceLoc][jfaceLoc] * potDif[jfaceLoc][ip];
        dVolFlux_dPres[ifaceLoc][ip] += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dPres[jfaceLoc][ip];
        dVolFlux_dSat[ifaceLoc][ip]  += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dSat[jfaceLoc][ip];
        dVolFlux_dFacePres[ifaceLoc][jfaceLoc][ip] += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dFacePres[jfaceLoc][ip];
      }
    }
  }
}


void FluxKernelHelper::UpdateUpwindedCoefficients( array2d< localIndex > const & elemRegionList,
                                                   array2d< localIndex > const & elemSubRegionList,
                                                   array2d< localIndex > const & elemList,
                                                   SortedArray< localIndex > const & regionFilter,
                                                   arraySlice1d< localIndex const > const elemToFaces,
                                                   ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & mob,
                                                   ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & dMob_dPres,
                                                   ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & dMob_dSat,
                                                   ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst const & elemDofNumber,
                                                   stackArray1d< localIndex, 3 > const & elemIds,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & volFlux,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & upwMobility,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & dUpwMobility_dPres,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & dUpwMobility_dSat,
                                                   stackArray2d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > & upwDofNumber,
                                                   stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > & neighborDofNumber )
{

  localIndex constexpr numPhases   = TwoPhaseBase::NUM_PHASES;
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;

  localIndex const numFacesInElem = elemToFaces.size();

  stackArray2d< localIndex, 3 *maxNumFaces > neighborIds( numFacesInElem, 3 );
  neighborIds = -1;

  localIndex const er  = elemIds[0];
  localIndex const esr = elemIds[1];
  localIndex const ei  = elemIds[2];

  // On the fly, find the neighbors at all the interfaces of this element
  // We can decide later to precompute and save these indices
  FluxKernelHelper::FindNeighborsInTarget( elemRegionList,
                                           elemSubRegionList,
                                           elemList,
                                           regionFilter,
                                           elemToFaces,
                                           elemIds,
                                           elemDofNumber,
                                           neighborIds,
                                           neighborDofNumber );

  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {

    // we initialize these upw quantities with the values of the local elem
    // they will be updated below if there is a neighbor
    bool const foundNeighborInTarget = (neighborIds[ifaceLoc][0] != -1 &&
                                        neighborIds[ifaceLoc][1] != -1 &&
                                        neighborIds[ifaceLoc][2] != -1);

    for( localIndex ip = 0; ip < numPhases; ++ip )
    {

      // if true, the neighbor element is upwind
      // if false, the local element is uwwind
      localIndex const neighborIsUpw = foundNeighborInTarget && volFlux[ifaceLoc][ip] < 0;

      // we evaluate the mobilities with the upwind element wrt the volumetrix flux
      localIndex const erUpw  = neighborIsUpw * neighborIds[ifaceLoc][0] + (1-neighborIsUpw) * er;
      localIndex const esrUpw = neighborIsUpw * neighborIds[ifaceLoc][1] + (1-neighborIsUpw) * esr;
      localIndex const eiUpw  = neighborIsUpw * neighborIds[ifaceLoc][2] + (1-neighborIsUpw) * ei;

      // get the dof number in the upwind element
      globalIndex const dofNumber = neighborIsUpw * neighborDofNumber[ifaceLoc]
                                    + (1-neighborIsUpw) * elemDofNumber[er][esr][ei];

      upwMobility[ifaceLoc][ip]        = mob[erUpw][esrUpw][eiUpw][ip];
      dUpwMobility_dPres[ifaceLoc][ip] = dMob_dPres[erUpw][esrUpw][eiUpw][ip];
      dUpwMobility_dSat[ifaceLoc][ip]  = dMob_dSat[erUpw][esrUpw][eiUpw][ip];
      upwDofNumber[ifaceLoc][ip]       = dofNumber;
    }
  }
}

void FluxKernelHelper::AssembleOneSidedMassFluxes( real64 const & dt,
                                                   arrayView1d< globalIndex const > const & faceDofNumber,
                                                   arraySlice1d< localIndex const > const elemToFaces,
                                                   globalIndex const elemDofNumber,
                                                   array1d< localIndex > const & phaseToRow,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & volFlux,
                                                   stackArray2d< real64,
                                                                 HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & dVolFlux_dPres,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & dVolFlux_dSat,
                                                   stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *HybridFVMInnerProduct::MAX_NUM_FACES
                                                                 *TwoPhaseBase::NUM_PHASES > const & dVolFlux_dFacePres,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & upwMobility,
                                                   stackArray2d< real64,
                                                                 HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & dUpwMobility_dPres,
                                                   stackArray2d< real64,
                                                                 HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & dUpwMobility_dSat,
                                                   stackArray2d< globalIndex,
                                                                 HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & upwDofNumber,
                                                   stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > const & neighborDofNumber,
                                                   ParallelMatrix * const matrix,
                                                   ParallelVector * const rhs )
{
  localIndex constexpr numDof      = TwoPhaseBase::NUM_DOF;
  localIndex constexpr numPhases   = TwoPhaseBase::NUM_PHASES;
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;

  localIndex constexpr dPres = TwoPhaseBase::ColOffset::DPRES;
  localIndex constexpr dSat = TwoPhaseBase::ColOffset::DSAT;

  localIndex const numFacesInElem = elemToFaces.size();

  // dof numbers
  stackArray1d< globalIndex, numPhases >               eqnRowIndices( numPhases );
  stackArray1d< globalIndex, numDof *(1+maxNumFaces) > elemDofColIndices( numDof * (1+numFacesInElem) );
  stackArray1d< globalIndex, numPhases * maxNumFaces > faceDofColIndices( numPhases * numFacesInElem );
  eqnRowIndices[0] = elemDofNumber + phaseToRow[0];
  eqnRowIndices[1] = elemDofNumber + phaseToRow[1];
  elemDofColIndices[dPres] = elemDofNumber + dPres;
  elemDofColIndices[dSat]  = elemDofNumber + dSat;

  // fluxes
  stackArray1d< real64, numPhases >                          sumLocalMassFluxes( numPhases );
  stackArray2d< real64, numPhases *numDof *(1+maxNumFaces) > dSumLocalMassFluxes_dElemVars( numPhases, numDof * (1+numFacesInElem) );
  stackArray2d< real64, numPhases * numPhases *maxNumFaces > dSumLocalMassFluxes_dFaceVars( numPhases, numPhases * numFacesInElem );
  sumLocalMassFluxes = 0;
  dSumLocalMassFluxes_dElemVars = 0;
  dSumLocalMassFluxes_dFaceVars = 0;

  // for each element, loop over the one-sided faces
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {
    localIndex const elemVarsOffset = numDof*(ifaceLoc+1);

    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      // 1) residual
      sumLocalMassFluxes[ip] += dt * upwMobility[ifaceLoc][ip] * volFlux[ifaceLoc][ip];

      // 2) derivatives wrt the elem centered vars of the local elem
      dSumLocalMassFluxes_dElemVars[ip][dPres] += dt * upwMobility[ifaceLoc][ip] * dVolFlux_dPres[ifaceLoc][ip];
      dSumLocalMassFluxes_dElemVars[ip][dSat]  += dt * upwMobility[ifaceLoc][ip] * dVolFlux_dSat[ifaceLoc][ip];

      // 3) derivatives of the upwinded mobilities
      if( upwDofNumber[ifaceLoc][ip] == elemDofNumber )
      {
        dSumLocalMassFluxes_dElemVars[ip][dPres] += dt * dUpwMobility_dPres[ifaceLoc][ip] * volFlux[ifaceLoc][ip];
        dSumLocalMassFluxes_dElemVars[ip][dSat]  += dt * dUpwMobility_dSat[ifaceLoc][ip] * volFlux[ifaceLoc][ip];
      }
      else
      {
        dSumLocalMassFluxes_dElemVars[ip][elemVarsOffset+dPres] = dt * dUpwMobility_dPres[ifaceLoc][ip] * volFlux[ifaceLoc][ip];
        dSumLocalMassFluxes_dElemVars[ip][elemVarsOffset+dSat]  = dt * dUpwMobility_dSat[ifaceLoc][ip] * volFlux[ifaceLoc][ip];
      }

      // 4) derivatives wrt the face centered var
      for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
      {
        localIndex const faceVarsOffset = numPhases*jfaceLoc;
        dSumLocalMassFluxes_dFaceVars[ip][faceVarsOffset + ip] += dt * upwMobility[ifaceLoc][ip] * dVolFlux_dFacePres[ifaceLoc][jfaceLoc][ip];
      }
    }

    // collect the relevant dof numbers
    elemDofColIndices[elemVarsOffset+dPres] = neighborDofNumber[ifaceLoc] + dPres;
    elemDofColIndices[elemVarsOffset+dSat] = neighborDofNumber[ifaceLoc] + dSat;

    localIndex const faceVarsOffset = numPhases*ifaceLoc;
    faceDofColIndices[faceVarsOffset] = faceDofNumber[elemToFaces[ifaceLoc]];
    faceDofColIndices[faceVarsOffset+1] = faceDofNumber[elemToFaces[ifaceLoc]] + 1;
  }

  // we are ready to assemble the local flux and its derivatives

  // residual
  rhs->add( eqnRowIndices,
            sumLocalMassFluxes );

  // jacobian -- derivative wrt elem centered vars
  matrix->add( eqnRowIndices,
               elemDofColIndices,
               dSumLocalMassFluxes_dElemVars );

  // jacobian -- derivatives wrt face centered vars
  matrix->add( eqnRowIndices,
               faceDofColIndices,
               dSumLocalMassFluxes_dFaceVars );
}

void FluxKernelHelper::AssembleConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                                            arraySlice1d< localIndex const > const elemToFaces,
                                            globalIndex const elemDofNumber,
                                            stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & volFlux,
                                            stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & dVolFlux_dPres,
                                            stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *TwoPhaseBase::NUM_PHASES > const & dVolFlux_dSat,
                                            stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES *HybridFVMInnerProduct::MAX_NUM_FACES
                                                          *TwoPhaseBase::NUM_PHASES > const & dVolFlux_dFacePres,
                                            ParallelMatrix * const matrix,
                                            ParallelVector * const rhs )
{
  localIndex constexpr numDof      = TwoPhaseBase::NUM_DOF;
  localIndex constexpr numPhases   = TwoPhaseBase::NUM_PHASES;
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;

  localIndex const numFacesInElem = elemToFaces.size();

  localIndex const dPres = TwoPhaseBase::ColOffset::DPRES;
  localIndex const dSat = TwoPhaseBase::ColOffset::DSAT;

  // dof numbers
  stackArray1d< globalIndex, numPhases > eqnRowIndices( numPhases );
  stackArray1d< globalIndex, numPhases *numDof > elemDofColIndices( numDof );
  stackArray1d< globalIndex, numPhases *maxNumFaces > faceDofColIndices( numPhases*numFacesInElem );

  elemDofColIndices[dPres] = elemDofNumber + dPres;
  elemDofColIndices[dSat] = elemDofNumber + dSat;

  // fluxes
  stackArray1d< real64, numPhases > flux( numPhases );
  stackArray2d< real64, numPhases *numDof > dFlux_dElemVars( numPhases, numDof );
  stackArray2d< real64, numPhases *numPhases *maxNumFaces > dFlux_dFaceVars( numPhases, numPhases*numFacesInElem );

  // for each element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {

    dFlux_dFaceVars = 0.;
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {

      // dof numbers
      eqnRowIndices[ip] = faceDofNumber[elemToFaces[ifaceLoc]] + ip;

      // fluxes
      flux[ip] = volFlux[ifaceLoc][ip];
      dFlux_dElemVars[ip][dPres] = dVolFlux_dPres[ifaceLoc][ip];
      dFlux_dElemVars[ip][dSat] = dVolFlux_dSat[ifaceLoc][ip];

      for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
      {
        localIndex const faceVarsOffset = numPhases*jfaceLoc;
        faceDofColIndices[faceVarsOffset+ip] = faceDofNumber[elemToFaces[jfaceLoc]] + ip;
        dFlux_dFaceVars[ip][faceVarsOffset+ip] = dVolFlux_dFacePres[ifaceLoc][jfaceLoc][ip];
      }
    }

    // residual
    rhs->add( eqnRowIndices,
              flux );

    // jacobian -- derivative wrt local elem centered vars
    matrix->add( eqnRowIndices,
                 elemDofColIndices,
                 dFlux_dElemVars );

    // jacobian -- derivatives wrt face pressure terms
    matrix->add( eqnRowIndices,
                 faceDofColIndices,
                 dFlux_dFaceVars );

  }
}


void FluxKernelHelper::FindNeighborsInTarget( array2d< localIndex > const & elemRegionList,
                                              array2d< localIndex > const & elemSubRegionList,
                                              array2d< localIndex > const & elemList,
                                              SortedArray< localIndex > const & regionFilter,
                                              arraySlice1d< localIndex const > const elemToFaces,
                                              stackArray1d< localIndex, 3 > const & elemIds,
                                              ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst const & elemDofNumber,
                                              stackArray2d< localIndex, 3*HybridFVMInnerProduct::MAX_NUM_FACES > & neighborIds,
                                              stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > & neighborDofNumber )
{
  localIndex const numFacesInElem = elemToFaces.size();

  localIndex const er  = elemIds[0];
  localIndex const esr = elemIds[1];
  localIndex const ei  = elemIds[2];

  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {

    localIndex const faceId = elemToFaces[ifaceLoc];
    bool foundNeighborInTarget = false;
    neighborDofNumber[ifaceLoc] = elemDofNumber[er][esr][ei];

    // the face has at most two adjacent elements
    // one of these two elements is the current element indexed by er, esr, ei
    // but here we are interested in the indices of the other element
    // this other element is "the neighbor" for this one-sided face
    for( localIndex k=0; k < elemRegionList.size( 1 ) && !foundNeighborInTarget; ++k )
    {

      // this element is not the current element
      // we have found the neighbor or we are at the boundary
      if( elemRegionList[faceId][k]    != er  ||
          elemSubRegionList[faceId][k] != esr ||
          elemList[faceId][k]          != ei )
      {

        bool const onBoundary = (elemRegionList[faceId][k]    == -1 ||
                                 elemSubRegionList[faceId][k] == -1 ||
                                 elemList[faceId][k]          == -1);
        bool const neighborInTarget = regionFilter.contains( elemRegionList[faceId][k] );

        // if not on boundary, save the mobility and the upwDofNumber
        if( !onBoundary && neighborInTarget )
        {
          foundNeighborInTarget = true;
          localIndex const erNeighbor  = elemRegionList[faceId][k];
          localIndex const esrNeighbor = elemSubRegionList[faceId][k];
          localIndex const eiNeighbor  = elemList[faceId][k];

          neighborIds[ifaceLoc][0] = erNeighbor;
          neighborIds[ifaceLoc][1] = esrNeighbor;
          neighborIds[ifaceLoc][2] = eiNeighbor;

          // always save the indices of the neighbor for the grav flux
          neighborDofNumber[ifaceLoc] = elemDofNumber[erNeighbor][esrNeighbor][eiNeighbor];
        }
      }
    }
  }
}


} // namespace TwoPhaseHybridFVMKernels

} // namespace geosx
