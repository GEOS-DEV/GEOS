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
 * @file SinglePhaseHybridFVMKernels.cpp
 */

#include "SinglePhaseHybridFVMKernels.hpp"

namespace geosx
{

namespace SinglePhaseHybridFVMKernels
{

/******************************** FluxKernelHelper ********************************/

void FluxKernelHelper::ComputeOneSidedVolFluxes( arrayView1d< real64 const > const & facePres,
                                                 arrayView1d< real64 const > const & dFacePres,
                                                 arrayView1d< real64 const > const & faceGravCoef,
                                                 arraySlice1d< localIndex const > const elemToFaces,
                                                 real64 const & elemPres,
                                                 real64 const & dElemPres,
                                                 real64 const & elemGravCoef,
                                                 real64 const & elemDens,
                                                 real64 const & dElemDens_dp,
                                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                               *HybridFVMInnerProduct::MAX_NUM_FACES > const & transMatrix,
                                                 stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > & oneSidedVolFlux,
                                                 stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > & dOneSidedVolFlux_dp,
                                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                               *HybridFVMInnerProduct::MAX_NUM_FACES > & dOneSidedVolFlux_dfp )
{
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;

  localIndex const numFacesInElem = elemToFaces.size();

  oneSidedVolFlux      = 0;
  dOneSidedVolFlux_dp  = 0;
  dOneSidedVolFlux_dfp = 0;

  stackArray1d< real64, maxNumFaces > potDif( numFacesInElem );
  stackArray1d< real64, maxNumFaces > dPotDif_dp( numFacesInElem );
  stackArray1d< real64, maxNumFaces > dPotDif_dfp( numFacesInElem );
  potDif      = 0;
  dPotDif_dp  = 0;
  dPotDif_dfp = 0;

  // 1) precompute the potential difference at each one-sided face
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {

    // 1) compute the potential diff between the cell center and the face center
    real64 const ccPres = elemPres + dElemPres;
    real64 const fPres  = facePres[elemToFaces[ifaceLoc]] + dFacePres[elemToFaces[ifaceLoc]];

    real64 const ccGravCoef = elemGravCoef;
    real64 const fGravCoef  = faceGravCoef[elemToFaces[ifaceLoc]];

    real64 const ccDens     = elemDens;
    real64 const dCcDens_dp = dElemDens_dp;
    // no density evaluated at the face center

    // pressure difference
    real64 const presDif      = ccPres - fPres;
    real64 const dPresDif_dp  =  1;
    real64 const dPresDif_dfp = -1;

    // gravity term
    real64 const gravCoefDif  = ccGravCoef - fGravCoef;
    real64 const gravTerm     = ccDens     * gravCoefDif;
    real64 const dGravTerm_dp = dCcDens_dp * gravCoefDif;

    // potential difference
    potDif[ifaceLoc]      = presDif     - gravTerm;
    dPotDif_dp[ifaceLoc]  = dPresDif_dp - dGravTerm_dp;
    dPotDif_dfp[ifaceLoc] = dPresDif_dfp;

  }

  // 2) multiply the potential difference by the transmissibility
  //    to obtain the total volumetrix flux
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {
    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
    for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
    {
      // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
      oneSidedVolFlux[ifaceLoc]                += transMatrix[ifaceLoc][jfaceLoc] * potDif[jfaceLoc];
      dOneSidedVolFlux_dp[ifaceLoc]            += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dp[jfaceLoc];
      dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc] += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dfp[jfaceLoc];
    }
  }
}

void FluxKernelHelper::UpdateUpwindedCoefficients( arrayView2d< localIndex const > const & elemRegionList,
                                                   arrayView2d< localIndex const > const & elemSubRegionList,
                                                   arrayView2d< localIndex const > const & elemList,
                                                   SortedArrayView< localIndex const > const & regionFilter,
                                                   arraySlice1d< localIndex const > const elemToFaces,
                                                   ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >::ViewTypeConst const & mob,
                                                   ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >::ViewTypeConst const & dMob_dp,
                                                   ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst const & elemDofNumber,
                                                   localIndex const er,
                                                   localIndex const esr,
                                                   localIndex const ei,
                                                   stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > const & oneSidedVolFlux,
                                                   stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > & upwMobility,
                                                   stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > & dUpwMobility_dp,
                                                   stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > & upwDofNumber )
{
  localIndex const numFacesInElem = elemToFaces.size();

  // for this element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {

    // we initialize these upw quantities with the values of the local elem
    upwMobility[ifaceLoc]     = mob[er][esr][ei];
    dUpwMobility_dp[ifaceLoc] = dMob_dp[er][esr][ei];
    upwDofNumber[ifaceLoc]    = elemDofNumber[er][esr][ei];

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

          // if not on boundary, save the mobility and the upwDofNumber
          if( !onBoundary && neighborInTarget )
          {
            upwMobility[ifaceLoc]     = mob[erNeighbor][esrNeighbor][eiNeighbor];
            dUpwMobility_dp[ifaceLoc] = dMob_dp[erNeighbor][esrNeighbor][eiNeighbor];
            upwDofNumber[ifaceLoc]    = elemDofNumber[erNeighbor][esrNeighbor][eiNeighbor];
          }
          // if the face is on the boundary, use the properties of the local elem
        }
      }
    }
  }
}


void FluxKernelHelper::AssembleOneSidedMassFluxes( real64 const & dt,
                                                   arrayView1d< globalIndex const > const & faceDofNumber,
                                                   arraySlice1d< localIndex const > const elemToFaces,
                                                   globalIndex const elemDofNumber,
                                                   stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > const & oneSidedVolFlux,
                                                   stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > const & dOneSidedVolFlux_dp,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                 *HybridFVMInnerProduct::MAX_NUM_FACES > const & dOneSidedVolFlux_dfp,
                                                   stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > const & upwMobility,
                                                   stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > const & dUpwMobility_dp,
                                                   stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > const & upwDofNumber,
                                                   ParallelMatrix * const matrix,
                                                   ParallelVector * const rhs )
{
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;

  localIndex const numFacesInElem = elemToFaces.size();

  // fluxes
  real64 sumLocalMassFluxes     = 0;
  stackArray1d< real64, 1+maxNumFaces > dSumLocalMassFluxes_dElemVars( 1+numFacesInElem );
  stackArray1d< real64, maxNumFaces >   dSumLocalMassFluxes_dFaceVars( numFacesInElem );
  sumLocalMassFluxes = 0;
  dSumLocalMassFluxes_dElemVars = 0;
  dSumLocalMassFluxes_dFaceVars = 0;

  // dof numbers
  globalIndex const eqnRowIndex = elemDofNumber;
  stackArray1d< globalIndex, 1+maxNumFaces > elemDofColIndices( 1+numFacesInElem );
  stackArray1d< globalIndex, maxNumFaces >   faceDofColIndices( numFacesInElem );
  elemDofColIndices[0] = elemDofNumber;

  // for each element, loop over the one-sided faces
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {

    // compute the mass flux at the one-sided face plus its derivatives
    // add the newly computed flux to the sum
    sumLocalMassFluxes                        += dt * upwMobility[ifaceLoc] * oneSidedVolFlux[ifaceLoc];
    dSumLocalMassFluxes_dElemVars[0]          += dt * upwMobility[ifaceLoc] * dOneSidedVolFlux_dp[ifaceLoc];
    dSumLocalMassFluxes_dElemVars[ifaceLoc+1]  = dt * dUpwMobility_dp[ifaceLoc] * oneSidedVolFlux[ifaceLoc];
    for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
    {
      dSumLocalMassFluxes_dFaceVars[jfaceLoc] += dt * upwMobility[ifaceLoc] * dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc];
    }

    // collect the relevant dof numbers
    elemDofColIndices[ifaceLoc+1] = upwDofNumber[ifaceLoc]; // if upwDofNumber == elemDofNumber, the derivative is zero
    faceDofColIndices[ifaceLoc] = faceDofNumber[elemToFaces[ifaceLoc]];
  }

  // we are ready to assemble the local flux and its derivatives

  // residual
  rhs->add( eqnRowIndex,
            sumLocalMassFluxes );

  // jacobian -- derivative wrt elem centered vars
  matrix->add( eqnRowIndex,
               elemDofColIndices,
               dSumLocalMassFluxes_dElemVars );

  // jacobian -- derivatives wrt face centered vars
  matrix->add( eqnRowIndex,
               faceDofColIndices,
               dSumLocalMassFluxes_dFaceVars );

}


void FluxKernelHelper::AssembleConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                                            arraySlice1d< localIndex const > const elemToFaces,
                                            globalIndex const elemDofNumber,
                                            stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > const & oneSidedVolFlux,
                                            stackArray1d< real64, HybridFVMInnerProduct::MAX_NUM_FACES > const & dOneSidedVolFlux_dp,
                                            stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                          *HybridFVMInnerProduct::MAX_NUM_FACES > const & dOneSidedVolFlux_dfp,
                                            ParallelMatrix * const matrix,
                                            ParallelVector * const rhs )
{
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;

  localIndex const numFacesInElem = elemToFaces.size();

  // fluxes
  stackArray1d< real64, maxNumFaces > dFlux_dfp( numFacesInElem );

  // dof numbers
  stackArray1d< globalIndex, maxNumFaces > dofColIndicesFacePres( numFacesInElem );
  globalIndex const dofColIndexElemPres = elemDofNumber;


  // for each element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {
    // flux at this face
    real64 const flux      = oneSidedVolFlux[ifaceLoc];
    real64 const dFlux_dp  = dOneSidedVolFlux_dp[ifaceLoc];

    // dof number of this face constraint
    globalIndex const eqnRowIndex = faceDofNumber[elemToFaces[ifaceLoc]];

    dFlux_dfp = 0.0;
    for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
    {
      dFlux_dfp[jfaceLoc] = dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc];
      dofColIndicesFacePres[jfaceLoc] = faceDofNumber[elemToFaces[jfaceLoc]];
    }

    // residual
    rhs->add( eqnRowIndex,
              flux );

    // jacobian -- derivative wrt local cell centered pressure term
    matrix->add( eqnRowIndex,
                 dofColIndexElemPres,
                 dFlux_dp );

    // jacobian -- derivatives wrt face pressure terms
    matrix->add( eqnRowIndex,
                 dofColIndicesFacePres,
                 dFlux_dfp );
  }
}

} // namespace SinglePhaseHybridFVMKernels

} // namespace geosx
