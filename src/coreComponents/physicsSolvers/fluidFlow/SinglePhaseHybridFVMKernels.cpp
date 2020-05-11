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

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "linearAlgebra/DofManager.hpp"


namespace geosx
{

namespace SinglePhaseHybridFVMKernels
{

/******************************** AssemblerKernelHelper ********************************/

template< localIndex NF >
void
AssemblerKernelHelper::ComputeOneSidedVolFluxes( arrayView1d< real64 const > const & facePres,
                                                 arrayView1d< real64 const > const & dFacePres,
                                                 arrayView1d< real64 const > const & faceGravCoef,
                                                 arraySlice1d< localIndex const > const elemToFaces,
                                                 real64 const & elemPres,
                                                 real64 const & dElemPres,
                                                 real64 const & elemGravCoef,
                                                 real64 const & elemDens,
                                                 real64 const & dElemDens_dp,
                                                 arraySlice2d< real64 const > const & transMatrix,
                                                 arraySlice1d< real64 > const & oneSidedVolFlux,
                                                 arraySlice1d< real64 > const & dOneSidedVolFlux_dp,
                                                 arraySlice2d< real64 > const & dOneSidedVolFlux_dfp )
{
  stackArray1d< real64, NF > potDif( NF );
  stackArray1d< real64, NF > dPotDif_dp( NF );
  stackArray1d< real64, NF > dPotDif_dfp( NF );

  // 1) precompute the potential difference at each one-sided face
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
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
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
      oneSidedVolFlux[ifaceLoc]                += transMatrix[ifaceLoc][jfaceLoc] * potDif[jfaceLoc];
      dOneSidedVolFlux_dp[ifaceLoc]            += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dp[jfaceLoc];
      dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc] += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dfp[jfaceLoc];
    }
  }
}

template< localIndex NF >
void
AssemblerKernelHelper::UpdateUpwindedCoefficients( arrayView2d< localIndex const > const & elemRegionList,
                                                   arrayView2d< localIndex const > const & elemSubRegionList,
                                                   arrayView2d< localIndex const > const & elemList,
                                                   SortedArrayView< localIndex const > const & regionFilter,
                                                   arraySlice1d< localIndex const > const elemToFaces,
                                                   ElementView< arrayView1d< real64 const > > const & mob,
                                                   ElementView< arrayView1d< real64 const > > const & dMob_dp,
                                                   ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
                                                   localIndex const er,
                                                   localIndex const esr,
                                                   localIndex const ei,
                                                   arraySlice1d< real64 const > const & oneSidedVolFlux,
                                                   arraySlice1d< real64 > const & upwMobility,
                                                   arraySlice1d< real64 > const & dUpwMobility_dp,
                                                   arraySlice1d< globalIndex > const & upwDofNumber )
{
  // for this element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
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

template< localIndex NF >
void
AssemblerKernelHelper::AssembleOneSidedMassFluxes( real64 const & dt,
                                                   globalIndex const rankOffset,
                                                   arrayView1d< globalIndex const > const & faceDofNumber,
                                                   arraySlice1d< localIndex const > const elemToFaces,
                                                   globalIndex const elemDofNumber,
                                                   arraySlice1d< real64 const > const & oneSidedVolFlux,
                                                   arraySlice1d< real64 const > const & dOneSidedVolFlux_dp,
                                                   arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp,
                                                   arraySlice1d< real64 const > const & upwMobility,
                                                   arraySlice1d< real64 const > const & dUpwMobility_dp,
                                                   arraySlice1d< globalIndex const > const & upwDofNumber,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  // fluxes
  real64 sumLocalMassFluxes = 0;
  stackArray1d< real64, 1+NF > dSumLocalMassFluxes_dElemVars( 1+NF );
  stackArray1d< real64, NF >   dSumLocalMassFluxes_dFaceVars( NF );
  for( localIndex i = 0; i < NF+1; ++i )
  {
    dSumLocalMassFluxes_dElemVars( i ) = 0.;
  }
  for( localIndex i = 0; i < NF; ++i )
  {
    dSumLocalMassFluxes_dFaceVars( i ) = 0.;
  }

  // dof numbers
  globalIndex const eqnRowIndex = elemDofNumber - rankOffset;
  stackArray1d< globalIndex, 1+NF > elemDofColIndices( 1+NF );
  stackArray1d< globalIndex, NF >   faceDofColIndices( NF );
  elemDofColIndices[0] = elemDofNumber;

  // for each element, loop over the one-sided faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {

    // compute the mass flux at the one-sided face plus its derivatives
    // add the newly computed flux to the sum
    sumLocalMassFluxes                        += dt * upwMobility[ifaceLoc] * oneSidedVolFlux[ifaceLoc];
    dSumLocalMassFluxes_dElemVars[0]          += dt * upwMobility[ifaceLoc] * dOneSidedVolFlux_dp[ifaceLoc];
    dSumLocalMassFluxes_dElemVars[ifaceLoc+1]  = dt * dUpwMobility_dp[ifaceLoc] * oneSidedVolFlux[ifaceLoc];
    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      dSumLocalMassFluxes_dFaceVars[jfaceLoc] += dt * upwMobility[ifaceLoc] * dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc];
    }

    // collect the relevant dof numbers
    elemDofColIndices[ifaceLoc+1] = upwDofNumber[ifaceLoc]; // if upwDofNumber == elemDofNumber, the derivative is zero
    faceDofColIndices[ifaceLoc] = faceDofNumber[elemToFaces[ifaceLoc]];
  }

  // we are ready to assemble the local flux and its derivatives
  // no need for atomic adds - each row is assembled by a single thread

  // residual
  localRhs[eqnRowLocalIndex] += sumLocalMassFluxes;

  // jacobian -- derivative wrt elem centered vars
  localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                            elemDofColIndices.data(),
                                                            dSumLocalMassFluxes_dElemVars,
                                                            elemDofColIndices.size() );

  // jacobian -- derivatives wrt face centered vars
  localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                            faceDofColIndices.data(),
                                                            dSumLocalMassFluxes_dFaceVars,
                                                            faceDofColIndices.size() );
}


template< localIndex NF >
void
AssemblerKernelHelper::AssembleConstraints( globalIndex const rankOffset,
                                            arrayView1d< globalIndex const > const & faceDofNumber,
                                            arraySlice1d< localIndex const > const elemToFaces,
                                            globalIndex const elemDofNumber,
                                            arraySlice1d< real64 const > const & oneSidedVolFlux,
                                            arraySlice1d< real64 const > const & dOneSidedVolFlux_dp,
                                            arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
{
  // fluxes
  stackArray1d< real64, NF > dFlux_dfp( NF );

  // dof numbers
  stackArray1d< globalIndex, NF > dofColIndicesFacePres( NF );
  //globalIndex const dofColIndexElemPres = elemDofNumber;

  // for each element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    // flux at this face
    real64 const flux      = oneSidedVolFlux[ifaceLoc];
    real64 const dFlux_dp  = dOneSidedVolFlux_dp[ifaceLoc];

    // dof number of this face constraint
    localIndex const eqnLocalRowIndex = faceDofNumber[elemToFaces[ifaceLoc]] - rankOffset;

    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      dFlux_dfp[jfaceLoc] = dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc];
      dofColIndicesFacePres[jfaceLoc] = faceDofNumber[elemToFaces[jfaceLoc]];
    }

    // residual
    atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnLocalRowIndex], flux );

    // jacobian -- derivative wrt local cell centered pressure term
    localMatrix.addToRow< parallelDeviceAtomic >( eqnLocalRowIndex, &dofColIndexElemPres, &dFlux_dp, 1 );

    // jacobian -- derivatives wrt face pressure terms
    localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                      dofColIndicesFacePres.data(),
                                                                      dFlux_dfp.data(),
                                                                      numFacesInElem );
  }
}

/******************************** AssemblerKernel ********************************/

template< localIndex NF >
void
AssemblerKernel::Compute( localIndex const er,
                          localIndex const esr,
                          localIndex const ei,
                          SortedArrayView< localIndex const > const & regionFilter,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          arrayView1d< globalIndex const > const & faceDofNumber,
                          arrayView1d< real64 const > const & facePres,
                          arrayView1d< real64 const > const & dFacePres,
                          arrayView1d< real64 const > const & faceGravCoef,
                          arraySlice1d< localIndex const > const elemToFaces,
                          real64 const & elemPres,
                          real64 const & dElemPres,
                          real64 const & elemGravCoef,
                          real64 const & elemDens,
                          real64 const & dElemDens_dp,
                          ElementView< arrayView1d< real64 const > > const & mobility,
                          ElementView< arrayView1d< real64 const > > const & dMobility_dp,
                          ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
                          arraySlice2d< real64 const > const & transMatrix,
                          real64 const & dt,
                          ParallelMatrix * const matrix,
                          ParallelVector * const rhs )
{

  // one sided flux
  stackArray1d< real64, NF > oneSidedVolFlux( NF );
  stackArray1d< real64, NF > dOneSidedVolFlux_dp( NF );
  stackArray2d< real64, NF *NF > dOneSidedVolFlux_dfp( NF, NF );
  for( localIndex i = 0; i < NF; ++i )
  {
    oneSidedVolFlux( i ) = 0.;
    dOneSidedVolFlux_dp( i ) = 0.;
    for( localIndex j = 0; j < NF; ++j )
    {
      dOneSidedVolFlux_dfp( i, j ) = 0.; // assume row major
    }
  }

  // upwinded mobility
  stackArray1d< real64, NF > upwMobility( NF );
  stackArray1d< real64, NF > dUpwMobility_dp( NF );
  stackArray1d< globalIndex, NF > upwDofNumber( NF );

  /*
   * compute auxiliary quantities at the one sided faces of this element:
   * 1) One-sided volumetric fluxes
   * 2) Upwinded mobilities
   */

  // for each one-sided face of the elem,
  // compute the volumetric flux using transMatrix
  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
                                                         dFacePres,
                                                         faceGravCoef,
                                                         elemToFaces,
                                                         elemPres,
                                                         dElemPres,
                                                         elemGravCoef,
                                                         elemDens,
                                                         dElemDens_dp,
                                                         transMatrix,
                                                         oneSidedVolFlux,
                                                         dOneSidedVolFlux_dp,
                                                         dOneSidedVolFlux_dfp );

  // at this point, we know the local flow direction in the element
  // so we can upwind the transport coefficients (mobilities) at the one sided faces
  // ** this function needs non-local information **
  AssemblerKernelHelper::UpdateUpwindedCoefficients< NF >( elemRegionList,
                                                           elemSubRegionList,
                                                           elemList,
                                                           regionFilter,
                                                           elemToFaces,
                                                           mobility,
                                                           dMobility_dp,
                                                           elemDofNumber,
                                                           er, esr, ei,
                                                           oneSidedVolFlux,
                                                           upwMobility,
                                                           dUpwMobility_dp,
                                                           upwDofNumber );

  /*
   * perform assembly in this element in two steps:
   * 1) mass conservation equations
   * 2) face constraints
   */

  // use the computed one sided vol fluxes and the upwinded mobilities
  // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
  AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >( dt,
                                                           faceDofNumber,
                                                           elemToFaces,
                                                           elemDofNumber[er][esr][ei],
                                                           oneSidedVolFlux,
                                                           dOneSidedVolFlux_dp,
                                                           dOneSidedVolFlux_dfp,
                                                           upwMobility,
                                                           dUpwMobility_dp,
                                                           upwDofNumber,
                                                           matrix,
                                                           rhs );

  // use the computed one sided vol fluxes to assemble the constraints
  // enforcing flux continuity at this element's faces
  AssemblerKernelHelper::AssembleConstraints< NF >( faceDofNumber,
                                                    elemToFaces,
                                                    elemDofNumber[er][esr][ei],
                                                    oneSidedVolFlux,
                                                    dOneSidedVolFlux_dp,
                                                    dOneSidedVolFlux_dfp,
                                                    matrix,
                                                    rhs );
}


/******************************** FluxKernel ********************************/

template< localIndex NF >
void
FluxKernel::Launch( localIndex er,
                    localIndex esr,
                    CellElementSubRegion const & subRegion,
                    SortedArrayView< localIndex const > const & regionFilter,
                    MeshLevel const & mesh,
                    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                    arrayView2d< localIndex const > const & elemRegionList,
                    arrayView2d< localIndex const > const & elemSubRegionList,
                    arrayView2d< localIndex const > const & elemList,
                    ArrayOfArraysView< localIndex const > const & faceToNodes,
                    arrayView1d< globalIndex const > const & faceDofNumber,
                    arrayView1d< real64 const > const & facePres,
                    arrayView1d< real64 const > const & dFacePres,
                    arrayView1d< real64 const > const & faceGravCoef,
                    arrayView2d< real64 const > const & elemDens,
                    arrayView2d< real64 const > const & dElemDens_dp,
                    ElementView< arrayView1d< real64 const > > const & mobility,
                    ElementView< arrayView1d< real64 const > > const & dMobility_dp,
                    real64 const lengthTolerance,
                    real64 const dt,
                    DofManager const * const dofManager,
                    ParallelMatrix * const matrix,
                    ParallelVector * const rhs )
{

  // get the cell-centered DOF numbers and ghost rank for the assembly
  string const elemDofKey = dofManager->getKey( SinglePhaseBase::viewKeyStruct::pressureString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  dofNumberAccessor = mesh.getElemManager()->ConstructViewAccessor< array1d< globalIndex >,
                                                                    arrayView1d< globalIndex const > >( elemDofKey );

  ElementView< arrayView1d< globalIndex const > > const & elemDofNumber = dofNumberAccessor.toViewConst();

  arrayView1d< integer const > const & elemGhostRank =
    subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

  // get the map from elem to faces
  arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

  // get the cell-centered pressures
  arrayView1d< real64 const > const & elemPres  =
    subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dElemPres =
    subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaPressureString );

  // get the element data needed for transmissibility computation
  arrayView1d< R1Tensor const > const & elemCenter =
    subRegion.getReference< array1d< R1Tensor > >( CellBlock::viewKeyStruct::elementCenterString );
  arrayView1d< real64 const > const & elemVolume =
    subRegion.getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );
  arrayView1d< R1Tensor const > const & elemPerm =
    subRegion.getReference< array1d< R1Tensor > >( SinglePhaseBase::viewKeyStruct::permeabilityString );

  // get the cell-centered depth
  arrayView1d< real64 const > const & elemGravCoef =
    subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::gravityCoefString );

  // assemble the residual and Jacobian element by element
  // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
  using KERNEL_POLICY = parallelDevicePolicy< 32 >;
  forAll< KERNEL_POLICY >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const ei )
  {

    if( elemGhostRank[ei] < 0 )
    {

      // transmissibility matrix
      stackArray2d< real64, NF *NF > transMatrix( NF, NF );

      // recompute the local transmissibility matrix at each iteration
      // we can decide later to precompute transMatrix if needed
      HybridFVMInnerProduct::QTPFACellInnerProductKernel::Compute< NF >( nodePosition,
                                                                         faceToNodes,
                                                                         elemToFaces[ei],
                                                                         elemCenter[ei],
                                                                         elemVolume[ei],
                                                                         elemPerm[ei],
                                                                         2,
                                                                         lengthTolerance,
                                                                         transMatrix );

      // perform flux assembly in this element
      SinglePhaseHybridFVMKernels::AssemblerKernel::Compute< NF >( er, esr, ei,
                                                                   regionFilter,
                                                                   elemRegionList,
                                                                   elemSubRegionList,
                                                                   elemList,
                                                                   faceDofNumber,
                                                                   facePres,
                                                                   dFacePres,
                                                                   faceGravCoef,
                                                                   elemToFaces[ei],
                                                                   elemPres[ei],
                                                                   dElemPres[ei],
                                                                   elemGravCoef[ei],
                                                                   elemDens[ei][0],
                                                                   dElemDens_dp[ei][0],
                                                                   mobility,
                                                                   dMobility_dp,
                                                                   elemDofNumber.toViewConst(),
                                                                   transMatrix,
                                                                   dt,
                                                                   matrix,
                                                                   rhs );

    }
  } );
}


#define INST_AssembleKernelHelper( NF ) \
  template \
  void \
  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( arrayView1d< real64 const > const & facePres, \
                                                         arrayView1d< real64 const > const & dFacePres, \
                                                         arrayView1d< real64 const > const & faceGravCoef, \
                                                         arraySlice1d< localIndex const > const elemToFaces, \
                                                         real64 const & elemPres, \
                                                         real64 const & dElemPres, \
                                                         real64 const & elemGravDepth, \
                                                         real64 const & elemDens, \
                                                         real64 const & dElemDens_dp, \
                                                         arraySlice2d< real64 const > const & transMatrix, \
                                                         arraySlice1d< real64 > const & oneSidedVolFlux, \
                                                         arraySlice1d< real64 > const & dOneSidedVolFlux_dp, \
                                                         arraySlice2d< real64 > const & dOneSidedVolFlux_dfp ); \
  template \
  void \
  AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >( real64 const & dt, \
                                                           arrayView1d< globalIndex const > const & faceDofNumber, \
                                                           arraySlice1d< localIndex const > const elemToFaces, \
                                                           globalIndex const elemDofNumber, \
                                                           arraySlice1d< real64 const > const & oneSidedVolFlux, \
                                                           arraySlice1d< real64 const > const & dOneSidedVolFlux_dp, \
                                                           arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp, \
                                                           arraySlice1d< real64 const > const & upwMobility, \
                                                           arraySlice1d< real64 const > const & dUpwMobility_dp, \
                                                           arraySlice1d< globalIndex const > const & upwDofNumber, \
                                                           ParallelMatrix * const matrix, \
                                                           ParallelVector * const rhs ); \
  template \
  void \
  AssemblerKernelHelper::AssembleConstraints< NF >( arrayView1d< globalIndex const > const & faceDofNumber, \
                                                    arraySlice1d< localIndex const > const elemToFaces, \
                                                    globalIndex const elemDofNumber, \
                                                    arraySlice1d< real64 const > const & oneSidedVolFlux, \
                                                    arraySlice1d< real64 const > const & dOneSidedVolFlux_dp, \
                                                    arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp, \
                                                    ParallelMatrix * const matrix, \
                                                    ParallelVector * const rhs )


INST_AssembleKernelHelper( 4 );
INST_AssembleKernelHelper( 5 );
INST_AssembleKernelHelper( 6 );

#undef INST_AssembleKernelHelper


#define INST_FluxKernel( NF ) \
  template \
  void FluxKernel::Launch< NF >( localIndex er, \
                                 localIndex esr, \
                                 CellElementSubRegion const & subRegion, \
                                 SortedArrayView< localIndex const > const & regionFilter, \
                                 MeshLevel const & mesh, \
                                 arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
                                 arrayView2d< localIndex const > const & elemRegionList, \
                                 arrayView2d< localIndex const > const & elemSubRegionList, \
                                 arrayView2d< localIndex const > const & elemList, \
                                 ArrayOfArraysView< localIndex const > const & faceToNodes, \
                                 arrayView1d< globalIndex const > const & faceDofNumber, \
                                 arrayView1d< real64 const > const & facePres, \
                                 arrayView1d< real64 const > const & dFacePres, \
                                 arrayView1d< real64 const > const & faceGravCoef, \
                                 arrayView2d< real64 const > const & elemDens, \
                                 arrayView2d< real64 const > const & dElemDens_dp, \
                                 ElementView< arrayView1d< real64 const > > const & mobility, \
                                 ElementView< arrayView1d< real64 const > > const & dMobility_dp, \
                                 real64 const lengthTolerance, \
                                 real64 const dt, \
                                 DofManager const * const dofManager, \
                                 ParallelMatrix * const matrix, \
                                 ParallelVector * const rhs )

INST_FluxKernel( 4 );
INST_FluxKernel( 5 );
INST_FluxKernel( 6 );

#undef INST_FluxKernel


} // namespace SinglePhaseHybridFVMKernels

} // namespace geosx
