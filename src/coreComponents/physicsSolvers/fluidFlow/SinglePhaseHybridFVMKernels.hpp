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
 * @file SinglePhaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiRTInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geosx
{
namespace SinglePhaseHybridFVMKernels
{

/******************************** AssemblerKernelHelper ********************************/

struct AssemblerKernelHelper
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
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh facesb
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] elemDens the density at this elenent's center
   * @param[in] dElemDens_dp the derivative of the density at this element's center
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[out] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[out] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  ComputeOneSidedVolFluxes( arrayView1d< real64 const > const & facePres,
                            arrayView1d< real64 const > const & dFacePres,
                            arrayView1d< real64 const > const & faceGravCoef,
                            arraySlice1d< localIndex const > const & elemToFaces,
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
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {
      // now in the following nested loop,
      // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
      for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        // 1) compute the potential diff between the cell center and the face center
        real64 const ccPres = elemPres + dElemPres;
        real64 const fPres  = facePres[elemToFaces[jfaceLoc]] + dFacePres[elemToFaces[jfaceLoc]];

        real64 const ccGravCoef = elemGravCoef;
        real64 const fGravCoef  = faceGravCoef[elemToFaces[jfaceLoc]];

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
        real64 const potDif      = presDif     - gravTerm;
        real64 const dPotDif_dp  = dPresDif_dp - dGravTerm_dp;
        real64 const dPotDif_dfp = dPresDif_dfp;

        // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
        oneSidedVolFlux[ifaceLoc]                = oneSidedVolFlux[ifaceLoc]
                                                   + transMatrix[ifaceLoc][jfaceLoc] * potDif;
        dOneSidedVolFlux_dp[ifaceLoc]            = dOneSidedVolFlux_dp[ifaceLoc]
                                                   + transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dp;
        dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc] = dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc]
                                                   + transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dfp;
      }
    }
  }

  /**
   * @brief In a given element, collect the upwinded mobilities at this element's faces
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemToFaces elem-to-faces maps
   * @param[in] mobility the mobilities in the domain (non-local)
   * @param[in] dMobility_dp the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] elemDofNumber the dof numbers of all the cell centered pressures (non-local)
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[inout] upwMobility the upwinded mobilities at this element's faces
   * @param[inout] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or
   * neighbor)
   * @param[inout] upwDofNumber  the dof number of the upwind pressure
   *
   * Note: because of the upwinding, this function requires non-local information
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  UpdateUpwindedCoefficients( localIndex const er,
                              localIndex const esr,
                              localIndex const ei,
                              arrayView2d< localIndex const > const & elemRegionList,
                              arrayView2d< localIndex const > const & elemSubRegionList,
                              arrayView2d< localIndex const > const & elemList,
                              SortedArrayView< localIndex const > const & regionFilter,
                              arraySlice1d< localIndex const > const & elemToFaces,
                              ElementViewConst< arrayView1d< real64 const > > const & mobility,
                              ElementViewConst< arrayView1d< real64 const > > const & dMobility_dp,
                              ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                              arraySlice1d< real64 const > const & oneSidedVolFlux,
                              arraySlice1d< real64 > const & upwMobility,
                              arraySlice1d< real64 > const & dUpwMobility_dp,
                              arraySlice1d< globalIndex > const & upwDofNumber )
  {
    // for this element, loop over the local (one-sided) faces
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {

      // we initialize these upw quantities with the values of the local elem
      upwMobility[ifaceLoc]     = mobility[er][esr][ei];
      dUpwMobility_dp[ifaceLoc] = dMobility_dp[er][esr][ei];
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
              upwMobility[ifaceLoc]     = mobility[erNeighbor][esrNeighbor][eiNeighbor];
              dUpwMobility_dp[ifaceLoc] = dMobility_dp[erNeighbor][esrNeighbor][eiNeighbor];
              upwDofNumber[ifaceLoc]    = elemDofNumber[erNeighbor][esrNeighbor][eiNeighbor];
            }
            // if the face is on the boundary, use the properties of the local elem
          }
        }
      }
    }
  }


  /**
   * @brief In a given element, assemble the mass conservation equation
   * @param[in] dt the time step size
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] upwMobility the upwinded mobilities at this element's faces
   * @param[in] dUpwMobility_dp the derivatives of the upwinded mobilities wrt the cell-centered pressures (local or
   * neighbor)
   * @param[in] upwDofNumber  the dof number of the upwind pressure
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  AssembleOneSidedMassFluxes( arrayView1d< globalIndex const > const & faceDofNumber,
                              arraySlice1d< localIndex const > const & elemToFaces,
                              globalIndex const elemDofNumber,
                              globalIndex const rankOffset,
                              arraySlice1d< real64 const > const & oneSidedVolFlux,
                              arraySlice1d< real64 const > const & dOneSidedVolFlux_dp,
                              arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp,
                              arraySlice1d< real64 const > const & upwMobility,
                              arraySlice1d< real64 const > const & dUpwMobility_dp,
                              arraySlice1d< globalIndex const > const & upwDofNumber,
                              real64 const & dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
  {
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
      globalIndex const eqnRowLocalIndex = elemDofNumber - rankOffset;
      stackArray1d< globalIndex, 1+NF > elemDofColIndices( 1+NF );
      stackArray1d< globalIndex, NF >   faceDofColIndices( NF );
      elemDofColIndices[0] = elemDofNumber;

      // for each element, loop over the one-sided faces
      for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
      {

        real64 const dt_upwMob = dt * upwMobility[ifaceLoc];
        // compute the mass flux at the one-sided face plus its derivatives
        // add the newly computed flux to the sum
        sumLocalMassFluxes                        = sumLocalMassFluxes + dt_upwMob * oneSidedVolFlux[ifaceLoc];
        dSumLocalMassFluxes_dElemVars[0]          = dSumLocalMassFluxes_dElemVars[0] + dt * upwMobility[ifaceLoc] * dOneSidedVolFlux_dp[ifaceLoc];
        dSumLocalMassFluxes_dElemVars[ifaceLoc+1] = dt * dUpwMobility_dp[ifaceLoc] * oneSidedVolFlux[ifaceLoc];
        for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
        {
          dSumLocalMassFluxes_dFaceVars[jfaceLoc] = dSumLocalMassFluxes_dFaceVars[jfaceLoc]
                                                    + dt_upwMob * dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc];
        }

        // collect the relevant dof numbers
        elemDofColIndices[ifaceLoc+1] = upwDofNumber[ifaceLoc]; // if upwDofNumber == elemDofNumber, the derivative is zero
        faceDofColIndices[ifaceLoc] = faceDofNumber[elemToFaces[ifaceLoc]];
      }

      // we are ready to assemble the local flux and its derivatives
      // no need for atomic adds - each row is assembled by a single thread

      // residual
      localRhs[eqnRowLocalIndex] = localRhs[eqnRowLocalIndex] + sumLocalMassFluxes;

      // jacobian -- derivative wrt elem centered vars
      localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                                elemDofColIndices.data(),
                                                                dSumLocalMassFluxes_dElemVars.data(),
                                                                elemDofColIndices.size() );

      // jacobian -- derivatives wrt face centered vars
      localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                                faceDofColIndices.data(),
                                                                dSumLocalMassFluxes_dFaceVars.data(),
                                                                faceDofColIndices.size() );
    }
  }


  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  AssembleConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                       arrayView1d< integer const > const & faceGhostRank,
                       arraySlice1d< localIndex const > const & elemToFaces,
                       globalIndex const elemDofNumber,
                       globalIndex const rankOffset,
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
    globalIndex const dofColIndexElemPres = elemDofNumber;

    // for each element, loop over the local (one-sided) faces
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {
      if( faceGhostRank[elemToFaces[ifaceLoc]] >= 0 )
      {
        continue;
      }

      // flux at this face
      real64 const flux      = oneSidedVolFlux[ifaceLoc];
      real64 const dFlux_dp  = dOneSidedVolFlux_dp[ifaceLoc];

      // dof number of this face constraint
      localIndex const eqnLocalRowIndex = LvArray::integerConversion< localIndex >( faceDofNumber[elemToFaces[ifaceLoc]] - rankOffset );

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
                                                                        NF );
    }
  }

};

/******************************** AssemblerKernel ********************************/

struct AssemblerKernel
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
   * @brief In a given element, assemble the mass conservation equation and the contribution of this element to the face
   * constraints
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] faceGhostRank ghost rank of each face
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] eleomGravCoef the depth at this element's center
   * @param[in] elemDens the density at this elenent's center
   * @param[in] dElemDens_dp the derivative of the density at this element's center
   * @param[in] mobility the mobilities in the domain (non-local)
   * @param[in] dMobility_dPres the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[in] dt time step size
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  Compute( localIndex const er,
           localIndex const esr,
           localIndex const ei,
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
           real64 const & elemDens,
           real64 const & dElemDens_dp,
           ElementViewConst< arrayView1d< real64 const > > const & mobility,
           ElementViewConst< arrayView1d< real64 const > > const & dMobility_dp,
           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
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
    if( elemGhostRank < 0 )
    {
      AssemblerKernelHelper::UpdateUpwindedCoefficients< NF >( er, esr, ei,
                                                               elemRegionList,
                                                               elemSubRegionList,
                                                               elemList,
                                                               regionFilter,
                                                               elemToFaces,
                                                               mobility,
                                                               dMobility_dp,
                                                               elemDofNumber,
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
      AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >( faceDofNumber,
                                                               elemToFaces,
                                                               elemDofNumber[er][esr][ei],
                                                               rankOffset,
                                                               oneSidedVolFlux,
                                                               dOneSidedVolFlux_dp,
                                                               dOneSidedVolFlux_dfp,
                                                               upwMobility,
                                                               dUpwMobility_dp,
                                                               upwDofNumber,
                                                               dt,
                                                               localMatrix,
                                                               localRhs );
    }

    // use the computed one sided vol fluxes to assemble the constraints
    // enforcing flux continuity at this element's faces
    AssemblerKernelHelper::AssembleConstraints< NF >( faceDofNumber,
                                                      faceGhostRank,
                                                      elemToFaces,
                                                      elemDofNumber[er][esr][ei],
                                                      rankOffset,
                                                      oneSidedVolFlux,
                                                      dOneSidedVolFlux_dp,
                                                      dOneSidedVolFlux_dfp,
                                                      localMatrix,
                                                      localRhs );
  }


};

/******************************** FluxKernel ********************************/

struct FluxKernel
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
   * @brief Assemble the mass conservation equations and face constraints in the cell subregion
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion pointer to the cell element subregion
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] mesh the mesh object (single level only)
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] elemDens the density in the elements of the subregion
   * @param[in] dElemDens_dp the derivative of the density in the elements of the subregion
   * @param[in] mobility the mobilities in the domain (non-local)
   * @param[in] dMobility_dPres the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] lengthTolerance tolerance used in the transmissibility calculations
   * @param[in] dt time step size
   * @param[in] dofManager the dof manager
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< typename IP_TYPE, localIndex NF >
  static void
  Launch( localIndex er,
          localIndex esr,
          CellElementSubRegion const & subRegion,
          constitutive::SingleFluidBase const & fluid,
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
          ElementViewConst< arrayView1d< real64 const > > const & mobility,
          ElementViewConst< arrayView1d< real64 const > > const & dMobility_dp,
          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
          localIndex const rankOffset,
          real64 const lengthTolerance,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    // get the cell-centered DOF numbers and ghost rank for the assembly
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    // get the map from elem to faces
    arrayView2d< localIndex const > const elemToFaces = subRegion.faceList().toViewConst();

    // get the cell-centered pressures
    arrayView1d< real64 const > const elemPres  =
      subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::pressureString );
    arrayView1d< real64 const > const dElemPres =
      subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaPressureString );

    // get the element data needed for transmissibility computation
    arrayView2d< real64 const > const elemCenter =
      subRegion.getReference< array2d< real64 > >( CellBlock::viewKeyStruct::elementCenterString );
    arrayView1d< real64 const > const elemVolume =
      subRegion.getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );
    arrayView1d< R1Tensor const > const elemPerm =
      subRegion.getReference< array1d< R1Tensor > >( SinglePhaseBase::viewKeyStruct::permeabilityString );

    // get the cell-centered depth
    arrayView1d< real64 const > const elemGravCoef =
      subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::gravityCoefString );

    // get the fluid data
    arrayView2d< real64 const > const elemDens = fluid.density();
    arrayView2d< real64 const > const dElemDens_dp = fluid.dDensity_dPressure();

    // assemble the residual and Jacobian element by element
    // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
    using KERNEL_POLICY = parallelDevicePolicy< 32 >;
    forAll< KERNEL_POLICY >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const ei )
    {

      // transmissibility matrix
      stackArray2d< real64, NF *NF > transMatrix( NF, NF );

      real64 const perm[ 3 ] = { elemPerm[ei][0], elemPerm[ei][1], elemPerm[ei][2] };

      // recompute the local transmissibility matrix at each iteration
      // we can decide later to precompute transMatrix if needed
      IP_TYPE::template Compute< NF >( nodePosition,
                                       faceToNodes,
                                       elemToFaces[ei],
                                       elemCenter[ei],
                                       elemVolume[ei],
                                       perm,
                                       lengthTolerance,
                                       transMatrix );

      // perform flux assembly in this element
      SinglePhaseHybridFVMKernels::AssemblerKernel::Compute< NF >( er, esr, ei,
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
                                                                   elemDens[ei][0],
                                                                   dElemDens_dp[ei][0],
                                                                   mobility,
                                                                   dMobility_dp,
                                                                   elemDofNumber,
                                                                   elemGhostRank[ei],
                                                                   rankOffset,
                                                                   dt,
                                                                   transMatrix,
                                                                   localMatrix,
                                                                   localRhs );

    } );
  }


};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{
  template< typename VIEWTYPE >
  using ElementViewConst = typename ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static void
  Launch( LOCAL_VECTOR const localResidual,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & facePresDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ElementViewConst< arrayView1d< real64 const > > const & elemVolume,
          real64 const & defaultViscosity,
          real64 * localResidualNorm )
  {

    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( facePresDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      // if not ghost face and if adjacent to target region
      if( faceGhostRank[iface] < 0 && facePresDofNumber[iface] >= 0 )
      {
        real64 normalizer = 0;
        localIndex elemCounter = 0;
        for( localIndex k=0; k<elemRegionList.size( 1 ); ++k )
        {
          localIndex const er  = elemRegionList[iface][k];
          localIndex const esr = elemSubRegionList[iface][k];
          localIndex const ei  = elemList[iface][k];

          bool const onBoundary = (er == -1 || esr == -1 || ei == -1);

          // if not on boundary, save the mobility and the upwDofNumber
          if( !onBoundary )
          {
            normalizer += elemVolume[er][esr][ei];
            elemCounter++;
          }
        }
        normalizer /= elemCounter;
        normalizer /= defaultViscosity;

        localIndex const lid = facePresDofNumber[iface] - rankOffset;
        real64 const val = localResidual[lid] / normalizer; // to get something dimensionless
        sumScaled += val * val;
      }
    } );

    *localResidualNorm = *localResidualNorm + sumScaled.get();
  }

};

/******************************** Kernel switches ********************************/

namespace helpers
{

template< typename T, typename LAMBDA >
void KernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernelLaunchSelectorFaceSwitch: type should be integral" );

  switch( value )
  {
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return;}
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return;}
    case 6:
    { lambda( std::integral_constant< T, 6 >() ); return;}
    default: GEOSX_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

} // namespace helpers

template< typename IP_TYPE, typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector( localIndex numFacesInElem, ARGS && ... args )
{
  helpers::KernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
  {
    KERNELWRAPPER::template Launch< IP_TYPE, NF() >( std::forward< ARGS >( args )... );
  } );
}


} // namespace SinglePhaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
