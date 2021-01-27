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
 * @file CompositionalMultiphaseHybridFVMAssemblerHelperKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMASSEMBLERHELPERKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMASSEMBLERHELPERKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVMUpwindingHelperKernels.hpp"
#include "physicsSolvers/fluidFlow/HybridFVMHelperKernels.hpp"

namespace geosx
{

namespace CompositionalMultiphaseHybridFVMKernels
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
   * @brief In a given element, compute the transmissibility-weighted pressure gradients in the cell
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh facesb
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] phaseDens the phase densities in the element
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the element wrt pressure
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the element wrt component fraction
   * @param[in] phaseMob the phase mobilities in the element
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the element wrt pressure
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the element wrt component fraction
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the element wrt component density
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[out] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[out] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[out] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  applyGradient( arrayView1d< real64 const > const & facePres,
                 arrayView1d< real64 const > const & dFacePres,
                 arrayView1d< real64 const > const & faceGravCoef,
                 arraySlice1d< localIndex const > const & elemToFaces,
                 real64 const & elemPres,
                 real64 const & dElemPres,
                 real64 const & elemGravCoef,
                 arraySlice1d< real64 const > const & elemPhaseMassDens,
                 arraySlice1d< real64 const > const & dElemPhaseMassDens_dPres,
                 arraySlice2d< real64 const > const & dElemPhaseMassDens_dCompFrac,
                 arraySlice1d< real64 const > const & elemPhaseMob,
                 arraySlice1d< real64 const > const & dElemPhaseMob_dPres,
                 arraySlice2d< real64 const > const & dElemPhaseMob_dCompDens,
                 arraySlice2d< real64 const > const & dElemCompFrac_dCompDens,
                 arraySlice2d< real64 const > const & transMatrix,
                 real64 ( & oneSidedVolFlux )[ NF ],
                 real64 ( & dOneSidedVolFlux_dPres )[ NF ],
                 real64 ( & dOneSidedVolFlux_dFacePres )[ NF ][ NF ],
                 real64 ( & dOneSidedVolFlux_dCompDens )[ NF ][ NC ] )
  {
    real64 dPhaseMassDens_dC[ NP ][ NC ]{};
    real64 dPresDif_dCompDens[ NC ]{};
    real64 dPhaseGravDif_dCompDens[ NC ]{};
    real64 dPhaseMobPotDif_dCompDens[ NC ]{};

    // 0) precompute dPhaseDens_dC since it is always computed at the element center
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      applyChainRule( NC,
                      dElemCompFrac_dCompDens,
                      dElemPhaseMassDens_dCompFrac[ip],
                      dPhaseMassDens_dC[ip] );
    }

    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {
      // now in the following nested loop,
      // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
      for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {

        // depth difference between element center and face center
        real64 const ccGravCoef = elemGravCoef;
        real64 const fGravCoef = faceGravCoef[elemToFaces[jfaceLoc]];
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
          real64 const phaseGravDif = elemPhaseMassDens[ip] * gravCoefDif;
          real64 const dPhaseGravDif_dPres = dElemPhaseMassDens_dPres[ip] * gravCoefDif;
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            dPhaseGravDif_dCompDens[ic] = dPhaseMassDens_dC[ip][ic] * gravCoefDif;
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


  /**
   * @brief In a given element, compute the flux divergence, i.e, sum the fluxes at this element's faces
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] rankOffset offset of this rank
   * @param[in] faceDofNumber the dof numbers of the face pressures in the domain
   * @param[in] elemRegionList map from face to element region index
   * @param[in] elemSubRegionList map from face to element subRegion index
   * @param[in] elemList map from face to element index
   * @param[in] regionFilter set of target regions of the solver
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] phaseDens phaseDens the phase densities in the domain
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities wrt pressure
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities wrt component fraction
   * @param[in] phaseMob the phase mobilities in the domain
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the element wrt pressure
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the element wrt component density
   * @param[in] dCompFrac_dCompDens the derivatives of the component fractions wrt component density
   * @param[in] phaseCompFrac the phase component fractions in the domain
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions wrt pressure
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions wrt component fractions
   * @param[in] elemDofNumber the dof numbers of the element in the domain
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[in] dt the time step
   * @param[inout] localMatrix the Jacobian matrix
   * @param[inout] localRhs the residual
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  assembleFluxDivergence( localIndex const (&localIds)[ 3 ],
                          globalIndex const rankOffset,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          SortedArrayView< localIndex const > const & regionFilter,
                          arrayView1d< globalIndex const > const & faceDofNumber,
                          arrayView1d< real64 const > const & mimFaceGravCoef,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          real64 const & elemGravCoef,
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
                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                          arraySlice2d< real64 const > const & transMatrixGrav,
                          real64 const (&oneSidedVolFlux)[ NF ],
                          real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
  {
    localIndex constexpr NDOF = NC+1;

    // dof numbers
    globalIndex dofColIndicesElemVars[ NDOF*(NF+1) ]{};
    globalIndex dofColIndicesFaceVars[ NF ]{};
    for( localIndex idof = 0; idof < NDOF; ++idof )
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
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {

      // 1) Find if there is a neighbor, and if there is, grab the indices of the neighbor element

      localIndex neighborIds[ 3 ] = { localIds[0], localIds[1], localIds[2] };
      hybridFVMKernels::CellConnectivity::findNeighbor( localIds,
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
                                         faceDofNumber[elemToFaces[ifaceLoc]],
                                         dt,
                                         divMassFluxes,
                                         dDivMassFluxes_dElemVars,
                                         dDivMassFluxes_dFaceVars,
                                         dofColIndicesElemVars,
                                         dofColIndicesFaceVars );

      // 3) *************** Assemble buoyancy terms ******************

      real64 const transGravCoef = (localIds[0] != neighborIds[0] || localIds[1] != neighborIds[1] || localIds[2] != neighborIds[2])
                                   * transMatrixGrav[ifaceLoc][ifaceLoc] * (elemGravCoef - mimFaceGravCoef[elemToFaces[ifaceLoc]]);

      // 3.a) Compute the upwinded x_{c, \ell} \rho_{\ell} \frac{\lambda_{\ell}\lambda_m}{\lambda_T}
      //      and (\rho_{\ell} - \rho_m) g \Delta z for each phase at this face
      UpwindingHelper::upwindBuoyancyCoefficient< NC, NP >( localIds,
                                                            neighborIds,
                                                            transGravCoef,
                                                            phaseDens,
                                                            dPhaseDens_dPres,
                                                            dPhaseDens_dCompFrac,
                                                            phaseMassDens,
                                                            dPhaseMassDens_dPres,
                                                            dPhaseMassDens_dCompFrac,
                                                            phaseMob,
                                                            dPhaseMob_dPres,
                                                            dPhaseMob_dCompDens,
                                                            dCompFrac_dCompDens,
                                                            phaseCompFrac,
                                                            dPhaseCompFrac_dPres,
                                                            dPhaseCompFrac_dCompFrac,
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

    // we are ready to assemble the local flux and its derivatives
    // no need for atomic adds - each row is assembled by a single thread

    for( localIndex ic = 0; ic < NC; ++ic )
    {
      localIndex const eqnRowLocalIndex =
        LvArray::integerConversion< localIndex >( elemDofNumber[localIds[0]][localIds[1]][localIds[2]] + ic - rankOffset );

      GEOSX_ASSERT_GE( eqnRowLocalIndex, 0 );
      GEOSX_ASSERT_GT( localMatrix.numRows(), eqnRowLocalIndex );

      // residual
      localRhs[eqnRowLocalIndex] = localRhs[eqnRowLocalIndex] + divMassFluxes[ic];

      // jacobian -- derivative wrt elem centered vars
      localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                                &dofColIndicesElemVars[0],
                                                                &dDivMassFluxes_dElemVars[0][0] + ic * NDOF * (NF+1),
                                                                NDOF * (NF+1) );

      // jacobian -- derivatives wrt face centered vars
      localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                                &dofColIndicesFaceVars[0],
                                                                &dDivMassFluxes_dFaceVars[0][0] + ic * NF,
                                                                NF );
    }
  }


  /**
   * @brief In a given element, compute the viscous flux divergence, i.e, sum the viscous fluxes at this element's faces
   * @param[in] ifaceLoc the local index of the face
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face centered pressure
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[in] upwPhaseViscCoef the upwinded viscous transport coefficient at this face (ifaceLoc)
   * @param[in] dUpwPhaseViscCoef_dPres the derivative of the upwinded viscous transport coefficient wrt pressure
   * @param[in] dUpwPhaseViscCoef_dCompDens the derivative of the upwinded viscous transport coefficient wrt component density
   * @param[in] upwViscDofNumber degree of freedom number of the upwind element
   * @param[in] faceDofNumber degree of freedom number of the face
   * @param[in] dt the time step
   * @param[inout] divMassFluxes the divergence of the fluxes in the element
   * @param[inout] dDivMassFluxes_dElemVars the derivatives of the flux divergence wrt the element centered vars (pres and comp dens)
   * @param[inout] dDivMassFluxes_dFaceVars the derivatives of the flux divergence wrt the face centered vars
   * @param[inout] dofColIndicesElemVars degrees of freedom of the cells involved in the flux divergence
   * @param[inout] dofColIndicesFaceVars degrees of freedom of the faces involved in the flux divergence
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
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
                       globalIndex const faceDofNumber,
                       real64 const & dt,
                       real64 ( & divMassFluxes )[ NC ],
                       real64 ( & dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ],
                       real64 ( & dDivMassFluxes_dFaceVars )[ NC ][ NF ],
                       globalIndex ( & dofColIndicesElemVars )[ (NC+1)*(NF+1) ],
                       globalIndex ( & dofColIndicesFaceVars )[ NF ] )
  {
    localIndex constexpr NDOF = NC+1;
    localIndex const elemVarsOffset = NDOF*(ifaceLoc+1);

    for( localIndex ip = 0; ip < NP; ++ip )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
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
        for( localIndex jc = 0; jc < NC; ++jc )
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
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                              + ( elemDofNumber != upwViscDofNumber )
                                                              * dt * dUpwPhaseViscCoef_dCompDens[ip][ic][jc] * oneSidedVolFlux[ifaceLoc];
        }

        for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
        {
          dDivMassFluxes_dFaceVars[ic][jfaceLoc] = dDivMassFluxes_dFaceVars[ic][jfaceLoc]
                                                   + dt_upwPhaseViscCoef * dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc];
        }
      }
    }

    // collect the relevant dof numbers
    for( localIndex idof = 0; idof < NDOF; ++idof )
    {
      dofColIndicesElemVars[elemVarsOffset+idof] = neighborDofNumber + idof;
    }
    dofColIndicesFaceVars[ifaceLoc] = faceDofNumber;
  }

  /**
   * @brief In a given element, compute the buoyancy flux divergence, i.e, sum the buoyancy fluxes at this element's faces
   * @param[in] ifaceLoc the local index of the face
   * @param[in] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[in] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[in] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   * @param[in] upwPhaseGravCoef the upwinded buoyancy transport coefficient at this face (ifaceLoc)
   * @param[in] dUpwPhaseGravCoef_dPres the derivative of the upwinded buoyancy transport coefficient wrt pressure
   * @param[in] dUpwPhaseGravCoef_dCompDens the derivative of the upwinded buoyancy transport coefficient wrt component density
   * @param[in] neighborDofNumber the degree of freedom number of the neighbor element
   * @param[in] dt the time step
   * @param[inout] divMassFluxes the divergence of the fluxes in the element
   * @param[inout] dDivMassFluxes_dElemVars the derivatives of the flux divergence wrt the element centered vars (pres and comp dens)
   * @param[inout] dofColIndicesElemVars degrees of freedom of the cells involved in the flux divergence
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
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
    localIndex constexpr NDOF = NC+1;
    localIndex const elemVarsOffset = NDOF*(ifaceLoc+1);

    for( localIndex ip = 0; ip < NP; ++ip )
    {
      for( localIndex jp = 0; jp < NP - 1; ++jp )
      {
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          real64 const dt_upwPhaseGravCoef = dt * upwPhaseGravCoef[ip][jp][ic];

          // residual
          divMassFluxes[ic] = divMassFluxes[ic] + dt_upwPhaseGravCoef * phaseGravTerm[ip][jp];

          // local derivatives
          dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                            + dt * dUpwPhaseGravCoef_dPres[ip][jp][ic][Pos::LOCAL] * phaseGravTerm[ip][jp];
          dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                            + dt_upwPhaseGravCoef * dPhaseGravTerm_dPres[ip][jp][Pos::LOCAL];

          for( localIndex jc = 0; jc < NC; ++jc )
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

          for( localIndex jc = 0; jc < NC; ++jc )
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

  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures in the subRegion
   * @param[in] faceGhostRank the ghost ranks of the face pressures in the subRegion
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] rankOffset the offset of this rank
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  assembleFaceConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
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
    real64 dFlux_dElemVars[ NDOF ]{};
    real64 dFlux_dFaceVars[ NF ]{};

    // dof numbers
    globalIndex dofColIndicesElemVars[ NDOF ]{};
    globalIndex dofColIndicesFaceVars[ NF ]{};
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

};

} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMASSEMBLERHELPERKERNELS_HPP
