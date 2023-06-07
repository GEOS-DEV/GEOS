/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file AssemblerKernelHelper.cpp
 */

#include "AssemblerKernelHelper.hpp"

#include "KernelUtilities.hpp"
#include "CellConnectivity.hpp"
#include "UpwindingHelper.hpp"

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"


namespace geos
{

namespace compositionalMultiphaseHybridFVMKernels
{

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

template< integer NF, integer NC, integer NP >
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
                          arrayView1d< real64 const > const & mimFaceGravCoef,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          real64 const & elemGravCoef,
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
                          arraySlice2d< real64 const > const & transMatrixGrav,
                          real64 const (&oneSidedVolFlux)[ NF ],
                          real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
{
  using namespace compositionalMultiphaseUtilities;
  integer constexpr NDOF = NC+1;

  // dof numbers
  globalIndex dofColIndicesElemVars[ NDOF*(NF+1) ]{};
  globalIndex dofColIndicesFaceVars[ NF ]{};
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

  // Apply equation/variable change transformation(s)
  real64 work[NDOF*(NF+1)];
  shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NDOF * ( NF + 1 ), dDivMassFluxes_dElemVars, work );
  shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NF, dDivMassFluxes_dFaceVars, work );
  shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, divMassFluxes );

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
                       globalIndex const faceDofNumber,
                       real64 const & dt,
                       real64 ( & divMassFluxes )[ NC ],
                       real64 ( & dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ],
                       real64 ( & dDivMassFluxes_dFaceVars )[ NC ][ NF ],
                       globalIndex ( & dofColIndicesElemVars )[ (NC+1)*(NF+1) ],
                       globalIndex ( & dofColIndicesFaceVars )[ NF ] )
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
  dofColIndicesFaceVars[ifaceLoc] = faceDofNumber;
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

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
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
    if( faceGhostRank[elemToFaces[ifaceLoc]] >= 0 )
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


#define INST_AssemblerKernelHelper( NF, NC, NP ) \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    applyGradient< NF, NC, NP >( arrayView1d< real64 const > const & facePres, \
                                 arrayView1d< real64 const > const & faceGravCoef, \
                                 arraySlice1d< localIndex const > const & elemToFaces, \
                                 real64 const & elemPres, \
                                 real64 const & elemGravCoef, \
                                 arraySlice1d< real64 const, multifluid::USD_PHASE-2 > const & elemPhaseMassDens, \
                                 arraySlice2d< real64 const, multifluid::USD_PHASE_DC-2 > const & dElemPhaseMassDens_dCompFrac, \
                                 arraySlice1d< real64 const, compflow::USD_PHASE-1 > const & elemPhaseMob, \
                                 arraySlice2d< real64 const, compflow::USD_PHASE_DC-1 > const & dElemPhaseMob, \
                                 arraySlice2d< real64 const, compflow::USD_COMP_DC-1 > const & dElemCompFrac_dCompDens, \
                                 arraySlice2d< real64 const > const & transMatrix, \
                                 real64 ( &oneSidedVolFlux )[ NF ], \
                                 real64 ( &dOneSidedVolFlux_dPres )[ NF ], \
                                 real64 ( &dOneSidedVolFlux_dFacePres )[ NF ][ NF ], \
                                 real64 ( &dOneSidedVolFlux_dCompDens )[ NF ][ NC ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    assembleFluxDivergence< NF, NC, NP >( localIndex const (&localIds)[ 3 ], \
                                          globalIndex const rankOffset, \
                                          arrayView2d< localIndex const > const & elemRegionList, \
                                          arrayView2d< localIndex const > const & elemSubRegionList, \
                                          arrayView2d< localIndex const > const & elemList, \
                                          SortedArrayView< localIndex const > const & regionFilter, \
                                          arrayView1d< globalIndex const > const & faceDofNumber, \
                                          arrayView1d< real64 const > const & mimFaceGravCoef, \
                                          arraySlice1d< localIndex const > const & elemToFaces, \
                                          real64 const & elemGravCoef, \
                                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                          arraySlice2d< real64 const > const & transMatrixGrav, \
                                          real64 const (&oneSidedVolFlux)[ NF ], \
                                          real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                          real64 const & dt, \
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                          arrayView1d< real64 > const & localRhs ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    assembleViscousFlux< NF, NC, NP >( localIndex const ifaceLoc, \
                                       real64 const (&oneSidedVolFlux)[ NF ], \
                                       real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                       real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                       real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                       real64 const (&upwPhaseViscCoef)[ NP ][ NC ], \
                                       real64 const (&dUpwPhaseViscCoef_dPres)[ NP ][ NC ], \
                                       real64 const (&dUpwPhaseViscCoef_dCompDens)[ NP ][ NC ][ NC ], \
                                       globalIndex const elemDofNumber, \
                                       globalIndex const neighborDofNumber, \
                                       globalIndex const upwViscDofNumber, \
                                       globalIndex const faceDofNumber, \
                                       real64 const & dt, \
                                       real64 ( &divMassFluxes )[ NC ], \
                                       real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ], \
                                       real64 ( &dDivMassFluxes_dFaceVars )[ NC ][ NF ], \
                                       globalIndex ( &dofColIndicesElemVars )[ (NC+1)*(NF+1) ], \
                                       globalIndex ( &dofColIndicesFaceVars )[ NF ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    assembleBuoyancyFlux< NF, NC, NP >( localIndex const ifaceLoc, \
                                        real64 const (&phaseGravTerm)[ NP ][ NP-1 ], \
                                        real64 const (&dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ], \
                                        real64 const (&dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ], \
                                        real64 const (&upwPhaseGravCoef)[ NP ][ NP-1 ][ NC ], \
                                        real64 const (&dUpwPhaseGravCoef_dPres)[ NP ][ NP-1 ][ NC ][ 2 ], \
                                        real64 const (&dUpwPhaseGravCoef_dCompDens)[ NP ][ NP-1 ][ NC ][ 2 ][ NC ], \
                                        real64 const & dt, \
                                        real64 ( &divMassFluxes )[ NC ], \
                                        real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    assembleFaceConstraints< NF, NC, NP >( arrayView1d< globalIndex const > const & faceDofNumber, \
                                           arrayView1d< integer const > const & faceGhostRank, \
                                           arraySlice1d< localIndex const > const & elemToFaces, \
                                           globalIndex const elemDofNumber, \
                                           globalIndex const rankOffset,  \
                                           real64 const (&oneSidedVolFlux)[ NF ], \
                                           real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                           real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                           real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                           arrayView1d< real64 > const & localRhs )

INST_AssemblerKernelHelper( 4, 1, 2 );
INST_AssemblerKernelHelper( 4, 2, 2 );
INST_AssemblerKernelHelper( 4, 3, 2 );
INST_AssemblerKernelHelper( 4, 4, 2 );
INST_AssemblerKernelHelper( 4, 5, 2 );

INST_AssemblerKernelHelper( 4, 1, 3 );
INST_AssemblerKernelHelper( 4, 2, 3 );
INST_AssemblerKernelHelper( 4, 3, 3 );
INST_AssemblerKernelHelper( 4, 4, 3 );
INST_AssemblerKernelHelper( 4, 5, 3 );

INST_AssemblerKernelHelper( 5, 1, 2 );
INST_AssemblerKernelHelper( 5, 2, 2 );
INST_AssemblerKernelHelper( 5, 3, 2 );
INST_AssemblerKernelHelper( 5, 4, 2 );
INST_AssemblerKernelHelper( 5, 5, 2 );

INST_AssemblerKernelHelper( 5, 1, 3 );
INST_AssemblerKernelHelper( 5, 2, 3 );
INST_AssemblerKernelHelper( 5, 3, 3 );
INST_AssemblerKernelHelper( 5, 4, 3 );
INST_AssemblerKernelHelper( 5, 5, 3 );

INST_AssemblerKernelHelper( 6, 1, 2 );
INST_AssemblerKernelHelper( 6, 2, 2 );
INST_AssemblerKernelHelper( 6, 3, 2 );
INST_AssemblerKernelHelper( 6, 4, 2 );
INST_AssemblerKernelHelper( 6, 5, 2 );

INST_AssemblerKernelHelper( 6, 1, 3 );
INST_AssemblerKernelHelper( 6, 2, 3 );
INST_AssemblerKernelHelper( 6, 3, 3 );
INST_AssemblerKernelHelper( 6, 4, 3 );
INST_AssemblerKernelHelper( 6, 5, 3 );

#undef INST_AssemblerKernelHelper

}

}
