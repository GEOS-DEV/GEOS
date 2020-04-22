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
 * @file CompositionalMultiphaseHybridFVMKernels.cpp
 */

#include "CompositionalMultiphaseHybridFVMKernels.hpp"
#include "CompositionalMultiphaseBase.hpp"

namespace geosx
{

namespace CompositionalMultiphaseHybridFVMKernels
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
                                                 arraySlice2d< real64 const > const dElemDens_dComp,
                                                 arraySlice2d< real64 const > const dElemCompFrac_dCompDens,
                                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                      *HybridFVMInnerProduct::MAX_NUM_FACES > const & transMatrix,
                                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                      *constitutive::MultiFluidBase::MAX_NUM_PHASES > & volFlux,
                                                 stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                      *constitutive::MultiFluidBase::MAX_NUM_PHASES > & dVolFlux_dPres,
                                                 stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                      *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                      *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > & dVolFlux_dComp,
                                                 stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                      *HybridFVMInnerProduct::MAX_NUM_FACES
                                                                      *constitutive::MultiFluidBase::MAX_NUM_PHASES > & dVolFlux_dFacePotential )
{
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;
  localIndex constexpr maxNumPhases = constitutive::MultiFluidBase::MAX_NUM_PHASES;
  localIndex constexpr maxNumComponents = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  
  localIndex const numFacesInElem = elemToFaces.size();
  localIndex const numPhases = dVolFlux_dComp.size( 1 );
  localIndex const numComponents = dVolFlux_dComp.size( 2 );
  
  volFlux = 0;
  dVolFlux_dPres = 0;
  dVolFlux_dComp = 0;
  dVolFlux_dFacePotential = 0;

  // this is going to store \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
  stackArray2d< real64, maxNumFaces *maxNumPhases > potDif( numFacesInElem, numPhases );
  stackArray2d< real64, maxNumFaces *maxNumPhases > dPotDif_dPres( numFacesInElem, numPhases );
  stackArray3d< real64, maxNumFaces *maxNumPhases *maxNumComponents > dPotDif_dComp( numFacesInElem, numPhases, numComponents );
  stackArray2d< real64, maxNumFaces *maxNumPhases > dPotDif_dFacePotential( numFacesInElem, numPhases );
  potDif = 0;
  dPotDif_dPres = 0;
  dPotDif_dComp = 0;
  dPotDif_dFacePotential = 0;

  stackArray1d< real64, maxNumComponents > dPhasePres_dComp( numComponents );
  stackArray1d< real64, maxNumComponents > dPhaseGravTerm_dComp( numComponents );    
  stackArray2d< real64, maxNumPhases* maxNumComponents > dPhaseDens_dC( numPhases, numComponents );

  // 0) precompute dPhaseDens_dC since it is always computed at the element center  
  dPhaseDens_dC = 0;
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    applyChainRule( numComponents,
                    dElemCompFrac_dCompDens,
                    dElemDens_dComp[ip],
                    dPhaseDens_dC[ip] );
  }
  
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
      dPhasePres_dComp = 0; // for capillary pressure
      
      real64 const phaseGravTerm = elemDens[ip] * depthDif;
      real64 const dPhaseGravTerm_dPres = dElemDens_dPres[ip] * depthDif;
      for( localIndex ic = 0; ic < numComponents; ++ic )
      {
        dPhaseGravTerm_dComp[ic] = dPhaseDens_dC[ip][ic] * depthDif;    
      }

      // face phase potential
      real64 const facePhasePotential = facePotential[elemToFaces[ifaceLoc]][ip]
                                      + dFacePotential[elemToFaces[ifaceLoc]][ip];

      // phase potential difference
      potDif[ifaceLoc][ip] = phasePres - phaseGravTerm - facePhasePotential;
      dPotDif_dPres[ifaceLoc][ip] = dPhasePres_dPres - dPhaseGravTerm_dPres;
      dPotDif_dFacePotential[ifaceLoc][ip] = -1;
      for( localIndex ic = 0; ic < numComponents; ++ic )
      {
        dPotDif_dComp[ifaceLoc][ip][ic] = dPhasePres_dComp[ic] - dPhaseGravTerm_dComp[ic];
      }
    }
  }

  // 2) multiply the potential difference by the transmissibility
  //    to obtain the volumetric flux
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {

    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided volumetric flux at face iface
    for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
    {

      for( localIndex ip = 0; ip < numPhases; ++ip )
      {
        // this is going to store T (\nabla p_{\ell} - \rho_p_{\ell} g \nabla d)
        volFlux[ifaceLoc][ip]                      += transMatrix[ifaceLoc][jfaceLoc] * potDif[jfaceLoc][ip];
        dVolFlux_dPres[ifaceLoc][ip]               += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dPres[jfaceLoc][ip];
        dVolFlux_dFacePotential[ifaceLoc][jfaceLoc][ip] += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dFacePotential[jfaceLoc][ip];
         
        for( localIndex ic = 0; ic < numComponents; ++ic )
        {         
          dVolFlux_dComp[ifaceLoc][ip][ic] += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dComp[jfaceLoc][ip][ic];
        }
      }
    }
  }  
}

  
void FluxKernelHelper::UpdateUpwindedCoefficients( array2d< localIndex > const & elemRegionList,
                                                   array2d< localIndex > const & elemSubRegionList,
                                                   array2d< localIndex > const & elemList,
                                                   SortedArray< localIndex > const & regionFilter,
                                                   arraySlice1d< localIndex const > const elemToFaces,
                                                   localIndex const fluidIndex,                                            
                                                   ElementRegionManager::MaterialViewAccessor< arrayView4d< real64 const > >::ViewTypeConst const & phaseCompFrac,
                                                   ElementRegionManager::MaterialViewAccessor< arrayView4d< real64 const > >::ViewTypeConst const & dPhaseCompFrac_dPres,
                                                   ElementRegionManager::MaterialViewAccessor< arrayView5d< real64 const > >::ViewTypeConst const & dPhaseCompFrac_dComp,
                                                   ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >::ViewTypeConst const & dCompFrac_dCompDens,                                            
                                                   ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & phaseMob,
                                                   ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >::ViewTypeConst const & dPhaseMob_dPres,
                                                   ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >::ViewTypeConst const & dPhaseMob_dComp,
                                                   ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst const & elemDofNumber,
                                                   stackArray1d< localIndex, 3 > const & elemIds,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & volFlux,
                                                   stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > & upwPhaseCompFrac,
                                                   stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > & dUpwPhaseCompFrac_dPres,
                                                   stackArray4d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > & dUpwPhaseCompFrac_dComp,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES > & upwPhaseMobility,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES > & dUpwPhaseMobility_dPres,
                                                   stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > & dUpwPhaseMobility_dComp,
                                                   stackArray2d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                             *constitutive::MultiFluidBase::MAX_NUM_PHASES > & upwDofNumber,
                                                   stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > & neighborDofNumber )
{
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;
  localIndex constexpr maxNumComponents = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  
  localIndex const numFacesInElem = elemToFaces.size();
  localIndex const numPhases = dUpwPhaseMobility_dComp.size( 1 );
  localIndex const numComponents = dUpwPhaseMobility_dComp.size( 2 );

  stackArray2d< localIndex, 3 *maxNumFaces > neighborIds( numFacesInElem, 3 );
  neighborIds = -1;

  localIndex const er  = elemIds[0];
  localIndex const esr = elemIds[1];
  localIndex const ei  = elemIds[2];

  stackArray1d< real64, maxNumComponents > dPhaseCompFrac_dC( numComponents );  
  
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
      // if false, the local element is upwind
      localIndex const neighborIsUpw = foundNeighborInTarget && volFlux[ifaceLoc][ip] < 0;

      // we evaluate the mobilities with the upwind element wrt the volumetrix flux
      localIndex const erUpw  = neighborIsUpw * neighborIds[ifaceLoc][0] + (1-neighborIsUpw) * er;
      localIndex const esrUpw = neighborIsUpw * neighborIds[ifaceLoc][1] + (1-neighborIsUpw) * esr;
      localIndex const eiUpw  = neighborIsUpw * neighborIds[ifaceLoc][2] + (1-neighborIsUpw) * ei;
     
      // get the dof number in the upwind element
      globalIndex const dofNumber = neighborIsUpw * neighborDofNumber[ifaceLoc]
                                    + (1-neighborIsUpw) * elemDofNumber[er][esr][ei];

      // save upwinded phase mobilities
      upwPhaseMobility[ifaceLoc][ip] = phaseMob[erUpw][esrUpw][eiUpw][ip];
      dUpwPhaseMobility_dPres[ifaceLoc][ip] = dPhaseMob_dPres[erUpw][esrUpw][eiUpw][ip];
      for( localIndex ic = 0; ic < numComponents; ++ic )
      {
        dUpwPhaseMobility_dComp[ifaceLoc][ip][ic] = dPhaseMob_dComp[erUpw][esrUpw][eiUpw][ip][ic];
      }

      // save upwinded phase component fractions
      for( localIndex ic = 0; ic < numComponents; ++ic )
      {
        upwPhaseCompFrac[ifaceLoc][ip][ic] = phaseCompFrac[erUpw][esrUpw][fluidIndex][eiUpw][0][ip][ic];
        dUpwPhaseCompFrac_dPres[ifaceLoc][ip][ic] = dPhaseCompFrac_dPres[erUpw][esrUpw][fluidIndex][eiUpw][0][ip][ic];

        dPhaseCompFrac_dC = 0;
        applyChainRule( numComponents,
                        dCompFrac_dCompDens[erUpw][esrUpw][eiUpw],
                        dPhaseCompFrac_dComp[erUpw][esrUpw][fluidIndex][eiUpw][0][ip][ic],
                        dPhaseCompFrac_dC );
        
        for( localIndex jc = 0; jc < numComponents; ++jc )
        {
          dUpwPhaseCompFrac_dComp[ifaceLoc][ip][ic][jc] = dPhaseCompFrac_dC[jc]; 
        }
      }
      upwDofNumber[ifaceLoc][ip] = dofNumber;
    }
  }
}

  
void FluxKernelHelper::AssembleOneSidedMassFluxes( real64 const & dt,
                                                   arrayView1d< globalIndex const > const & faceDofNumber,
                                                   arraySlice1d< localIndex const > const elemToFaces,
                                                   globalIndex const elemDofNumber,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & volFlux,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & dVolFlux_dPres,
                                                   stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > const & dVolFlux_dComp,
                                                   stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & dVolFlux_dFacePotential,
                                                   stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > const & upwPhaseCompFrac,
                                                   stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > const & dUpwPhaseCompFrac_dPres,
                                                   stackArray4d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > const & dUpwPhaseCompFrac_dComp,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & upwPhaseMobility,
                                                   stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & dUpwPhaseMobility_dPres,
                                                   stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                        *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > const & dUpwPhaseMobility_dComp,
                                                   stackArray2d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                             *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & upwDofNumber,
                                                   stackArray1d< globalIndex, HybridFVMInnerProduct::MAX_NUM_FACES > const & neighborDofNumber,
                                                   ParallelMatrix * const matrix,
                                                   ParallelVector * const rhs )
{
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;
  localIndex constexpr maxNumPhases = constitutive::MultiFluidBase::MAX_NUM_PHASES;
  localIndex constexpr maxNumComponents = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDofs = maxNumComponents + 1;
  
  localIndex const numFacesInElem = elemToFaces.size();
  localIndex const numPhases = dVolFlux_dComp.size( 1 );
  localIndex const numComponents = dVolFlux_dComp.size( 2 );
  localIndex const numDofs = numComponents + 1;

  // dof numbers
  stackArray1d< globalIndex, maxNumComponents >            eqnRowIndices( numComponents );
  stackArray1d< globalIndex, maxNumDofs *(1+maxNumFaces) > elemDofColIndices( numDofs * (1+numFacesInElem) );
  stackArray1d< globalIndex, maxNumPhases * maxNumFaces >  faceDofColIndices( numPhases * numFacesInElem );
  for( localIndex ic = 0; ic < numComponents; ++ic )
  {
    eqnRowIndices[ic] = elemDofNumber + ic;
  }
  for( localIndex idof = 0; idof < numDofs; ++idof )
  {    
    elemDofColIndices[idof] = elemDofNumber + idof;
  }
  
  // fluxes
  stackArray1d< real64, maxNumComponents >                              sumFluxes( numComponents );
  stackArray2d< real64, maxNumComponents *maxNumDofs *(1+maxNumFaces) > dSumFluxes_dElemVars( numComponents, numDofs * (1+numFacesInElem) );
  stackArray2d< real64, maxNumComponents *maxNumPhases *maxNumFaces >   dSumFluxes_dFaceVars( numComponents, numPhases * numFacesInElem );
  sumFluxes = 0;
  dSumFluxes_dElemVars = 0;
  dSumFluxes_dFaceVars = 0;

  // for each element, loop over the one-sided faces
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {
    localIndex const elemVarsOffset = numDofs*(ifaceLoc+1);

    for( localIndex ip = 0; ip < numPhases; ++ip )
    {

      for( localIndex ic = 0; ic < numComponents; ++ic )
      { 

        real64 const phaseCompCoef = dt * upwPhaseCompFrac[ifaceLoc][ip][ic] * upwPhaseMobility[ifaceLoc][ip];
        real64 const dPhaseCompCoef_dPres = dt * ( dUpwPhaseCompFrac_dPres[ifaceLoc][ip][ic] * upwPhaseMobility[ifaceLoc][ip]
                                                 + upwPhaseCompFrac[ifaceLoc][ip][ic] * dUpwPhaseMobility_dPres[ifaceLoc][ip] );
         
        // 1) residual
        sumFluxes[ic] += phaseCompCoef * volFlux[ifaceLoc][ip];

        // 2) derivatives wrt the elem centered vars of the local elem
        dSumFluxes_dElemVars[ic][0] += phaseCompCoef * dVolFlux_dPres[ifaceLoc][ip];
        for( localIndex jc = 0; jc < numComponents; ++jc )
        {  
          dSumFluxes_dElemVars[ic][jc+1] += phaseCompCoef * dVolFlux_dComp[ifaceLoc][ip][jc];
        }

        // 3) derivatives of the upwinded mobilities
        if( upwDofNumber[ifaceLoc][ip] == elemDofNumber )
        {
          dSumFluxes_dElemVars[ic][0] += dPhaseCompCoef_dPres * volFlux[ifaceLoc][ip];
          
          for( localIndex jc = 0; jc < numComponents; ++jc )
          {
            real64 const dPhaseCompCoef_dComp = dt * ( dUpwPhaseCompFrac_dComp[ifaceLoc][ip][ic][jc] * upwPhaseMobility[ifaceLoc][ip]
                                                     + upwPhaseCompFrac[ifaceLoc][ip][ic] * dUpwPhaseMobility_dComp[ifaceLoc][ip][jc] );
            dSumFluxes_dElemVars[ic][jc+1] += dPhaseCompCoef_dComp * volFlux[ifaceLoc][ip];
          }
        }
        else
        {
          dSumFluxes_dElemVars[ic][elemVarsOffset] = dPhaseCompCoef_dPres * volFlux[ifaceLoc][ip];
          
          for( localIndex jc = 0; jc < numComponents; ++jc )
          {
            real64 const dPhaseCompCoef_dComp = dt * ( dUpwPhaseCompFrac_dComp[ifaceLoc][ip][ic][jc] * upwPhaseMobility[ifaceLoc][ip]
                                                     + upwPhaseCompFrac[ifaceLoc][ip][ic] * dUpwPhaseMobility_dComp[ifaceLoc][ip][jc] );
            dSumFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dPhaseCompCoef_dComp * volFlux[ifaceLoc][ip];
          }
        }

        // 4) derivatives wrt the face centered var
        for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
        {
          localIndex const faceVarsOffset = numPhases*jfaceLoc;
          dSumFluxes_dFaceVars[ic][faceVarsOffset + ip] += phaseCompCoef * dVolFlux_dFacePotential[ifaceLoc][jfaceLoc][ip];
        }
      }
    }
      
    // collect the relevant dof numbers
    for( localIndex idof = 0; idof < numDofs; ++idof )
    {    
      elemDofColIndices[elemVarsOffset+idof] = neighborDofNumber[ifaceLoc] + idof;
    }

    localIndex const faceVarsOffset = numPhases*ifaceLoc;
    for( localIndex ip = 0; ip < numPhases; ++ip )
    {
      faceDofColIndices[faceVarsOffset+ip] = faceDofNumber[elemToFaces[ifaceLoc]] + ip;
    }
  }

  // we are ready to assemble the local flux and its derivatives

  // residual
  rhs->add( eqnRowIndices,
            sumFluxes );

  // jacobian -- derivative wrt elem centered vars
  matrix->add( eqnRowIndices,
               elemDofColIndices,
               dSumFluxes_dElemVars );

  // jacobian -- derivatives wrt face centered vars
  matrix->add( eqnRowIndices,
               faceDofColIndices,
               dSumFluxes_dFaceVars );  
}
  
void FluxKernelHelper::AssembleConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                                            arraySlice1d< localIndex const > const elemToFaces,
                                            globalIndex const elemDofNumber,
                                            stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                 *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & volFlux,
                                            stackArray2d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                 *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & dVolFlux_dPres,
                                            stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                 *constitutive::MultiFluidBase::MAX_NUM_PHASES
                                                                 *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS > const & dVolFlux_dComp,
                                            stackArray3d< real64, HybridFVMInnerProduct::MAX_NUM_FACES
                                                                 *HybridFVMInnerProduct::MAX_NUM_FACES
                                                                 *constitutive::MultiFluidBase::MAX_NUM_PHASES > const & dVolFlux_dFacePotential,
                                            ParallelMatrix * const matrix,
                                            ParallelVector * const rhs )
{
  localIndex constexpr maxNumFaces = HybridFVMInnerProduct::MAX_NUM_FACES;
  localIndex constexpr maxNumPhases = constitutive::MultiFluidBase::MAX_NUM_PHASES;
  localIndex constexpr maxNumComponents = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDofs = maxNumComponents + 1;
  
  localIndex const numFacesInElem = elemToFaces.size();
  localIndex const numPhases = dVolFlux_dComp.size( 1 );
  localIndex const numComponents = dVolFlux_dComp.size( 2 );
  localIndex const numDofs = numComponents + 1;

  // dof numbers
  stackArray1d< globalIndex, maxNumPhases > eqnRowIndices( numPhases );
  stackArray1d< globalIndex, maxNumDofs > elemDofColIndices( numDofs );
  stackArray1d< globalIndex, maxNumPhases *maxNumFaces > faceDofColIndices( numPhases*numFacesInElem );

  for( localIndex idof = 0; idof < numDofs; ++idof )
  {    
    elemDofColIndices[idof] = elemDofNumber + idof;
  }

  // fluxes
  stackArray1d< real64, maxNumPhases > flux( numPhases );
  stackArray2d< real64, maxNumPhases *maxNumDofs > dFlux_dElemVars( numPhases, numDofs );
  stackArray2d< real64, maxNumPhases *maxNumPhases *maxNumFaces > dFlux_dFaceVars( numPhases, numPhases*numFacesInElem );

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
      dFlux_dElemVars[ip][0] = dVolFlux_dPres[ifaceLoc][ip];
      for( localIndex ic = 0; ic < numComponents; ++ic )
      {
        dFlux_dElemVars[ip][ic+1] = dVolFlux_dComp[ifaceLoc][ip][ic];
      }
      
      for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
      {
        localIndex const faceVarsOffset = numPhases*jfaceLoc;
        faceDofColIndices[faceVarsOffset+ip] = faceDofNumber[elemToFaces[jfaceLoc]] + ip;
        dFlux_dFaceVars[ip][faceVarsOffset+ip] = dVolFlux_dFacePotential[ifaceLoc][jfaceLoc][ip];
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


} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx
