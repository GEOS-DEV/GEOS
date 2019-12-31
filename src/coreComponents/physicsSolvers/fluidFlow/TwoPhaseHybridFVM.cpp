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
 * @file TwoPhaseHybridFVM.cpp
 */

#include "TwoPhaseHybridFVM.hpp"

#include "common/TimingMacros.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;


TwoPhaseHybridFVM::TwoPhaseHybridFVM( const std::string& name,
                                      Group * const parent ):
  TwoPhaseBase(name, parent),
  m_faceDofKey(""),
  m_areaRelTol(1e-8),
  m_minTotalMob(1e-8)
{
  // two cell-centered dof per cell
  m_numDofPerCell = 2;
}

  
void TwoPhaseHybridFVM::RegisterDataOnMesh(Group * const MeshBodies)
{
  // 1) Register the cell-centered data
  TwoPhaseBase::RegisterDataOnMesh(MeshBodies);

  // 2) Register the face data 
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    FaceManager * const faceManager = meshLevel->getFaceManager();

    // primary variables: face pressures
    faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::facePressureString )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the pressures at the faces.");

    faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::deltaFacePressureString )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the accumulated pressure updates at the faces.");
   
  }
}
  
void TwoPhaseHybridFVM::ImplicitStepSetup( real64 const & time_n,
                                         real64 const & dt,
                                         DomainPartition * const domain,
                                         DofManager & dofManager,
                                         ParallelMatrix & matrix,
                                         ParallelVector & rhs,
                                         ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  // setup the cell-centered fields
  TwoPhaseBase::ImplicitStepSetup( time_n, dt, domain, dofManager, matrix, rhs, solution );

  // setup the face fields
  MeshLevel * const meshLevel     = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = meshLevel->getFaceManager();

  // get the accumulated pressure updates
  arrayView1d<real64> & dFacePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);
 
  forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  {
    dFacePres[iface] = 0;
  });
  
}

void TwoPhaseHybridFVM::ImplicitStepComplete( real64 const & time_n,
                                              real64 const & dt,
                                              DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  // increment the cell-centered fields
  TwoPhaseBase::ImplicitStepComplete( time_n, dt, domain );

  // increment the face fields
  MeshLevel * const meshLevel     = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = meshLevel->getFaceManager();

  // get the face-based DOF numbers 
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference< array1d<globalIndex> >( m_faceDofKey );

  // get the face-based pressures 
  arrayView1d<real64> const & facePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::facePressureString);
  arrayView1d<real64 const> const & dFacePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

  forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  {
    // update if face is in target region
    if (faceDofNumber[iface] >= 0)
    {
      facePres[iface] += dFacePres[iface];
    }
  });
}

void TwoPhaseHybridFVM::SetupDofs( DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                                   DofManager & dofManager ) const
{
  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleUpwindedOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::elemDofFieldString,
                       DofManager::Location::Elem,
                       DofManager::Connectivity::Face, 
                       m_targetRegions );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::faceDofFieldString,
                       DofManager::Location::Face,
                       DofManager::Connectivity::Elem,
                       m_targetRegions );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::faceDofFieldString,
                          viewKeyStruct::elemDofFieldString,
                          DofManager::Connectivity::Elem,
                          true );

 
}

void TwoPhaseHybridFVM::AssembleFluxTerms( real64 const GEOSX_UNUSED_ARG( time_n ),
                                           real64 const dt,
                                           DomainPartition const * const domain,
                                           DofManager const * const dofManager,
                                           ParallelMatrix * const matrix,
                                           ParallelVector * const rhs )
{
  MeshLevel const * const mesh                   = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  NodeManager const * const nodeManager          = mesh->getNodeManager();
  FaceManager const * const faceManager          = mesh->getFaceManager();   

  // in this function we need to make sure that we act only on the target regions
  // for that, we need the following region filter
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }
  
  // node data (for transmissibility computation)

  arrayView1d<R1Tensor const> const & nodePosition = nodeManager->referencePosition();

  
  // face data

  // get the face-based DOF numbers for the assembly
  string const faceDofKey = dofManager->getKey( viewKeyStruct::faceDofFieldString );
  // save the face Dof key for use in two functions
  // that do not have acces to the coupled solver dofManager
  // namely ResetStateToBeginningOfStep and ImplicitStepComplete
  m_faceDofKey = faceDofKey;
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference< array1d<globalIndex> >( faceDofKey );  

  // get the face-centered pressures
  arrayView1d<real64 const> const & facePres =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::facePressureString );
  arrayView1d<real64 const> const & dFacePres =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::deltaFacePressureString );

  // get the face-centered depth 
  arrayView1d<real64> const & faceGravDepth =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::gravityDepthString);
  
  // get the face-to-nodes connectivity for the transmissibility calculation
  ArrayOfArraysView<localIndex const> const & faceToNodes = faceManager->nodeList();

  array2d<localIndex> const & elemRegionList    = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList          = faceManager->elementList(); 

  
  // max number of faces allowed in an element 
  localIndex constexpr maxNumFaces = MAX_NUM_FACES_IN_ELEM;
  localIndex constexpr numPhases   = NUM_PHASES;
  
  // tolerance for transmissibility calculation
  real64 const lengthTolerance = domain->getMeshBody( 0 )->getGlobalLengthScale() * m_areaRelTol; 


  elemManager->
    forElementSubRegionsComplete<CellElementSubRegion,
                                 FaceElementSubRegion>( m_targetRegions,
                                                      [&]( localIndex const er,
                                                           localIndex const esr,
                                                           ElementRegionBase const * const,
                                                           auto const * const subRegion )
  {

    // elem data
    
    // get the cell-centered DOF numbers and ghost rank for the assembly
    string const elemDofKey = dofManager->getKey( viewKeyStruct::elemDofFieldString );
    arrayView1d<globalIndex const> const & elemDofNumber =
      subRegion->template getReference< array1d<globalIndex> >( elemDofKey ); 
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];   

    // get the map from elem to faces
    arrayView2d< localIndex const > const & elemToFaces = subRegion->faceList();
    
    // get the cell-centered pressures
    arrayView1d<real64 const> const & elemPres  = m_pressure[er][esr];
    arrayView1d<real64 const> const & dElemPres = m_deltaPressure[er][esr];

    // get the cell centered densities and mobilities
    arrayView3d<real64 const> const & elemPhaseDens     = m_phaseDens[er][esr][m_fluidIndex];
    arrayView3d<real64 const> const & dElemPhaseDens_dp = m_dPhaseDens_dPres[er][esr][m_fluidIndex];

    arrayView2d<real64 const> const & elemPhaseMob     = m_phaseMob[er][esr];
    arrayView2d<real64 const> const & dElemPhaseMob_dp = m_dPhaseMob_dPres[er][esr];
    arrayView2d<real64 const> const & dElemPhaseMob_dS = m_dPhaseMob_dSat[er][esr];    
    
    // get the element data needed for transmissibility computation
    arrayView1d<R1Tensor const> const & elemCenter =
     subRegion->template getReference< array1d<R1Tensor> >( CellBlock::viewKeyStruct::elementCenterString );
    arrayView1d<real64 const> const & elemVolume =
     subRegion->template getReference< array1d<real64> >( CellBlock::viewKeyStruct::elementVolumeString ); 
    arrayView1d<R1Tensor const> const & elemPerm =
     subRegion->template getReference< array1d<R1Tensor> >( viewKeyStruct::permeabilityString ); 

    // get the cell-centered depth 
    arrayView1d<real64> const & elemGravDepth =
      subRegion->template getReference<array1d<real64>>(viewKeyStruct::gravityDepthString);

    
    // assemble the residual and Jacobian element by element
    // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {

      if (elemGhostRank[ei] < 0)
      {
        
        localIndex const numFacesInElem = elemToFaces.size(1);

        // transmissibility matrix
        stackArray2d<real64, maxNumFaces*maxNumFaces> transMatrix( numFacesInElem,
                                                                   numFacesInElem );
      
        // one sided flux
        stackArray1d<real64, maxNumFaces> oneSidedTotalVolFlux( numFacesInElem );
        stackArray1d<real64, maxNumFaces> dOneSidedTotalVolFlux_dp( numFacesInElem );
        stackArray1d<real64, maxNumFaces> dOneSidedTotalVolFlux_dS( numFacesInElem );
        stackArray1d<real64, maxNumFaces> dOneSidedTotalVolFlux_dfp( numFacesInElem );
      
        // upwinded coefficient for the viscous term
        stackArray2d<real64, numPhases*maxNumFaces> viscousCoef( numFacesInElem, numPhases  );
        stackArray2d<real64, numPhases*maxNumFaces> dViscousCoef_dp( numFacesInElem, numPhases );
        stackArray2d<real64, numPhases*maxNumFaces> dViscousCoef_dS( numFacesInElem, numPhases );
        stackArray1d<globalIndex, maxNumFaces> viscousDofNumber( numFacesInElem ); 

        // upwinded coefficient for the buoyancy term
        stackArray2d<real64, numPhases*maxNumFaces> buoyancyCoef( numFacesInElem, numPhases );
        stackArray3d<real64, 2*numPhases*maxNumFaces> dBuoyancyCoef_dp( numFacesInElem, numPhases, 2 );
        stackArray3d<real64, 2*numPhases*maxNumFaces> dBuoyancyCoef_dS( numFacesInElem, numPhases, 2 );
        stackArray2d<globalIndex, 2*maxNumFaces> buoyancyDofNumber( numFacesInElem, 2 );      
        
        // recompute the local transmissibility matrix at each iteration
        // we can decide later to precompute transMatrix if needed
        ComputeTransmissibilityMatrix( nodePosition,     
                                       faceToNodes,        
                                       elemToFaces[ei],
                                       elemCenter[ei],
                                       elemVolume[ei],
                                       elemPerm[ei],
                                       lengthTolerance,
                                       transMatrix );   

        
        /*
         * compute auxiliary quantities at the one sided faces of this element:
         * 1) One-sided volumetric fluxes
         * 2) Upwinded mobilities
         */     
 
        // for each one-sided face of the elem,
        // compute the volumetric flux using transMatrix
        ComputeOneSidedVolFluxes( facePres,
                                  dFacePres,
                                  faceGravDepth,
                                  elemToFaces[ei],          
                                  elemPres[ei],
                                  dElemPres[ei],
                                  elemGravDepth[ei],
                                  elemPhaseMob[ei],
                                  dElemPhaseMob_dp[ei],
                                  dElemPhaseMob_dS[ei],
                                  elemPhaseDens[ei][0],
                                  dElemPhaseDens_dp[ei][0],
                                  transMatrix,
                                  oneSidedTotalVolFlux,             
                                  dOneSidedTotalVolFlux_dp,
                                  dOneSidedTotalVolFlux_dS,
                                  dOneSidedTotalVolFlux_dfp );

        // at this point, we know the local flow direction in the element
        // so we can upwind the transport coefficients (mobilities) at the one sided faces 
        // ** this function needs non-local information **
        UpdateUpwindedCoefficients( mesh,
                                    elemRegionList, 
                                    elemSubRegionList,
                                    elemList,
                                    regionFilter,
                                    elemToFaces[ei],
                                    m_phaseMob, 
                                    m_dPhaseMob_dPres,
                                    m_dPhaseMob_dSat,
                                    m_phaseDens, 
                                    m_dPhaseDens_dPres,
                                    er, esr, ei,
                                    elemDofNumber[ei],
                                    elemDofKey,
                                    oneSidedTotalVolFlux,
                                    viscousCoef,      
                                    dViscousCoef_dp,
                                    dViscousCoef_dS,
                                    viscousDofNumber,
                                    buoyancyCoef,      
                                    dBuoyancyCoef_dp,
                                    dBuoyancyCoef_dS,
                                    buoyancyDofNumber );
        
        /*
         * perform assembly in this element in two steps:
         * 1) mass conservation equations
         * 2) face constraints
         */
        
        // use the computed one sided vol fluxes and the upwinded mobilities
        // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
        AssembleOneSidedMassFluxes( dt,
                                    faceDofNumber,
                                    elemToFaces[ei],
                                    elemDofNumber[ei],
                                    oneSidedTotalVolFlux,
                                    dOneSidedTotalVolFlux_dp,
                                    dOneSidedTotalVolFlux_dS,
                                    dOneSidedTotalVolFlux_dfp,
                                    viscousCoef,      
                                    dViscousCoef_dp,
                                    dViscousCoef_dS,
                                    viscousDofNumber,
                                    buoyancyCoef,
                                    dBuoyancyCoef_dp,
                                    dBuoyancyCoef_dS,
                                    buoyancyDofNumber,
                                    matrix,
                                    rhs );

        // use the computed one sided vol fluxes to assemble the constraints
        // enforcing flux continuity at this element's faces
        AssembleConstraints( faceDofNumber,
                             elemToFaces[ei],
                             elemDofNumber[ei],
                             oneSidedTotalVolFlux,
                             dOneSidedTotalVolFlux_dp,
                             dOneSidedTotalVolFlux_dS,                       
                             dOneSidedTotalVolFlux_dfp,
                             matrix,
                             rhs );

      }
    });
  });
}


void TwoPhaseHybridFVM::ComputeOneSidedVolFluxes( arrayView1d<real64 const> const & facePres,
                                                  arrayView1d<real64 const> const & dFacePres,
                                                  arrayView1d<real64 const> const & faceGravDepth,
                                                  arraySlice1d<localIndex const> const elemToFaces,
                                                  real64 const & elemPres,
                                                  real64 const & dElemPres,
                                                  real64 const & elemGravDepth,
                                                  arraySlice1d<real64 const> const elemPhaseMob,
                                                  arraySlice1d<real64 const> const dElemPhaseMob_dp,
                                                  arraySlice1d<real64 const> const dElemPhaseMob_dS,
                                                  arraySlice1d<real64 const> const elemPhaseDens,
                                                  arraySlice1d<real64 const> const dElemPhaseDens_dp,
                                                  stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                                      *MAX_NUM_FACES_IN_ELEM> const & transMatrix,
                                                  stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & oneSidedVolFlux,
                                                  stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dOneSidedVolFlux_dp,
                                                  stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dOneSidedVolFlux_dS,
                                                  stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dOneSidedVolFlux_dfp ) const
{
  localIndex constexpr numPhases  = NUM_PHASES;  
  localIndex const numFacesInElem = elemToFaces.size();
  
  // for this element, loop over the local (one-sided) faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    oneSidedVolFlux[ifaceLoc]      = 0;
    dOneSidedVolFlux_dp[ifaceLoc]  = 0;
    dOneSidedVolFlux_dfp[ifaceLoc] = 0;
         
    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided flux at face iface
    for (localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc)
    {

      // pressure and gravity terms
      // TODO: account for capillary pressure here
      real64 const ccPres = elemPres + dElemPres;
      real64 const fPres  = facePres[elemToFaces[jfaceLoc]] + dFacePres[elemToFaces[jfaceLoc]];

      real64 const ccGravDepth = elemGravDepth;
      real64 const fGravDepth  = faceGravDepth[elemToFaces[jfaceLoc]];

      // pressure difference
      real64 const presDif      = ccPres - fPres;
      real64 const dPresDif_dp  =  1;
      real64 const dPresDif_dfp = -1;

      // gravity term
      real64 const depthDif = ccGravDepth - fGravDepth;

      real64 sumWeightedPotDif      = 0;
      real64 dSumWeightedPotDif_dp  = 0;
      real64 dSumWeightedPotDif_dS  = 0;
      real64 dSumWeightedPotDif_dfp = 0;
      
      for (localIndex ip = 0; ip < numPhases; ++ip)
      {
        real64 const phaseGravTerm     = elemPhaseDens[ip]     * depthDif;
        real64 const dPhaseGravTerm_dp = dElemPhaseDens_dp[ip] * depthDif;

        // 1) compute the potential diff between the cell center and the face center
        real64 const phasePotDif       = presDif     - phaseGravTerm;
        real64 const dPhasePotDif_dp   = dPresDif_dp - dPhaseGravTerm_dp;
        real64 const dPhasePotDif_dfp  = dPresDif_dfp;

        // sum of potential differences weighted by phase mobility
        sumWeightedPotDif       += elemPhaseMob[ip]     * phasePotDif;
        dSumWeightedPotDif_dp   += dElemPhaseMob_dp[ip] * phasePotDif
                                 + elemPhaseMob[ip]     * dPhasePotDif_dp;
        dSumWeightedPotDif_dS   += dElemPhaseMob_dS[ip] * phasePotDif;
        dSumWeightedPotDif_dfp  += elemPhaseMob[ip]     * dPhasePotDif_dfp;
      }
      
      // 2) compute the contribution of this face to the volumetric fluxes in the cell
      oneSidedVolFlux[ifaceLoc]      += transMatrix[ifaceLoc][jfaceLoc] * sumWeightedPotDif;
      dOneSidedVolFlux_dp[ifaceLoc]  += transMatrix[ifaceLoc][jfaceLoc] * dSumWeightedPotDif_dp;
      dOneSidedVolFlux_dS[ifaceLoc]  += transMatrix[ifaceLoc][jfaceLoc] * dSumWeightedPotDif_dS;      
      dOneSidedVolFlux_dfp[ifaceLoc] += transMatrix[ifaceLoc][jfaceLoc] * dSumWeightedPotDif_dfp;
    }
  }    
}


void TwoPhaseHybridFVM::UpdateUpwindedCoefficients( MeshLevel const * const  mesh,
                                                    array2d<localIndex> const & elemRegionList,
                                                    array2d<localIndex> const & elemSubRegionList,
                                                    array2d<localIndex> const & elemList,
                                                    set<localIndex> const & regionFilter,
                                                    arraySlice1d<localIndex const> const elemToFaces,
                                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & domainMob,
                                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dDomainMob_dp,
                                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dDomainMob_dS,
                                                    ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & domainDens,
                                                    ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dDomainDens_dp,
                                                    localIndex const er,
                                                    localIndex const esr,
                                                    localIndex const ei,
                                                    globalIndex const elemDofNumber,
                                                    string const elemDofKey,
                                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES_IN_ELEM> & viscousCoef,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES_IN_ELEM> & dViscousCoef_dp,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES_IN_ELEM> & dViscousCoef_dS,
                                                    stackArray1d<globalIndex, MAX_NUM_FACES_IN_ELEM> & viscousDofNumber,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES_IN_ELEM>   & buoyancyCoef,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES_IN_ELEM> & dBuoyancyCoef_dp,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES_IN_ELEM> & dBuoyancyCoef_dS,
                                                    stackArray2d<globalIndex, 2*MAX_NUM_FACES_IN_ELEM> & buoyancyDofNumber ) const
{
  localIndex const numFacesInElem = elemToFaces.size();
  localIndex constexpr numPhases  = NUM_PHASES;
  
  // for this element, loop over the local (one-sided) faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {

    // we initialize these upw quantities with the values of the local elem    
    localIndex erUpw  = er;
    localIndex esrUpw = esr;
    localIndex eiUpw  = ei;
    viscousDofNumber[ifaceLoc] = elemDofNumber;

    // if the local elem if upstream, we are done, we can proceed to the next one-sided face
    // otherwise, we have to access the properties of the neighbor element
    // this is done on the fly below
    if (oneSidedVolFlux[ifaceLoc] < 0)
    {
    
      // the face has at most two adjacent elements
      // one of these two elements is the current element indexed by er, esr, ei
      // but here we are interested in the indices of the other element
      // this other element is "the neighbor" for this one-sided face 
      for (localIndex k=0; k<elemRegionList.size(1); ++k)
      {
        localIndex const erNeighbor  = elemRegionList[elemToFaces[ifaceLoc]][k];
        localIndex const esrNeighbor = elemSubRegionList[elemToFaces[ifaceLoc]][k];
        localIndex const eiNeighbor  = elemList[elemToFaces[ifaceLoc]][k];
                                                    
        // this element is not the current element
        // we have found the neighbor or we are at the boundary 
        if ( erNeighbor != er || esrNeighbor != esr || eiNeighbor != ei ) 
        {

          bool const onBoundary       = (erNeighbor == -1 || esrNeighbor == -1 || eiNeighbor != -1);
          bool const neighborInTarget = regionFilter.contains(erNeighbor); 

          // if not on boundary, save the mobility and the upwDofNumber  
          if ( !onBoundary && neighborInTarget )
          {
            ElementRegionBase const * const neighborRegion =
              Group::group_cast<ElementRegionBase const *>(mesh->getElemManager()->GetRegion(erNeighbor));
            ElementSubRegionBase const * const neighborSubRegion =
              Group::group_cast<ElementSubRegionBase const *>(neighborRegion->GetSubRegion(esrNeighbor));

            arrayView1d<globalIndex const> const & neighborDofNumber =
              neighborSubRegion->getReference<array1d<globalIndex>>( elemDofKey );  

            erUpw  = erNeighbor;
            esrUpw = esrNeighbor;
            eiUpw  = eiNeighbor;
            viscousDofNumber[ifaceLoc] = neighborDofNumber[eiNeighbor];
            
          }
          // if the face is on the boundary, use the properties of the local elem
        }
      }
    }

    // upwinded coefficient for the viscous terms

    // total mobility
    real64 totalMob     = 0;
    real64 dTotalMob_dp = 0;
    real64 dTotalMob_dS = 0;

    for (localIndex ip = 0; ip < numPhases; ++ip)
    {
      totalMob     += domainMob[erUpw][esrUpw][eiUpw][ip];
      dTotalMob_dp += dDomainMob_dp[erUpw][esrUpw][eiUpw][ip];
      dTotalMob_dS += dDomainMob_dS[erUpw][esrUpw][eiUpw][ip];
    }
    if (totalMob < m_minTotalMob)
    {
      totalMob     = m_minTotalMob;
      dTotalMob_dp = 0;
      dTotalMob_dS = 0;
    }
    
    for (localIndex ip = 0; ip < numPhases; ++ip)
    {
      // upwinded mobilities
      real64 const mob     = domainMob[erUpw][esrUpw][eiUpw][ip];
      real64 const dMob_dp = dDomainMob_dp[erUpw][esrUpw][eiUpw][ip];
      real64 const dMob_dS = dDomainMob_dS[erUpw][esrUpw][eiUpw][ip];        

      // upwinded mobility ratios
      real64 const mobRatio     = mob / totalMob;
      real64 const dMobRatio_dp = (dMob_dp * totalMob - mob * dTotalMob_dp)
                                / (totalMob*totalMob);
      real64 const dMobRatio_dS = (dMob_dS * totalMob - mob * dTotalMob_dS)
                                / (totalMob*totalMob);
      
      // upwinded densities
      real64 const dens     = domainDens[erUpw][esrUpw][m_fluidIndex][eiUpw][0][ip];
      real64 const dDens_dp = dDomainDens_dp[erUpw][esrUpw][m_fluidIndex][eiUpw][0][ip];

      viscousCoef[ifaceLoc][ip]     = dens * mobRatio;
      dViscousCoef_dp[ifaceLoc][ip] = dDens_dp * mobRatio + dens * dMobRatio_dp;
      dViscousCoef_dS[ifaceLoc][ip] = dens * dMobRatio_dS;
    }

    // TODO: add the buoyancy coefficients
    
    // upwinded coefficient for the buoyancy terms
    for (localIndex ip = 0; ip < numPhases; ++ip)
    {
      buoyancyCoef[ifaceLoc][ip] = 0.0;
      for (localIndex nei = 0; nei < 2; ++nei)
      {
        dBuoyancyCoef_dp[ifaceLoc][ip][nei] = 0.0;
        dBuoyancyCoef_dS[ifaceLoc][ip][nei] = 0.0;
        buoyancyDofNumber[ifaceLoc][nei]    = -1;
      }
    }
  }
}


void TwoPhaseHybridFVM::AssembleOneSidedMassFluxes( real64 const & dt,
                                                    arrayView1d<globalIndex const> const & faceDofNumber,
                                                    arraySlice1d<localIndex const> const elemToFaces,
                                                    globalIndex const elemDofNumber,
                                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dp,
                                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dS,
                                                    stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dfp,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES_IN_ELEM> const & viscousCoef,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES_IN_ELEM> const & dViscousCoef_dp,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES_IN_ELEM> const & dViscousCoef_dS,
                                                    stackArray1d<globalIndex, MAX_NUM_FACES_IN_ELEM> const & viscousDofNumber,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES_IN_ELEM>   const & GEOSX_UNUSED_ARG( buoyancyCoef ),
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES_IN_ELEM> const & GEOSX_UNUSED_ARG( dBuoyancyCoef_dp ),
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES_IN_ELEM> const & GEOSX_UNUSED_ARG( dBuoyancyCoef_dS ),
                                                    stackArray2d<globalIndex, 2*MAX_NUM_FACES_IN_ELEM> const & GEOSX_UNUSED_ARG( buoyancyDofNumber ),
                                                    ParallelMatrix * const matrix,
                                                    ParallelVector * const rhs ) const
{
  localIndex constexpr maxNumFaces = MAX_NUM_FACES_IN_ELEM;
  localIndex constexpr numPhases   = NUM_PHASES;
  localIndex constexpr numDof      = NUM_DOF;

  localIndex const numFacesInElem = elemToFaces.size();

  localIndex const w  = TwoPhaseBase::RowOffset::WETTING;
  localIndex const nw = TwoPhaseBase::RowOffset::NONWETTING;
  localIndex const dp = TwoPhaseBase::ColOffset::DPRES;
  localIndex const dS = TwoPhaseBase::ColOffset::DSAT;

  // dof numbers
  stackArray1d<globalIndex, numPhases>              eqnRowIndices( numPhases );
  stackArray1d<globalIndex, numDof*(1+maxNumFaces)> elemDofColIndices( numDof * (1+numFacesInElem) );
  stackArray1d<globalIndex, maxNumFaces>            faceDofColIndices( numFacesInElem );
  eqnRowIndices[w]      = elemDofNumber + w;
  eqnRowIndices[nw]     = elemDofNumber + nw;
  elemDofColIndices[dp] = elemDofNumber + dp;
  elemDofColIndices[dS] = elemDofNumber + dS;

  // fluxes
  stackArray1d<real64, numPhases> sumLocalMassFluxes( numPhases );
  stackArray2d<real64, numPhases*numDof*(1+maxNumFaces)> dSumLocalMassFluxes_dElemVars( numPhases, numDof * (1+numFacesInElem) );
  stackArray2d<real64, numPhases*maxNumFaces> dSumLocalMassFluxes_dFaceVars( numPhases, numFacesInElem );
  sumLocalMassFluxes = 0;
  dSumLocalMassFluxes_dElemVars = 0;
  dSumLocalMassFluxes_dFaceVars = 0;
  
  // TODO: add buoyancy terms
  
  // for each element, loop over the one-sided faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    localIndex const fOffset = numDof*(ifaceLoc+1);
    for (localIndex ip = 0; ip < numPhases; ++ip)
    {
      localIndex const rowId = m_phaseToRow[ip];
      
      // compute the mass flux at the one-sided face plus its derivatives
      // add the newly computed flux to the sum
      real64 const coef     = dt * viscousCoef[ifaceLoc][ip];
      real64 const dCoef_dp = dt * dViscousCoef_dp[ifaceLoc][ip];
      real64 const dCoef_dS = dt * dViscousCoef_dS[ifaceLoc][ip];

      // residual 
      sumLocalMassFluxes[rowId]  += coef * oneSidedVolFlux[ifaceLoc];
 
      // derivatives wrt the cell centered vars of the local elem
      dSumLocalMassFluxes_dElemVars[rowId][dp] += coef * dOneSidedVolFlux_dp[ifaceLoc];
      dSumLocalMassFluxes_dElemVars[rowId][dS] += coef * dOneSidedVolFlux_dS[ifaceLoc];

      // derivatives wrt the cell centered vars of the neighbor
      dSumLocalMassFluxes_dElemVars[rowId][fOffset+dp] = dCoef_dp * oneSidedVolFlux[ifaceLoc];
      dSumLocalMassFluxes_dElemVars[rowId][fOffset+dS] = dCoef_dS * oneSidedVolFlux[ifaceLoc];

      // derivatives wrt the face centered var 
      dSumLocalMassFluxes_dFaceVars[rowId][ifaceLoc] = coef * dOneSidedVolFlux_dfp[ifaceLoc];
    }

    // collect the relevant dof numbers
    elemDofColIndices[fOffset+dp] = viscousDofNumber[ifaceLoc] + dp;
    elemDofColIndices[fOffset+dS] = viscousDofNumber[ifaceLoc] + dS;
    faceDofColIndices[ifaceLoc]   = faceDofNumber[elemToFaces[ifaceLoc]];
    
  }
  
  // we are ready to assemble the local flux and its derivatives

  // residual
  rhs->add( eqnRowIndices.data(),
            sumLocalMassFluxes.data(),
            numPhases );

  // jacobian -- derivative wrt elem centered vars
  matrix->add( eqnRowIndices.data(),
               elemDofColIndices.data(),
               dSumLocalMassFluxes_dElemVars.data(),
               numPhases,
               numDof * (1+numFacesInElem) );

  // jacobian -- derivatives wrt face centered vars
  matrix->add( eqnRowIndices.data(),
               faceDofColIndices.data(),
               dSumLocalMassFluxes_dFaceVars.data(),
               numPhases,
               numFacesInElem );
}


void TwoPhaseHybridFVM::AssembleConstraints( arrayView1d<globalIndex const> const & faceDofNumber,
                                             arraySlice1d<localIndex const> const elemToFaces,
                                             globalIndex const elemDofNumber,
                                             stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                             stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dp,
                                             stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dS,
                                             stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dfp,
                                             ParallelMatrix * const matrix,
                                             ParallelVector * const rhs ) const 
{
  localIndex constexpr numDof     = NUM_DOF;
  localIndex const numFacesInElem = elemToFaces.size();

  localIndex const dp = TwoPhaseBase::ColOffset::DPRES;
  localIndex const dS = TwoPhaseBase::ColOffset::DSAT;

  // dof numbers
  stackArray1d<globalIndex, numDof> elemDofColIndices( numDof );
  elemDofColIndices[dp] = elemDofNumber + dp;
  elemDofColIndices[dS] = elemDofNumber + dS;

  // fluxes
  stackArray1d<real64, numDof> dFlux_dElemVars( numDof );
  
  // for each element, loop over the local (one-sided) faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    // dof numbers 
    globalIndex const eqnRowIndex     = faceDofNumber[elemToFaces[ifaceLoc]];
    globalIndex const faceDofColIndex = faceDofNumber[elemToFaces[ifaceLoc]];

    // fluxes
    real64 const flux            = oneSidedVolFlux[ifaceLoc];
    real64 const dFlux_dFaceVars = dOneSidedVolFlux_dfp[ifaceLoc];
    dFlux_dElemVars[dp]          = dOneSidedVolFlux_dp[ifaceLoc];
    dFlux_dElemVars[dS]          = dOneSidedVolFlux_dS[ifaceLoc];

    // residual
    rhs->add( &eqnRowIndex,
              &flux,
              1 );

    // jacobian -- derivative wrt local cell centered vars
    matrix->add( &eqnRowIndex,
                 elemDofColIndices.data(),
                 dFlux_dElemVars.data(),
                 1, numDof );
   
    // jacobian -- derivatives wrt face pressure terms
    matrix->add( &eqnRowIndex,
                 &faceDofColIndex,
                 &dFlux_dFaceVars,
                 1, 1 );

  }      
}
  

void
TwoPhaseHybridFVM::ApplyBoundaryConditions( real64 const GEOSX_UNUSED_ARG( time_n ),
                                            real64 const GEOSX_UNUSED_ARG( dt ),
                                            DomainPartition * const GEOSX_UNUSED_ARG( domain ),
                                            DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                                            ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                            ParallelVector & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  // TODO: implement boundary conditions the mimetic way

}

real64 TwoPhaseHybridFVM::CalculateResidualNorm( DomainPartition const * const domain,
                                                 DofManager const & dofManager,
                                                 ParallelVector const & rhs )
{
  MeshLevel const * const mesh          = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = mesh->getFaceManager();

  localIndex constexpr numPhases = NUM_PHASES;
  
  // here we compute the cell-centered residual norm in the derived class
  // to avoid duplicating a synchronization point 
  
  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString );
  
  // local residual
  array1d<real64> localResidualNorm( numPhases+1 ), globalResidualNorm( numPhases+1 );
  localResidualNorm  = 0;
  globalResidualNorm = 0;
  
  // 1. Compute the residual for the mass conservation equations

  // compute the norm of local residual scaled by cell pore volume
  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<globalIndex const> const & elemDofNumber =
      subRegion->getReference< array1d<globalIndex> >( elemDofKey );

    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<real64 const> const & volume         = m_volume[er][esr];
    arrayView1d<real64 const> const & porosityOld    = m_porosityOld[er][esr];    
    arrayView1d<real64 const> const & satOld         = m_wettingPhaseSat[er][esr];
    arrayView2d<real64> const & phaseDensOld         = m_phaseDensOld[er][esr];

    localIndex const subRegionSize = subRegion->size();
    for ( localIndex a = 0; a < subRegionSize; ++a )
    {
      if (elemGhostRank[a] < 0)
      {
        // for the normalization, compute a saturation-weighted total density
        real64 const totalDensOld = satOld[a]     * phaseDensOld[a][m_wettingPh]
                                  + (1-satOld[a]) * phaseDensOld[a][m_nonWettingPh];
        real64 const normalizer   = porosityOld[a] * totalDensOld * volume[a];
        
        for (localIndex ip = 0; ip < numPhases; ++ip)
        {  
          localIndex const lid = rhs.getLocalRowID( elemDofNumber[a] + m_phaseToRow[ip] );
          real64 const val = localResidual[lid] / normalizer;
          localResidualNorm[m_phaseToRow[ip]] += val * val;
        }
      }
    }
  });

  
  // 2. Compute the residual for the face-based constraints

  arrayView1d<integer const> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference< array1d<globalIndex> >( faceDofKey );

  for( localIndex iface = 0 ; iface < faceManager->size(); ++iface )
  {
    // if not ghost face and if adjacent to target region
    if (faceGhostRank[iface] < 0 && faceDofNumber[iface] >= 0)
    {
      localIndex const lid    = rhs.getLocalRowID( faceDofNumber[iface] );
      real64 const normalizer = 1; // TODO: compute the normalizer here
      real64 const val        = localResidual[lid] / normalizer;
      localResidualNorm[2] += val * val;
    }
  }

  
  // 3. Combine the two norms

  MpiWrapper::allReduce(localResidualNorm.data(), globalResidualNorm.data(), numPhases+1, MPI_SUM, MPI_COMM_GEOSX);
 
  real64 maxNorm = 0;
  for (localIndex rowId = 0; rowId < numPhases + 1; ++rowId)
  {
    if (globalResidualNorm[rowId] > maxNorm)
    {
      maxNorm = globalResidualNorm[rowId];
    }
  }
  return sqrt( maxNorm );
}


void TwoPhaseHybridFVM::ApplySystemSolution( DofManager const & dofManager,
                                             ParallelVector const & solution,
                                             real64 const scalingFactor,
                                             DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager * const faceManager = mesh->getFaceManager();
  
  // 1. apply the elem-based update

  applyToSubRegions( mesh, [&] ( localIndex const GEOSX_UNUSED_ARG( er ),
                                 localIndex const GEOSX_UNUSED_ARG( esr ),
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    dofManager.addVectorToField( solution,
                                 viewKeyStruct::elemDofFieldString,
                                 scalingFactor,
                                 subRegion,
                                 viewKeyStruct::deltaPressureString,
                                 0, 1 );

    dofManager.addVectorToField( solution,
                                 viewKeyStruct::elemDofFieldString,
                                 scalingFactor,
                                 subRegion,
                                 viewKeyStruct::deltaWettingPhaseSatString,
                                 1, 2 );
  } );

  // 2. apply the face-based update

  dofManager.addVectorToField( solution,
                               viewKeyStruct::faceDofFieldString,
                               scalingFactor,
                               faceManager,
                               viewKeyStruct::deltaFacePressureString );

  // 3. synchronize
  
  std::map<string, string_array> fieldNames;
  fieldNames["face"].push_back( viewKeyStruct::deltaFacePressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaWettingPhaseSatString );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         mesh,
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    UpdateState( subRegion );
  } );
}


void TwoPhaseHybridFVM::ResetStateToBeginningOfStep( DomainPartition * const GEOSX_UNUSED_ARG( domain ) )
{
  
}


// this function is obviously redundant with computeCellStencil in the TwoPointFluxApproximation class
// this is here for now, but I will have to find a better place for this type of function at some point 
void TwoPhaseHybridFVM::ComputeTransmissibilityMatrix( arrayView1d<R1Tensor const> const & nodePosition, 
                                                       ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                                       arraySlice1d<localIndex const> const elemToFaces,
                                                       R1Tensor const & elemCenter,
                                                       real64   const & GEOSX_UNUSED_ARG( elemVolume ),
                                                       R1Tensor const & elemPerm,
                                                       real64   const & lengthTolerance,
                                                       stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                                           *MAX_NUM_FACES_IN_ELEM> & transMatrix ) const 
{
  ComputeTPFAInnerProduct( nodePosition,
                           faceToNodes,
                           elemToFaces,
                           elemCenter,
                           elemPerm,
                           lengthTolerance,
                           transMatrix );
}

namespace
{

// for now, I just copy-pasted this function from TwoPointFluxApproximation
// this will go away at some point 
void makeFullTensor(R1Tensor const & values, R2SymTensor & result)
{
  result = 0.0;
  R1Tensor axis;
  R2SymTensor temp;

  // assemble full tensor from eigen-decomposition
  for (unsigned icoord = 0; icoord < 3; ++icoord)
  {
    // assume principal axis aligned with global coordinate system
    axis = 0.0;
    axis(icoord) = 1.0;

    // XXX: is there a more elegant way to do this?
    temp.dyadic_aa(axis);
    temp *= values(icoord);
    result += temp;
  }
}

}

// this function is obviously redundant with computeCellStencil in the TwoPointFluxApproximation class
// this is here for now, but I will have to find a better place for this type of function at some point 
void TwoPhaseHybridFVM::ComputeTPFAInnerProduct( arrayView1d<R1Tensor const> const & nodePosition, 
                                                 ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                                 arraySlice1d<localIndex const> const elemToFaces,
                                                 R1Tensor const & elemCenter,
                                                 R1Tensor const & elemPerm,
                                                 real64   const & lengthTolerance,
                                                 stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                                        *MAX_NUM_FACES_IN_ELEM> const & transMatrix ) const
{
  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor permeabilityTensor;

  real64 const areaTolerance   = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  localIndex const numFacesInElem = elemToFaces.size();
  
  // we are ready to compute the transmissibility matrix
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    for (localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc)
    {
      // for now, TPFA trans
      if (ifaceLoc == jfaceLoc)
      {
        
        // 1) compute the face geometry data: center, normal, vector from cell center to face center
        real64 const faceArea =
          computationalGeometry::Centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                     faceToNodes.sizeOfArray(elemToFaces[ifaceLoc]),
                                                     nodePosition,
                                                     faceCenter,
                                                     faceNormal,
                                                     areaTolerance );
        
        cellToFaceVec  = faceCenter;
        cellToFaceVec -= elemCenter;

        if (Dot(cellToFaceVec, faceNormal) < 0.0)
        {
          faceNormal *= -1;
        }

        real64 const c2fDistance = cellToFaceVec.Normalize();

        // 2) assemble full coefficient tensor from principal axis/components
        makeFullTensor(elemPerm, permeabilityTensor);

        faceConormal.AijBj(permeabilityTensor, faceNormal);

        // 3) compute the one-sided face transmissibility
        transMatrix[ifaceLoc][jfaceLoc]  = Dot(cellToFaceVec, faceConormal);
        transMatrix[ifaceLoc][jfaceLoc] *= faceArea / c2fDistance;
        transMatrix[ifaceLoc][jfaceLoc]  = std::max( transMatrix[ifaceLoc][jfaceLoc],
                                                     weightTolerance );
      }
      else
      {
        transMatrix[ifaceLoc][jfaceLoc] = 0;
      }
    }
  }
}
 
REGISTER_CATALOG_ENTRY( SolverBase, TwoPhaseHybridFVM, std::string const &, Group * const )
} /* namespace geosx */
