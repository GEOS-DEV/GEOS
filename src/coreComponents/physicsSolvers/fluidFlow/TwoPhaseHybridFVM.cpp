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
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"

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
  m_elemDofKey(""),  
  m_areaRelTol(1e-8),
  m_minTotalMob(1e-8),
  m_ipType(InnerProductType::QUASI_TPFA),
  m_orthonormalizeWithSVD(false)
{
  // two elem-centered dof per elem
  m_numDofPerCell = 2;
}

  
void TwoPhaseHybridFVM::RegisterDataOnMesh(Group * const MeshBodies)
{
  // 1) Register the elem-centered data
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

    faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::sumTransmissibilityString );    
    faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::weightedGravityCoefString );
    
  }
}

void TwoPhaseHybridFVM::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  TwoPhaseBase::InitializePostInitialConditions_PreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  PrecomputeData( domain );

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

  // setup the elem-centered fields
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

  // increment the elem-centered fields
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

void TwoPhaseHybridFVM::SetupDofs( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                   DofManager & dofManager ) const
{
  localIndex constexpr numDof = NUM_DOF;
  
  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::elemDofFieldString,
                       DofManager::Location::Elem,
                       numDof,
                       m_targetRegions );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString,
                          viewKeyStruct::elemDofFieldString,
                          DofManager::Connectivity::Face );
  
  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::faceDofFieldString,
                       DofManager::Location::Face,
                       m_targetRegions );

  dofManager.addCoupling( viewKeyStruct::faceDofFieldString,
                          viewKeyStruct::faceDofFieldString,
                          DofManager::Connectivity::Elem );
  
  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::faceDofFieldString,
                          viewKeyStruct::elemDofFieldString,
                          DofManager::Connectivity::Elem,
                          true );

 
}

void TwoPhaseHybridFVM::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
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

  string const faceDofKey = dofManager->getKey( viewKeyStruct::faceDofFieldString );  
  string const elemDofKey = dofManager->getKey( viewKeyStruct::elemDofFieldString );
  // save the face Dof key for use in two functions
  // that do not have acces to the coupled solver dofManager
  // namely ResetStateToBeginningOfStep and ImplicitStepComplete
  m_faceDofKey = faceDofKey;
  // save the elem Dof key for use in UpdateUpwindedCoefficients  
  m_elemDofKey = elemDofKey;
  
  // in this function we need to make sure that we act only on the target regions
  // for that, we need the following region filter
  SortedArray<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }
  
  // node data (for transmissibility computation)

  arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition =
    nodeManager->referencePosition();

  
  // face data

  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference< array1d<globalIndex> >( faceDofKey );  

  // get the face-centered pressures
  arrayView1d<real64 const> const & facePres =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::facePressureString );
  arrayView1d<real64 const> const & dFacePres =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::deltaFacePressureString );
  
  // get the face-centered depth 
  arrayView1d<real64> const & faceGravCoef =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::gravityCoefString);
  arrayView1d<real64> const & weightedFaceGravCoef =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::weightedGravityCoefString);
  
  // get the face-to-nodes connectivity for the transmissibility calculation
  ArrayOfArraysView<localIndex const> const & faceToNodes = faceManager->nodeList();

  array2d<localIndex> const & elemRegionList    = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList          = faceManager->elementList(); 

  
  // max number of faces allowed in an element 
  localIndex constexpr maxNumFaces = MAX_NUM_FACES;
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
    
    // get the elem-centered DOF numbers and ghost rank for the assembly
    arrayView1d<globalIndex const> const & elemDofNumber =
      subRegion->template getReference< array1d<globalIndex> >( elemDofKey ); 
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];   

    // get the map from elem to faces
    arrayView2d< localIndex const > const & elemToFaces = subRegion->faceList();
    
    // get the elem-centered pressures
    arrayView1d<real64 const> const & elemPres  = m_pressure[er][esr];
    arrayView1d<real64 const> const & dElemPres = m_deltaPressure[er][esr];

    // get the elem centered densities and mobilities
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

    // get the elem-centered depth 
    arrayView1d<real64 const> const & elemGravCoef =
      subRegion->template getReference<array1d<real64>>(viewKeyStruct::gravityCoefString);

    
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

        stackArray1d<localIndex,  3>             elemIds( 3 );
        stackArray2d<localIndex,  3*maxNumFaces> neighborIds( numFacesInElem, 3 );
        stackArray1d<globalIndex, maxNumFaces>   neighborDofNumbers( numFacesInElem );
        elemIds[0]  = er;
        elemIds[1]  = esr;
        elemIds[2]  = ei;
        neighborIds = -1;
        neighborDofNumbers = -1;
        
        // one sided total flux: \sum_p T \lambda_p ( \nabla p - \rho_p g \nabla d)
        stackArray1d<real64, maxNumFaces> totalVolFlux( numFacesInElem );
        stackArray1d<real64, maxNumFaces> dTotalVolFlux_dp( numFacesInElem );
        stackArray1d<real64, maxNumFaces> dTotalVolFlux_dS( numFacesInElem );
        stackArray2d<real64, maxNumFaces*maxNumFaces> dTotalVolFlux_dfp( numFacesInElem, numFacesInElem );

        // one sided grav term: T (\rho_l - \rho_m) g \nabla d 
        stackArray2d<real64, numPhases*maxNumFaces>   difGravHead( numFacesInElem, numPhases );
        stackArray3d<real64, 2*numPhases*maxNumFaces> dDifGravHead_dp( numFacesInElem, numPhases, 2 );
        
        // upwinded coefficient for the viscous term: \rho_p \lambda_p / \lambda_T    
        stackArray2d<real64, numPhases*maxNumFaces>   viscousCoef( numFacesInElem, numPhases  );
        stackArray3d<real64, 2*numPhases*maxNumFaces> dViscousCoef_dp( numFacesInElem, numPhases, 2 );
        stackArray3d<real64, 2*numPhases*maxNumFaces> dViscousCoef_dS( numFacesInElem, numPhases, 2 );

        // upwinded coefficient for the grav term: \rho_p \lambda_p \lambda_m / \lambda_T  
        stackArray2d<real64, numPhases*maxNumFaces>   gravCoef( numFacesInElem, numPhases );
        stackArray3d<real64, 2*numPhases*maxNumFaces> dGravCoef_dp( numFacesInElem, numPhases, 2 );
        stackArray3d<real64, 2*numPhases*maxNumFaces> dGravCoef_dS( numFacesInElem, numPhases, 2 );

        // On the fly, find the neighbors at all the interfaces of this element
        // We can decide later to precompute and save these indices
        FindAllNeighborsInTarget( mesh,
                                  elemRegionList,
                                  elemSubRegionList,
                                  elemList,
                                  regionFilter,
                                  elemToFaces[ei],
                                  elemIds,
                                  elemDofNumber[ei],
                                  neighborIds,
                                  neighborDofNumbers );

        
        // Recompute the local transmissibility matrix at each iteration
        // We can decide later to precompute transMatrix if needed
        ComputeTransmissibilityMatrix( nodePosition,     
                                       faceToNodes,        
                                       elemToFaces[ei],
                                       elemCenter[ei],
                                       elemVolume[ei],
                                       elemPerm[ei],
                                       lengthTolerance,
                                       transMatrix );   

        
        /*
         * In this function, we want to assemble the flux in a total flux formulation
         * To do that, we first compute at all the faces of this element: 
         * 1) One-sided total volumetric fluxes
         * 2) One-sided gravity term
         * Once we know these quantities, we can upwind the mobility ratios
         * that will multiply the total flux and the gravity term
         */     
 
        // For each one-sided face of the elem,
        // Compute the total volumetric flux using transMatrix
        ComputeTotalVolFluxes( facePres,
                               dFacePres,
                               faceGravCoef,
                               elemToFaces[ei],          
                               elemPres[ei],
                               dElemPres[ei],
                               elemGravCoef[ei],
                               elemPhaseMob[ei],
                               dElemPhaseMob_dp[ei],
                               dElemPhaseMob_dS[ei],
                               elemPhaseDens[ei][0],
                               dElemPhaseDens_dp[ei][0],
                               transMatrix,
                               totalVolFlux,             
                               dTotalVolFlux_dp,
                               dTotalVolFlux_dS,
                               dTotalVolFlux_dfp );

        
        // For each one-sided face of the elem,
        // Compute the grav term T (\rho_p - \rho_m) g \nabla g using transMatrix
        // ** this function needs non-local information because (\rho_p - \rho_m) is averaged across the faces 
        ComputeGravityHead( weightedFaceGravCoef,
                            elemToFaces[ei],          
                            elemGravCoef[ei],
                            m_phaseDens,
                            m_dPhaseDens_dPres,
                            elemIds,
                            neighborIds,
                            transMatrix,
                            difGravHead,
                            dDifGravHead_dp );

        
        // At this point, we know the local flow direction (for the viscous and gravity contributtions) in the element
        // So we can upwind the transport coefficients (mobility ratios) at the one sided faces 
        // ** this function needs non-local information because of the upwinding **
        UpdateUpwindedCoefficients( elemToFaces[ei],
                                    m_phaseMob, 
                                    m_dPhaseMob_dPres,
                                    m_dPhaseMob_dSat,
                                    m_phaseDens, 
                                    m_dPhaseDens_dPres,
                                    elemIds,
                                    neighborIds, 
                                    totalVolFlux,
                                    difGravHead,
                                    viscousCoef,      
                                    dViscousCoef_dp,
                                    dViscousCoef_dS,
                                    gravCoef,      
                                    dGravCoef_dp,
                                    dGravCoef_dS );
        
        /*
         * At this point, we have computed all the fluxes (viscous and gravity contributions) in the element
         * We perform the assembly in this element in two steps:
         * 1) mass conservation equations
         * 2) face constraints
         */
        
        // Use the computed one sided total vol fluxes (+ grav terms) and the upwinded mobilities
        // to assemble the upwinded mass fluxes in the mass conservation eqns of the elem
        AssembleOneSidedMassFluxes( dt,
                                    faceDofNumber,
                                    elemToFaces[ei],
                                    elemDofNumber[ei],
                                    neighborDofNumbers,
                                    totalVolFlux,
                                    dTotalVolFlux_dp,
                                    dTotalVolFlux_dS,
                                    dTotalVolFlux_dfp,
                                    difGravHead,
                                    dDifGravHead_dp,
                                    viscousCoef,      
                                    dViscousCoef_dp,
                                    dViscousCoef_dS,
                                    gravCoef,
                                    dGravCoef_dp,
                                    dGravCoef_dS,
                                    matrix,
                                    rhs );

        // Use the computed one sided vol fluxes to assemble the constraints
        // enforcing total vol flux continuity at this element's faces
        AssembleConstraints( faceDofNumber,
                             elemToFaces[ei],
                             elemDofNumber[ei],
                             totalVolFlux,
                             dTotalVolFlux_dp,
                             dTotalVolFlux_dS,                       
                             dTotalVolFlux_dfp,
                             matrix,
                             rhs );

      }
    });
  });
}


void TwoPhaseHybridFVM::ComputeTotalVolFluxes( arrayView1d<real64 const> const & facePres,
                                               arrayView1d<real64 const> const & dFacePres,
                                               arrayView1d<real64 const> const & faceGravCoef,
                                               arraySlice1d<localIndex const> const elemToFaces,
                                               real64 const & elemPres,
                                               real64 const & dElemPres,
                                               real64 const & elemGravCoef,
                                               arraySlice1d<real64 const> const elemMob,
                                               arraySlice1d<real64 const> const dElemMob_dp,
                                               arraySlice1d<real64 const> const dElemMob_dS,
                                               arraySlice1d<real64 const> const elemDens,
                                               arraySlice1d<real64 const> const dElemDens_dp,
                                               stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> const & transMatrix,
                                               stackArray1d<real64, MAX_NUM_FACES> & totalVolFlux,
                                               stackArray1d<real64, MAX_NUM_FACES> & dTotalVolFlux_dp,
                                               stackArray1d<real64, MAX_NUM_FACES> & dTotalVolFlux_dS,
                                               stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> & dTotalVolFlux_dfp ) const
{
  localIndex constexpr numPhases   = NUM_PHASES;
  localIndex constexpr maxNumFaces = MAX_NUM_FACES;
  
  localIndex const numFacesInElem = elemToFaces.size();

  totalVolFlux      = 0;
  dTotalVolFlux_dp  = 0;
  dTotalVolFlux_dS  = 0;    
  dTotalVolFlux_dfp = 0;
  
  // this is going to store \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
  stackArray1d<real64, maxNumFaces> sumWeightedPotDif( numFacesInElem );
  stackArray1d<real64, maxNumFaces> dSumWeightedPotDif_dp( numFacesInElem );
  stackArray1d<real64, maxNumFaces> dSumWeightedPotDif_dS( numFacesInElem );
  stackArray1d<real64, maxNumFaces> dSumWeightedPotDif_dfp( numFacesInElem );
  sumWeightedPotDif      = 0;
  dSumWeightedPotDif_dp  = 0;
  dSumWeightedPotDif_dS  = 0;
  dSumWeightedPotDif_dfp = 0;

  // 1) precompute the potential difference at each one-sided face
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    // pressure and gravity terms
    // TODO: account for capillary pressure here
    real64 const ccPres = elemPres + dElemPres;
    real64 const fPres  = facePres[elemToFaces[ifaceLoc]] + dFacePres[elemToFaces[ifaceLoc]];
      
    real64 const ccGravCoef = elemGravCoef;
    real64 const fGravCoef  = faceGravCoef[elemToFaces[ifaceLoc]];

    // pressure difference
    real64 const presDif      = ccPres - fPres;
    real64 const dPresDif_dp  =  1;
    real64 const dPresDif_dfp = -1;

    // depth difference between element center and face center
    real64 const depthDif = ccGravCoef - fGravCoef;

    for (localIndex ip = 0; ip < numPhases; ++ip)
    {
      // TODO: check if we can use the faceDensity here
      real64 const phaseGravTerm     = elemDens[ip]     * depthDif;
      real64 const dPhaseGravTerm_dp = dElemDens_dp[ip] * depthDif;

      // obtain the phase potential difference
      real64 const phasePotDif       = presDif     - phaseGravTerm;
      real64 const dPhasePotDif_dp   = dPresDif_dp - dPhaseGravTerm_dp;
      real64 const dPhasePotDif_dfp  = dPresDif_dfp;

      // sum of potential differences weighted by phase mobility
      sumWeightedPotDif[ifaceLoc]      += elemMob[ip]     * phasePotDif;
      dSumWeightedPotDif_dp[ifaceLoc]  += dElemMob_dp[ip] * phasePotDif
                                        + elemMob[ip]     * dPhasePotDif_dp;
      dSumWeightedPotDif_dS[ifaceLoc]  += dElemMob_dS[ip] * phasePotDif;
      dSumWeightedPotDif_dfp[ifaceLoc] += elemMob[ip]     * dPhasePotDif_dfp;
    }
  }
  
  // 2) multiply the potential difference by the transmissibility
  //    to obtain the total volumetrix flux
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
       
    // now in the following nested loop, 
    // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
    for (localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc)
    {
      // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
      totalVolFlux[ifaceLoc]                += transMatrix[ifaceLoc][jfaceLoc] * sumWeightedPotDif[jfaceLoc];
      dTotalVolFlux_dp[ifaceLoc]            += transMatrix[ifaceLoc][jfaceLoc] * dSumWeightedPotDif_dp[jfaceLoc];
      dTotalVolFlux_dS[ifaceLoc]            += transMatrix[ifaceLoc][jfaceLoc] * dSumWeightedPotDif_dS[jfaceLoc];      
      dTotalVolFlux_dfp[ifaceLoc][jfaceLoc] += transMatrix[ifaceLoc][jfaceLoc] * dSumWeightedPotDif_dfp[jfaceLoc];
    }
  }
}


void TwoPhaseHybridFVM::ComputeGravityHead( arrayView1d<real64 const> const & weightedFaceGravCoef,
                                            arraySlice1d<localIndex const> const elemToFaces,
                                            real64 const & elemGravCoef,
                                            ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dens,
                                            ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dDens_dp,
                                            stackArray1d<localIndex, 3>                       const & elemIds,
                                            stackArray2d<localIndex, 3*MAX_NUM_FACES>         const & neighborIds,
                                            stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> const & transMatrix,
                                            stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>          & difGravHead,
                                            stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>        & dDifGravHead_dp ) const
{
  localIndex const fi = m_fluidIndex;  
  localIndex const numFacesInElem = elemToFaces.size();

  localIndex const er  = elemIds[0];
  localIndex const esr = elemIds[1];
  localIndex const ei  = elemIds[2];
  
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    difGravHead[ifaceLoc][m_ipw] = 0;
    localIndex const ern  = neighborIds[ifaceLoc][0];
    localIndex const esrn = neighborIds[ifaceLoc][1];
    localIndex const ein  = neighborIds[ifaceLoc][2];
    
    for (localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc)
    {  
      localIndex const faceId = elemToFaces[jfaceLoc];
    
      // 1) multiply the depth by the transmissibility
      // note: flux continuity at the face is guaranteed
      // because we use this weightedFaceGravCoef (a static Lagrange multiplier)
      // precomputed such that the two gravity fluxes match.
      // TODO: investigate the impact of this (so far not tested with a consistent inner product)
      real64 depthDif = elemGravCoef - weightedFaceGravCoef[faceId];
      difGravHead[ifaceLoc][m_ipw] += transMatrix[ifaceLoc][jfaceLoc] * depthDif;
    }

    dDifGravHead_dp[ifaceLoc][m_ipw][0] = difGravHead[ifaceLoc][m_ipw];
    dDifGravHead_dp[ifaceLoc][m_ipw][1] = difGravHead[ifaceLoc][m_ipw];    

    // 2) multiply the gravity term by the density difference evaluated at the center of the elem
    //    we want to obtain: T (\rho_p - \rho_m) g \nabla d

    // TODO: check if we can use a face density here
    
    // first we have to compute the average density at the interface
    real64 difAvgDens             = dens[er][esr][fi][ei][0][m_ipw]     - dens[er][esr][fi][ei][0][m_ipnw]; 
    real64 dDifAvgDens_dpLoc      = dDens_dp[er][esr][fi][ei][0][m_ipw] - dDens_dp[er][esr][fi][ei][0][m_ipnw] ;
    real64 dDifAvgDens_dpNeighbor = 0;
    if (ern != -1 && esrn != -1 && ein != -1)
    {
      difAvgDens            *= 0.5;
      difAvgDens            += 0.5 * ( dens[ern][esrn][fi][ein][0][m_ipw] - dens[ern][esrn][fi][ein][0][m_ipnw] );
      dDifAvgDens_dpLoc     *= 0.5;  
      dDifAvgDens_dpNeighbor = 0.5 * ( dDens_dp[ern][esrn][fi][ein][0][m_ipw] - dDens_dp[ern][esrn][fi][ein][0][m_ipnw] );
    }

    // then we multiply the depth term by the difference between the phase avg densities
    // wetting phase 
    difGravHead[ifaceLoc][m_ipw]        *= difAvgDens;
    dDifGravHead_dp[ifaceLoc][m_ipw][0] *= dDifAvgDens_dpLoc;
    dDifGravHead_dp[ifaceLoc][m_ipw][1] *= dDifAvgDens_dpNeighbor;
    // non wetting phase (stored for convenience)
    difGravHead[ifaceLoc][m_ipnw]        = - difGravHead[ifaceLoc][m_ipw]; 
    dDifGravHead_dp[ifaceLoc][m_ipnw][0] = - dDifGravHead_dp[ifaceLoc][m_ipw][0];
    dDifGravHead_dp[ifaceLoc][m_ipnw][1] = - dDifGravHead_dp[ifaceLoc][m_ipw][1];    
  }
}

  
void TwoPhaseHybridFVM::UpdateUpwindedCoefficients( arraySlice1d<localIndex const> const elemToFaces,
                                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & mob,
                                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dp,
                                                    ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dS,
                                                    ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dens,
                                                    ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dDens_dp,
                                                    stackArray1d<localIndex, 3>                    const & elemIds,
                                                    stackArray2d<localIndex, 3*MAX_NUM_FACES>      const & neighborIds,
                                                    stackArray1d<real64, MAX_NUM_FACES>            const & totalVolFlux,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES> const & difGravHead,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>       & viscousCoef,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>     & dViscousCoef_dp,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>     & dViscousCoef_dS,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>       & gravCoef,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>     & dGravCoef_dp,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>     & dGravCoef_dS ) const
{
  localIndex const numFacesInElem = elemToFaces.size();

  viscousCoef     = 0;
  dViscousCoef_dp = 0;
  dViscousCoef_dS = 0;
  
  gravCoef     = 0;
  dGravCoef_dp = 0;
  dGravCoef_dS = 0;
  
  // for this element, loop over the local (one-sided) faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    // TODO: upwind the density separately
    // TODO: implement PPU as an option here
     
    // upwinded coefficient for the viscous terms
    // compute \rho_p \lambda_p / \lambda_T
    UpdateLocalViscousCoefficients( mob, 
                                    dMob_dp,
                                    dMob_dS,
                                    dens,
                                    dDens_dp,
                                    elemIds,
                                    neighborIds[ifaceLoc],
                                    totalVolFlux[ifaceLoc],
                                    viscousCoef[ifaceLoc],
                                    dViscousCoef_dp[ifaceLoc],
                                    dViscousCoef_dS[ifaceLoc] );

    // upwinded coefficient for the grav terms
    // compute \rho_p \lambda_p \lambda_m / \lambda_T
    if (neighborIds[ifaceLoc][0] != -1 &&
        neighborIds[ifaceLoc][1] != -1 &&
        neighborIds[ifaceLoc][2] != -1)
    {
      UpdateLocalGravCoefficients( mob,
                                   dMob_dp,
                                   dMob_dS,
                                   dens,
                                   dDens_dp,
                                   elemIds,
                                   neighborIds[ifaceLoc],
                                   difGravHead[ifaceLoc],                                  
                                   gravCoef[ifaceLoc],
                                   dGravCoef_dp[ifaceLoc],
                                   dGravCoef_dS[ifaceLoc] );
    }
  }                             
}

void TwoPhaseHybridFVM::UpdateLocalViscousCoefficients( ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const  & mob,
                                                        ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const  & dMob_dp,
                                                        ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const  & dMob_dS,
                                                        ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dens,
                                                        ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dDens_dp,
                                                        stackArray1d<localIndex, 3>    const & elemIds,
                                                        arraySlice1d<localIndex const> const neighborIds,
                                                        real64 const & totalVolFlux,
                                                        arraySlice1d<real64> const viscousCoef,
                                                        arraySlice2d<real64> const dViscousCoef_dp,
                                                        arraySlice2d<real64> const dViscousCoef_dS ) const
{
  localIndex constexpr numPhases = NUM_PHASES;

  // we initialize these upw quantities with the values of the local elem
  // they will be updated below if there is a neighbor
  bool const foundNeighborInTarget = (neighborIds[0] != -1 &&
                                      neighborIds[1] != -1 &&
                                      neighborIds[2] != -1);
  localIndex const neighborIsUpw = foundNeighborInTarget && totalVolFlux < 0;
  
  // we evaluate the mobilities with the upwind element wrt the total volumetrix flux
  // it is my understanding that I shoud avoid if statements so I am doing that  
  localIndex const erUpw  = neighborIsUpw * neighborIds[0] + (1-neighborIsUpw) * elemIds[0];
  localIndex const esrUpw = neighborIsUpw * neighborIds[1] + (1-neighborIsUpw) * elemIds[1];
  localIndex const eiUpw  = neighborIsUpw * neighborIds[2] + (1-neighborIsUpw) * elemIds[2];

  // 1) compute the total mobility
  real64 totalMob     = 0;
  real64 dTotalMob_dp = 0;
  real64 dTotalMob_dS = 0;

  for (localIndex ip = 0; ip < numPhases; ++ip)
  {
    totalMob     += mob[erUpw][esrUpw][eiUpw][ip];
    dTotalMob_dp += dMob_dp[erUpw][esrUpw][eiUpw][ip];
    dTotalMob_dS += dMob_dS[erUpw][esrUpw][eiUpw][ip];
  }
  // physically the total mobility cannot be zero so I don't check its value here
   
  // 2) for each phase, compute the upwinded viscous coefficients (+derivatives) 
  //    using the upwinded phase mobility ratioss
  for (localIndex ip = 0; ip < numPhases; ++ip)
  {
    // upwinded mobility ratios
    real64 const mobRatio     = mob[erUpw][esrUpw][eiUpw][ip] / totalMob;
    real64 const dMobRatio_dp = ( dMob_dp[erUpw][esrUpw][eiUpw][ip] * totalMob
                                - mob[erUpw][esrUpw][eiUpw][ip]     * dTotalMob_dp)
                              / (totalMob*totalMob);
    real64 const dMobRatio_dS = ( dMob_dS[erUpw][esrUpw][eiUpw][ip] * totalMob
                                - mob[erUpw][esrUpw][eiUpw][ip]     * dTotalMob_dS)
                              / (totalMob*totalMob);

    localIndex const dupw  = (erUpw  == neighborIds[0] &&
                              esrUpw == neighborIds[1] &&
                              eiUpw  == neighborIds[2]);
    localIndex const fi = m_fluidIndex;
    
    // upwinded viscous coef
    // we want to obtain \rho_p \lambda_p / \lambda_T
    viscousCoef[ip]           = dens[erUpw][esrUpw][fi][eiUpw][0][ip]     * mobRatio;
    dViscousCoef_dp[ip][dupw] = dDens_dp[erUpw][esrUpw][fi][eiUpw][0][ip] * mobRatio
                              + dens[erUpw][esrUpw][fi][eiUpw][0][ip]     * dMobRatio_dp;
    dViscousCoef_dS[ip][dupw] = dens[erUpw][esrUpw][fi][eiUpw][0][ip]     * dMobRatio_dS;
  }
}


void TwoPhaseHybridFVM::UpdateLocalGravCoefficients( ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & mob,
                                                     ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dp,
                                                     ElementRegionManager::ElementViewAccessor<arrayView2d<real64>>  const & dMob_dS,
                                                     ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dens,
                                                     ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dDens_dp,
                                                     stackArray1d<localIndex, 3>    const & elemIds,
                                                     arraySlice1d<localIndex const> const neighborIds,
                                                     arraySlice1d<real64 const> const difGravCoef,
                                                     arraySlice1d<real64> const gravCoef,
                                                     arraySlice2d<real64> const dGravCoef_dp,
                                                     arraySlice2d<real64> const dGravCoef_dS ) const
{
  // here, we evaluate:
  // - the mobility of the heavy phase with the props of the bottom elem
  // - the mobility of the light phase with the props of the top elem
  localIndex const neighborIsUpw = difGravCoef[m_ipw] < 0;
  // it is my understanding that I shoud avoid if statements so I am doing that
  localIndex const erw   = neighborIsUpw * neighborIds[0] + (1-neighborIsUpw) * elemIds[0];
  localIndex const esrw  = neighborIsUpw * neighborIds[1] + (1-neighborIsUpw) * elemIds[1];
  localIndex const eiw   = neighborIsUpw * neighborIds[2] + (1-neighborIsUpw) * elemIds[2];
  localIndex const ernw  = neighborIsUpw * elemIds[0]     + (1-neighborIsUpw) * neighborIds[0];
  localIndex const esrnw = neighborIsUpw * elemIds[1]     + (1-neighborIsUpw) * neighborIds[1];
  localIndex const einw  = neighborIsUpw * elemIds[2]     + (1-neighborIsUpw) * neighborIds[2];
  
  // upwinded total mobility
  real64 totalMob          = mob[erw][esrw][eiw][m_ipw] + mob[ernw][esrnw][einw][m_ipnw];
  real64 dTotalMob_dp_eiw  = dMob_dp[erw][esrw][eiw][m_ipw];
  real64 dTotalMob_dp_einw = dMob_dp[ernw][esrnw][einw][m_ipnw];
  real64 dTotalMob_dS_eiw  = dMob_dS[erw][esrw][eiw][m_ipw];
  real64 dTotalMob_dS_einw = dMob_dS[ernw][esrnw][einw][m_ipnw];
  if (totalMob < m_minTotalMob)
  {
    totalMob = m_minTotalMob;
    dTotalMob_dp_eiw  = 0;
    dTotalMob_dp_einw = 0;
    dTotalMob_dS_eiw  = 0;
    dTotalMob_dS_einw = 0;    
  }

  // upwinded phase mobility ratio
  real64 const mobProd            = mob[erw][esrw][eiw][m_ipw]     * mob[ernw][esrnw][einw][m_ipnw];
  real64 const dMobProd_dp_eiw    = dMob_dp[erw][esrw][eiw][m_ipw] * mob[ernw][esrnw][einw][m_ipnw];
  real64 const dMobProd_dp_einw   = mob[erw][esrw][eiw][m_ipw]     * dMob_dp[ernw][esrnw][einw][m_ipnw];
  real64 const dMobProd_dS_eiw    = dMob_dS[erw][esrw][eiw][m_ipw] * mob[ernw][esrnw][einw][m_ipnw];
  real64 const dMobProd_dS_einw   = mob[erw][esrw][eiw][m_ipw]     * dMob_dS[ernw][esrnw][einw][m_ipnw];
  
  real64 const mobRatio           = mobProd / totalMob;  
  real64 const dMobRatio_dp_eiw   = (dMobProd_dp_eiw * totalMob  - mobProd * dTotalMob_dp_eiw)
                                  / (totalMob * totalMob);
  real64 const dMobRatio_dp_einw  = (dMobProd_dp_einw * totalMob - mobProd * dTotalMob_dp_einw)
                                  / (totalMob * totalMob);
  real64 const dMobRatio_dS_eiw   = (dMobProd_dS_eiw * totalMob  - mobProd * dTotalMob_dS_eiw)
                                  / (totalMob * totalMob);  
  real64 const dMobRatio_dS_einw  = (dMobProd_dS_einw * totalMob - mobProd * dTotalMob_dS_einw)
                                  / (totalMob * totalMob);  

  localIndex const deiw  = (erw  == neighborIds[0] && esrw  == neighborIds[1] && eiw  == neighborIds[2]);
  localIndex const deinw = (ernw == neighborIds[0] && esrnw == neighborIds[1] && einw == neighborIds[2]);
  localIndex const fi = m_fluidIndex;
    
  // upwinded grav coef
  // we want to obtain \rho_p \lambda_p \lambda_m / \lambda_T  
  gravCoef[m_ipw]  = dens[erw][esrw][fi][eiw][0][m_ipw]     * mobRatio;
  gravCoef[m_ipnw] = dens[ernw][esrnw][fi][einw][0][m_ipnw] * mobRatio;

  dGravCoef_dp[m_ipw][deiw]   = dDens_dp[erw][esrw][fi][eiw][0][m_ipw]     * mobRatio
                              + dens[erw][esrw][fi][eiw][0][m_ipw]         * dMobRatio_dp_eiw;
  dGravCoef_dp[m_ipw][deinw]  = dens[erw][esrw][fi][eiw][0][m_ipw]         * dMobRatio_dp_einw;
  dGravCoef_dS[m_ipw][deiw]   = dens[erw][esrw][fi][eiw][0][m_ipw]         * dMobRatio_dS_eiw;
  dGravCoef_dS[m_ipw][deinw]  = dens[erw][esrw][fi][eiw][0][m_ipw]         * dMobRatio_dS_einw;

  dGravCoef_dp[m_ipnw][deiw]  = dens[ernw][esrnw][fi][einw][0][m_ipnw]     * dMobRatio_dp_eiw;
  dGravCoef_dp[m_ipnw][deinw] = dDens_dp[ernw][esrnw][fi][einw][0][m_ipnw] * mobRatio
                              + dens[ernw][esrnw][fi][einw][0][m_ipnw]     * dMobRatio_dp_einw;
  dGravCoef_dS[m_ipnw][deiw]  = dens[ernw][esrnw][fi][einw][0][m_ipnw]     * dMobRatio_dS_eiw;
  dGravCoef_dS[m_ipnw][deinw] = dens[ernw][esrnw][fi][einw][0][m_ipnw]     * dMobRatio_dS_einw;
}

void TwoPhaseHybridFVM::AssembleOneSidedMassFluxes( real64 const & dt,
                                                    arrayView1d<globalIndex const> const & faceDofNumber,
                                                    arraySlice1d<localIndex const> const elemToFaces,
                                                    globalIndex const elemDofNumber,
                                                    stackArray1d<globalIndex, MAX_NUM_FACES>          const & neighborDofNumbers,
                                                    stackArray1d<real64, MAX_NUM_FACES>               const & totalVolFlux,
                                                    stackArray1d<real64, MAX_NUM_FACES>               const & dTotalVolFlux_dp,
                                                    stackArray1d<real64, MAX_NUM_FACES>               const & dTotalVolFlux_dS,
                                                    stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> const & dTotalVolFlux_dfp,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>    const & difGravHead,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>  const & dDifGravHead_dp,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>    const & viscousCoef,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>  const & dViscousCoef_dp,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>  const & dViscousCoef_dS,
                                                    stackArray2d<real64, NUM_PHASES*MAX_NUM_FACES>    const & gravCoef,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>  const & dGravCoef_dp,
                                                    stackArray3d<real64, 2*NUM_PHASES*MAX_NUM_FACES>  const & dGravCoef_dS,
                                                    ParallelMatrix * const matrix,
                                                    ParallelVector * const rhs ) const
{
  localIndex constexpr maxNumFaces = MAX_NUM_FACES;
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
  stackArray1d<real64, numPhases>                        sumLocalMassFluxes( numPhases );
  stackArray2d<real64, numPhases*numDof*(1+maxNumFaces)> dSumLocalMassFluxes_dElemVars( numPhases, numDof * (1+numFacesInElem) );
  stackArray2d<real64, numPhases*maxNumFaces>            dSumLocalMassFluxes_dFaceVars( numPhases, numFacesInElem );
  sumLocalMassFluxes = 0;
  dSumLocalMassFluxes_dElemVars = 0;
  dSumLocalMassFluxes_dFaceVars = 0;
  
  // for each element, loop over the one-sided faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    localIndex const fOffset = numDof*(ifaceLoc+1);
    
    for (localIndex ip = 0; ip < numPhases; ++ip)
    {
      localIndex const rowId = m_phaseToRow[ip];
      
      // compute the mass flux at the one-sided face plus its derivatives
      // add the newly computed flux to the sum
      real64 const localViscCoef             = dt * viscousCoef[ifaceLoc][ip];
      real64 const dLocalViscCoef_dpLoc      = dt * dViscousCoef_dp[ifaceLoc][ip][0];
      real64 const dLocalViscCoef_dpNeighbor = dt * dViscousCoef_dp[ifaceLoc][ip][1];      
      real64 const dLocalViscCoef_dSLoc      = dt * dViscousCoef_dS[ifaceLoc][ip][0];
      real64 const dLocalViscCoef_dSNeighbor = dt * dViscousCoef_dS[ifaceLoc][ip][1];      
      
      real64 const localGravCoef             = dt * gravCoef[ifaceLoc][ip];
      real64 const dLocalGravCoef_dpLoc      = dt * dGravCoef_dp[ifaceLoc][ip][0];
      real64 const dLocalGravCoef_dpNeighbor = dt * dGravCoef_dp[ifaceLoc][ip][1];
      real64 const dLocalGravCoef_dSLoc      = dt * dGravCoef_dS[ifaceLoc][ip][0];
      real64 const dLocalGravCoef_dSNeighbor = dt * dGravCoef_dS[ifaceLoc][ip][1];

      // 1) residual
      // the flux of phase p in total velocity formulation is:
      //    T ( \rho_p \lambda_p / \lambda_T q_T
      //      + \rho_p \lambda_p * \lambda_m / \lambda_T (\rho_p - \rho_m) g \nabla z )
      sumLocalMassFluxes[rowId] += localViscCoef * totalVolFlux[ifaceLoc]
                                 + localGravCoef * difGravHead[ifaceLoc][ip];

      // note: dDifGravHead_dS = 0
      
      // 2) derivatives wrt the elem centered vars of the local elem
      dSumLocalMassFluxes_dElemVars[rowId][dp] += localViscCoef        * dTotalVolFlux_dp[ifaceLoc]
                                                + dLocalViscCoef_dpLoc * totalVolFlux[ifaceLoc]
                                                + localGravCoef        * dDifGravHead_dp[ifaceLoc][ip][0]
                                                + dLocalGravCoef_dpLoc * difGravHead[ifaceLoc][ip];
      dSumLocalMassFluxes_dElemVars[rowId][dS] += localViscCoef        * dTotalVolFlux_dS[ifaceLoc]
                                                + dLocalViscCoef_dSLoc * totalVolFlux[ifaceLoc]
                                                + dLocalGravCoef_dSLoc * difGravHead[ifaceLoc][ip];

      // 3) derivatices wrt the elem centered vars of the neighbor
      dSumLocalMassFluxes_dElemVars[rowId][fOffset+dp] += dLocalViscCoef_dpNeighbor * totalVolFlux[ifaceLoc]
                                                        + localGravCoef             * dDifGravHead_dp[ifaceLoc][ip][1]
                                                        + dLocalGravCoef_dpNeighbor * difGravHead[ifaceLoc][ip];
      dSumLocalMassFluxes_dElemVars[rowId][fOffset+dS] += dLocalViscCoef_dSNeighbor * totalVolFlux[ifaceLoc]
                                                        + dLocalGravCoef_dSNeighbor * difGravHead[ifaceLoc][ip];

      // 4) derivatives wrt the face centered var
      for (localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc)
      { 	
        dSumLocalMassFluxes_dFaceVars[rowId][jfaceLoc] += localViscCoef * dTotalVolFlux_dfp[ifaceLoc][jfaceLoc];
      }
    }

    // collect the relevant dof numbers
    elemDofColIndices[fOffset+dp] = neighborDofNumbers[ifaceLoc] + dp;
    elemDofColIndices[fOffset+dS] = neighborDofNumbers[ifaceLoc] + dS;
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
                                             stackArray1d<real64, MAX_NUM_FACES>               const & totalVolFlux,
                                             stackArray1d<real64, MAX_NUM_FACES>               const & dTotalVolFlux_dp,
                                             stackArray1d<real64, MAX_NUM_FACES>               const & dTotalVolFlux_dS,
                                             stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> const & dTotalVolFlux_dfp,
                                             ParallelMatrix * const matrix,
                                             ParallelVector * const rhs ) const 
{
  localIndex constexpr numDof     = NUM_DOF;
  localIndex constexpr maxNumFaces = MAX_NUM_FACES; 
  
  localIndex const numFacesInElem = elemToFaces.size();

  localIndex const dp = TwoPhaseBase::ColOffset::DPRES;
  localIndex const dS = TwoPhaseBase::ColOffset::DSAT;

  // dof numbers
  stackArray1d<globalIndex, numDof>      elemDofColIndices( numDof );
  stackArray1d<globalIndex, maxNumFaces> faceDofColIndices( numFacesInElem );
  
  elemDofColIndices[dp] = elemDofNumber + dp;
  elemDofColIndices[dS] = elemDofNumber + dS;

  // fluxes
  stackArray1d<real64, numDof>      dFlux_dElemVars( numDof );
  stackArray1d<real64, maxNumFaces> dFlux_dFaceVars( numFacesInElem );  
  
  // for each element, loop over the local (one-sided) faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    // dof numbers 
    globalIndex const eqnRowIndex = faceDofNumber[elemToFaces[ifaceLoc]];

    // fluxes
    real64 const flux   = totalVolFlux[ifaceLoc];
    dFlux_dElemVars[dp] = dTotalVolFlux_dp[ifaceLoc];
    dFlux_dElemVars[dS] = dTotalVolFlux_dS[ifaceLoc];

    dFlux_dFaceVars = 0.0;
    for (localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc)
    {
      faceDofColIndices[jfaceLoc] = faceDofNumber[elemToFaces[jfaceLoc]];      
      dFlux_dFaceVars[jfaceLoc]   = dTotalVolFlux_dfp[ifaceLoc][jfaceLoc]; 
    }
    
    // residual
    rhs->add( &eqnRowIndex,
              &flux,
              1 );

    // jacobian -- derivative wrt local elem centered vars
    matrix->add( &eqnRowIndex,
                 elemDofColIndices.data(),
                 dFlux_dElemVars.data(),
                 1, numDof );
   
    // jacobian -- derivatives wrt face pressure terms
    matrix->add( &eqnRowIndex,
                 faceDofColIndices.data(),
                 dFlux_dFaceVars.data(),
                 1, numFacesInElem );

  }      
}

void
TwoPhaseHybridFVM::ApplyBoundaryConditions( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                            real64 const GEOSX_UNUSED_PARAM( dt ),
                                            DomainPartition * const GEOSX_UNUSED_PARAM( domain ),
                                            DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                            ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                            ParallelVector & GEOSX_UNUSED_PARAM( rhs ) )
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
  
  // here we compute the elem-centered residual norm in the derived class
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

  // compute the norm of local residual scaled by elem pore volume
  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_PARAM( region ),
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
        real64 const totalDensOld = satOld[a]     * phaseDensOld[a][m_ipw]
                                  + (1-satOld[a]) * phaseDensOld[a][m_ipnw];
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
  
  // 1. apply the elem-based update

  dofManager.addVectorToField( solution,
                               viewKeyStruct::elemDofFieldString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( solution,
                               viewKeyStruct::elemDofFieldString,
                               viewKeyStruct::deltaWettingPhaseSatString,				 
                               scalingFactor,
                               1, 2 );

  // 2. apply the face-based update

  dofManager.addVectorToField( solution,
                               viewKeyStruct::faceDofFieldString,
                               viewKeyStruct::deltaFacePressureString,
                               scalingFactor );

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

void TwoPhaseHybridFVM::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // 1. Reset the cell-centered fields
  TwoPhaseBase::ResetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  MeshLevel * const mesh          = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();

  // get the accumulated face pressure updates
  arrayView1d<real64> & dFacePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

  forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  {
    dFacePres[iface] = 0;
  });  
}

void TwoPhaseHybridFVM::FindAllNeighborsInTarget( MeshLevel const * const mesh,
                                                  array2d<localIndex> const & elemRegionList,
                                                  array2d<localIndex> const & elemSubRegionList,
                                                  array2d<localIndex> const & elemList,
                                                  SortedArray<localIndex> const & regionFilter,
                                                  arraySlice1d<localIndex const> const elemToFaces,
                                                  stackArray1d<localIndex, 3> const & elemIds,
                                                  globalIndex const elemDofNumber,   
                                                  stackArray2d<localIndex, 3*MAX_NUM_FACES> & neighborIds,
                                                  stackArray1d<globalIndex, MAX_NUM_FACES>  & neighborDofNumber ) const
{
  localIndex const numFacesInElem = elemToFaces.size(1);
  
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    neighborDofNumber[ifaceLoc] = elemDofNumber;
    bool foundNeighborInTarget = false;
    
    // the face has at most two adjacent elements
    // one of these two elements is the current element indexed by er, esr, ei
    // but here we are interested in the indices of the other element
    // this other element is "the neighbor" for this one-sided face 
    for (localIndex k=0; k < elemRegionList.size(1) && !foundNeighborInTarget; ++k)
    {
      localIndex const faceId = elemToFaces[ifaceLoc];
      
      // this element is not the current element
      // we have found the neighbor or we are at the boundary 
      if ( elemRegionList[faceId][k]    != elemIds[0] ||
           elemSubRegionList[faceId][k] != elemIds[1] ||
           elemList[faceId][k]          != elemIds[2] ) 
      {

        bool const onBoundary = (elemRegionList[faceId][k]    == -1 ||
                                 elemSubRegionList[faceId][k] == -1 ||
                                 elemList[faceId][k]          == -1);
        bool const neighborInTarget = regionFilter.contains(elemRegionList[faceId][k]); 

        // if not on boundary, save the mobility and the upwDofNumber  
        if ( !onBoundary && neighborInTarget )
        {
          foundNeighborInTarget = true;
          neighborIds[ifaceLoc][0] = elemRegionList[faceId][k];
          neighborIds[ifaceLoc][1] = elemSubRegionList[faceId][k];
          neighborIds[ifaceLoc][2] = elemList[faceId][k];

          ElementRegionBase const * const neighborRegion =
            Group::group_cast<ElementRegionBase const *>(mesh->getElemManager()->GetRegion( neighborIds[ifaceLoc][0] ));
          ElementSubRegionBase const * const neighborSubRegion =
            Group::group_cast<ElementSubRegionBase const *>(neighborRegion->GetSubRegion( neighborIds[ifaceLoc][1] ));

          arrayView1d<globalIndex const> const & dofNumber =
            neighborSubRegion->getReference<array1d<globalIndex>>(m_elemDofKey);  

          // always save the indices of the neighbor for the grav flux
          neighborDofNumber[ifaceLoc] = dofNumber[neighborIds[ifaceLoc][2]];
        }
      }
    }
  }
}


// for the gravity flux to be mass conservative at a mesh face, we
// need to make sure that the (one-sided) gravity flux coming from
// an element is the same as the (one-sided) gravity flux coming
// from the neighbor. Given an inner product, this can be done by
// solving a linear system once at the beginning of the simulation
// to obtain a transmissibility-weighted depth at each mesh face
// (then used in ComputeGravityHead). For the TPFA, these
// trans-weighted depths can be computed by hand, which is what
// is done here. The general (non-TPFA) case is not supported yet. 
void TwoPhaseHybridFVM::PrecomputeData( DomainPartition * const domain )
{
  localIndex constexpr maxNumFaces = MAX_NUM_FACES;
  
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  NodeManager const * const nodeManager          = mesh->getNodeManager();
  FaceManager const * const faceManager          = mesh->getFaceManager();

  // node data (for transmissibility computation)

  arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition =
    nodeManager->referencePosition();

  // get the face-to-nodes connectivity for the transmissibility calculation
  ArrayOfArraysView<localIndex const> const & faceToNodes = faceManager->nodeList();
  
  // get the face-centered depth
  arrayView1d<real64> const & weightedFaceGravCoef =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::weightedGravityCoefString);
  arrayView1d<real64> const & sumTrans =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::sumTransmissibilityString);
  weightedFaceGravCoef = 0;
  sumTrans = 0;
  
  // tolerance for transmissibility calculation
  real64 const lengthTolerance = domain->getMeshBody( 0 )->getGlobalLengthScale() * m_areaRelTol; 
  
  elemManager->
    forElementSubRegionsComplete<CellElementSubRegion,
                                 FaceElementSubRegion>( m_targetRegions,
                                                      [&]( localIndex const GEOSX_UNUSED_PARAM( er ),
                                                           localIndex const GEOSX_UNUSED_PARAM( esr ),
                                                           ElementRegionBase const * const,
                                                           auto const * const subRegion )
  {
    // get the map from elem to faces
    arrayView2d< localIndex const > const & elemToFaces = subRegion->faceList();

    // get the elem-centered depth 
    arrayView1d<real64 const> const & elemGravCoef =
      subRegion->template getReference<array1d<real64>>(viewKeyStruct::gravityCoefString);
    
    // get the element data needed for transmissibility computation
    arrayView1d<R1Tensor const> const & elemCenter =
     subRegion->template getReference< array1d<R1Tensor> >( CellBlock::viewKeyStruct::elementCenterString );
    arrayView1d<real64 const> const & elemVolume =
     subRegion->template getReference< array1d<real64> >( CellBlock::viewKeyStruct::elementVolumeString ); 
    arrayView1d<R1Tensor const> const & elemPerm =
     subRegion->template getReference< array1d<R1Tensor> >( viewKeyStruct::permeabilityString ); 
   
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      localIndex const numFacesInElem = elemToFaces.size(1);

      // transmissibility matrix
      stackArray2d<real64, maxNumFaces*maxNumFaces> transMatrix( numFacesInElem,
                                                                 numFacesInElem );
      
      ComputeTransmissibilityMatrix( nodePosition,     
                                     faceToNodes,        
                                     elemToFaces[ei],
                                     elemCenter[ei],
                                     elemVolume[ei],
                                     elemPerm[ei],
                                     lengthTolerance,
                                     transMatrix );

      for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
      {
        localIndex faceId = elemToFaces[ei][ifaceLoc];
        sumTrans[faceId] += transMatrix[ifaceLoc][ifaceLoc];
        weightedFaceGravCoef[faceId] += transMatrix[ifaceLoc][ifaceLoc] * elemGravCoef[ei];
      }
    });
  });

  forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  {
    weightedFaceGravCoef[iface] /= sumTrans[iface];
  });  

}



// this function is obviously redundant with computeCellStencil in the TwoPointFluxApproximation class
// this is here for now, but I will have to find a better place for this type of function at some point 
void TwoPhaseHybridFVM::ComputeTransmissibilityMatrix( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition, 
                                                       ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                                       arraySlice1d<localIndex const> const elemToFaces,
                                                       R1Tensor const & elemCenter,
                                                       real64   const & elemVolume,
                                                       R1Tensor const & elemPerm,
                                                       real64   const & lengthTolerance,
                                                       stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> & transMatrix ) const 
{
  switch (m_ipType)
  {
    case InnerProductType::TPFA:
    {
      ComputeTPFAInnerProduct( nodePosition,
                               faceToNodes,
                               elemToFaces,
                               elemCenter,
                               elemPerm,
                               lengthTolerance,
                               transMatrix );
      break;
    }

    case InnerProductType::QUASI_TPFA:
    {
      // for now, Q-TPFA is useful for debugging the IP computation
      // since it reduces to TPFA on orthogonal meshes...
      ComputeQFamilyInnerProduct( nodePosition,
                                  faceToNodes,
                                  elemToFaces,
                                  elemCenter,
                                  elemVolume,
                                  elemPerm,
                                  2,   
                                  lengthTolerance,
                                  transMatrix );
      break;
    }

    default:
    {
      GEOSX_LOG_RANK("Unknown inner product in SinglePhaseHybridFVM");
      break;
    }
  }
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
void TwoPhaseHybridFVM::ComputeTPFAInnerProduct( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition, 
                                                 ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                                 arraySlice1d<localIndex const> const elemToFaces,
                                                 R1Tensor const & elemCenter,
                                                 R1Tensor const & elemPerm,
                                                 real64   const & lengthTolerance,
                                                 stackArray2d<real64, MAX_NUM_FACES*MAX_NUM_FACES> const & transMatrix ) const
{
  R1Tensor faceCenter, faceNormal, faceConormal, elemToFaceVec;
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
        
        // 1) compute the face geometry data: center, normal, vector from elem center to face center
        real64 const faceArea =
          computationalGeometry::Centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                     faceToNodes.sizeOfArray(elemToFaces[ifaceLoc]),
                                                     nodePosition,
                                                     faceCenter,
                                                     faceNormal,
                                                     areaTolerance );
        
        elemToFaceVec  = faceCenter;
        elemToFaceVec -= elemCenter;

        if (Dot(elemToFaceVec, faceNormal) < 0.0)
        {
          faceNormal *= -1;
        }

        real64 const c2fDistance = elemToFaceVec.Normalize();

        // 2) assemble full coefficient tensor from principal axis/components
        makeFullTensor(elemPerm, permeabilityTensor);

        faceConormal.AijBj(permeabilityTensor, faceNormal);

        // 3) compute the one-sided face transmissibility
        transMatrix[ifaceLoc][jfaceLoc]  = Dot(elemToFaceVec, faceConormal);
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


void TwoPhaseHybridFVM::ComputeQFamilyInnerProduct( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition,
                                                    ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                                    arraySlice1d<localIndex const> const elemToFaces,
                                                    R1Tensor const & elemCenter,
                                                    real64   const & elemVolume,
                                                    R1Tensor const & elemPerm,
                                                    real64   const & tParam, 
                                                    real64   const & lengthTolerance,
                                                    stackArray2d<real64, MAX_NUM_FACES
                                                                        *MAX_NUM_FACES> const & transMatrix ) const
{
  R1Tensor faceCenter, faceNormal, cellToFaceVec;
  real64 const areaTolerance = lengthTolerance * lengthTolerance;

  localIndex const numFacesInElem = elemToFaces.size();
  localIndex const dim = 3;
  
  // TODO: remove all the array2ds of this function once the BlasLapackLA calls have been removed
  //       work with stackArray2ds instead 
  array2d<real64> cellToFaceMat;
  array2d<real64> normalsMat( numFacesInElem, dim );
  array2d<real64> permMat( dim, dim );
  array2d<real64> transMat( numFacesInElem, numFacesInElem );
  
  // TODO: figure out if it is possible/beneficial to preallocate these arrays
  array1d<real64> work_dim;
  array2d<real64> work_dimByDim; 
  array2d<real64> work_dimByNumFaces( dim, numFacesInElem );
  array2d<real64> work_numFacesByDim( numFacesInElem, dim );
  array2d<real64> worka_numFacesByNumFaces( numFacesInElem, numFacesInElem );
  array2d<real64> workb_numFacesByNumFaces( numFacesInElem, numFacesInElem );
  array2d<real64> workc_numFacesByNumFaces( numFacesInElem, numFacesInElem );

  array1d<real64> q0; 
  array1d<real64> q1; 
  array1d<real64> q2; 

  if (m_orthonormalizeWithSVD)
  {
    cellToFaceMat.resizeDimension<0>( numFacesInElem );
    cellToFaceMat.resizeDimension<1>( dim );
    work_dim.resize( dim );
    work_dimByDim.resize( dim, dim );
  }
  else
  {
    q0.resize( numFacesInElem );
    q1.resize( numFacesInElem );
    q2.resize( numFacesInElem );
  }
  
  // 1) fill the matrices cellToFaceMat and normalsMat row by row 
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {

    // compute the face geometry data: center, normal, vector from cell center to face center
    real64 const faceArea =
      computationalGeometry::Centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                 faceToNodes.sizeOfArray(elemToFaces[ifaceLoc]),
                                                 nodePosition,
                                                 faceCenter,
                                                 faceNormal,
                                                 areaTolerance );
        
    cellToFaceVec  = faceCenter;
    cellToFaceVec -= elemCenter;

    if (m_orthonormalizeWithSVD)
    {  
      cellToFaceMat(ifaceLoc,0) = cellToFaceVec(0);
      cellToFaceMat(ifaceLoc,1) = cellToFaceVec(1);
      cellToFaceMat(ifaceLoc,2) = cellToFaceVec(2);
    }
    else
    {
      q0(ifaceLoc) = cellToFaceVec(0);
      q1(ifaceLoc) = cellToFaceVec(1);
      q2(ifaceLoc) = cellToFaceVec(2);
    }
    
    if (Dot(cellToFaceVec, faceNormal) < 0.0)
    {
      faceNormal *= -1;
    }

    faceNormal.Normalize();
    
    normalsMat(ifaceLoc,0) = faceArea*faceNormal(0);
    normalsMat(ifaceLoc,1) = faceArea*faceNormal(1);
    normalsMat(ifaceLoc,2) = faceArea*faceNormal(2);
    
  }

  // 2) assemble full coefficient tensor from principal axis/components
  // TODO: figure out if there is a better way to that 
  R2SymTensor permeabilityTensor;
  makeFullTensor(elemPerm, permeabilityTensor);
  for (localIndex i = 0; i < dim; ++i)
  {
    for (localIndex j = 0; j < dim; ++j)
    {
      permMat(i,j) = permeabilityTensor(i,j);
    }
  }

  // TODO: replace the BlasLapack calls below with explicitly for loops
  //       this should be easy if MGS orthonormalization is as robust as SVD
  
  // 3) compute N K N' 
  BlasLapackLA::matrixMatrixTMultiply( permMat,
                                       normalsMat,
                                       work_dimByNumFaces );
  BlasLapackLA::matrixMatrixMultiply( normalsMat,
                                      work_dimByNumFaces,
                                      transMat );

  // 4) compute the orthonormalization of the matrix cellToFaceVec 
  //    This is done either with SVD or MGS
  //    If we find that MGS is stable, I will remove SVD

  if (m_orthonormalizeWithSVD)
  {
    // calling SVD seems to be an overkill to orthonormalize the 3 columns of cellToFaceMat...
    BlasLapackLA::matrixSVD( cellToFaceMat,
                             work_numFacesByDim,
                             work_dim,
                             work_dimByDim );
  }
  else
  {
    // q0
    BlasLapackLA::vectorScale( 1.0/BlasLapackLA::vectorNorm2( q0 ), q0 );

    // q1
    real64 const q0Dotq1 = BlasLapackLA::vectorDot( q0, q1 );
    BlasLapackLA::vectorVectorAdd( q0, q1, -q0Dotq1 );
    BlasLapackLA::vectorScale( 1.0/BlasLapackLA::vectorNorm2( q1 ), q1 );
  
    // q2
    real64 const q0Dotq2 = BlasLapackLA::vectorDot( q0, q2 );
    BlasLapackLA::vectorVectorAdd( q0, q2, -q0Dotq2 );
    real64 const q1Dotq2 = BlasLapackLA::vectorDot( q1, q2 );
    BlasLapackLA::vectorVectorAdd( q1, q2, -q1Dotq2 );
    BlasLapackLA::vectorScale( 1.0/BlasLapackLA::vectorNorm2( q2 ), q2 );

    // TODO: remove the copies once the BlasLapackLA calls have been removed
    for (int i = 0; i < numFacesInElem; ++i)
    {  
      work_numFacesByDim(i,0) = q0(i);
      work_numFacesByDim(i,1) = q1(i);
      work_numFacesByDim(i,2) = q2(i);
    }
  }
  
  // 5) compute P_Q = I - QQ'
  BlasLapackLA::matrixMatrixTMultiply( work_numFacesByDim,
                                       work_numFacesByDim,
                                       worka_numFacesByNumFaces );
  BlasLapackLA::matrixScale( -1, worka_numFacesByNumFaces );
  for (localIndex i = 0; i < numFacesInElem; ++i)
  {
    worka_numFacesByNumFaces(i,i)++;
  }

  // 6) compute P_Q D P_Q where D = diag(diag(N K N')  
  for (localIndex i = 0; i < numFacesInElem; ++i)
  {
    for (localIndex j = 0; j < numFacesInElem; ++j)
    {
      workb_numFacesByNumFaces(i,j) = (i == j ) ?  transMat(i,j) : 0.0;
    }
  }
  BlasLapackLA::matrixMatrixMultiply( workb_numFacesByNumFaces,
                                      worka_numFacesByNumFaces,
                                      workc_numFacesByNumFaces );
  BlasLapackLA::matrixMatrixMultiply( worka_numFacesByNumFaces,
                                      workc_numFacesByNumFaces,
                                      workb_numFacesByNumFaces );
  
  // 7) compute T = ( N K N' + t U diag(diag(N K N')) U ) / elemVolume
  BlasLapackLA::matrixScale( tParam, workb_numFacesByNumFaces );
  BlasLapackLA::matrixMatrixAdd( workb_numFacesByNumFaces, transMat );
  BlasLapackLA::matrixScale( 1/elemVolume, transMat );
  
  // for now, I have this copy to transfer the data from the array2d to the stackArray2d
  // I need the array2d to call the BlasLapackLA functions
  // if and when everything works with explicit for loops in this kernel function,
  // I will be able to do remove all the array2ds and then I won't need this copy anymore
  for (localIndex i = 0; i < numFacesInElem; ++i)
  {
    for (localIndex j = 0; j < numFacesInElem; ++j)
    {
      transMatrix(i,j) = transMat(i,j);
    }
  }

}

REGISTER_CATALOG_ENTRY( SolverBase, TwoPhaseHybridFVM, std::string const &, Group * const )
} /* namespace geosx */
