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
 * @file SinglePhaseMimetic.cpp
 */

#include "SinglePhaseMimetic.hpp"

#include "common/TimingMacros.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFlowKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseFlowKernels;

SinglePhaseMimetic::SinglePhaseMimetic( const std::string& name,
                                        Group * const parent ):
  SinglePhaseFlowBase(name, parent),
  m_faceDofKey(""),
  m_areaRelTol(1e-8)
{
  // one cell-centered dof per cell
  m_numDofPerCell = 1;
}

  
void SinglePhaseMimetic::RegisterDataOnMesh(Group * const MeshBodies)
{

  // 1) Register the cell-centered data
  SinglePhaseFlowBase::RegisterDataOnMesh(MeshBodies);

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

  // one-sided face-based auxiliary data
  RegisterOneSidedFaceData(MeshBodies);
  
}

void SinglePhaseMimetic::RegisterOneSidedFaceData( Group * const MeshBodies )
{

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase * const subRegion )
    {
      // one-sided face-based auxiliary data
      subRegion->registerWrapper< array2d<globalIndex> >( viewKeyStruct::neighborDofNumberString );
      subRegion->registerWrapper< array2d<localIndex> >( viewKeyStruct::neighborRegionIdString );
      subRegion->registerWrapper< array2d<localIndex> >( viewKeyStruct::neighborSubRegionIdString );
      subRegion->registerWrapper< array2d<localIndex> >( viewKeyStruct::neighborElemIdString );
    });   
  }
}
  

void SinglePhaseMimetic::ResizeOneSidedFaceFields( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  
  // resize the one-sided face-based arrays using the number of faces per element in the subregion
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    localIndex const numFacesPerElem = subRegion->numFacesPerElement();
    
    subRegion->getReference< array2d<globalIndex> >( viewKeyStruct::neighborDofNumberString ).resizeDimension<1>(numFacesPerElem);
    subRegion->getReference< array2d<localIndex> >( viewKeyStruct::neighborRegionIdString ).resizeDimension<1>(numFacesPerElem);
    subRegion->getReference< array2d<localIndex> >( viewKeyStruct::neighborSubRegionIdString ).resizeDimension<1>(numFacesPerElem);
    subRegion->getReference< array2d<localIndex> >( viewKeyStruct::neighborElemIdString ).resizeDimension<1>(numFacesPerElem);
  });
}


void SinglePhaseMimetic::ConstructOneSidedFaceMaps( DomainPartition * const domain,
                                                    DofManager const & dofManager )
{
  // the cell-centered solvers have the stencil object, maybe the mimetic method can have its own object to store these maps  
  
  MeshLevel * const mesh                   = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager const * const faceManager    = mesh->getFaceManager();

  array2d<localIndex> const & elemRegionList    = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList          = faceManager->elementList(); 

  // in this function we need to make sure that we act only on the target regions
  // for that, we need the following region filter
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }
  
  // in this loop we collect the one-sided faces corresponding to an element of the target regions
  // TODO: figure out what to do with fracture elements (can we do everything in one loop?) 
  elemManager->
    forElementSubRegionsComplete<CellElementSubRegion>( m_targetRegions,
                                                      [&]( localIndex const er,
                                                           localIndex const esr,
                                                           ElementRegionBase * const,
                                                           CellElementSubRegion * const subRegion )
  {
    arrayView2d< localIndex const > const & elemToFaces = subRegion->faceList();

    // get the "non-local" maps that relate each one-sided face to the neighbor data
    arrayView2d<globalIndex> const & neighborDofNumber =
      subRegion->getReference<array2d<globalIndex>>( viewKeyStruct::neighborDofNumberString );
    arrayView2d<localIndex> const & neighborRegionId =
      subRegion->getReference<array2d<localIndex>>( viewKeyStruct::neighborRegionIdString );
    arrayView2d<localIndex> const & neighborSubRegionId =
      subRegion->getReference<array2d<localIndex>>( viewKeyStruct::neighborSubRegionIdString );
    arrayView2d<localIndex> const & neighborElemId =
      subRegion->getReference<array2d<localIndex>>( viewKeyStruct::neighborElemIdString );

    // cell-centered dof numbers   
    string const elemDofKey = dofManager.getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & elemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( elemDofKey );  

    for (localIndex ei = 0; ei < subRegion->size(); ++ei)
    {
      
      // for each element, visit all the one-sided faces
      for (localIndex ifaceLoc=0; ifaceLoc<elemToFaces.size(1); ++ifaceLoc)
      {
        
        // get the ID of the face
        localIndex const iface = elemToFaces[ei][ifaceLoc];

        // One-sided face-based maps (used to assemble the mass one-sided fluxes)

        // the face has at most two adjacent elements
        // one of these two elements is the current element indexed by er, esr, ei
        // but here, we are interested in saving the indices of the other element
        // this other element is "the neighbor" for this one-sided face 
        for (localIndex k=0; k<elemRegionList.size(1); ++k)
        {
          
          localIndex const erNeighbor  = elemRegionList[iface][k];
          localIndex const esrNeighbor = elemSubRegionList[iface][k];
          localIndex const eiNeighbor  = elemList[iface][k];

          // this element is not the current element
          // we have found the neighbor or we are at the boundary 
          if ( erNeighbor != er || esrNeighbor != esr || eiNeighbor != ei ) 
          {

            // save the neighbor info
            neighborRegionId[ei][ifaceLoc]    = erNeighbor;
            neighborSubRegionId[ei][ifaceLoc] = esrNeighbor;
            neighborElemId[ei][ifaceLoc]      = eiNeighbor;

            bool const onBoundary       = (erNeighbor == -1 || esrNeighbor == -1 || eiNeighbor != -1);
            bool const neighborInTarget = regionFilter.contains(erNeighbor); 
      
            // if not on boundary, save the dof number of the neighbor element (if in target)
            if ( !onBoundary && neighborInTarget )
            {
              ElementRegionBase const * const neighborRegion =
                Group::group_cast<ElementRegionBase const *>(mesh->getElemManager()->GetRegion(erNeighbor));
              ElementSubRegionBase const * const neighborSubRegion =
                Group::group_cast<ElementSubRegionBase const *>(neighborRegion->GetSubRegion(esrNeighbor));

              arrayView1d<globalIndex const> const & neighborElemDofNumber =
                neighborSubRegion->getReference<array1d<globalIndex>>( elemDofKey );  
              
              neighborDofNumber[ei][ifaceLoc] = neighborElemDofNumber[eiNeighbor];
            }
            // if on boundary
            else
            {
              // for the case of a boundary face, we assign the "local" dof number
              // to neighborDofNumber. Then, in the assembly, we add a zero derivative
              // in this slot of the Jacobian matrix
              neighborDofNumber[ei][ifaceLoc] = elemDofNumber[ei];      
            }
          }

        } // end of loop on elems adjacent to the face

      } // end of loop on one-sided faces of the elem
      
    } // end of loop on elements of the subregion

  }); // end of loop on subregions
  
}

  
void SinglePhaseMimetic::ImplicitStepSetup( real64 const & time_n,
                                            real64 const & dt,
                                            DomainPartition * const domain,
                                            DofManager & dofManager,
                                            ParallelMatrix & matrix,
                                            ParallelVector & rhs,
                                            ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  // setup the cell-centered fields
  SinglePhaseFlowBase::ImplicitStepSetup( time_n, dt, domain, dofManager, matrix, rhs, solution );
  
  // computes the number of one-sided faces and resizes the one-sided face fields accordingly
  // if no topological change in the mesh, this does not have to be recomputed 
  ResizeOneSidedFaceFields( domain );
  
  // build all the maps that are necessary to iterate over one-sided faces
  // if no topological change in the mesh, this does not have to be recomputed 
  ConstructOneSidedFaceMaps( domain, dofManager );

  // setup the face fields
  MeshLevel * const meshLevel     = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = meshLevel->getFaceManager();

  // get the face-based DOF numbers 
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString );
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference< array1d<globalIndex> >( faceDofKey );  

  // save the face Dof key for use in two functions
  // that do not have acces to the coupled solver dofManager
  // namely ResetStateToBeginningOfStep and ImplicitStepComplete
  m_faceDofKey = faceDofKey;
  
  // get the accumulated pressure updates
  arrayView1d<real64> & dFacePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);
 
  forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  {
    // zero out if face is in target region
    if (faceDofNumber[iface] >= 0)
    {
      dFacePres[iface] = 0;
    }
  });

}

void SinglePhaseMimetic::ImplicitStepComplete( real64 const & time_n,
                                               real64 const & dt,
                                               DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  // increment the cell-centered fields
  SinglePhaseFlowBase::ImplicitStepComplete( time_n, dt, domain );

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

void SinglePhaseMimetic::SetupDofs( DomainPartition const * const GEOSX_UNUSED_ARG( domain ),
                                    DofManager & dofManager ) const
{
  
  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleUpwindedOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::pressureString,
                       DofManager::Location::Elem,
                       DofManager::Connectivity::Face, 
                       m_targetRegions );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::facePressureString,
                       DofManager::Location::Face,
                       DofManager::Connectivity::Elem,
                       m_targetRegions );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::facePressureString,
                          viewKeyStruct::pressureString,
                          DofManager::Connectivity::Elem,
                          true );

}

void SinglePhaseMimetic::AssembleFluxTerms( real64 const GEOSX_UNUSED_ARG( time_n ),
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

  // node data (for transmissibility computation)

  arrayView1d<R1Tensor const> const & nodePosition = nodeManager->referencePosition();

  
  // face data
  
  // get the face-based DOF numbers for the assembly
  string const faceDofKey = dofManager->getKey( viewKeyStruct::facePressureString );
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
 
  // max number of faces allowed in an element 
  localIndex constexpr maxNumFaces = MAX_NUM_FACES_IN_ELEM;

  // tolerance for transmissibility calculation
  real64 const lengthTolerance = domain->getMeshBody( 0 )->getGlobalLengthScale() * m_areaRelTol; 


  elemManager->
    forElementSubRegionsComplete<CellElementSubRegion>( m_targetRegions,
                                                      [&]( localIndex const er,
                                                           localIndex const esr,
                                                           ElementRegionBase const * const,
                                                           CellElementSubRegion const * const subRegion )
  {

    // elem data
    
    // get the cell-centered DOF numbers and ghost rank for the assembly
    string const elemDofKey = dofManager->getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & elemDofNumber =
      subRegion->getReference< array1d<globalIndex> >( elemDofKey ); 
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];   

    // get the indices of the neighbor
    // this is needed for the upwinding only
    // this could also be computed on the fly in each element   
    arrayView2d<localIndex const> const & neighborRegionId =
      subRegion->getReference< array2d<localIndex> >( viewKeyStruct::neighborRegionIdString );
    arrayView2d<localIndex const> const & neighborSubRegionId =
      subRegion->getReference< array2d<localIndex> >( viewKeyStruct::neighborSubRegionIdString );
    arrayView2d<localIndex const> const & neighborElemId =
      subRegion->getReference< array2d<localIndex> >( viewKeyStruct::neighborElemIdString );
    arrayView2d<globalIndex const> const & neighborDofNumber = 
      subRegion->getReference< array2d<globalIndex> >( viewKeyStruct::neighborDofNumberString );

    // get the map from elem to faces
    arrayView2d< localIndex const > const & elemToFaces = subRegion->faceList();
    
    // get the cell-centered pressures
    arrayView1d<real64 const> const & elemPres  = m_pressure[er][esr];
    arrayView1d<real64 const> const & dElemPres = m_deltaPressure[er][esr];
    
    // get the cell centered densities
    arrayView2d<real64 const> const & elemDens     = m_density[er][esr][m_fluidIndex];
    arrayView2d<real64 const> const & dElemDens_dp = m_dDens_dPres[er][esr][m_fluidIndex];
    
    // get the element data needed for transmissibility computation
    arrayView1d<R1Tensor const> const & elemCenter =
     subRegion->getReference< array1d<R1Tensor> >( CellBlock::viewKeyStruct::elementCenterString ); 
    arrayView1d<R1Tensor const> const & elemPerm =
     subRegion->getReference< array1d<R1Tensor> >( viewKeyStruct::permeabilityString ); 

    // get the cell-centered depth 
    arrayView1d<real64> const & elemGravDepth =
      subRegion->getReference<array1d<real64>>(viewKeyStruct::gravityDepthString);

    
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
        stackArray1d<real64, maxNumFaces> oneSidedVolFlux( numFacesInElem );
        stackArray1d<real64, maxNumFaces> dOneSidedVolFlux_dp( numFacesInElem );
        stackArray1d<real64, maxNumFaces> dOneSidedVolFlux_dfp( numFacesInElem );
      
        // upwinded mobility
        stackArray1d<real64, maxNumFaces> upwMobility( numFacesInElem );
        stackArray1d<real64, maxNumFaces> dUpwMobility_dp( numFacesInElem );
        stackArray1d<globalIndex, maxNumFaces> upwDofNumber( numFacesInElem ); 
      
        // recompute the local transmissibility matrix at each iteration
        // we can decide later to precompute transMatrix if needed
        ComputeTransmissibilityMatrix( nodePosition,     
                                       faceToNodes,        
                                       elemToFaces[ei],
                                       elemCenter[ei],     
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
                                  elemDens[ei][0],
                                  dElemDens_dp[ei][0],
                                  transMatrix,
                                  oneSidedVolFlux,             
                                  dOneSidedVolFlux_dp,
                                  dOneSidedVolFlux_dfp );

      
        // at this point, we know the local flow direction in the element
        // so we can upwind the transport coefficients (mobilities) at the one sided faces 
        // ** this function needs non-local information **
        UpdateUpwindedCoefficients( neighborRegionId[ei], 
                                    neighborSubRegionId[ei],
                                    neighborElemId[ei],
                                    neighborDofNumber[ei],
                                    m_mobility, 
                                    m_dMobility_dPres,
                                    m_mobility[er][esr][ei],
                                    m_dMobility_dPres[er][esr][ei],
                                    elemDofNumber[ei],
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
        AssembleOneSidedMassFluxes( dt,
                                    faceDofNumber,
                                    elemToFaces[ei],
                                    elemDofNumber[ei],
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
        AssembleConstraints( faceDofNumber,
                             elemToFaces[ei],
                             elemDofNumber[ei],
                             oneSidedVolFlux,
                             dOneSidedVolFlux_dp,
                             dOneSidedVolFlux_dfp,
                             matrix,
                             rhs );

      }
    });
  });
}

void SinglePhaseMimetic::ComputeOneSidedVolFluxes( arrayView1d<real64 const> const & facePres,
                                                   arrayView1d<real64 const> const & dFacePres,
                                                   arrayView1d<real64 const> const & faceGravDepth,
                                                   arraySlice1d<localIndex const> const elemToFaces,
                                                   real64 const & elemPres,
                                                   real64 const & dElemPres,
                                                   real64 const & elemGravDepth,
                                                   real64 const & elemDens,
                                                   real64 const & dElemDens_dp,
                                                   stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                                       *MAX_NUM_FACES_IN_ELEM> const & transMatrix,
                                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & oneSidedVolFlux,
                                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dOneSidedVolFlux_dp,
                                                   stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dOneSidedVolFlux_dfp ) const
{
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
      // 1) compute the potential diff between the cell center and the face center
      real64 const ccPres = elemPres + dElemPres;
      real64 const fPres  = facePres[elemToFaces[jfaceLoc]] + dFacePres[elemToFaces[jfaceLoc]];

      real64 const ccGravDepth = elemGravDepth;
      real64 const fGravDepth  = faceGravDepth[elemToFaces[jfaceLoc]];

      real64 const ccDens     = elemDens;
      real64 const dCcDens_dp = dElemDens_dp; 
      // no density evaluated at the face center
          
      // pressure difference
      real64 const presDif      = ccPres - fPres;
      real64 const dPresDif_dp  =  1;
      real64 const dPresDif_dfp = -1;

      // gravity term
      real64 const depthDif     = ccGravDepth - fGravDepth; 
      real64 const gravTerm     = ccDens     * depthDif;
      real64 const dGravTerm_dp = dCcDens_dp * depthDif;

      // potential difference
      real64 const potDif      = presDif     - gravTerm;
      real64 const dPotDif_dp  = dPresDif_dp - dGravTerm_dp;    
      real64 const dPotDif_dfp = dPresDif_dfp; 

      // 2) compute the contribution of this face to the volumetric fluxes in the cell
      oneSidedVolFlux[ifaceLoc]      += transMatrix[ifaceLoc][jfaceLoc] * potDif;
      dOneSidedVolFlux_dp[ifaceLoc]  += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dp;
      dOneSidedVolFlux_dfp[ifaceLoc] += transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dfp;
    }
  }  
}

void SinglePhaseMimetic::UpdateUpwindedCoefficients( arraySlice1d<localIndex const> const neighborRegionId,
                                                     arraySlice1d<localIndex const> const neighborSubRegionId,
                                                     arraySlice1d<localIndex const> const neighborElemId,
                                                     arraySlice1d<globalIndex const> const neighborDofNumber,
                                                     ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & domainMobility,
                                                     ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dDomainMobility_dp,
                                                     real64 const & elemMobility,
                                                     real64 const & dElemMobility_dp,
                                                     globalIndex const elemDofNumber,
                                                     stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                                     stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & upwMobility,
                                                     stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> & dUpwMobility_dp,
                                                     stackArray1d<globalIndex, MAX_NUM_FACES_IN_ELEM> & upwDofNumber ) const
{
  localIndex const numFacesInElem = neighborElemId.size();
  
  // for this element, loop over the local (one-sided) faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    bool const isBoundaryFace = (neighborElemId[ifaceLoc] < 0);
   
    // if the local element is upwind
    if (oneSidedVolFlux[ifaceLoc] >= 0 || isBoundaryFace)
    {
      upwMobility[ifaceLoc]     = elemMobility;
      dUpwMobility_dp[ifaceLoc] = dElemMobility_dp;
      upwDofNumber[ifaceLoc]    = elemDofNumber;
    }
    // else the neighbor is upwind
    else
    {
      // here, instead of using the precomputed maps,
      // we could recompute the indices on the fly
      localIndex const erNeighbor  = neighborRegionId[ifaceLoc];
      localIndex const esrNeighbor = neighborSubRegionId[ifaceLoc];
      localIndex const eiNeighbor  = neighborElemId[ifaceLoc];  
          
      upwMobility[ifaceLoc]     = domainMobility[erNeighbor][esrNeighbor][eiNeighbor];
      dUpwMobility_dp[ifaceLoc] = dDomainMobility_dp[erNeighbor][esrNeighbor][eiNeighbor];
      upwDofNumber[ifaceLoc]    = neighborDofNumber[ifaceLoc];
    }
  }   
}


void SinglePhaseMimetic::AssembleOneSidedMassFluxes( real64 const & dt,
                                                     arrayView1d<globalIndex const> const & faceDofNumber,
                                                     arraySlice1d<localIndex const> const elemToFaces,
                                                     globalIndex const elemDofNumber,
                                                     stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                                     stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dp,
                                                     stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dfp,
                                                     stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & upwMobility,
                                                     stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dUpwMobility_dp,
                                                     stackArray1d<globalIndex, MAX_NUM_FACES_IN_ELEM> const & upwDofNumber,
                                                     ParallelMatrix * const matrix,
                                                     ParallelVector * const rhs ) const
{
  localIndex constexpr maxNumFaces = MAX_NUM_FACES_IN_ELEM;

  localIndex const numFacesInElem = elemToFaces.size();
  
  // fluxes
  real64 sumLocalMassFluxes     = 0;
  real64 dSumLocalMassFluxes_dp = 0;
  stackArray1d<real64, maxNumFaces> dSumLocalMassFluxes_dp_neighbor( numFacesInElem );
  stackArray1d<real64, maxNumFaces> dSumLocalMassFluxes_dfp( numFacesInElem );
  dSumLocalMassFluxes_dp_neighbor = 0;
  dSumLocalMassFluxes_dfp = 0;
          
  // dof numbers 
  globalIndex const eqnRowIndex         = elemDofNumber;
  globalIndex const dofColIndexElemPres = elemDofNumber;
  stackArray1d<globalIndex, maxNumFaces> dofColIndicesNeighborPres( numFacesInElem );
  stackArray1d<globalIndex, maxNumFaces> dofColIndicesFacePres( numFacesInElem );

  // for each element, loop over the one-sided faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    // compute the mass flux at the one-sided face plus its derivatives
    // add the newly computed flux to the sum
    sumLocalMassFluxes      += dt * upwMobility[ifaceLoc] * oneSidedVolFlux[ifaceLoc];
    dSumLocalMassFluxes_dp  += dt * upwMobility[ifaceLoc] * dOneSidedVolFlux_dp[ifaceLoc];
    if (elemDofNumber == upwDofNumber[ifaceLoc])
    {
      dSumLocalMassFluxes_dp  += dt * dUpwMobility_dp[ifaceLoc] * oneSidedVolFlux[ifaceLoc];
    }
    else
    {
      dSumLocalMassFluxes_dp_neighbor[ifaceLoc] = dt * dUpwMobility_dp[ifaceLoc] * oneSidedVolFlux[ifaceLoc];
    }
    dSumLocalMassFluxes_dfp[ifaceLoc] = dt * upwMobility[ifaceLoc] * dOneSidedVolFlux_dfp[ifaceLoc];

    // collect the relevant dof numbers
    dofColIndicesNeighborPres[ifaceLoc] = upwDofNumber[ifaceLoc]; // if upwDofNumber == elemDofNumber, the derivative is zero 
    dofColIndicesFacePres[ifaceLoc] = faceDofNumber[elemToFaces[ifaceLoc]];

  }
  
  // we are ready to assemble the local flux and its derivatives

  // residual
  rhs->add( &eqnRowIndex,
            &sumLocalMassFluxes,
            1 );
                  
  // jacobian -- derivative wrt local cell centered pressure term
  matrix->add( &eqnRowIndex,
               &dofColIndexElemPres,
               &dSumLocalMassFluxes_dp,
               1,
               1 );

  // jacobian -- derivatives wrt neighbor cell centered pressure terms
  matrix->add( &eqnRowIndex,
               dofColIndicesNeighborPres.data(),
               dSumLocalMassFluxes_dp_neighbor.data(),
               1,
               numFacesInElem );
        
  // jacobian -- derivatives wrt face pressure terms
  matrix->add( &eqnRowIndex,
               dofColIndicesFacePres.data(),
               dSumLocalMassFluxes_dfp.data(),
               1,
               numFacesInElem );

}


void SinglePhaseMimetic::AssembleConstraints( arrayView1d<globalIndex const> const & faceDofNumber,
                                              arraySlice1d<localIndex const> const elemToFaces,
                                              globalIndex const elemDofNumber,
                                              stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & oneSidedVolFlux,
                                              stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dp,
                                              stackArray1d<real64, MAX_NUM_FACES_IN_ELEM> const & dOneSidedVolFlux_dfp,
                                              ParallelMatrix * const matrix,
                                              ParallelVector * const rhs ) const 
{
  localIndex const numFacesInElem = elemToFaces.size();
  
  globalIndex const dofColIndexElemPres = elemDofNumber; 
  
  // for each element, loop over the local (one-sided) faces
  for (localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc)
  {
    // dof numbers 
    globalIndex const eqnRowIndex         = faceDofNumber[elemToFaces[ifaceLoc]];
    globalIndex const dofColIndexFacePres = faceDofNumber[elemToFaces[ifaceLoc]];

    // fluxes
    real64 const flux      = oneSidedVolFlux[ifaceLoc];
    real64 const dFlux_dp  = dOneSidedVolFlux_dp[ifaceLoc];
    real64 const dFlux_dfp = dOneSidedVolFlux_dfp[ifaceLoc];
    
    // residual
    rhs->add( &eqnRowIndex,
              &flux,
              1 );

    // jacobian -- derivative wrt local cell centered pressure term
    matrix->add( &eqnRowIndex,
                 &dofColIndexElemPres,
                 &dFlux_dp,
                 1, 1 );

    // jacobian -- derivatives wrt face pressure terms
    matrix->add( &eqnRowIndex,
                 &dofColIndexFacePres,
                 &dFlux_dfp,
                 1, 1 );

  }      
}
  

void
SinglePhaseMimetic::ApplyBoundaryConditions( real64 const GEOSX_UNUSED_ARG( time_n ),
                                             real64 const GEOSX_UNUSED_ARG( dt ),
                                             DomainPartition * const GEOSX_UNUSED_ARG( domain ),
                                             DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                                             ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                             ParallelVector & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  // TODO: implement boundary conditions the mimetic way

}

real64 SinglePhaseMimetic::CalculateResidualNorm( DomainPartition const * const domain,
                                                  DofManager const & dofManager,
                                                  ParallelVector const & rhs )
{
  MeshLevel const * const mesh          = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = mesh->getFaceManager();

  // here we compute the cell-centered residual norm in the derived class
  // to avoid duplicating a synchronization point 
  
  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  string const elemDofKey = dofManager.getKey( viewKeyStruct::pressureString );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString );
  
  // local residual
  array1d<real64> localResidualNorm( 2 ), globalResidualNorm( 2 );
  globalResidualNorm[EqType::MASS_CONS]  = 0; 
  globalResidualNorm[EqType::CONSTRAINT] = 0; 
  localResidualNorm[EqType::MASS_CONS]   = 0; 
  localResidualNorm[EqType::CONSTRAINT]  = 0; 

  
  // 1. Compute the residual for the mass conservation equations

  // compute the norm of local residual scaled by cell pore volume
  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<globalIndex const> const & elemDofNumber =
      subRegion->getReference< array1d<globalIndex> >( elemDofKey );

    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<real64 const> const & refPoro        = m_porosityRef[er][esr];
    arrayView1d<real64 const> const & volume         = m_volume[er][esr];
    arrayView1d<real64 const> const & dVol           = m_deltaVolume[er][esr];
    arrayView2d<real64 const> const & dens           = m_density[er][esr][m_fluidIndex];

    localIndex const subRegionSize = subRegion->size();
    for ( localIndex a = 0; a < subRegionSize; ++a )
    {
      if (elemGhostRank[a] < 0)
      {
        localIndex const lid = rhs.getLocalRowID( elemDofNumber[a] );
        // TODO: implement a normalization that matches the normalization used in SinglePhaseCellCentered
        real64 const val = localResidual[lid] / (refPoro[a] * dens[a][0] * ( volume[a] + dVol[a]));
        localResidualNorm[EqType::MASS_CONS] += val * val;
      }
    }
  });

  
  // 2. Compute the residual for the face-based constraints

  arrayView1d<integer const> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference< array1d<globalIndex> >( faceDofKey );

  // TODO: figure out how to use the forall loop
  //forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  for( localIndex iface = 0 ; iface < faceManager->size(); ++iface )
  {
    // if not ghost face and if adjacent to target region
    if (faceGhostRank[iface] < 0 && faceDofNumber[iface] >= 0)
    {
      localIndex const lid    = rhs.getLocalRowID( faceDofNumber[iface] );
      real64 const normalizer = 1; // TODO: compute the normalizer here
      real64 const val        = localResidual[lid] / normalizer;
      localResidualNorm[EqType::CONSTRAINT] += val * val;
    }
  }

  
  // 3. Combine the two norms

  MpiWrapper::allReduce(localResidualNorm.data(), globalResidualNorm.data(), 2, MPI_SUM, MPI_COMM_GEOSX);
 
  real64 maxNorm = (globalResidualNorm[EqType::MASS_CONS] < globalResidualNorm[EqType::CONSTRAINT])
                 ? globalResidualNorm[EqType::MASS_CONS]
                 : globalResidualNorm[EqType::CONSTRAINT];
  
  return sqrt( maxNorm );
}

void SinglePhaseMimetic::ApplySystemSolution( DofManager const & dofManager,
                                              ParallelVector const & solution,
                                              real64 const scalingFactor,
                                              DomainPartition * const domain )
{
  MeshLevel * const mesh          = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();

  // here we apply the cell-centered update in the derived class
  // to avoid duplicating a synchronization point 
  
  // 1. apply the cell-centered update

  applyToSubRegions( mesh, [&] ( localIndex const GEOSX_UNUSED_ARG( er ),
                                 localIndex const GEOSX_UNUSED_ARG( esr ),
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    dofManager.addVectorToField( solution,
                                 viewKeyStruct::pressureString,
                                 scalingFactor,
                                 subRegion,
                                 viewKeyStruct::deltaPressureString );
  } );

  // 2. apply the face-based update

  dofManager.addVectorToField( solution,
                               viewKeyStruct::facePressureString,
                               scalingFactor,
                               faceManager,
                               viewKeyStruct::deltaFacePressureString );

  // 3. synchronize
  
  std::map<string, string_array> fieldNames;
  fieldNames["faces"].push_back( viewKeyStruct::deltaFacePressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion )
  {
    UpdateState( subRegion );
  } );

}


void SinglePhaseMimetic::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // 1. Reset the cell-centered fields
  SinglePhaseFlowBase::ResetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  MeshLevel * const mesh          = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();

  // get the face-based DOF numbers 
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference< array1d<globalIndex> >( m_faceDofKey );  

  // get the accumulated face pressure updates
  arrayView1d<real64> & dFacePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

  forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  {
    if (faceDofNumber[iface] >= 0)
    {  
      dFacePres[iface] = 0;
    }
  });
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
void SinglePhaseMimetic::ComputeTransmissibilityMatrix( arrayView1d<R1Tensor const> const & nodePosition, 
                                                        ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                                        arraySlice1d<localIndex const> const elemToFaces,
                                                        R1Tensor const & elemCenter, 
                                                        R1Tensor const & elemPerm,
                                                        real64   const & lengthTolerance,
                                                        stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                                            *MAX_NUM_FACES_IN_ELEM> & transMatrix ) const 
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
        localIndex const iface = elemToFaces[ifaceLoc];
        
        // 1) compute the face geometry data: center, normal, vector from cell center to face center
        real64 const faceArea =
          computationalGeometry::Centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                     faceToNodes.sizeOfArray( iface ),
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

 
REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseMimetic, std::string const &, Group * const )
} /* namespace geosx */
