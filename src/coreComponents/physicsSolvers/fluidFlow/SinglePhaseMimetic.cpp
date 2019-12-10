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
  m_numOneSidedFaces(0),
  m_areaRelTol(1e-8)
{
  // TODO: decide what to do with m_numDofsPerCell here
}


void SinglePhaseMimetic::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  SinglePhaseFlowBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

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
 
    // face-based auxiliary data
    faceManager->registerWrapper<array1d<localIndex> >( viewKeyStruct::numAdjacentElementsString );
    faceManager->registerWrapper<array2d<localIndex> >( viewKeyStruct::faceToOneSidedFaceString );
    faceManager->registerWrapper<array2d<globalIndex> >( viewKeyStruct::faceToElemDofNumberString );    
    
  }

  // one-sided face-based auxiliary data
  RegisterOneSidedFaceData(MeshBodies);
  
}

void SinglePhaseMimetic::RegisterOneSidedFaceData( Group * const MeshBodies )
{

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    FaceManager * const faceManager = meshLevel->getFaceManager();
    
    // one-sided face-based auxiliary data
    faceManager->registerWrapper< array1d<localIndex> >( viewKeyStruct::oneSidedFaceToFaceString )->
      setSizedFromParent(0);

    faceManager->registerWrapper< array1d<globalIndex> >( viewKeyStruct::neighborDofNumberString )->
      setSizedFromParent(0);
    faceManager->registerWrapper< array1d<globalIndex> >( viewKeyStruct::neighborRegionIdString )->
      setSizedFromParent(0);
    faceManager->registerWrapper< array1d<globalIndex> >( viewKeyStruct::neighborSubRegionIdString )->
      setSizedFromParent(0);
    faceManager->registerWrapper< array1d<globalIndex> >( viewKeyStruct::neighborElemIdString )->
      setSizedFromParent(0);
    
    faceManager->registerWrapper< array1d<real64> >( viewKeyStruct::oneSidedVolFluxString )->
      setSizedFromParent(0);
    faceManager->registerWrapper< array1d<real64> >( viewKeyStruct::dOneSidedVolFlux_dPressureString )->
      setSizedFromParent(0);
    faceManager->registerWrapper< array1d<real64> >( viewKeyStruct::dOneSidedVolFlux_dFacePressureString )->
      setSizedFromParent(0);

    faceManager->registerWrapper< array1d<real64> >( viewKeyStruct::upwMobilityString )->
      setSizedFromParent(0);
    faceManager->registerWrapper< array2d<real64> >( viewKeyStruct::dUpwMobility_dPressureString )->
      setSizedFromParent(0);

    // elem auxiliary data used to access the one-sided face arrays
    applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase * const subRegion )
    {
      subRegion->registerWrapper< array1d<localIndex> >( viewKeyStruct::elemOffsetString )->
        setSizedFromParent(0);
    });   
  }
}
  
void SinglePhaseMimetic::ResizeFaceFields( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();

  // resize the face-based data 
  faceManager->getReference<array2d<localIndex> >( viewKeyStruct::faceToOneSidedFaceString ).resizeDimension<1>(2);
  faceManager->getReference<array2d<globalIndex> >( viewKeyStruct::faceToElemDofNumberString ).resizeDimension<1>(2);    

  // resize the one-sided face-based data
  ResizeOneSidedFaceFields( domain );
  
}

void SinglePhaseMimetic::ResizeOneSidedFaceFields( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();
  
  // first count the number of half-faces in the target region
  m_numOneSidedFaces = 0;
  elemManager->
  forElementSubRegionsComplete<CellElementSubRegion,
			       FaceElementSubRegion>( m_targetRegions,
                                                      [&]( localIndex const GEOSX_UNUSED_ARG( er ),
                                                           localIndex const GEOSX_UNUSED_ARG( esr ),
                                                           ElementRegionBase * const,
                                                           auto const * const subRegion )
  {
    arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();
   
    for (localIndex ei = 0; ei < subRegion->size(); ++ei)                                  
    {

      for(localIndex iface=0; iface<elemsToFaces.size(1); ++iface)
      {
        ++m_numOneSidedFaces;
      }
    }
    
  });

  // then resize the one-sided face-based arrays 
  faceManager->getReference< array1d<localIndex> >( viewKeyStruct::oneSidedFaceToFaceString ).resizeDimension<0>(m_numOneSidedFaces);

  faceManager->getReference< array1d<globalIndex> >( viewKeyStruct::neighborDofNumberString ).resizeDimension<0>(m_numOneSidedFaces);
  faceManager->getReference< array1d<globalIndex> >( viewKeyStruct::neighborRegionIdString ).resizeDimension<0>(m_numOneSidedFaces);
  faceManager->getReference< array1d<globalIndex> >( viewKeyStruct::neighborSubRegionIdString ).resizeDimension<0>(m_numOneSidedFaces);
  faceManager->getReference< array1d<globalIndex> >( viewKeyStruct::neighborElemIdString ).resizeDimension<0>(m_numOneSidedFaces);

  faceManager->getReference< array1d<real64> >( viewKeyStruct::oneSidedVolFluxString ).resizeDimension<0>(m_numOneSidedFaces);
  faceManager->getReference< array1d<real64> >( viewKeyStruct::dOneSidedVolFlux_dPressureString ).resizeDimension<0>(m_numOneSidedFaces);
  faceManager->getReference< array1d<real64> >( viewKeyStruct::dOneSidedVolFlux_dFacePressureString ).resizeDimension<0>(m_numOneSidedFaces);

  faceManager->getReference< array1d<real64> >( viewKeyStruct::upwMobilityString ).resizeDimension<0>(m_numOneSidedFaces);
  faceManager->getReference< array2d<real64> >( viewKeyStruct::dUpwMobility_dPressureString ).resizeDimension<0>(m_numOneSidedFaces);
  faceManager->getReference< array2d<real64> >( viewKeyStruct::dUpwMobility_dPressureString ).resizeDimension<1>(2);

  // finaly resize the element-based array used to access the one-sided face arrays  
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    subRegion->getReference< array1d<localIndex> >( viewKeyStruct::elemOffsetString ).resizeDimension<0>(subRegion->size()+1);
  });
}

void SinglePhaseMimetic::ResetFaceMaps( DomainPartition * const domain )
{
  MeshLevel * const mesh                = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = mesh->getFaceManager();

  // get the reverse maps from face to one-sided face data
  arrayView1d<localIndex> const & numAdjacentElems =
    faceManager->getReference<array1d<localIndex> >( viewKeyStruct::numAdjacentElementsString );
  arrayView2d<localIndex> const & faceToOneSidedFace =
    faceManager->getReference<array2d<localIndex> >( viewKeyStruct::faceToOneSidedFaceString );
  arrayView2d<globalIndex> const & faceToElemDofNumber = 
    faceManager->getReference<array2d<globalIndex> >( viewKeyStruct::faceToElemDofNumberString );

  for (localIndex iface = 0; iface < faceManager->size(); ++iface)
  {
    numAdjacentElems[iface]       = 0;
    faceToOneSidedFace[iface][0]  = -1;
    faceToOneSidedFace[iface][1]  = -1;
    faceToElemDofNumber[iface][0] = -1;
    faceToElemDofNumber[iface][1] = -1;
  }
}

void SinglePhaseMimetic::ConstructFaceMaps( DomainPartition * const domain )
{
  ResetFaceMaps( domain );
  
  MeshLevel * const mesh                   = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager const * const faceManager    = mesh->getFaceManager();
  
  DofManager const & dofManager = this->getDofManager();

  // TODO: is there a way to constify this?
  array2d<localIndex> const & elemRegionList    = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList          = faceManager->elementList(); 
 
  // get the map from one-sided face to face
  arrayView1d<localIndex> const & oneSidedFaceToFace =
    faceManager->getReference<array1d<localIndex>>( viewKeyStruct::oneSidedFaceToFaceString );  

  // get the reverse maps from face to one-sided face data
  arrayView1d<localIndex> const & numAdjacentElems =
    faceManager->getReference<array1d<localIndex> >( viewKeyStruct::numAdjacentElementsString );
  arrayView2d<localIndex> const & faceToOneSidedFace =
    faceManager->getReference<array2d<localIndex> >( viewKeyStruct::faceToOneSidedFaceString );
  arrayView2d<globalIndex> const & faceToElemDofNumber = 
    faceManager->getReference<array2d<globalIndex> >( viewKeyStruct::faceToElemDofNumberString );
  
  // get the "non-local" maps that relate each one-sided face to the neighbor data
  arrayView1d<globalIndex> const & neighborDofNumber =
    faceManager->getReference<array1d<globalIndex>>( viewKeyStruct::neighborDofNumberString );
  arrayView1d<globalIndex> const & neighborRegionId =
    faceManager->getReference<array1d<globalIndex>>( viewKeyStruct::neighborRegionIdString );
  arrayView1d<globalIndex> const & neighborSubRegionId =
    faceManager->getReference<array1d<globalIndex>>( viewKeyStruct::neighborSubRegionIdString );
  arrayView1d<globalIndex> const & neighborElemId =
    faceManager->getReference<array1d<globalIndex>>( viewKeyStruct::neighborElemIdString );

  localIndex oneSidedFaceCounter = 0;
  
  // in this loop we collect the one-sided faces corresponding to an element
  // TODO: figure out why it does not compile with <CellElementSubRegion,FaceElementSubRegion>
  elemManager->
    forElementSubRegionsComplete<CellElementSubRegion>( m_targetRegions,
                                                      [&]( localIndex const er,
                                                           localIndex const esr,
                                                           ElementRegionBase * const,
                                                           CellElementSubRegion * const subRegion )
  {
    arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();

    // get the offsets to access the local one-sided faces of an element
    arrayView1d<localIndex> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );

    // cell-centered dof numbers   
    string const elemCenteredDofKey = dofManager.getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & elemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( elemCenteredDofKey );  


    
    // for each element, visit the one-sided faces
    for (localIndex ei = 0; ei < subRegion->size(); ++ei)
    {
      // save the position of the first one-sided face of the element
      elemOffset[ei] = oneSidedFaceCounter;

      // for each one-sided face, save the data of the neighbor element
      for (localIndex iosface=0; iosface<elemsToFaces.size(1); ++iosface)
      {
	
        // get the ID of the face
        localIndex const currentFaceId = elemsToFaces[ei][iosface];


        // 1) One-sided face-based maps (used to assemble the mass one-sided fluxes)

	
        // save the index of the face corresponding to this one-sided face
        // the array "oneSidedFaceToFace" is redundant, I may get rid of it later
        oneSidedFaceToFace[oneSidedFaceCounter] = currentFaceId;
 	
        // the face has at most two adjacent elements
        // one of these two elements is the current element indexed by er, esr, ei
        // but here, we are interested in saving the indices of the other element
        // this other element is "the neighbor" for this one-sided face 
        for ( localIndex k=0 ; k<elemRegionList.size(1) ; ++k )
        {
          localIndex const erNeighbor  = elemRegionList[currentFaceId][k];
          localIndex const esrNeighbor = elemSubRegionList[currentFaceId][k];
          localIndex const eiNeighbor  = elemList[currentFaceId][k];
	  
          // if we have found the neighbor or a boundary face
          if ( erNeighbor != er || esrNeighbor != esr || eiNeighbor != ei ) 
          {

            // save the neighbor info
            neighborRegionId[oneSidedFaceCounter]    = erNeighbor;
            neighborSubRegionId[oneSidedFaceCounter] = esrNeighbor;
            neighborElemId[oneSidedFaceCounter]      = eiNeighbor;
	    
            // if not on boundary, save the dof number of the element
            if ( erNeighbor != -1 && esrNeighbor != -1 && eiNeighbor != -1 )
            {
              ElementRegionBase const * const neighborRegion =
                Group::group_cast<ElementRegionBase const *>(mesh->getElemManager()->GetRegion(erNeighbor));
              ElementSubRegionBase const * const neighborSubRegion =
                Group::group_cast<ElementSubRegionBase const *>(neighborRegion->GetSubRegion(esrNeighbor));

              arrayView1d<globalIndex const> const & dofNumber =
                neighborSubRegion->getReference<array1d<globalIndex>>( elemCenteredDofKey );  
              
              neighborDofNumber[oneSidedFaceCounter] = dofNumber[eiNeighbor];
            }
            // if on boundary
            else
            {
              // for the case of a boundary face, we assign the "local" dof number
              // to neighborDofNumber. Then, in the assembly, we add a zero derivative
              // in this slot of the Jacobian matrix
              neighborDofNumber[oneSidedFaceCounter] = elemDofNumber[ei];      
            }
          }

        } // end of loop on elems adjacent to the face


	
	// 2) Reverse maps: Face-based maps used to assemble the constraints 

	
        // now save index of the one-sided face corresponding to this face
        if (faceToOneSidedFace[currentFaceId][0] == -1)
	{ 
	  faceToOneSidedFace[currentFaceId][0]  = oneSidedFaceCounter;
	  faceToElemDofNumber[currentFaceId][0] = elemDofNumber[ei];
	}
        else
	{
	  faceToOneSidedFace[currentFaceId][1]  = oneSidedFaceCounter;
	  faceToElemDofNumber[currentFaceId][1] = elemDofNumber[ei];
	}
	++numAdjacentElems[currentFaceId];
	
        ++oneSidedFaceCounter;
	
      } // end of loop on one-sided faces of the elem
      
      
    } // end of loop on elements of the subregion

    elemOffset[subRegion->size()] = oneSidedFaceCounter;


    
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
  // reset the views into cell-centered fields
  ResetViews( domain );
  
  // computes the number of one-sided faces and resizes the one-sided face fields accordingly
  // if no topological change in the mesh, this does not have to be recomputed 
  ResizeFaceFields( domain );
  
  // build all the maps that are necessary to iterate over one-sided faces
  // if no topological change in the mesh, this does not have to be recomputed 
  ConstructFaceMaps( domain );

  // setup the cell-centered fields
  SinglePhaseFlowBase::ImplicitStepSetup( time_n, dt, domain, dofManager, matrix, rhs, solution );

  // setup the face fields
  MeshLevel * const meshLevel     = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = meshLevel->getFaceManager();

  arrayView1d<real64> & dFacePres = faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

  localIndex const numFaces = faceManager->size();
  for( localIndex iface = 0 ; iface < numFaces ; ++iface )
  {
    dFacePres[iface] = 0;
  }

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
  
  arrayView1d<real64> const & facePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::facePressureString);
  arrayView1d<real64 const> const & dFacePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

  localIndex const numFaces = faceManager->size();
  for( localIndex iface = 0 ; iface < numFaces ; ++iface )
  {
    facePres[iface] += dFacePres[iface];
  }
}

void SinglePhaseMimetic::SetupDofs( DomainPartition const * const domain,
                                    DofManager & dofManager ) const
{
  // setup the connectivity of cell-centered fields
  SinglePhaseFlowBase::SetupDofs( domain, dofManager );
  
  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::facePressureString,
                       DofManager::Location::Face,
                       DofManager::Connectivity::Elem,
                       m_targetRegions );
  
  // TODO: figure out what to do with target regions for faces
  dofManager.addCoupling( viewKeyStruct::facePressureString,
                          viewKeyStruct::pressureString,
                          DofManager::Connectivity::Elem,
   			  true );


}


void SinglePhaseMimetic::AssembleFluxTerms( real64 const time_n,
                                            real64 const dt,
                                            DomainPartition const * const domain,
                                            DofManager const * const dofManager,
                                            ParallelMatrix * const matrix,
                                            ParallelVector * const rhs )
{
  // TODO: figure out if it is better to:
  //  1. Compute all the vol fluxes for all the elements,
  //     then upwind all the mobs at all the faces,
  //	 then compute all the mass fluxes for all the elements
  //  or 
  //  2. In one element, compute the vol fluxes, upwind the mobs, then compute the mass fluxes
  //     then move to the next element
  // Currently I use 1.

  
  // for each one-sided face, compute the corresponding volumetric flux
  ComputeOneSidedVolFluxes( domain );

  // we know the local flow direction => we can upwind the transport coefficients
  UpdateUpwindedTransportCoefficients( domain );
  
  // use the computed one sided vol fluxes to obtain the upwinded mass fluxes 
  AssembleUpwindedOneSidedMassFluxes( time_n, dt, domain, dofManager, matrix, rhs );

  // assemble the constraints stating the onesided fluxes should be equal on both sides
  AssembleConstraints( time_n, dt, domain, dofManager, matrix, rhs );
}


void SinglePhaseMimetic::ComputeOneSidedVolFluxes( DomainPartition const * const domain ) 
{
  GEOSX_MARK_FUNCTION;

  real64 const lengthTolerance = domain->getMeshBody( 0 )->getGlobalLengthScale() * this->m_areaRelTol;
  
  MeshLevel const * const mesh          = domain->getMeshBody(0)->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();   

  // node data
  
  arrayView1d<R1Tensor const> const & X = nodeManager->referencePosition();
  
  // face data

  // get the face-to-nodes connectivity
  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager->nodeList();

  // get the map from one-sided face to face
  arrayView1d<localIndex const> const & oneSidedFaceToFace =
    faceManager->getReference<array1d<localIndex>>( viewKeyStruct::oneSidedFaceToFaceString );  
  
  // get the face-centered depth
  arrayView1d<real64> const & faceGravDepth =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::gravityDepthString);

  // get the face-centered pressures
  arrayView1d<real64 const> const & facePres =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::facePressureString );
  arrayView1d<real64 const> const & dFacePres =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::deltaFacePressureString );

  // get the one-sided volumetric fluxes
  arrayView1d<real64> const & oneSidedVolFlux =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::oneSidedVolFluxString );
  arrayView1d<real64> const & dOneSidedVolFlux_dp =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dPressureString );
  arrayView1d<real64> const & dOneSidedVolFlux_dfp =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dFacePressureString );

  // elem data
  
  // get the cell-centered props
  FluxKernel::ElementView < arrayView1d<real64 const> > const & elemCenteredGravDepth = m_gravDepth.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & elemCenteredDens      = m_density.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dElemCenteredDens_dp  = m_dDens_dPres.toViewConst();
  
  // max number of faces allowed in an element 
  localIndex constexpr maxNumFacesInElem = SinglePhaseMimetic::MAX_NUM_FACES_IN_ELEM;

  // compute the one-sided volumetric fluxes element by element

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    // get the offsets to access the data retrieved below  
    arrayView1d<localIndex const> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );

    // get the element data needed for transmissibility computation
    arrayView1d<R1Tensor const> const & elemCenter =
     subRegion->getReference<array1d<R1Tensor>>( CellBlock::viewKeyStruct::elementCenterString ); 
    arrayView1d<R1Tensor const> const & permeability =
     subRegion->getReference<array1d<R1Tensor>>( viewKeyStruct::permeabilityString ); 
   
    // get the cell-centered pressures
    arrayView1d<real64 const> const & elemCenteredPres =
      subRegion->getReference< array1d<real64> >( viewKeyStruct::pressureString );
    arrayView1d<real64 const> const & dElemCenteredPres =
      subRegion->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );   

    // assemble the fluxes element by element
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      localIndex const eOffset        = elemOffset[ei];
      localIndex const numFacesInElem = elemOffset[ei+1] - elemOffset[ei];
      
      stackArray2d<real64, maxNumFacesInElem*maxNumFacesInElem> oneSidedTrans( numFacesInElem,
									       numFacesInElem );
      
      // we currently recompute the transmissibilities at each iteration
      RecomputeOneSidedTransmissibilities( X, // node data
 			                   faceToNodes, // face data
					   oneSidedFaceToFace, // element data
					   elemCenter[ei],
					   permeability[ei],
					   eOffset,
					   numFacesInElem,
					   lengthTolerance,
					   oneSidedTrans) ;
      
      // for each element, loop over the local (one-sided) faces
      for (localIndex iosface = 0; iosface < numFacesInElem; ++iosface)
      {
        localIndex const ifOffset = eOffset + iosface;
	
        oneSidedVolFlux[ifOffset]      = 0;
        dOneSidedVolFlux_dp[ifOffset]  = 0;
        dOneSidedVolFlux_dfp[ifOffset] = 0;
         
        // now in the following nested loop,
        // we compute the contribution of face josface to the one sided flux at face iosface
        for (localIndex josface = 0; josface < numFacesInElem; ++josface)
        {
          localIndex const jfOffset = eOffset + josface;
          
          // 1) compute the potential diff between the cell center and the face pressure
          real64 const ccPres = elemCenteredPres[ei] + dElemCenteredPres[ei];
          real64 const fPres  = facePres[oneSidedFaceToFace[jfOffset]] + dFacePres[oneSidedFaceToFace[jfOffset]];

          real64 const ccGravDepth = elemCenteredGravDepth[er][esr][ei];
	  real64 const fGravDepth  = faceGravDepth[oneSidedFaceToFace[jfOffset]];

          real64 const ccDens     = elemCenteredDens[er][esr][m_fluidIndex][ei][0];
	  real64 const dCcDens_dp = dElemCenteredDens_dp[er][esr][m_fluidIndex][ei][0]; 
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
          oneSidedVolFlux[ifOffset]      += oneSidedTrans[iosface][josface] * potDif;
          dOneSidedVolFlux_dp[ifOffset]  += oneSidedTrans[iosface][josface] * dPotDif_dp;
          dOneSidedVolFlux_dfp[ifOffset] += oneSidedTrans[iosface][josface] * dPotDif_dfp;
        }
	
        // TODO: decide if we want to upwind here instead of in a separate function
      }  
    });
  });
}


void SinglePhaseMimetic::UpdateUpwindedTransportCoefficients( DomainPartition const * const domain )
{
  GEOSX_MARK_FUNCTION;
  
  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = mesh->getFaceManager();
  
  // face data
  
  // get the one-sided volumetric fluxes
  arrayView1d<real64> const & oneSidedVolFlux =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::oneSidedVolFluxString );

  // get the upwinded mobilities
  arrayView1d<real64> const & upwMobility =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::upwMobilityString );
  arrayView2d<real64> const & dUpwMobility_dp =
    faceManager->getReference<array2d<real64>>( viewKeyStruct::dUpwMobility_dPressureString );
  
  // get the indices of the neighbor
  arrayView1d<globalIndex const> const & neighborRegionId =
    faceManager->getReference<array1d<globalIndex>>( viewKeyStruct::neighborRegionIdString );
  arrayView1d<globalIndex const> const & neighborSubRegionId =
    faceManager->getReference<array1d<globalIndex>>( viewKeyStruct::neighborSubRegionIdString );
  arrayView1d<globalIndex const> const & neighborElemId =
    faceManager->getReference<array1d<globalIndex>>( viewKeyStruct::neighborElemIdString );
  
  // elem data
  
  // get the elem-centered mobilities
  FluxKernel::ElementView < arrayView1d<real64 const> > const & elemCenteredMobility     = m_mobility.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & dElemCenteredMobility_dp = m_dMobility_dPres.toViewConst();
  
  // compute the one-sided volumetric fluxes element by element
  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {

    // get the offsets to access the data retrieved below  
    arrayView1d<localIndex const> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );

    // assemble the fluxes element by element
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      localIndex const eOffset        = elemOffset[ei];
      localIndex const numFacesInElem = elemOffset[ei+1] - elemOffset[ei];

      // for each element, loop over the local (one-sided) faces
      for (localIndex iface = 0; iface < numFacesInElem; ++iface)
      {
        localIndex const fOffset = eOffset + iface;

        bool const isBoundaryFace = (neighborElemId[fOffset] < 0);
        
        // if the local element is upwind
        if (oneSidedVolFlux[fOffset] > 0 || isBoundaryFace)
        {
          upwMobility[fOffset] = elemCenteredMobility[er][esr][ei];
          dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::LOCAL] = dElemCenteredMobility_dp[er][esr][ei];
          dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::NEIGHBOR] = 0;
        }
        // else the neighbor is upwind
        else
        {
          localIndex const erNeighbor  = neighborRegionId[fOffset];
          localIndex const esrNeighbor = neighborSubRegionId[fOffset];
          localIndex const eiNeighbor  = neighborElemId[fOffset];  
          
          upwMobility[fOffset] = elemCenteredMobility[erNeighbor][esrNeighbor][eiNeighbor];
          dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::LOCAL] = 0;
          dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::NEIGHBOR] = dElemCenteredMobility_dp[erNeighbor][esrNeighbor][eiNeighbor];
        }
      }
    });
  });
}
                                                          

void SinglePhaseMimetic::AssembleUpwindedOneSidedMassFluxes( real64 const GEOSX_UNUSED_ARG( time_n ),
                                                             real64 const dt,
                                                             DomainPartition const * const domain,
                                                             DofManager const * const dofManager,
                                                             ParallelMatrix * const matrix,
                                                             ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh          = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = mesh->getFaceManager();

  // get the map from one-sided face to face
  arrayView1d<localIndex const> const & oneSidedFaceToFace =
    faceManager->getReference<array1d<localIndex>>( viewKeyStruct::oneSidedFaceToFaceString );  
  
  // get the face-based DOF numbers for the assembly
  string const faceDofKey = dofManager->getKey( viewKeyStruct::facePressureString );
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference<array1d<globalIndex>>( faceDofKey );  
  arrayView1d<globalIndex const> const & neighborDofNumber =
    faceManager->getReference<array1d<globalIndex>>( viewKeyStruct::neighborDofNumberString );

  // get the upwinded mobilities
  arrayView1d<real64 const> const & upwMobility =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::upwMobilityString );
  arrayView2d<real64 const> const & dUpwMobility_dp =
    faceManager->getReference<array2d<real64>>( viewKeyStruct::dUpwMobility_dPressureString );

  // get the one-sided volumetric fluxes
  arrayView1d<real64 const> const & oneSidedVolFlux =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::oneSidedVolFluxString );
  arrayView1d<real64 const> const & dOneSidedVolFlux_dp =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dPressureString );
  arrayView1d<real64 const> const & dOneSidedVolFlux_dfp =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dFacePressureString );
  
  // max number of faces allowed in an element 
  localIndex constexpr maxNumFacesInElem = SinglePhaseMimetic::MAX_NUM_FACES_IN_ELEM;

  // compute the one-sided volumetric fluxes element by element
  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    // get the offsets to access the data retrieved below  
    arrayView1d<localIndex const> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );

    // get the cell-centered DOF numbers and ghost rank for the assembly
    string const elemCenteredDofKey = dofManager->getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & elemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( elemCenteredDofKey );  
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];   

    // assemble the fluxes element by element
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {

      if (elemGhostRank[ei] < 0)
      {
 
        localIndex const eOffset        = elemOffset[ei];
        localIndex const numFacesInElem = elemOffset[ei+1] - elemOffset[ei];

        // fluxes
        real64 sumLocalMassFluxes           = 0;
        real64 dSumLocalMassFluxes_dp_local = 0;
        stackArray1d<real64, maxNumFacesInElem> dSumLocalMassFluxes_dp_neighbor( numFacesInElem );
        stackArray1d<real64, maxNumFacesInElem> dSumLocalMassFluxes_dfp( numFacesInElem );
        dSumLocalMassFluxes_dp_neighbor = 0;
	dSumLocalMassFluxes_dfp = 0;
	  
        // dof numbers
        globalIndex const eqnRowIndex          = elemDofNumber[ei];
        globalIndex const dofColIndexPresLocal = elemDofNumber[ei];
        stackArray1d<globalIndex, maxNumFacesInElem > dofColIndicesPresNeighbors( numFacesInElem );
        stackArray1d<globalIndex, maxNumFacesInElem > dofColIndicesFacePres( numFacesInElem );

        // for each element, loop over the local (one-sided) faces
        for (localIndex iface = 0; iface < numFacesInElem; ++iface)
        {
          localIndex const fOffset = eOffset + iface;

          // compute the mass flux at the one-sided face plus its derivatives
          // add the newly computed flux to the sum 
          IncrementLocalMassFluxSum( dt,
                                     upwMobility[fOffset],
                                     dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::LOCAL],
                                     dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::NEIGHBOR],
                                     oneSidedVolFlux[fOffset],
                                     dOneSidedVolFlux_dp[fOffset],
                                     dOneSidedVolFlux_dfp[fOffset],
                                     sumLocalMassFluxes,
                                     dSumLocalMassFluxes_dp_local,
                                     dSumLocalMassFluxes_dp_neighbor[iface],
                                     dSumLocalMassFluxes_dfp[iface] );

          // collect the relevant dof numbers
          dofColIndicesPresNeighbors[iface] = neighborDofNumber[fOffset];
          dofColIndicesFacePres[iface] = faceDofNumber[oneSidedFaceToFace[fOffset]];
  
        }

        // we are ready to assemble the local flux and its derivatives

        // residual
        rhs->add( &eqnRowIndex,
                  &sumLocalMassFluxes,
                  1 );

		  
        // jacobian -- derivative wrt local cell centered pressure term
        matrix->add( &eqnRowIndex,
                     &dofColIndexPresLocal,
                     &dSumLocalMassFluxes_dp_local,
                     1,
                     1 );

        // jacobian -- derivatives wrt neighbor cell centered pressure terms
        matrix->add( &eqnRowIndex,
                     dofColIndicesPresNeighbors.data(),
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
    });
  });
}

                                 
void SinglePhaseMimetic::AssembleConstraints( real64 const GEOSX_UNUSED_ARG( time_n ),
                                              real64 const GEOSX_UNUSED_ARG( dt ),
                                              DomainPartition const * const domain,
                                              DofManager const * const dofManager,
                                              ParallelMatrix * const matrix,
                                              ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh          = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = mesh->getFaceManager();

  // get the face-based DOF numbers and ghost ranks for the assembly
  string const faceDofKey = dofManager->getKey( viewKeyStruct::facePressureString );
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference<array1d<globalIndex>>( faceDofKey );

  arrayView1d<integer const> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  
  // get the number of adjacent elements of a given face (to handle domain boundaries)
  arrayView1d<localIndex const> const & numAdjacentElements =
    faceManager->getReference<array1d<localIndex>>( viewKeyStruct::numAdjacentElementsString );

  // get the reverse maps (face-to-sided face)
  arrayView2d<localIndex const> const & faceToOneSidedFace =
    faceManager->getReference<array2d<localIndex>>( viewKeyStruct::faceToOneSidedFaceString );
  arrayView2d<globalIndex const> const & faceToElemDofNumber = 
    faceManager->getReference<array2d<globalIndex>>( viewKeyStruct::faceToElemDofNumberString );

  
  // get the one-sided volumetric fluxes
  arrayView1d<real64 const> const & oneSidedVolFlux =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::oneSidedVolFluxString );
  arrayView1d<real64 const> const & dOneSidedVolFlux_dp =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dPressureString );
  arrayView1d<real64 const> const & dOneSidedVolFlux_dfp =
    faceManager->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dFacePressureString );

  
  // assemble the flux continuity constraints face by face
  forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  {
    // number of elements adjacent to this face
    // == 0 if the face is not in a target region
    localIndex const numAdjElems = numAdjacentElements[iface]; 
    
    if (faceGhostRank[iface] < 0 && numAdjElems > 0 )
    {
      
      // dof numbers
      globalIndex const eqnRowIndex         = faceDofNumber[iface];
      globalIndex const dofColIndexFacePres = faceDofNumber[iface];
      stackArray1d<globalIndex, 2> const dofColIndicesElemCenteredPres( numAdjElems );

      // sum of the fluxes for this face 
      real64 sumVolFluxes      = 0;
      real64 dSumVolFluxes_dfp = 0;
      stackArray1d<real64, 2> dSumVolFluxes_dp( numAdjElems );
      dSumVolFluxes_dp = 0;

      // collect the fluxes coming from the adjacent elements
      for (localIndex iAdjElem = 0; iAdjElem < numAdjElems; ++iAdjElem)
      {
	// use the reverse maps to access the one-sided face info
        sumVolFluxes               += oneSidedVolFlux[faceToOneSidedFace[iface][iAdjElem]];
        dSumVolFluxes_dfp          += dOneSidedVolFlux_dfp[faceToOneSidedFace[iface][iAdjElem]];
        dSumVolFluxes_dp[iAdjElem]  = dOneSidedVolFlux_dp[faceToOneSidedFace[iface][iAdjElem]];

	dofColIndicesElemCenteredPres[iAdjElem] = faceToElemDofNumber[iface][iAdjElem];
      }

      // we are ready to assemble the sum of fluxes and its derivatives
      
      // residual
      rhs->add( &eqnRowIndex,
                &sumVolFluxes, 
                1 );
      
      // jacobian -- derivative wrt face centered pressure term
      matrix->add( &eqnRowIndex,
                   &dofColIndexFacePres,
                   &dSumVolFluxes_dfp,
                   1,
                   1 );

      // jacobian -- derivatives wrt element centered pressure terms
      matrix->add( &eqnRowIndex,
                   dofColIndicesElemCenteredPres.data(),
                   dSumVolFluxes_dp.data(),
                   1,
                   numAdjElems );

    }
  });  
}

void SinglePhaseMimetic::IncrementLocalMassFluxSum( real64 const & dt,
                                                    real64 const & upwMobility,
                                                    real64 const & dUpwMobility_dp,
                                                    real64 const & dUpwMobility_dp_neighbor,
                                                    real64 const & oneSidedVolFlux,
                                                    real64 const & dOneSidedVolFlux_dp,
                                                    real64 const & dOneSidedVolFlux_dfp,
                                                    real64       & sumOneSidedMassFluxes,
                                                    real64       & dSumOneSidedMassFluxes_dp,
                                                    real64       & dSumOneSidedMassFluxes_dp_neighbor,
                                                    real64       & dSumOneSidedMassFluxes_dfp ) const
{ 
  sumOneSidedMassFluxes             += dt * upwMobility              * oneSidedVolFlux;
  dSumOneSidedMassFluxes_dp         += dt * ( dUpwMobility_dp        * oneSidedVolFlux 
                                            + upwMobility            * dOneSidedVolFlux_dp );
  dSumOneSidedMassFluxes_dp_neighbor = dt * dUpwMobility_dp_neighbor * oneSidedVolFlux;
  dSumOneSidedMassFluxes_dfp         = dt * upwMobility              * dOneSidedVolFlux_dfp;
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
  // 1. Compute the residual for the mass conservation equations
  real64 const massConservationResNorm = SinglePhaseFlowBase::CalculateResidualNorm( domain, 
                                                                                     dofManager,
                                                                                     rhs );

  // 2. Compute the residual for the face-based constraints
  MeshLevel const * const mesh          = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = mesh->getFaceManager();

  real64 const * localResidual = rhs.extractLocalVector();
  string const faceDofKey      = dofManager.getKey( viewKeyStruct::facePressureString );

  arrayView1d<integer const> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d<globalIndex const> const & dofNumber =
    faceManager->getReference< array1d<globalIndex> >( faceDofKey );

  localIndex const numFaces = faceManager->size();
  real64 localResidualNorm  = 0.0;
  for( localIndex iface = 0 ; iface < numFaces ; ++iface )
  {
    if (faceGhostRank[iface] < 0)
    {
      localIndex const lid    = rhs.getLocalRowID( dofNumber[iface] );
      real64 const normalizer = 1; // TODO: compute the normalizer here
      real64 const val        = localResidual[lid] / normalizer;
      localResidualNorm       += val * val;
    }
  }

  real64 globalResidualNorm;
  MpiWrapper::allReduce(&localResidualNorm, &globalResidualNorm, 1, MPI_SUM, MPI_COMM_GEOSX);

  // 3. Combine the two norms
  return sqrt( massConservationResNorm*massConservationResNorm + globalResidualNorm );
}

void SinglePhaseMimetic::ApplySystemSolution( DofManager const & dofManager,
                                              ParallelVector const & solution,
                                              real64 const scalingFactor,
                                              DomainPartition * const domain )
{
  // 1. apply the cell-centered update
  SinglePhaseFlowBase::ApplySystemSolution( dofManager, solution, scalingFactor, domain );

  // 2. apply the face-based update
  MeshLevel * const mesh          = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();

  dofManager.addVectorToField( solution,
                               viewKeyStruct::facePressureString,
                               scalingFactor,
                               faceManager,
                               viewKeyStruct::deltaFacePressureString );

  std::map<string, string_array> fieldNames;
  fieldNames["faces"].push_back( viewKeyStruct::deltaFacePressureString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );
}


void SinglePhaseMimetic::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // 1. Reset the cell-centered fields
  SinglePhaseFlowBase::ResetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  MeshLevel * const mesh          = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();
  
  arrayView1d<real64> & dFacePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

  localIndex const numFaces = faceManager->size();
  for( localIndex iface = 0 ; iface < numFaces ; ++iface )
  {
    dFacePres[iface] = 0;
  }
}

void SinglePhaseMimetic::ResetViews( DomainPartition * const domain )
{
  SinglePhaseFlowBase::ResetViews( domain );

  // so far we do not have face-based views to reset here
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
void SinglePhaseMimetic::RecomputeOneSidedTransmissibilities( arrayView1d<R1Tensor const> const & X, // node data
							      ArrayOfArraysView<localIndex const> const & faceToNodes, // face data
							      arrayView1d<localIndex const> const & oneSidedFaceToFace, // element data
							      R1Tensor const & elemCenter, 
 					                      R1Tensor const & permeability,
							      real64   const & elemOffset,
							      real64   const & numFacesInElem,
							      real64   const & lengthTolerance,
 					                      stackArray2d<real64, MAX_NUM_FACES_IN_ELEM
                                                                                  *MAX_NUM_FACES_IN_ELEM> & oneSidedTrans ) const 
{
  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor permeabilityTensor;

  real64 const areaTolerance   = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?
  
  // we are ready to compute the transmissibility matrix
  for (localIndex iosface = 0; iosface < numFacesInElem; ++iosface)
  {
    for (localIndex josface = 0; josface < numFacesInElem; ++josface)
    {
      localIndex const ifOffset = elemOffset + iosface;
	
      // for now, TPFA trans
      if (iosface == josface)
      {
	localIndex const currentFaceId = oneSidedFaceToFace[ifOffset];
	
        // 1) compute the face geometry data: center, normal, vector from cell center to face center
        real64 const faceArea =
	  computationalGeometry::Centroid_3DPolygon( faceToNodes[currentFaceId],
	    	 			             faceToNodes.sizeOfArray( currentFaceId ),
						     X,
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
        makeFullTensor(permeability, permeabilityTensor);

        faceConormal.AijBj(permeabilityTensor, faceNormal);

	// 3) compute the one-sided face transmissibility
	oneSidedTrans[iosface][josface]  = Dot(cellToFaceVec, faceConormal);
        oneSidedTrans[iosface][josface] *= faceArea / c2fDistance;
        oneSidedTrans[iosface][josface]  = std::max( oneSidedTrans[iosface][josface],
				  		     weightTolerance );
      }
      else
      {
        oneSidedTrans[iosface][josface] = 0;
      }
    }
  }
}

 
REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseMimetic, std::string const &, Group * const )
} /* namespace geosx */
