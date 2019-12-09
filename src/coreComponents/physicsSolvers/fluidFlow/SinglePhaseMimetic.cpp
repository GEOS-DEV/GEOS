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
  SinglePhaseFlowBase(name, parent)
{
  // TODO: decide what to do with m_numDofsPerCell here
}

void SinglePhaseMimetic::RegisterDataOnMesh(Group * const MeshBodies)
{

  // 1. Register the cell-centered data
  SinglePhaseFlowBase::RegisterDataOnMesh(MeshBodies);

  // 2. Register the face data 
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    FaceManager * const faceManager = meshLevel->getFaceManager();

    faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::facePressureString )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the pressures at the faces.");

    faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::deltaFacePressureString )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the accumulated pressure updates at the faces.");

    faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::faceGravDepthString )->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the depth at the face centers.");
    
    RegisterOneSidedFaceData( meshLevel );

    ResizeOneSidedFaceFields( meshLevel );
  }
}

  
void SinglePhaseMimetic::RegisterOneSidedFaceData( MeshLevel * const mesh )
{
  // in this loop we register the one-face data for each subregion
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    // this first map is a little bit redundant with subRegion->faceList() so I may get rid of it
    subRegion->registerWrapper< array1d<localIndex> >( viewKeyStruct::oneSidedFaceToFaceString )->
      setSizedFromParent(0);
    subRegion->registerWrapper< array1d<localIndex> >( viewKeyStruct::elemOffsetString )->
      setSizedFromParent(0);

    subRegion->registerWrapper< array1d<globalIndex> >( viewKeyStruct::neighborDofNumberString )->
      setSizedFromParent(0);
    subRegion->registerWrapper< array1d<globalIndex> >( viewKeyStruct::neighborRegionIdString )->
      setSizedFromParent(0);
    subRegion->registerWrapper< array1d<globalIndex> >( viewKeyStruct::neighborSubRegionIdString )->
      setSizedFromParent(0);
    subRegion->registerWrapper< array1d<globalIndex> >( viewKeyStruct::neighborElemIdString )->
      setSizedFromParent(0);
    
    subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::oneSidedVolFluxString )->
      setSizedFromParent(0);
    subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::dOneSidedVolFlux_dPressureString )->
      setSizedFromParent(0);
    subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::dOneSidedVolFlux_dFacePressureString )->
      setSizedFromParent(0);

    subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::upwMobilityString )->
      setSizedFromParent(0);
    subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::dUpwMobility_dPressureString )->
      setSizedFromParent(0);
    
  });
}

void SinglePhaseMimetic::ResizeOneSidedFaceFields( MeshLevel * const mesh )
{
  ElementRegionManager * const elemManager = mesh->getElemManager();
    
  // first count the number of half-faces
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
   
    //forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex GEOSX_UNUSED_ARG( ei ) )
    for (localIndex ei = 0; ei < subRegion->size(); ++ei)                                  
    {

      for(localIndex iface=0; iface<elemsToFaces.size(1); ++iface)
      {
        ++m_numOneSidedFaces;
      }
    }
    
  });

  // then resize the arrays
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    subRegion->getReference< array1d<localIndex> >( viewKeyStruct::oneSidedFaceToFaceString ).resizeDimension<0>(m_numOneSidedFaces);
    subRegion->getReference< array1d<localIndex> >( viewKeyStruct::elemOffsetString ).resizeDimension<0>(subRegion->size()+1);

    subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::neighborDofNumberString ).resizeDimension<0>(m_numOneSidedFaces);
    subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::neighborRegionIdString ).resizeDimension<0>(m_numOneSidedFaces);
    subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::neighborSubRegionIdString ).resizeDimension<0>(m_numOneSidedFaces);
    subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::neighborElemIdString ).resizeDimension<0>(m_numOneSidedFaces);

    subRegion->getReference< array1d<real64> >( viewKeyStruct::oneSidedVolFluxString ).resizeDimension<0>(m_numOneSidedFaces);
    subRegion->getReference< array1d<real64> >( viewKeyStruct::dOneSidedVolFlux_dPressureString ).resizeDimension<0>(m_numOneSidedFaces);
    subRegion->getReference< array1d<real64> >( viewKeyStruct::dOneSidedVolFlux_dFacePressureString ).resizeDimension<0>(m_numOneSidedFaces);

    subRegion->getReference< array1d<real64> >( viewKeyStruct::upwMobilityString ).resizeDimension<0>(m_numOneSidedFaces);
    subRegion->getReference< array2d<real64> >( viewKeyStruct::dUpwMobility_dPressureString ).resizeDimension<0>(m_numOneSidedFaces);
    subRegion->getReference< array2d<real64> >( viewKeyStruct::dUpwMobility_dPressureString ).resizeDimension<1>(2);
    
  });
}

  
void SinglePhaseMimetic::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  SinglePhaseFlowBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

}

  
void SinglePhaseMimetic::ImplicitStepSetup( real64 const & time_n,
                                            real64 const & dt,
                                            DomainPartition * const domain,
                                            DofManager & dofManager,
                                            ParallelMatrix & matrix,
                                            ParallelVector & rhs,
                                            ParallelVector & solution )
{
  // build all the maps that are necessary to iterate over one-sided faces
  // if no topological change in the mesh, this does not have to be recomputed 
  ConstructOneSidedFaceMaps( domain );

  // reset the views into cell-centered fields
  ResetViews( domain );

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
}

  
void SinglePhaseMimetic::ConstructOneSidedFaceMaps( DomainPartition * const domain )
{
  MeshLevel * const mesh                   = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager const * const faceManager    = mesh->getFaceManager();
  
  DofManager const & dofManager = this->getDofManager();
  
  array2d<localIndex> const & elemRegionList    = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList          = faceManager->elementList(); 
 
  // in this loop we collect the one-sided faces corresponding to an element
  // TODO: figure out why it does not compile with <CellElementSubRegion,FaceElementSubRegion>
  elemManager->
  forElementSubRegionsComplete<CellElementSubRegion>( m_targetRegions,
                                                      [&]( localIndex const er,
                                                           localIndex const esr,
                                                           ElementRegionBase * const,
                                                           CellElementSubRegion const * const subRegion )
  {
    arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();

    // get the map from one-sided face to face
    arrayView1d<localIndex> const & oneSidedFaceToFace =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::oneSidedFaceToFaceString );  
    
    // get the offsets to access the local one-sided faces of an element
    arrayView1d<localIndex> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );

    // cell-centered dof numbers   
    string const elemCenteredDofKey = dofManager.getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & localDofNumber =
      subRegion->getReference<array1d<globalIndex>>( elemCenteredDofKey );  
    
    // get the "non-local" maps that relate each one-sided face to the neighbor data
    arrayView1d<globalIndex> const & neighborDofNumber =
      subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::neighborDofNumberString );
    arrayView1d<globalIndex> const & neighborRegionId =
      subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::neighborRegionIdString );
    arrayView1d<globalIndex> const & neighborSubRegionId =
      subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::neighborSubRegionIdString );
    arrayView1d<globalIndex> const & neighborElemId =
      subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::neighborElemIdString );
    
    localIndex oneSidedFaceCounter = 0;
    
    // collect the one-sided faces in each element
    //forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    for (localIndex ei = 0; ei < subRegion->size(); ++ei)
    {
      // save the position of the first one-sided faces of the element
      elemOffset[ei] = oneSidedFaceCounter;
       
      for ( localIndex iface=0; iface<elemsToFaces.size(1); ++iface )
      {
        // save the index of the face corresponding to this one-sided face
        // the array "oneSidedFaceToFace" is redundant, I may get rid of it later
        oneSidedFaceToFace[oneSidedFaceCounter] = elemsToFaces[ei][iface];

        // get the ID of the face
        localIndex const currentFaceId = elemsToFaces[ei][iface];
        
        // the face has two adjacent elements
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
	    
            // if not on boundary
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
              neighborDofNumber[oneSidedFaceCounter] = localDofNumber[ei];
            }
          }
        }
        ++oneSidedFaceCounter;
      }
    }
    elemOffset[subRegion->size()] = oneSidedFaceCounter;
  });
}
  
void SinglePhaseMimetic::AssembleFluxTerms( real64 const time_n,
                                            real64 const dt,
                                            DomainPartition const * const domain,
                                            DofManager const * const dofManager,
                                            ParallelMatrix * const matrix,
                                            ParallelVector * const rhs )
{
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
  
  MeshLevel const * const mesh          = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = mesh->getFaceManager();   
  
  // get the face-centered pressures
  arrayView1d<real64 const> const & facePres =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::facePressureString );
  arrayView1d<real64 const> const & dFacePres =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::deltaFacePressureString );
  // get the face-centered depth 
  arrayView1d<real64 const> const & faceGravDepth =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::faceGravDepthString );  
  
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

    // get the map from one-sided face to face
    arrayView1d<localIndex const> const & oneSidedFaceToFace =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::oneSidedFaceToFaceString );  
    
    // get the offsets to access the data retrieved below  
    arrayView1d<localIndex const> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );
    
    // get the cell-centered pressures
    arrayView1d<real64 const> const & elemCenteredPres =
      subRegion->getReference< array1d<real64> >( viewKeyStruct::pressureString );
    arrayView1d<real64 const> const & dElemCenteredPres =
      subRegion->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );
    
    // get the one-sided volumetric fluxes
    arrayView1d<real64> const & oneSidedVolFlux =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::oneSidedVolFluxString );
    arrayView1d<real64> const & dOneSidedVolFlux_dp =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dPressureString );
    arrayView1d<real64> const & dOneSidedVolFlux_dfp =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dFacePressureString );

    // assemble the fluxes element by element
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      localIndex const eOffset        = elemOffset[ei];
      localIndex const numFacesInElem = elemOffset[ei+1] - elemOffset[ei];
      
      stackArray2d<real64, maxNumFacesInElem*maxNumFacesInElem> oneSidedTrans( numFacesInElem,
									       numFacesInElem );
      
      // we currently recompute the transmissibilities at each iteration
      RecomputeOneSidedTransmissibilities( subRegion, ei, numFacesInElem,
                                           oneSidedTrans );
      
      // for each element, loop over the local (one-sided) faces
      for (localIndex iface = 0; iface < numFacesInElem; ++iface)
      {
        localIndex const ifOffset = eOffset + iface;
	
        oneSidedVolFlux[ifOffset]      = 0;
        dOneSidedVolFlux_dp[ifOffset]  = 0;
        dOneSidedVolFlux_dfp[ifOffset] = 0;
        
        // now in the following nested loop,
        // we compute the contribution of face jface to the one sided flux at face iface
        for (localIndex jface = 0; jface < numFacesInElem; ++jface)
        {
          localIndex const jfOffset = eOffset + jface;
          
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
          oneSidedVolFlux[ifOffset]      += oneSidedTrans[iface][jface] * potDif;
          dOneSidedVolFlux_dp[ifOffset]  += oneSidedTrans[iface][jface] * dPotDif_dp;
          dOneSidedVolFlux_dfp[ifOffset] += oneSidedTrans[iface][jface] * dPotDif_dfp;
        }
	
        // TODO: decide if we want to upwind here instead of in a separate function
      }  
    });
  });
}


void SinglePhaseMimetic::RecomputeOneSidedTransmissibilities( ElementSubRegionBase const * const GEOSX_UNUSED_ARG( subRegion ),
                                                              localIndex const GEOSX_UNUSED_ARG( ei ),
                                                              localIndex const numFacesInElem,
                                                              stackArray2d<real64, SinglePhaseMimetic::MAX_NUM_FACES_IN_ELEM
                                                                                  *SinglePhaseMimetic::MAX_NUM_FACES_IN_ELEM> & oneSidedTrans )
{
  // for each element, loop over the local (one-sided) faces
  for (localIndex iface = 0; iface < numFacesInElem; ++iface)
  {
    for (localIndex jface = 0; jface < numFacesInElem; ++jface)
    {
      oneSidedTrans[iface][jface] = 0; // for now. Will work on that now 
    }
  }
}


void SinglePhaseMimetic::UpdateUpwindedTransportCoefficients( DomainPartition const * const domain )
{
  GEOSX_MARK_FUNCTION;
  
  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  
  // get the cell-centered mobilities
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

    // get the indices of the neighbor
    arrayView1d<globalIndex const> const & neighborRegionId =
      subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::neighborRegionIdString );
    arrayView1d<globalIndex const> const & neighborSubRegionId =
      subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::neighborSubRegionIdString );
    arrayView1d<globalIndex const> const & neighborElemId =
      subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::neighborElemIdString );

    // get the upwinded mobilities
    arrayView1d<real64> const & upwMobility =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::upwMobilityString );
    arrayView2d<real64> const & dUpwMobility_dp =
      subRegion->getReference<array2d<real64>>( viewKeyStruct::dUpwMobility_dPressureString );
    
    // get the one-sided volumetric fluxes
    arrayView1d<real64 const> const & oneSidedVolFlux =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::oneSidedVolFluxString );

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
 
  // get the face-based DOF numbers for the assembly
  string const faceDofKey = dofManager->getKey( viewKeyStruct::facePressureString );
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference<array1d<globalIndex>>( faceDofKey );  
  
  // max number of faces allowed in an element 
  localIndex constexpr maxNumFacesInElem = SinglePhaseMimetic::MAX_NUM_FACES_IN_ELEM;

  // compute the one-sided volumetric fluxes element by element
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase const * const subRegion )
  {
    // get the map from one-sided face to face
    arrayView1d<localIndex const> const & oneSidedFaceToFace =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::oneSidedFaceToFaceString );  
     
    // get the offsets to access the data retrieved below  
    arrayView1d<localIndex const> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );

    // get the cell-centered DOF numbers for the assembly
    string const elemCenteredDofKey = dofManager->getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & localDofNumber =
      subRegion->getReference<array1d<globalIndex>>( elemCenteredDofKey );  
    arrayView1d<globalIndex const> const & neighborDofNumber =
      subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::neighborDofNumberString );
    
    // get the upwinded mobilities
    arrayView1d<real64 const> const & upwMobility =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::upwMobilityString );
    arrayView2d<real64 const> const & dUpwMobility_dp =
      subRegion->getReference<array2d<real64>>( viewKeyStruct::dUpwMobility_dPressureString );

    // get the one-sided volumetric fluxes
    arrayView1d<real64 const> const & oneSidedVolFlux =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::oneSidedVolFluxString );
    arrayView1d<real64 const> const & dOneSidedVolFlux_dp =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dPressureString );
    arrayView1d<real64 const> const & dOneSidedVolFlux_dfp =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dFacePressureString );

    // assemble the fluxes element by element
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      localIndex const eOffset        = elemOffset[ei];
      localIndex const numFacesInElem = elemOffset[ei+1] - elemOffset[ei];

      // fluxes
      real64 sumLocalMassFluxes           = 0;
      real64 dSumLocalMassFluxes_dp_local = 0;
      stackArray1d<real64, maxNumFacesInElem> dSumLocalMassFluxes_dp_neighbor( numFacesInElem );
      stackArray1d<real64, maxNumFacesInElem> dSumLocalMassFluxes_dfp( numFacesInElem );

      // dof numbers
      globalIndex const eqnRowIndex          = localDofNumber[ei];
      globalIndex const dofColIndexPresLocal = localDofNumber[ei];
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

  // get the face-based DOF numbers for the assembly
  string const faceDofKey = dofManager->getKey( viewKeyStruct::facePressureString );
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference<array1d<globalIndex>>( faceDofKey );  

  // max number of faces allowed in an element 
  localIndex constexpr maxNumFacesInElem = SinglePhaseMimetic::MAX_NUM_FACES_IN_ELEM;
  
  // compute the one-sided volumetric fluxes element by element
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase const * const subRegion )
  {
    // get the map from one-sided face to face
    arrayView1d<localIndex const> const & oneSidedFaceToFace =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::oneSidedFaceToFaceString );  

    // get the offsets to access the data retrieved below  
    arrayView1d<localIndex const> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );
    
    // get the cell-centered DOF numbers for the assembly
    string const elemCenteredDofKey = dofManager->getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & localDofNumber =
      subRegion->getReference<array1d<globalIndex>>( elemCenteredDofKey );  
    
    // get the one-sided volumetric fluxes
    arrayView1d<real64 const> const & oneSidedVolFlux =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::oneSidedVolFluxString );
    arrayView1d<real64 const> const & dOneSidedVolFlux_dp =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dPressureString );
    arrayView1d<real64 const> const & dOneSidedVolFlux_dfp =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dFacePressureString );

    // assemble the fluxes element by element
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      localIndex const eOffset        = elemOffset[ei];
      localIndex const numFacesInElem = elemOffset[ei+1] - elemOffset[ei];

      // fluxes and derivatives
      stackArray1d<real64, maxNumFacesInElem > localVolFluxes( numFacesInElem );
      stackArray1d<real64, maxNumFacesInElem > dLocalVolFluxes_dp( numFacesInElem );
      stackArray1d<real64, maxNumFacesInElem > dLocalVolFluxes_dfp( numFacesInElem );
      
      // dof numbers
      stackArray1d<globalIndex, maxNumFacesInElem > eqnRowIndices( numFacesInElem );
      stackArray1d<globalIndex, maxNumFacesInElem > dofColIndicesFacePres( numFacesInElem );
      globalIndex const dofColIndexPresLocal = localDofNumber[ei];

      // for each element, loop over the local (one-sided) faces
      for (localIndex iface = 0; iface < numFacesInElem; ++iface)
      {
        localIndex const fOffset = eOffset + iface;

        // collect fluxes
        localVolFluxes[iface]      = oneSidedVolFlux[fOffset];
        dLocalVolFluxes_dp[iface]  = dOneSidedVolFlux_dp[fOffset];
        dLocalVolFluxes_dfp[iface] = dOneSidedVolFlux_dfp[fOffset]; 

        // collect eqn numbers and dofs
        eqnRowIndices[iface]         = faceDofNumber[oneSidedFaceToFace[fOffset]];
        dofColIndicesFacePres[iface] = faceDofNumber[oneSidedFaceToFace[fOffset]];
      }      

      // we are ready to assemble the local flux and its derivatives

      // residual
      rhs->add( eqnRowIndices.data(),
                localVolFluxes.data(),
                numFacesInElem );
      
      // jacobian -- derivative wrt local cell centered pressure term
      matrix->add( eqnRowIndices.data(),
                   &dofColIndexPresLocal,
                   dLocalVolFluxes_dp.data(),
                   numFacesInElem,
                   1 );

      // jacobian -- derivatives wrt face pressure terms
      matrix->add( eqnRowIndices.data(),
                   dofColIndicesFacePres.data(),
                   dLocalVolFluxes_dfp.data(),
                   numFacesInElem,
                   numFacesInElem );
      
    });
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

  // so we do not have face-based views to reset here
}


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseMimetic, std::string const &, Group * const )
} /* namespace geosx */
