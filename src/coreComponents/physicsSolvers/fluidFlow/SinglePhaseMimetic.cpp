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

  //   faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::deltaFacePressureString )->
  //     setPlotLevel(PlotLevel::LEVEL_0)->
  //     setRegisteringObjects(this->getName())->
  //     setDescription( "An array that holds the accumulated pressure updates  at the faces.");
  }
}

void SinglePhaseMimetic::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  SinglePhaseFlowBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  // TODO: decide what to do here
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
  MeshLevel * const meshLevel                 = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager             = meshLevel->getFaceManager();
  arrayView1d<real64>       const & facePres  = faceManager->getReference<array1d<real64>>(viewKeyStruct::facePressureString);
  arrayView1d<real64 const> const & dFacePres = faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

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
    faceManager->getReference< array1d<real64> >( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dFacePres =
    faceManager->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );  

  // max number of faces allowed in an element 
  localIndex constexpr maxNumFacesInElem = SinglePhaseMimetic::MAX_NUM_FACES_IN_ELEM;

  // compute the one-sided volumetric fluxes element by element 
  applyToSubRegions( mesh, [&] ( localIndex const GEOSX_UNUSED_ARG( er ),
				 localIndex const GEOSX_UNUSED_ARG( esr ),
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
    arrayView1d<real64 const> const & cellCenteredPres =
      subRegion->getReference< array1d<real64> >( viewKeyStruct::pressureString );
    arrayView1d<real64 const> const & dCellCenteredPres =
      subRegion->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );
    
    // get the one-sided volumetric fluxes
    arrayView1d<real64> const & oneSidedVolFlux =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::oneSidedVolFluxString );
    arrayView1d<real64> const & dOneSidedVolFlux_dp =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dPressureString );
    arrayView1d<real64> const & dOneSidedVolFlux_dpi =
      subRegion->getReference<array1d<real64>>( viewKeyStruct::dOneSidedVolFlux_dFacePressureString );

    // assemble the fluxes element by element
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      localIndex const eOffset        = elemOffset[ei];
      localIndex const numFacesInElem = elemOffset[ei+1] - elemOffset[ei];

      stackArray2d<real64, maxNumFacesInElem*maxNumFacesInElem> halfTrans( numFacesInElem, numFacesInElem );

      // we currently recompute the transmissibilities at each iteration
      RecomputeHalfTransmissibilities( subRegion, ei, numFacesInElem,
				       halfTrans );
    
      // for each element, loop over the local (one-sided) faces
      for (localIndex iface = 0; iface < numFacesInElem; ++iface)
      {
        localIndex const ifOffset = eOffset + iface;

        oneSidedVolFlux[ifOffset]      = 0;
        dOneSidedVolFlux_dp[ifOffset]  = 0;
        dOneSidedVolFlux_dpi[ifOffset] = 0;
	
        // now in the following nested loop,
	// we compute the contribution of face jface to the one sided flux at face iface
        for (localIndex jface = 0; jface < numFacesInElem; ++jface)
        {
          localIndex const jfOffset = eOffset + jface;
	  
	  // 1) compute the pressure diff between the cell center and the face pressure
	  real64 const ccPres = cellCenteredPres[ei] + dCellCenteredPres[ei];
	  real64 const fPres  = facePres[oneSidedFaceToFace[jfOffset]] + dFacePres[oneSidedFaceToFace[jfOffset]];
	
          real64 const potDif      = ccPres - fPres;
          real64 const dPotDif_dp  =  1;
  	  real64 const dPotDif_dpi = -1;

          // 2) compute the contribution of this face to the volumetric fluxes in the cell
  	  oneSidedVolFlux[ifOffset]      += halfTrans[iface][jface] * potDif;
   	  dOneSidedVolFlux_dp[ifOffset]  += halfTrans[iface][jface] * dPotDif_dp;
    	  dOneSidedVolFlux_dpi[ifOffset] += halfTrans[iface][jface] * dPotDif_dpi;
	}

	// TODO: decide if we want to upwind here instead of in a separate function
      }  
    });
  });
}


void SinglePhaseMimetic::RecomputeHalfTransmissibilities( ElementSubRegionBase const * const GEOSX_UNUSED_ARG( subRegion ),
							  localIndex const GEOSX_UNUSED_ARG( ei ),
							  localIndex const numFacesInElem,
							  stackArray2d<real64, SinglePhaseMimetic::MAX_NUM_FACES_IN_ELEM
							                      *SinglePhaseMimetic::MAX_NUM_FACES_IN_ELEM> & halfTrans )
{
  // for each element, loop over the local (one-sided) faces
  for (localIndex iface = 0; iface < numFacesInElem; ++iface)
  {
    for (localIndex jface = 0; jface < numFacesInElem; ++jface)
    {
      halfTrans[iface][jface] = 0; // for now 
    }
  }
}


void SinglePhaseMimetic::UpdateUpwindedTransportCoefficients( DomainPartition const * const domain )
{
  GEOSX_MARK_FUNCTION;
  
  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  
  // get the cell-centered mobilities
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & cellCenteredMobility     = m_mobility;
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & dCellCenteredMobility_dp = m_dMobility_dPres;
  
  // compute the one-sided volumetric fluxes element by element 
  applyToSubRegions( mesh, [&] ( localIndex const er,
				 localIndex const esr,
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
	  upwMobility[fOffset] = cellCenteredMobility[er][esr][ei];
          dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::LOCAL] = dCellCenteredMobility_dp[er][esr][ei];
          dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::NEIGHBOR] = 0;
	}
	// else the neighbor is upwind
	else
	{
          localIndex const erNeighbor  = neighborRegionId[fOffset];
	  localIndex const esrNeighbor = neighborSubRegionId[fOffset];
	  localIndex const eiNeighbor  = neighborElemId[fOffset];  
	  
	  upwMobility[fOffset] = cellCenteredMobility[erNeighbor][esrNeighbor][eiNeighbor];
          dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::LOCAL] = 0;
          dUpwMobility_dp[fOffset][SinglePhaseMimetic::CellPos::NEIGHBOR] = dCellCenteredMobility_dp[erNeighbor][esrNeighbor][eiNeighbor];
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
  applyToSubRegions( mesh, [&] ( localIndex const GEOSX_UNUSED_ARG( er ),
				 localIndex const GEOSX_UNUSED_ARG( esr ),
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    // get the map from one-sided face to face
    arrayView1d<localIndex const> const & oneSidedFaceToFace =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::oneSidedFaceToFaceString );  
     
    // get the offsets to access the data retrieved below  
    arrayView1d<localIndex const> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );

    // get the cell-centered DOF numbers for the assembly
    string const cellCenteredDofKey = dofManager->getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & localDofNumber =
      subRegion->getReference<array1d<globalIndex>>( cellCenteredDofKey );  
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
                                   dSumLocalMassFluxes_dfp[iface]);

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
  applyToSubRegions( mesh, [&] ( localIndex const GEOSX_UNUSED_ARG( er ),
				 localIndex const GEOSX_UNUSED_ARG( esr ),
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    // get the map from one-sided face to face
    arrayView1d<localIndex const> const & oneSidedFaceToFace =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::oneSidedFaceToFaceString );  

    // get the offsets to access the data retrieved below  
    arrayView1d<localIndex const> const & elemOffset =
      subRegion->getReference<array1d<localIndex>>( viewKeyStruct::elemOffsetString );
    
    // get the cell-centered DOF numbers for the assembly
    string const cellCenteredDofKey = dofManager->getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & localDofNumber =
      subRegion->getReference<array1d<globalIndex>>( cellCenteredDofKey );  
    
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
  MeshLevel const * const mesh     = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = mesh->getFaceManager();

  real64 const * localResidual = rhs.extractLocalVector();
  string const dofKey          = dofManager.getKey( viewKeyStruct::facePressureString );

  arrayView1d<integer const> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d<globalIndex const> const & dofNumber =
    faceManager->getReference< array1d<globalIndex> >( dofKey );

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
  MeshLevel * const mesh     = domain->getMeshBody(0)->getMeshLevel(0);
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
  MeshLevel * const mesh     = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = mesh->getFaceManager();
  arrayView1d<real64> & dFacePres = faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

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
