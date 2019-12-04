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
#include "physicsSolvers/FiniteVolume/SinglePhaseFlowKernels.hpp"

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

  SinglePhaseMimetic::InitializePostInitialConditions_PreSubGroups( rootGroup );

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
  GEOSX_MARK_FUNCTION;

  MeshLevel * meshLevel                          = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();
  FaceManager const * const faceManager          = meshLevel->getFaceManager();

  // TODO: put this constant somewhere
  localIndex constexpr maxNumFaces     = 20;
  localIndex constexpr maxNumNeighbors = 20;

  elemManager->forElementSubRegionsComplete<CellElementSubRegion,
                                            FaceElementSubRegion>( this->m_targetRegions,
                                                                   [&] ( localIndex er,
                                                                         localIndex esr,
                                                                         ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                                                         auto const * const subRegion )
  {

    // get the connectivity information (elem to face)
    arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();
 
    // assemble the fluxes element by element
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {

      localIndex const numLocalFaces        = m_numLocalFaces[ei];
      localIndex const numInteriorNeighbors = m_numInteriorNeighbors[ei];

      // create local work arrays
      stackArray1d<real64, maxNumFaces>             volOneSidedVolFluxes( numLocalFaces );
      stackArray1d<real64, maxNumFaces>             dVolOneSidedVolFluxes_dp( numLocalFaces );
      stackArray2d<real64, maxNumFaces*maxNumFaces> dVolOneSidedVolFluxes_dpi( numLocalFaces, 
                                                                               numLocalFaces );
 
      stackArray1d<real64, maxNumFaces>               localFluxes( numLocalFaces );
      stackArray1d<real64, 2 * maxNumFaces>           localFluxesJacobianCellPres( numLocalFaces, 2 );
      stackArray2d<real64, maxNumFaces * maxNumFaces> localFluxesJacobianFacePres( numLocalFaces, 
                                                                                numLocalFaces );

      real64 sumLocalFluxes = 0;
      stackArray1d<real64, maxNumNeighbors + 1> sumLocalFluxesJacobianCellPres( numInteriorNeighbors + 1 );
      stackArray2d<real64, maxNumFaces>         sumLocalFluxesJacobianFacePres( numLocalFaces );

      stackArray1d<globalIndex, maxNumNeighbors+1> dofColIndicesMassConsCellPres( numInteriorNeighbors + 1 );
      stackArray1d<globalIndex, maxNumFaces>       dofColIndicesMassConsFacePres( numLocalFaces );

      stackArray1d<globalIndex, maxNumFaces>       eqnRowIndicesContinuity( numLocalFaces );
      stackArray1d<globalIndex, maxNumNeighbors+1> dofColIndicesContinuityCellPres( numInteriorNeighbors + 1 );
      stackArray1d<globalIndex, maxNumFaces>       dofColIndicesContinuityFacePres( numLocalFaces );

      // 1) For each local face, compute the one-sided potential gradient
      // TODO

      // 2) For each local face, multiply the potential gradient by upwinded coefficient
      ComputeLocalMassOneSidedFluxes( ei, numLocalFaces, numInteriorNeighbors, 
                                      oneSidedVolFluxes, 
                                      dOneSidedVolFluxes_dp, 
                                      dOneSidedVolFluxes_dpi,
                                      localFluxes, 
                                      localFluxesJacobianCellPres,
                                      localFluxesJacobianFacePres,
                                      sumLocalFluxes,
                                      sumLocalFluxesJacobianCellPres,
                                      sumLocalFluxesJacobianFacePres );
       
      // 3) Assemble the local contributions    
 

   
    });

  });

}

void SinglePhaseMimetic::AssembleFluxesForMassConservation()
{
  eqnRowIndex = m_localCellDofNumbers[ei][numInternalNeighbors];

  for (localIndex ineighbor = 0; ineighbor < numInteriorNeighbors; ++ineighbor)
  {
    dofColIndicesCellPres[ineighbor] = m_localCellDofNumbers[ei][ineighbor];    
  }

  for (localIndex iface = 0; iface < numFaces; ++iface)
  {
    dofColIndicesFacePres[iface] = m_localFaceDofNumbers[ei][iface];
  }
      
  // Add to global residual/jacobian
  rhs->add( eqnRowIndex,
            sumLocalFluxes,
            1 );

  matrix->add( eqnRowIndex,
               dofColIndicesCellPres.data(),
               localFluxJacobianCellPres.data(),
               1,
               numInternalNeighbors + 1 );

  matrix->add( eqnRowIndex,
               dofColIndicesFacePres.data(),
               localFluxJacobianFacePres.data(),
               1,
               numFaces );

}

void SinglePhaseMimetic::AssembleFluxesForConstraints()
{
  eqnRowIndicesConstraint = m_localCellDofNumbers[ei][numInternalNeighbors];

  for (localIndex ineighbor = 0; ineighbor < numInteriorNeighbors; ++ineighbor)
  {
    dofColIndicesCellPres[ineighbor] = m_localCellDofNumbers[ei][ineighbor];    
  }

  for (localIndex iface = 0; iface < numFaces; ++iface)
  {
    dofColIndicesMassConsFacePres[iface] = m_localFaceDofNumbers[ei][iface];
  }
      
  // Add to global residual/jacobian
  rhs->add( eqnRowIndexMassCons,
            sumLocalFluxes,
            1 );

  matrix->add( eqnRowIndex,
               dofColIndicesMassConsCellPres.data(),
               localFluxJacobianCellPres.data(),
               1,
               numInternalNeighbors + 1 );

  matrix->add( eqnRowIndex,
               dofColIndicesMassConsFacePres.data(),
               localFluxJacobianFacePres.data(),
               1,
               numFaces );

}


void SinglePhaseMimetic::ComputeLocalMassOneSidedFluxes( localIndex const ei,
                                                         localIndex const numLocalFaces, 
                                                         localIndex const numInteriorNeighbors,
                                                         arraySlice1d<real64 const> const & oneSidedVolFluxes, 
                                                         arraySlice1d<real64 const> const & dOneSidedVolFluxes_dp, 
                                                         arraySlice2d<real64 const> const & dOneSidedVolFluxes_dpi,
                                                         arraySlice1d<real64> const & localFluxes, 
                                                         arraySlice1d<real64> const & localFluxesJacobianCellPres, 
                                                         arraySlice2d<real64> const & localFluxesJacobianFacePres,
                                                         real64 & sumFluxes,
                                                         arraySlice1d<real64> const & sumLocalFluxesJacobianCellPres,
                                                         arraySlice2d<real64> const & sumLocalFluxesJacobianFacePres )
{
  // we need some secondary (cell-centered) variables 
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64 const>> const & mob        = m_mobility;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64 const>> const & dMob_dPres = m_dMobility_dPres;

  for (localIndex iface = 0; iface < numLocalFaces; ++iface)
  {

    // 1) upwind the mobility

    real64 upwMobility     = 0;
    real64 dUpwMobility_dp = 0;
    localIndex upwId       = -1;

    if (oneSidedVolFluxes[iface] > 0) // the local element is upwind
    {
      upwId = SinglePhaseMimetic::LocalElement;

      upwMobility     = mob[er][esr][ei];
      dUpwMobility_dp = dMob_dPres[er][esr][ei];
    }
    else // the neighbor element is upwind
    {
      upwId = SinglePhaseMimetic::NeighborElement;

      if (m_neighborElementIndex[ei][iface] == -1) 
      {
        // TODO
      }
      else 
      {
        localIndex const er_up  = m_neighborRegionIndex[ie][iface];
        localIndex const esr_up = m_neighborSubRegionIndex[ie][iface];
        localIndex const ei_up  = m_neighborElementIndex[ie][iface];

        upwMobility     = mob[er_up][esr_up][ei_up];
        dUpwMobility_dp = dMob_dPres[er_up][esr_up][ei_up];
      } 
    }

    // 2) convert the volumetric flux to a mass flux and compute derivatives

    localFluxes[iface]  = dt * upwMobility  * oneSidedVolFluxes[iface];
    sumLocalFluxes     += localFluxes[iface];
   
    // derivatives wrt cell-centered pressures
    localFluxesJacobianCellPres[iface][upwId]                             = dt * dUpwMobility_dp * oneSidedVolFluxes[iface];
    localFluxesJacobianCellPres[iface][SinglePhaseMimetic::LocalElement] += dt * upwMobility     * dOneSidedVolFluxes_dp[iface];
    if (m_localNeighborIndex[iface] >= 0)
    {
      sumLocalFluxesJacobianCellPres[m_localNeighborIndex[iface]] += localFluxesJacobianCellPres[iface][SinglePhaseMimetic::NeighborElement];
    }
    sumLocalFluxesJacobianCellPres[numInteriorNeighbors] += localFluxesJacobianCellPres[iface][SinglePhaseMimetic::LocalElement];
    
    // derivatives wrt interface pressures
    for (localIndex jface = 0; jface < numFaces; ++jface)
    {
      localFluxesJacobianFacePres[iface][jface]  = dt * upwMobility * dOneSidedVolFluxesFacePres[iface][jface];
      sumLocalFluxesJacobianFacePres[jface]     += localFluxesJacobianFacePres[iface][jface];
    }
  }

      /*
      // get the ID of the face
      localIndex const currentFaceId = elemsToFaces[ei][iface];
      // retrieve the indices of the neighbor element
      for( localIndex k=0 ; k<elemRegionList.size(1) ; ++k )
      {
        localIndex const er_k  = elemRegionList[currentFaceId][k];
        localIndex const esr_k = elemSubRegionList[currentFaceId][k];
        localIndex const ei_k  = elemIndexList[currentFaceId][k];

        if ( (er_k != -1 && esr_k !=-1   && ei_k !=-1)    // not on boundary 
                 (er_k != er && esr_k != esr && ei_k != ei) ) // not equal to the local element id
        {                
          er_up  = er_k;
          esr_up = esr_k;
          ei_up  = ei_k;
        }
        GEOS_ERROR_IF( (er_up == -1 && esr_up ==-1   && ei_up ==-1), 
                       "SinglePhaseMimetic: invalid upwind element id" );
      }
      */          
   
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
  MeshLevel const * const meshLevel     = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager const * const faceManager = meshLevel->getFaceManager();

  real64 const * localResidual = rhs.extractLocalVector();
  string const dofKey          = dofManager.getKey( viewKeyStruct::facePressureString );

  arrayView1d<integer const> const & faceGhostRank =
      faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d<globalIndex const> const & dofNumber = faceManager->getReference< array1d<globalIndex> >( dofKey );

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
  MeshLevel * const meshLevel     = domain->getMeshBody(0)->getMeshLevel(0);
  FaceManager * const faceManager = meshLevel->getFaceManager();

  dofManager.addVectorToField( solution,
                               viewKeyStruct::facePressureString,
                               scalingFactor,
                               faceManager,
                               viewKeyStruct::deltaFacePressureString );

  std::map<string, string_array> fieldNames;
  fieldNames["faces"].push_back( viewKeyStruct::deltaFacePressureString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, meshLevel, comms );
}


void SinglePhaseMimetic::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // 1. Reset the cell-centered fields
  SinglePhaseFlowBase::ResetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  MeshLevel * const meshLevel     = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = meshLevel->getFaceManager();
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
