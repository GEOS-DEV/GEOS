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
 * @file SinglePhaseHybridFVM.cpp
 */

#include "SinglePhaseHybridFVM.hpp"

#include "common/TimingMacros.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace SinglePhaseHybridFVMKernels;
using namespace HybridFVMInnerProduct;

SinglePhaseHybridFVM::SinglePhaseHybridFVM( const std::string & name,
                                            Group * const parent ):
  SinglePhaseBase( name, parent ),
  m_faceDofKey( "" ),
  m_areaRelTol( 1e-8 )
{

  // one cell-centered dof per cell
  m_numDofPerCell = 1;
}


void SinglePhaseHybridFVM::RegisterDataOnMesh( Group * const MeshBodies )
{

  // 1) Register the cell-centered data
  SinglePhaseBase::RegisterDataOnMesh( MeshBodies );

  // 2) Register the face data
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    FaceManager * const faceManager = meshLevel->getFaceManager();

    // primary variables: face pressures
    faceManager->registerWrapper< array1d< real64 > >( viewKeyStruct::facePressureString )->
      setPlotLevel( PlotLevel::LEVEL_0 )->
      setRegisteringObjects( this->getName())->
      setDescription( "An array that holds the pressures at the faces." );

    faceManager->registerWrapper< array1d< real64 > >( viewKeyStruct::deltaFacePressureString )->
      setPlotLevel( PlotLevel::LEVEL_0 )->
      setRegisteringObjects( this->getName())->
      setDescription( "An array that holds the accumulated pressure updates at the faces." );

  }
}


void SinglePhaseHybridFVM::ImplicitStepSetup( real64 const & time_n,
                                              real64 const & dt,
                                              DomainPartition * const domain,
                                              DofManager & dofManager,
                                              ParallelMatrix & matrix,
                                              ParallelVector & rhs,
                                              ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  // setup the cell-centered fields
  SinglePhaseBase::ImplicitStepSetup( time_n, dt, domain, dofManager, matrix, rhs, solution );

  // setup the face fields
  MeshLevel * const meshLevel     = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager * const faceManager = meshLevel->getFaceManager();

  // get the face-based DOF numbers
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString );

  // save the face Dof key for use in two functions
  // that do not have acces to the coupled solver dofManager
  // namely ResetStateToBeginningOfStep and ImplicitStepComplete
  m_faceDofKey = faceDofKey;

  // get the accumulated pressure updates
  arrayView1d< real64 > & dFacePres =
    faceManager->getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  forAll< serialPolicy >( faceManager->size(), [=] ( localIndex iface )
  {
    dFacePres[iface] = 0;
  } );

}

void SinglePhaseHybridFVM::ImplicitStepComplete( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  // increment the cell-centered fields
  SinglePhaseBase::ImplicitStepComplete( time_n, dt, domain );

  // increment the face fields
  MeshLevel * const meshLevel     = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager * const faceManager = meshLevel->getFaceManager();

  // get the face-based DOF numbers
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager->getReference< array1d< globalIndex > >( m_faceDofKey );

  // get the face-based pressures
  arrayView1d< real64 > const & facePres =
    faceManager->getReference< array1d< real64 > >( viewKeyStruct::facePressureString );
  arrayView1d< real64 > const & dFacePres =
    faceManager->getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  forAll< serialPolicy >( faceManager->size(), [=] ( localIndex iface )
  {
    // update if face is in target region
    if( faceDofNumber[iface] >= 0 )
    {
      facePres[iface] += dFacePres[iface];
    }
    dFacePres[iface] = 0;
  } );
}

void SinglePhaseHybridFVM::SetupDofs( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                      DofManager & dofManager ) const
{

  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::pressureString,
                       DofManager::Location::Elem,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::pressureString,
                          viewKeyStruct::pressureString,
                          DofManager::Connector::Face );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::facePressureString,
                       DofManager::Location::Face,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::facePressureString,
                          viewKeyStruct::facePressureString,
                          DofManager::Connector::Elem );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::facePressureString,
                          viewKeyStruct::pressureString,
                          DofManager::Connector::Elem,
                          true );

}

void SinglePhaseHybridFVM::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const dt,
                                              DomainPartition const * const domain,
                                              DofManager const * const dofManager,
                                              ParallelMatrix * const matrix,
                                              ParallelVector * const rhs )
{
  MeshLevel const & mesh                   = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *mesh.getElemManager();
  NodeManager const & nodeManager          = *mesh.getNodeManager();
  FaceManager const & faceManager          = *mesh.getFaceManager();

  // in this function we need to make sure that we act only on the target regions
  // for that, we need the following region filter
  SortedArray< localIndex > regionFilter;
  for( string const & regionName : targetRegionNames() )
  {
    regionFilter.insert( elemManager.GetRegions().getIndex( regionName ) );
  }

  // node data (for transmissibility computation)

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // face data

  // get the face-based DOF numbers for the assembly
  string const faceDofKey = dofManager->getKey( viewKeyStruct::facePressureString );
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );

  // get the element dof numbers for the assembly
  string const elemDofKey = dofManager->getKey( viewKeyStruct::pressureString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = mesh.getElemManager()->ConstructViewAccessor< array1d< globalIndex >,
                                                                arrayView1d< globalIndex const > >( elemDofKey );

  // get the face-centered pressures
  arrayView1d< real64 const > const & facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString );
  arrayView1d< real64 const > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  // get the face-centered depth
  arrayView1d< real64 const > const & faceGravCoef =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // get the face-to-nodes connectivity for the transmissibility calculation
  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList();

  array2d< localIndex > const & elemRegionList    = faceManager.elementRegionList();
  array2d< localIndex > const & elemSubRegionList = faceManager.elementSubRegionList();
  array2d< localIndex > const & elemList          = faceManager.elementList();

  // tolerance for transmissibility calculation
  real64 const lengthTolerance = domain->getMeshBody( 0 )->getGlobalLengthScale() * m_areaRelTol;

  forTargetSubRegionsComplete< CellElementSubRegion >( mesh,
                                                       [&]( localIndex const,
                                                            localIndex const er,
                                                            localIndex const esr,
                                                            ElementRegionBase const &,
                                                            auto const & subRegion )
  {
    KernelLaunchSelector< FluxKernel >( subRegion.numFacesPerElement(),
                                        er,
                                        esr,
                                        subRegion,
                                        regionFilter,
                                        mesh,
                                        nodePosition,
                                        elemRegionList,
                                        elemSubRegionList,
                                        elemList,
                                        faceToNodes,
                                        faceDofNumber,
                                        facePres,
                                        dFacePres,
                                        faceGravCoef,
                                        m_density[er][esr],
                                        m_dDens_dPres[er][esr],
                                        m_mobility.toViewConst(),
                                        m_dMobility_dPres.toViewConst(),
                                        lengthTolerance,
                                        dt,
                                        dofManager,
                                        matrix,
                                        rhs );
  } );
}

void
SinglePhaseHybridFVM::ApplyBoundaryConditions( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                               real64 const GEOSX_UNUSED_PARAM( dt ),
                                               DomainPartition * const GEOSX_UNUSED_PARAM( domain ),
                                               DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                               ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                               ParallelVector & GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  // TODO: implement boundary conditions the hybrid way

}

real64 SinglePhaseHybridFVM::CalculateResidualNorm( DomainPartition const * const domain,
                                                    DofManager const & dofManager,
                                                    ParallelVector const & rhs )
{
  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const * const faceManager = mesh.getFaceManager();

  // here we compute the cell-centered residual norm in the derived class
  // to avoid duplicating a synchronization point

  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  string const elemDofKey = dofManager.getKey( viewKeyStruct::pressureString );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString );

  // local residual
  array1d< real64 > localResidualNorm( 2 ), globalResidualNorm( 2 );
  globalResidualNorm[EquationType::MASS_CONS]  = 0;
  globalResidualNorm[EquationType::CONSTRAINT] = 0;
  localResidualNorm[EquationType::MASS_CONS]   = 0;
  localResidualNorm[EquationType::CONSTRAINT]  = 0;


  // 1. Compute the residual for the mass conservation equations

  // compute the norm of local residual scaled by cell pore volume
  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase const &,
                                    ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & elemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( elemDofKey );

    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d< real64 const > const & refPoro = m_porosityRef[er][esr];
    arrayView1d< real64 const > const & volume = m_volume[er][esr];
    arrayView1d< real64 const > const & dVol = m_deltaVolume[er][esr];
    arrayView2d< real64 const > const & dens = m_density[er][esr];

    localIndex const subRegionSize = subRegion.size();
    for( localIndex a = 0; a < subRegionSize; ++a )
    {
      if( elemGhostRank[a] < 0 )
      {
        localIndex const lid = rhs.getLocalRowID( elemDofNumber[a] );
        // TODO: implement a normalization that matches the normalization used in SinglePhaseFVM
        real64 const val = localResidual[lid] / ( refPoro[a] * dens[a][0] * ( volume[a] + dVol[a] ) );
        localResidualNorm[EquationType::MASS_CONS] += val * val;
      }
    }
  } );


  // 2. Compute the residual for the face-based constraints

  arrayView1d< integer const > const & faceGhostRank =
    faceManager->getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager->getReference< array1d< globalIndex > >( faceDofKey );

  array2d< localIndex > const & elemRegionList    = faceManager->elementRegionList();
  array2d< localIndex > const & elemSubRegionList = faceManager->elementSubRegionList();
  array2d< localIndex > const & elemList          = faceManager->elementList();

  for( localIndex iface = 0; iface < faceManager->size(); ++iface )
  {
    // if not ghost face and if adjacent to target region
    if( faceGhostRank[iface] < 0 && faceDofNumber[iface] >= 0 )
    {
      real64 normalizer = 0;
      localIndex elemCounter = 0;
      for( localIndex k=0; k<elemRegionList.size( 1 ); ++k )
      {
        localIndex const er  = elemRegionList[iface][k];
        localIndex const esr = elemSubRegionList[iface][k];
        localIndex const ei  = elemList[iface][k];

        bool const onBoundary = (er == -1 || esr == -1 || ei == -1);

        // if not on boundary, save the mobility and the upwDofNumber
        if( !onBoundary )
        {
          ElementRegionBase const * const region = mesh.getElemManager()->GetRegion( er );
          ElementSubRegionBase const * const subRegion = region->GetSubRegion( esr );

          arrayView1d< real64 const > const & elemVolume =
            subRegion->getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );

          normalizer += elemVolume[ei];
          elemCounter++;
        }
      }
      normalizer /= elemCounter;

      localIndex const lid    = rhs.getLocalRowID( faceDofNumber[iface] );
      real64 const val        = localResidual[lid] / normalizer;
      localResidualNorm[EquationType::CONSTRAINT] += val * val;
    }
  }


  // 3. Combine the two norms

  MpiWrapper::allReduce( localResidualNorm.data(), globalResidualNorm.data(), 2, MPI_SUM, MPI_COMM_GEOSX );

  real64 const maxNorm = (globalResidualNorm[EquationType::MASS_CONS] > globalResidualNorm[EquationType::CONSTRAINT])
                       ? globalResidualNorm[EquationType::MASS_CONS]
                       : globalResidualNorm[EquationType::CONSTRAINT];

  real64 const residual = sqrt( maxNorm );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output,
             "( Rfluid ) = (%4.2e) ; ",
             residual );
    std::cout<<output;
  }

  return residual;
}


bool
SinglePhaseHybridFVM::CheckSystemSolution( DomainPartition const * const domain,
                                           DofManager const & dofManager,
                                           ParallelVector const & solution,
                                           real64 const scalingFactor )
{
  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  real64 const * localSolution = solution.extractLocalVector();
  int localCheck = 1;

  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );

  forTargetSubRegionsComplete( mesh, [&]( localIndex const,
                                          localIndex const er,
                                          localIndex const esr,
                                          ElementRegionBase const &,
                                          ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d< real64 const > const & pres = m_pressure[er][esr];
    arrayView1d< real64 const > const & dPres = m_deltaPressure[er][esr];

    forAll< serialPolicy >( subRegion.size(), [&]( localIndex const ei )
    {
      if( elemGhostRank[ei] >= 0 )
      {
        return;
      }

      globalIndex const offset = dofNumber[ei];
      // extract solution and apply to dP
      {
        localIndex const lid = solution.getLocalRowID( offset );
        real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[lid];

        if( newPres < 0.0 )
        {
          localCheck = 0;
        }
      }

    } );
  } );

  int const globalCheck = MpiWrapper::Min( localCheck );

  return globalCheck;
}


void SinglePhaseHybridFVM::ApplySystemSolution( DofManager const & dofManager,
                                                ParallelVector const & solution,
                                                real64 const scalingFactor,
                                                DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // here we apply the cell-centered update in the derived class
  // to avoid duplicating a synchronization point

  // 1. apply the cell-centered update

  dofManager.addVectorToField( solution,
                               viewKeyStruct::pressureString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor );

  // 2. apply the face-based update

  dofManager.addVectorToField( solution,
                               viewKeyStruct::facePressureString,
                               viewKeyStruct::deltaFacePressureString,
                               scalingFactor );

  // 3. synchronize

  // the tags in fieldNames have to match the tags used in NeighborCommunicator.cpp
  std::map< string, string_array > fieldNames;
  fieldNames["face"].push_back( viewKeyStruct::deltaFacePressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain->getNeighbors() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateState( subRegion, targetIndex );
  } );
}


void SinglePhaseHybridFVM::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // 1. Reset the cell-centered fields
  SinglePhaseBase::ResetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  MeshLevel * const mesh          = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager * const faceManager = mesh->getFaceManager();

  // get the accumulated face pressure updates
  arrayView1d< real64 > & dFacePres =
    faceManager->getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString );

  forAll< serialPolicy >( faceManager->size(), [=] ( localIndex iface )
  {
    dFacePres[iface] = 0;
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseHybridFVM, std::string const &, Group * const )
} /* namespace geosx */
