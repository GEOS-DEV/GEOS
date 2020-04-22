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
 * @file CompositionalMultiphaseHybridFVM.cpp
 */

#include "CompositionalMultiphaseHybridFVM.hpp"

#include "common/TimingMacros.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVMKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace CompositionalMultiphaseHybridFVMKernels;
using namespace HybridFVMInnerProduct;

CompositionalMultiphaseHybridFVM::CompositionalMultiphaseHybridFVM( const std::string & name,
                                                                    Group * const parent ):
  CompositionalMultiphaseBase( name, parent ),
  m_faceDofKey( "" ),
  m_elemDofKey( "" ),
  m_areaRelTol( 1e-8 ),
  m_ipType( InnerProductType::QUASI_TPFA ), // hard-coded for now
  m_orthonormalizeWithSVD( false ) // hard-coded for now
{
}

void CompositionalMultiphaseHybridFVM::RegisterDataOnMesh( Group * const MeshBodies )
{
  // 1) Register the elem-centered data
  CompositionalMultiphaseBase::RegisterDataOnMesh( MeshBodies );

  // 2) Register the face data
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    FaceManager & faceManager = *meshLevel.getFaceManager();

    // primary variables: face potentials
    faceManager.registerWrapper< array2d< real64 > >( viewKeyStruct::facePhasePotentialString )->
      setPlotLevel( PlotLevel::LEVEL_0 )->
      setRegisteringObjects( this->getName())->
      setDescription( "An array that holds the phase potentials at the faces." );

    faceManager.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaFacePhasePotentialString )->
      setPlotLevel( PlotLevel::LEVEL_0 )->
      setRegisteringObjects( this->getName())->
      setDescription( "An array that holds the accumulated phase potential updates at the faces." );
  }
}

  
void CompositionalMultiphaseHybridFVM::ImplicitStepSetup( real64 const & time_n,
                                                          real64 const & dt,
                                                          DomainPartition * const domain,
                                                          DofManager & dofManager,
                                                          ParallelMatrix & matrix,
                                                          ParallelVector & rhs,
                                                          ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  // setup the elem-centered fields
  CompositionalMultiphaseBase::ImplicitStepSetup( time_n, dt, domain, dofManager, matrix, rhs, solution );

  // setup the face fields
  MeshLevel & meshLevel = *domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *meshLevel.getFaceManager();

  // get the accumulated pressure updates
  arrayView2d< real64 > & dFacePotential =
    faceManager.getReference< array2d< real64 > >( viewKeyStruct::deltaFacePhasePotentialString );

  forAll< serialPolicy >( faceManager.size(), [=] ( localIndex const iface )
  {
    for( localIndex ip = 0; ip < m_numPhases; ++ip )
    {  
      dFacePotential[iface][ip] = 0;
    }
  } );

}  


void CompositionalMultiphaseHybridFVM::ImplicitStepComplete( real64 const & time_n,
                                                             real64 const & dt,
                                                             DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  // increment the elem-centered fields
  CompositionalMultiphaseBase::ImplicitStepComplete( time_n, dt, domain );

  // increment the face fields
  MeshLevel & meshLevel = *domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *meshLevel.getFaceManager();

  // get the face-based DOF numbers
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( m_faceDofKey );

  // get the face-based potentials
  arrayView2d< real64 > const & facePotential =
    faceManager.getReference< array2d< real64 > >( viewKeyStruct::facePhasePotentialString );
  arrayView2d< real64 const > const & dFacePotential =
    faceManager.getReference< array2d< real64 > >( viewKeyStruct::deltaFacePhasePotentialString );

  forAll< serialPolicy >( faceManager.size(), [=] ( localIndex const iface )
  {
    // update if face is in target region
    if( faceDofNumber[iface] >= 0 )
    {
      for( localIndex ip = 0; ip < m_numPhases; ++ip )
      {
        facePotential[iface][ip] += dFacePotential[iface][ip];
      }
    }
  } );
}  


void CompositionalMultiphaseHybridFVM::SetupDofs( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                                  DofManager & dofManager ) const
{
  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::elemDofFieldString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString,
                          viewKeyStruct::elemDofFieldString,
                          DofManager::Connector::Face );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::faceDofFieldString,
                       DofManager::Location::Face,
                       m_numPhases,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::faceDofFieldString,
                          viewKeyStruct::faceDofFieldString,
                          DofManager::Connector::Elem );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::faceDofFieldString,
                          viewKeyStruct::elemDofFieldString,
                          DofManager::Connector::Elem,
                          true );

}  


void CompositionalMultiphaseHybridFVM::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                          real64 const dt,
                                                          DomainPartition const * const domain,
                                                          DofManager const * const dofManager,
                                                          ParallelMatrix * const matrix,
                                                          ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *mesh.getElemManager();
  NodeManager const & nodeManager = *mesh.getNodeManager();
  FaceManager const & faceManager = *mesh.getFaceManager();

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
  SortedArray< localIndex > regionFilter;
  for( string const & regionName : targetRegionNames() )
  {
    regionFilter.insert( elemManager.GetRegions().getIndex( regionName ) );
  }

  // node data (for transmissibility computation)

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition =
    nodeManager.referencePosition();


  // face data

  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );

  // get the face-centered pressures
  arrayView2d< real64 const > const & facePotential =
    faceManager.getReference< array2d< real64 > >( viewKeyStruct::facePhasePotentialString );
  arrayView2d< real64 const > const & dFacePotential =
    faceManager.getReference< array2d< real64 > >( viewKeyStruct::deltaFacePhasePotentialString );

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

  forTargetSubRegionsComplete< CellElementSubRegion,
                               FaceElementSubRegion >( mesh,
                                                       [&]( localIndex const,
                                                            localIndex const er,
                                                            localIndex const esr,
                                                            ElementRegionBase const &,
                                                            auto const & subRegion )
  {
    FluxLaunch( er,
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
                facePotential,
                dFacePotential,
                faceGravCoef,
                lengthTolerance,
                dt,
                dofManager,
                matrix,
                rhs );
  } );
}  

void CompositionalMultiphaseHybridFVM::FluxLaunch( localIndex GEOSX_UNUSED_PARAM( er ),
                                                   localIndex GEOSX_UNUSED_PARAM( esr ),
                                                   FaceElementSubRegion const & GEOSX_UNUSED_PARAM( subRegion ),
                                                   SortedArray< localIndex > GEOSX_UNUSED_PARAM( regionFilter ),
                                                   MeshLevel const & GEOSX_UNUSED_PARAM( mesh ),
                                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & GEOSX_UNUSED_PARAM( nodePosition ),
                                                   array2d< localIndex > const & GEOSX_UNUSED_PARAM( elemRegionList ),
                                                   array2d< localIndex > const & GEOSX_UNUSED_PARAM( elemSubRegionList ),
                                                   array2d< localIndex > const & GEOSX_UNUSED_PARAM( elemList ),
                                                   ArrayOfArraysView< localIndex const > const & GEOSX_UNUSED_PARAM( faceToNodes ),
                                                   arrayView1d< globalIndex const > const & GEOSX_UNUSED_PARAM( faceDofNumber ),
                                                   arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( facePotential ),
                                                   arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( dFacePotential ),
                                                   arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( faceGravCoef ),
                                                   real64 const GEOSX_UNUSED_PARAM( lengthTolerance ),
                                                   real64 const GEOSX_UNUSED_PARAM( dt ),
                                                   DofManager const * const GEOSX_UNUSED_PARAM( dofManager ),
                                                   ParallelMatrix * const GEOSX_UNUSED_PARAM( matrix ),
                                                   ParallelVector * const GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_LOG_RANK( "Support for FaceElementSubRegion is not implemented in the Hybrid FVM scheme" );
}

void CompositionalMultiphaseHybridFVM::FluxLaunch( localIndex er,
                                                   localIndex esr,
                                                   CellElementSubRegion const & subRegion,
                                                   SortedArray< localIndex > regionFilter,
                                                   MeshLevel const & mesh,
                                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                                   array2d< localIndex > const & elemRegionList,
                                                   array2d< localIndex > const & elemSubRegionList,
                                                   array2d< localIndex > const & elemList,
                                                   ArrayOfArraysView< localIndex const > const & faceToNodes,
                                                   arrayView1d< globalIndex const > const & faceDofNumber,
                                                   arrayView2d< real64 const > const & facePotential,
                                                   arrayView2d< real64 const > const & dFacePotential,
                                                   arrayView1d< real64 const > const & faceGravCoef,
                                                   real64 const lengthTolerance,
                                                   real64 const dt,
                                                   DofManager const * const dofManager,
                                                   ParallelMatrix * const matrix,
                                                   ParallelVector * const rhs )
{

  // max number of faces allowed in an element
  localIndex constexpr maxNumFaces = MAX_NUM_FACES;

  // elem data

  // get the elem-centered DOF numbers and ghost rank for the assembly
  string const elemDofKey = dofManager->getKey( viewKeyStruct::elemDofFieldString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = mesh.getElemManager()->ConstructViewAccessor< array1d< globalIndex >,
                                                                arrayView1d< globalIndex const > >( elemDofKey );
  arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

  // get the map from elem to faces
  arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

  // get the elem-centered pressures
  arrayView1d< real64 const > const & elemPres  = m_pressure[er][esr];
  arrayView1d< real64 const > const & dElemPres = m_deltaPressure[er][esr];

  // get the element data needed for transmissibility computation
  arrayView1d< R1Tensor const > const & elemCenter =
    subRegion.template getReference< array1d< R1Tensor > >( CellBlock::viewKeyStruct::elementCenterString );
  arrayView1d< real64 const > const & elemVolume =
    subRegion.template getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );
  arrayView1d< R1Tensor const > const & elemPerm =
    subRegion.template getReference< array1d< R1Tensor > >( viewKeyStruct::permeabilityString );

  // get the elem-centered depth
  arrayView1d< real64 const > const & elemGravCoef =
    subRegion.template getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // assemble the residual and Jacobian element by element
  // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces

  forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
  {
    if( elemGhostRank[ei] < 0 )
    {

      localIndex const numFacesInElem = elemToFaces.size( 1 );

      // transmissibility matrix
      stackArray2d< real64, maxNumFaces *maxNumFaces > transMatrix( numFacesInElem,
                                                                    numFacesInElem );

      // recompute the local transmissibility matrix at each iteration
      // we can decide later to precompute transMatrix if needed
      ComputeTransmissibilityMatrix( nodePosition,
                                     faceToNodes,
                                     elemToFaces[ei],
                                     elemCenter[ei],
                                     elemVolume[ei],
                                     elemPerm[ei],
                                     lengthTolerance,
                                     m_orthonormalizeWithSVD,
                                     transMatrix );

      // perform flux assembly in this element
      CompositionalMultiphaseHybridFVMKernels::FluxKernel< CellElementSubRegion >::Compute( m_numPhases,
                                                                                            m_numComponents,
                                                                                            er, esr, ei,
                                                                                            regionFilter,
                                                                                            elemRegionList,
                                                                                            elemSubRegionList,
                                                                                            elemList,
                                                                                            faceDofNumber,
                                                                                            facePotential,
                                                                                            dFacePotential,
                                                                                            faceGravCoef,
                                                                                            elemToFaces[ei],
                                                                                            elemPres[ei],
                                                                                            dElemPres[ei],
                                                                                            elemGravCoef[ei],
                                                                                            m_phaseDens.toViewConst(),
                                                                                            m_dPhaseDens_dPres.toViewConst(),
                                                                                            m_dPhaseDens_dComp.toViewConst(),
                                                                                            m_phaseCompFrac.toViewConst(),
                                                                                            m_dPhaseCompFrac_dPres.toViewConst(),
                                                                                            m_dPhaseCompFrac_dComp.toViewConst(),
                                                                                            m_dCompFrac_dCompDens.toViewConst(),
                                                                                            m_phaseMob.toViewConst(),
                                                                                            m_dPhaseMob_dPres.toViewConst(),
                                                                                            m_dPhaseMob_dCompDens.toViewConst(),
                                                                                            elemDofNumber.toViewConst(),
                                                                                            transMatrix,
                                                                                            dt,
                                                                                            matrix,
                                                                                            rhs );

    }
  } );
}

void
CompositionalMultiphaseHybridFVM::ApplyBoundaryConditions( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                           real64 const GEOSX_UNUSED_PARAM( dt ),
                                                           DomainPartition * const GEOSX_UNUSED_PARAM( domain ),
                                                           DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                           ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                                           ParallelVector & GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  // TODO: implement boundary conditions the mimetic way

}

  
void CompositionalMultiphaseHybridFVM::ResizeFields( MeshLevel & meshLevel )
{
  CompositionalMultiphaseBase::ResizeFields( meshLevel );

  FaceManager & faceManager = *meshLevel.getFaceManager();

  faceManager.getReference< array2d< real64 > >( viewKeyStruct::facePhasePotentialString ).resizeDimension< 1 >( m_numPhases );
  faceManager.getReference< array2d< real64 > >( viewKeyStruct::deltaFacePhasePotentialString ).resizeDimension< 1 >( m_numPhases );

}


real64 CompositionalMultiphaseHybridFVM::CalculateResidualNorm( DomainPartition const * const domain,
                                                                DofManager const & dofManager,
                                                                ParallelVector const & rhs )
{
  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const & faceManager = *mesh.getFaceManager();
  
  real64 const * localResidual = rhs.extractLocalVector();

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString );

  // local residual
  real64 localResidualNorm = 0;

  // 1. Compute the residual for the mass conservation equations 

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase const &,
                                    ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );

    arrayView1d< integer const >     const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d< real64 const >      const & refPoro       = m_porosityRef[er][esr];
    arrayView1d< real64 const >      const & volume        = m_volume[er][esr];
    arrayView2d< real64 const >      const & totalDens     = m_totalDens[er][esr];

    localIndex const subRegionSize = subRegion.size();
    for( localIndex ei = 0; ei < subRegionSize; ++ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        globalIndex const offset = elemDofNumber[ei];
        for( localIndex idof = 0; idof < m_numDofPerCell; ++idof )
        {
          localIndex const lid = rhs.getLocalRowID( offset + idof );
          real64 const val = localResidual[lid] / (totalDens[ei][0] * refPoro[ei] * volume[ei]);
          localResidualNorm += val * val;
        }
      }
    }
  } );


  // 2. Compute the residual for the face-based constraints

  arrayView1d< integer const > const & faceGhostRank =
    faceManager.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );

  for( localIndex iface = 0; iface < faceManager.size(); ++iface )
  {
    // if not ghost face and if adjacent to target region
    if( faceGhostRank[iface] < 0 && faceDofNumber[iface] >= 0 )
    {
      for( localIndex ip = 0; ip < m_numPhases; ++ip )
      {
        localIndex const lid = rhs.getLocalRowID( faceDofNumber[iface] + ip );
        real64 const normalizer = 1; // TODO: compute the normalizer here
        real64 const val = localResidual[lid] / normalizer;
        localResidualNorm += val * val;
      }
    }
  }
  
  // compute global residual norm
  realT globalResidualNorm;
  MpiWrapper::allReduce( &localResidualNorm, &globalResidualNorm, 1, MPI_SUM, MPI_COMM_GEOSX );

  real64 const residualNorm = sqrt( globalResidualNorm );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output,
             "( Rfluid ) = (%4.2e) ; ",
             residualNorm );
    std::cout<<output;
  }
  
  return residualNorm;
}


void CompositionalMultiphaseHybridFVM::ApplySystemSolution( DofManager const & dofManager, 
                                                            ParallelVector const & solution,
                                                            real64 const scalingFactor,
                                                            DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // 1. apply the elem-based update

  dofManager.addVectorToField( solution,
                               viewKeyStruct::elemDofFieldString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( solution,
                               viewKeyStruct::elemDofFieldString,
                               viewKeyStruct::deltaGlobalCompDensityString,
                               scalingFactor,
                               1, m_numDofPerCell );

  // 2. apply the face-based update

  dofManager.addVectorToField( solution,
                               viewKeyStruct::faceDofFieldString,
                               viewKeyStruct::deltaFacePhasePotentialString,
                               scalingFactor,
                               0, m_numPhases );

  // 3. synchronize

  std::map< string, string_array > fieldNames;
  fieldNames["face"].push_back( viewKeyStruct::deltaFacePhasePotentialString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaGlobalCompDensityString );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         &mesh,
                                         domain->getNeighbors() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateState( subRegion, targetIndex );
  } );
}

  
void CompositionalMultiphaseHybridFVM::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // 1. Reset the cell-centered fields
  CompositionalMultiphaseBase::ResetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  MeshLevel & mesh = *domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *mesh.getFaceManager();

  // get the accumulated face potential updates
  arrayView2d< real64 > & dFacePotential =
    faceManager.getReference< array2d< real64 > >( viewKeyStruct::deltaFacePhasePotentialString );

  forAll< serialPolicy >( faceManager.size(), [=] ( localIndex const iface )
  {
    for( localIndex ip = 0; ip < m_numPhases; ++ip )
    { 
      dFacePotential[iface][ip] = 0;
    }
  } );
}


void CompositionalMultiphaseHybridFVM::ComputeTransmissibilityMatrix( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                                                      ArrayOfArraysView< localIndex const > const & faceToNodes,
                                                                      arraySlice1d< localIndex const > const elemToFaces,
                                                                      R1Tensor const & elemCenter,
                                                                      real64 const & elemVolume,
                                                                      R1Tensor const & elemPerm,
                                                                      real64 const & lengthTolerance,
                                                                      bool const & orthonormalizeWithSVD,
                                                                      stackArray2d< real64, MAX_NUM_FACES
                                                                                           *MAX_NUM_FACES > & transMatrix ) const
{
  // TODO: remove the switchyard
  switch( m_ipType )
  {
    case InnerProductType::TPFA:
    {
      TPFACellInnerProductKernel::Compute( nodePosition,
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
      // since it reduces to TPFA on K-orthogonal meshes...
      QTPFACellInnerProductKernel::Compute( nodePosition,
                                            faceToNodes,
                                            elemToFaces,
                                            elemCenter,
                                            elemVolume,
                                            elemPerm,
                                            2,
                                            lengthTolerance,
                                            orthonormalizeWithSVD,
                                            transMatrix );
      break;
    }

    default:
    {
      GEOSX_LOG_RANK( "Unknown inner product in CompositionalMultiphaseHybridFVM" );
      break;
    }
  }
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseHybridFVM, std::string const &, Group * const )
} /* namespace geosx */
