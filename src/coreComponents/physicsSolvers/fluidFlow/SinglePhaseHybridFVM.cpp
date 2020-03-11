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
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVMKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace SinglePhaseHybridFVMKernels;
  
template< typename BASE >  
SinglePhaseHybridFVM<BASE>::SinglePhaseHybridFVM( const std::string& name,
                                                  Group * const parent ):
  BASE(name, parent),
  m_faceDofKey(""),
  m_areaRelTol(1e-8),
  m_ipType(InnerProductType::QUASI_TPFA),
  m_orthonormalizeWithSVD(false) 
{
 
  // one cell-centered dof per cell
  m_numDofPerCell = 1;
}

template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::RegisterDataOnMesh(Group * const MeshBodies)
{

  // 1) Register the cell-centered data
  BASE::RegisterDataOnMesh(MeshBodies);

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

  
template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::ImplicitStepSetup( real64 const & time_n,
                                                    real64 const & dt,
                                                    DomainPartition * const domain,
                                                    DofManager & dofManager,
                                                    ParallelMatrix & matrix,
                                                    ParallelVector & rhs,
                                                    ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  // setup the cell-centered fields
  BASE::ImplicitStepSetup( time_n, dt, domain, dofManager, matrix, rhs, solution );

  // setup the face fields
  MeshLevel * const meshLevel     = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = meshLevel->getFaceManager();

  // get the face-based DOF numbers 
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString );

  // save the face Dof key for use in two functions
  // that do not have acces to the coupled solver dofManager
  // namely ResetStateToBeginningOfStep and ImplicitStepComplete
  m_faceDofKey = faceDofKey;
  
  // get the accumulated pressure updates
  arrayView1d<real64> & dFacePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

  forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  {
    dFacePres[iface] = 0;
  });

}

template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::ImplicitStepComplete( real64 const & time_n,
                                                       real64 const & dt, 
                                                       DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  // increment the cell-centered fields
  BASE::ImplicitStepComplete( time_n, dt, domain );

  // increment the face fields
  MeshLevel * const meshLevel     = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  FaceManager * const faceManager = meshLevel->getFaceManager();

  // get the face-based DOF numbers 
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference< array1d<globalIndex> >( m_faceDofKey );  

  // get the face-based pressures 
  arrayView1d<real64> const & facePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::facePressureString);
  arrayView1d<real64> const & dFacePres =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::deltaFacePressureString);

  forall_in_range<serialPolicy>( 0, faceManager->size(), GEOSX_LAMBDA ( localIndex iface )
  {
    // update if face is in target region
    if (faceDofNumber[iface] >= 0)
    {
      facePres[iface] += dFacePres[iface];
    }
    dFacePres[iface] = 0;    
  });
}

template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::SetupDofs( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                            DofManager & dofManager ) const
{
  
  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::pressureString,
                       DofManager::Location::Elem,
                       m_targetRegions );

  dofManager.addCoupling( viewKeyStruct::pressureString,
                          viewKeyStruct::pressureString,
                          DofManager::Connectivity::Face );
  
  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::facePressureString,
                       DofManager::Location::Face,
                       m_targetRegions );

  dofManager.addCoupling( viewKeyStruct::facePressureString,
                          viewKeyStruct::facePressureString,
                          DofManager::Connectivity::Elem );
  
  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::facePressureString,
                          viewKeyStruct::pressureString,
                          DofManager::Connectivity::Elem,
                          true);

}

template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
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
  SortedArray<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }
  
  // node data (for transmissibility computation)

  arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition = nodeManager->referencePosition();

  
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
  arrayView1d<real64 const> const & faceGravCoef =
    faceManager->getReference<array1d<real64>>(viewKeyStruct::gravityCoefString);
  
  // get the face-to-nodes connectivity for the transmissibility calculation
  ArrayOfArraysView<localIndex const> const & faceToNodes = faceManager->nodeList();

  array2d<localIndex> const & elemRegionList    = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList          = faceManager->elementList(); 

  
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
                facePres,
                dFacePres,
                faceGravCoef,
                lengthTolerance,
                dt,             
                dofManager,
                matrix,
                rhs );
  });
}


template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::FluxLaunch( localIndex GEOSX_UNUSED_PARAM( er ),
                                             localIndex GEOSX_UNUSED_PARAM( esr ),
                                             FaceElementSubRegion const * const GEOSX_UNUSED_PARAM( subRegion ),
                                             SortedArray<localIndex> GEOSX_UNUSED_PARAM( regionFilter ),
                                             MeshLevel const * const GEOSX_UNUSED_PARAM( mesh ),
                                             arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & GEOSX_UNUSED_PARAM( nodePosition ),
                                             array2d<localIndex> const & GEOSX_UNUSED_PARAM( elemRegionList ),
                                             array2d<localIndex> const & GEOSX_UNUSED_PARAM( elemSubRegionList ),
                                             array2d<localIndex> const & GEOSX_UNUSED_PARAM( elemList ),
                                             ArrayOfArraysView<localIndex const> const & GEOSX_UNUSED_PARAM( faceToNodes ),
                                             arrayView1d<globalIndex const> const & GEOSX_UNUSED_PARAM( faceDofNumber ),
                                             arrayView1d<real64 const> const & GEOSX_UNUSED_PARAM( facePres ),
                                             arrayView1d<real64 const> const & GEOSX_UNUSED_PARAM( dFacePres ),
                                             arrayView1d<real64 const> const & GEOSX_UNUSED_PARAM( faceGravCoef ),
                                             real64 const GEOSX_UNUSED_PARAM( lengthTolerance ),
                                             real64 const GEOSX_UNUSED_PARAM( dt ),                                  
                                             DofManager const * const GEOSX_UNUSED_PARAM( dofManager ),
                                             ParallelMatrix * const GEOSX_UNUSED_PARAM( matrix ),
                                             ParallelVector * const GEOSX_UNUSED_PARAM( rhs ) )
{
  // not implemented yet
}


template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::FluxLaunch( localIndex er,
                                             localIndex esr,
                                             CellElementSubRegion const * const subRegion,
                                             SortedArray<localIndex> regionFilter,
                                             MeshLevel const * const mesh,
                                             arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition,
                                             array2d<localIndex> const & elemRegionList,
                                             array2d<localIndex> const & elemSubRegionList,
                                             array2d<localIndex> const & elemList, 
                                             ArrayOfArraysView<localIndex const> const & faceToNodes,                                
                                             arrayView1d<globalIndex const> const & faceDofNumber,
                                             arrayView1d<real64 const> const & facePres,
                                             arrayView1d<real64 const> const & dFacePres,
                                             arrayView1d<real64 const> const & faceGravCoef,
                                             real64 const lengthTolerance,
                                             real64 const dt,        
                                             DofManager const * const dofManager,
                                             ParallelMatrix * const matrix,
                                             ParallelVector * const rhs )
{
  // max number of faces allowed in an element 
  localIndex constexpr maxNumFaces = MAX_NUM_FACES;

  // elem data
    
  // get the cell-centered DOF numbers and ghost rank for the assembly
  string const elemDofKey = dofManager->getKey( viewKeyStruct::pressureString );
  arrayView1d<globalIndex const> const & elemDofNumber =
    subRegion->template getReference< array1d<globalIndex> >( elemDofKey ); 
  arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];   

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
   subRegion->template getReference< array1d<R1Tensor> >( CellBlock::viewKeyStruct::elementCenterString );
  arrayView1d<real64 const> const & elemVolume =
   subRegion->template getReference< array1d<real64> >( CellBlock::viewKeyStruct::elementVolumeString ); 
  arrayView1d<R1Tensor const> const & elemPerm =
   subRegion->template getReference< array1d<R1Tensor> >( viewKeyStruct::permeabilityString ); 

  // get the cell-centered depth 
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

      // perform flux assembly in this element in two steps:
      FluxKernel<CellElementSubRegion>::Compute( er, esr, ei,
                                                 regionFilter,
                                                 mesh,
                                                 elemRegionList,
                                                 elemSubRegionList,
                                                 elemList,
                                                 faceDofNumber,
                                                 facePres,
                                                 dFacePres,
                                                 faceGravCoef,
                                                 elemToFaces[ei],
                                                 elemDofNumber[ei],
                                                 elemDofKey,
                                                 elemPres[ei],
                                                 dElemPres[ei],
                                                 elemGravCoef[ei],
                                                 elemDens[ei][0],
                                                 dElemDens_dp[ei][0],
                                                 m_mobility, 
                                                 m_dMobility_dPres,
                                                 transMatrix,
                                                 dt,
                                                 matrix,
                                                 rhs );
                                                                           
    }
  });
}

template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::ApplyBoundaryConditions( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                          real64 const GEOSX_UNUSED_PARAM( dt ),
                                                          DomainPartition * const GEOSX_UNUSED_PARAM( domain ),
                                                          DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                          ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                                          ParallelVector & GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  // TODO: implement boundary conditions the mimetic way

}


template< typename BASE >  
real64 SinglePhaseHybridFVM<BASE>::CalculateResidualNorm( DomainPartition const * const domain,
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
  globalResidualNorm[EquationType::MASS_CONS]  = 0; 
  globalResidualNorm[EquationType::CONSTRAINT] = 0; 
  localResidualNorm[EquationType::MASS_CONS]   = 0; 
  localResidualNorm[EquationType::CONSTRAINT]  = 0; 

  
  // 1. Compute the residual for the mass conservation equations

  // compute the norm of local residual scaled by cell pore volume
  this->applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                       ElementRegionBase const * const GEOSX_UNUSED_PARAM( region ),
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
        // TODO: implement a normalization that matches the normalization used in SinglePhaseFVM
        real64 const val = localResidual[lid] / (refPoro[a] * dens[a][0] * ( volume[a] + dVol[a]));
        localResidualNorm[EquationType::MASS_CONS] += val * val;
      }
    }
  });

  
  // 2. Compute the residual for the face-based constraints

  arrayView1d<integer const> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d<globalIndex const> const & faceDofNumber =
    faceManager->getReference< array1d<globalIndex> >( faceDofKey );

  array2d<localIndex> const & elemRegionList    = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList          = faceManager->elementList(); 
  
  for( localIndex iface = 0 ; iface < faceManager->size(); ++iface )
  {
    // if not ghost face and if adjacent to target region
    if (faceGhostRank[iface] < 0 && faceDofNumber[iface] >= 0)
    {
      real64 normalizer = 0;
      localIndex elemCounter = 0;
      for (localIndex k=0; k<elemRegionList.size(1); ++k)
      {
        localIndex const er  = elemRegionList[iface][k];
        localIndex const esr = elemSubRegionList[iface][k];
        localIndex const ei  = elemList[iface][k];

        bool const onBoundary = (er == -1 || esr == -1 || ei == -1);

        // if not on boundary, save the mobility and the upwDofNumber  
        if ( !onBoundary )
        {
            ElementRegionBase const * const region =
              Group::group_cast<ElementRegionBase const *>(mesh->getElemManager()->GetRegion(er));
            ElementSubRegionBase const * const subRegion =
              Group::group_cast<ElementSubRegionBase const *>(region->GetSubRegion(esr));
  
            arrayView1d<real64 const> const & elemVolume =
              subRegion->getReference< array1d<real64> >( CellBlock::viewKeyStruct::elementVolumeString ); 

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

  MpiWrapper::allReduce(localResidualNorm.data(), globalResidualNorm.data(), 2, MPI_SUM, MPI_COMM_GEOSX);
 
  real64 maxNorm = (globalResidualNorm[EquationType::MASS_CONS] > globalResidualNorm[EquationType::CONSTRAINT])
                 ? globalResidualNorm[EquationType::MASS_CONS]
                 : globalResidualNorm[EquationType::CONSTRAINT];
  
  return sqrt( maxNorm );
}


template< typename BASE >  
bool SinglePhaseHybridFVM<BASE>::CheckSystemSolution( DomainPartition const * const domain,
                                                      DofManager const & dofManager,
                                                      ParallelVector const & solution,
                                                      real64 const scalingFactor )
{
  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  real64 const * localSolution = solution.extractLocalVector();
  int localCheck = 1;

  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );

  this->applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                       ElementRegionBase const * const GEOSX_UNUSED_PARAM( region ),
                                       ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64 const> const & pres      = m_pressure[er][esr];
    arrayView1d<real64 const> const & dPres     = m_deltaPressure[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), [&] ( localIndex ei )
    {
      if (elemGhostRank[ei] >= 0)
      {
        return;
      }

      globalIndex const offset = dofNumber[ei];
      // extract solution and apply to dP
      {
        localIndex const lid = solution.getLocalRowID( offset );
        real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[lid];

        if (newPres < 0.0)
        {
          localCheck = 0;
        }
      }

    });
  });
  int globalCheck;

  MpiWrapper::allReduce( &localCheck,
                         &globalCheck,
                         1,
                         MPI_MIN,
                         MPI_COMM_GEOSX );

  bool result = true;
  if (globalCheck == 0)
  {
    result = false;
  }
  return result;
}

template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::ApplySystemSolution( DofManager const & dofManager,
                                                      ParallelVector const & solution,
                                                      real64 const scalingFactor,
                                                      DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

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
  std::map<string, string_array> fieldNames;
  fieldNames["face"].push_back( viewKeyStruct::deltaFacePressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );
  
  this->applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion )
  {
    this->UpdateState( subRegion );
  } );

}


template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // 1. Reset the cell-centered fields
  BASE::ResetStateToBeginningOfStep( domain );

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
template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::ComputeTPFAInnerProduct( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition, 
                                                          ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                                          arraySlice1d<localIndex const> const elemToFaces,
                                                          R1Tensor const & elemCenter,
                                                          R1Tensor const & elemPerm,
                                                          real64   const & lengthTolerance,
                                                          stackArray2d<real64, SinglePhaseHybridFVMKernels::MAX_NUM_FACES
                                                                              *SinglePhaseHybridFVMKernels::MAX_NUM_FACES> const & transMatrix ) const
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


template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::ComputeQFamilyInnerProduct( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition, 
                                                             ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                                             arraySlice1d<localIndex const> const elemToFaces,
                                                             R1Tensor const & elemCenter,
                                                             real64   const & elemVolume,
                                                             R1Tensor const & elemPerm,
                                                             real64   const & tParam, 
                                                             real64   const & lengthTolerance,
                                                             stackArray2d<real64, SinglePhaseHybridFVMKernels::MAX_NUM_FACES
                                                                                 *SinglePhaseHybridFVMKernels::MAX_NUM_FACES> const & transMatrix ) const
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
   
    normalsMat(ifaceLoc,0) = faceNormal(0);
    normalsMat(ifaceLoc,1) = faceNormal(1);
    normalsMat(ifaceLoc,2) = faceNormal(2);
    
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
 

// TODO: template on the type of inner product
// TODO: template on the type of subRegion to have a special treatment for fractures
template< typename BASE >  
void SinglePhaseHybridFVM<BASE>::ComputeTransmissibilityMatrix( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & nodePosition, 
                                                                ArrayOfArraysView<localIndex const> const & faceToNodes, 
                                                                arraySlice1d<localIndex const> const elemToFaces,
                                                                R1Tensor const & elemCenter,
                                                                real64   const & elemVolume,
                                                                R1Tensor const & elemPerm,
                                                                real64   const & lengthTolerance,
                                                                stackArray2d<real64, SinglePhaseHybridFVMKernels::MAX_NUM_FACES
                                                                                    *SinglePhaseHybridFVMKernels::MAX_NUM_FACES> const & transMatrix ) const
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
typedef SinglePhaseHybridFVM<SinglePhaseBase> NoProppant;
typedef SinglePhaseHybridFVM<SinglePhaseProppantBase> Proppant;
REGISTER_CATALOG_ENTRY( SolverBase, NoProppant, std::string const &, Group * const )
REGISTER_CATALOG_ENTRY( SolverBase, Proppant, std::string const &, Group * const )
}

} /* namespace geosx */
