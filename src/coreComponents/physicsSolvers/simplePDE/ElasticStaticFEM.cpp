/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "ElasticStaticFEM.hpp"
#include "ElasticStaticFEMKernels.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "common/TimingMacros.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "managers/DomainPartition.hpp"

#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

ElasticStaticFEM::ElasticStaticFEM( std::string const & name,
                                    Group * const parent ):
  SolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::solidMaterialNamesString, &m_solidMaterialNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Material name for computing stress increment from strain increment" );
}

ElasticStaticFEM::~ElasticStaticFEM()
{}


void ElasticStaticFEM::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    NodeManager * const nodes = mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getNodeManager();

    nodes->registerWrapper< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( keys::TotalDisplacement )->
      setPlotLevel( PlotLevel::LEVEL_0 )->
      setDescription( "Nodal displacement." )->
      reference().resizeDimension< 1 >( 3 );

    nodes->registerWrapper< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( keys::IncrementalDisplacement )->
      setPlotLevel( PlotLevel::LEVEL_3 )->
      setDescription( "Nodal displacement increment at the current time step." )->
      reference().resizeDimension< 1 >( 3 );
  }
}

real64 ElasticStaticFEM::SolverStep( real64 const & time_n,
                                     real64 const & dt,
                                     int const cycleNumber,
                                     DomainPartition & domain )
{
  return LinearImplicitStep( time_n, dt, cycleNumber, domain );
}

void ElasticStaticFEM::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                          real64 const & GEOSX_UNUSED_PARAM( dt ),
                                          DomainPartition & domain )
{
  SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );
}


void ElasticStaticFEM::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                             real64 const & GEOSX_UNUSED_PARAM( dt ),
                                             DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void ElasticStaticFEM::SetupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                  DofManager & dofManager ) const
{
  dofManager.addField( keys::TotalDisplacement,
                       DofManager::Location::Node,
                       3);

  dofManager.addCoupling( keys::TotalDisplacement,
                          keys::TotalDisplacement,
                          DofManager::Connector::Elem );
}

void ElasticStaticFEM::SetupSystem( DomainPartition & domain,
                              DofManager & dofManager,
                              CRSMatrix< real64, globalIndex > & localMatrix,
                              array1d< real64 > & localRhs,
                              array1d< real64 > & localSolution,
                              bool const setSparsity )
{
  SolverBase::SetupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );

  MeshLevel * const mesh = domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  NodeManager const * const nodeManager = mesh->getNodeManager();
  arrayView1d< globalIndex const > const &
  dofIndex = nodeManager->getReference< globalIndex_array >( dofManager.getKey( keys::TotalDisplacement ) );

  SparsityPattern< globalIndex > sparsityPattern( dofManager.numLocalDofs(),
                                                  dofManager.numGlobalDofs(),
                                                  8*8*3 );

  finiteElement::fillSparsity< CellElementSubRegion,
                               ElasticStaticFEMKernels >( *mesh,
                                                          targetRegionNames(),
                                                          getDiscretizationName(),
                                                          dofIndex,
                                                          dofManager.rankOffset(),
                                                          sparsityPattern );

  sparsityPattern.compress();
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( sparsityPattern ) );
}

void ElasticStaticFEM::AssembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                       real64 const GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & domain,
                                       DofManager const & dofManager,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs )
{
  MeshLevel * const mesh = domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  NodeManager const & nodeManager = *( mesh->getNodeManager() );

  arrayView1d< globalIndex const > const & 
  dofIndex = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( dataRepository::keys::TotalDisplacement ) );

  finiteElement::regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                               constitutive::SolidBase,
                                               CellElementSubRegion,
                                               ElasticStaticFEMKernels >( *mesh,
                                                                          targetRegionNames(),
                                                                          getDiscretizationName(),
                                                                          m_solidMaterialNames,
                                                                          dofIndex,
                                                                          dofManager.rankOffset(),
                                                                          localMatrix,
                                                                          localRhs);
}

void ElasticStaticFEM::ApplySystemSolution( DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localSolution,
                                            real64 const scalingFactor,
                                            DomainPartition & domain )
{
  dofManager.addVectorToField( localSolution,
                               keys::TotalDisplacement,
                               keys::TotalDisplacement,
                               -scalingFactor );

  ResetStateToBeginningOfStep( domain );
  dofManager.addVectorToField( localSolution,
                               keys::TotalDisplacement,
                               keys::IncrementalDisplacement,
                               -scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( keys::IncrementalDisplacement );
  fieldNames["node"].emplace_back( keys::TotalDisplacement );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors(),
                                         true );
}

void ElasticStaticFEM::ApplyBoundaryConditions( real64 const time_n,
                                                real64 const dt,
                                                DomainPartition & domain,
                                                DofManager const & dofManager,
                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                arrayView1d< real64 > const & localRhs )
{
  ApplyDisplacementBC_implicit( time_n + dt, dofManager, domain, localMatrix, localRhs );
  CRSApplyTractionBC( time_n + dt, dofManager, domain, localRhs );
}

void ElasticStaticFEM::SolveSystem( DofManager const & dofManager,
                                    ParallelMatrix & matrix,
                                    ParallelVector & rhs,
                                    ParallelVector & solution )
{
  solution.zero();
  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}


void ElasticStaticFEM::ApplyDisplacementBC_implicit( real64 const time,
                                                     DofManager const & dofManager,
                                                     DomainPartition & domain,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();

  fsManager.Apply( time,
                   &domain,
                   "nodeManager",
                   keys::TotalDisplacement,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const targetGroup,
                        string const fieldName )
  {
    bc->ApplyBoundaryConditionToSystem< FieldSpecificationEqual,
                                        parallelDevicePolicy< 32 > >( targetSet,
                                                                      time,
                                                                      targetGroup,
                                                                      fieldName,
                                                                      dofManager.getKey( keys::TotalDisplacement ),
                                                                      dofManager.rankOffset(),
                                                                      localMatrix,
                                                                      localRhs );
  } );
}

void ElasticStaticFEM::CRSApplyTractionBC( real64 const time,
                                           DofManager const & dofManager,
                                           DomainPartition & domain,
                                           arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  FaceManager const & faceManager = *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getFaceManager();
  NodeManager const & nodeManager = *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getNodeManager();

  arrayView1d< real64 const > const faceArea  = faceManager.getReference< real64_array >( "faceArea" );
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  arrayView1d< globalIndex const > const 
  blockLocalDofNumber = nodeManager.getReference< globalIndex_array >( dofManager.getKey( keys::TotalDisplacement ) );

  globalIndex const dofRankOffset = dofManager.rankOffset();

  fsManager.Apply( time,
                   &domain,
                   "faceManager",
                   string( "Traction" ),
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const,
                        string const & )
  {
    integer const component = bc->GetComponent();
    real64 tractionStress = bc->GetScale();

    forAll< parallelDevicePolicy< 32 > >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex const kf = targetSet[ i ];
      int numNodes = faceToNodeMap.sizeOfArray( kf );
      real64 nodalForce = tractionStress * faceArea[ kf ] / numNodes;

      for( localIndex a=0; a<numNodes; ++a )
      {
        localIndex const dof = blockLocalDofNumber[ faceToNodeMap( kf, a ) ] + component - dofRankOffset;
        if( dof < 0 || dof >= localRhs.size() )
          continue;
        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[ dof ], nodalForce );
      }
    } );
  } );
}

void ElasticStaticFEM::ResetStateToBeginningOfStep( DomainPartition & domain )
{

  MeshLevel * const mesh = domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodeManager = *( mesh->getNodeManager() );

  arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & incdisp  = nodeManager.incrementalDisplacement();

  forAll< parallelDevicePolicy< 32 > >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    for( localIndex i = 0; i < 3; ++i )
    {
      incdisp( a, i ) = 0.0;
    }
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, ElasticStaticFEM, string const &, dataRepository::Group * const )
} /* namespace geosx */
