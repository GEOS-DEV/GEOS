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

#include "DiffusionFEM.hpp"
#include "DiffusionFEMKernels.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "common/TimingMacros.hpp"
#include "common/DataTypes.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{

using namespace dataRepository;

DiffusionFEM::DiffusionFEM( std::string const & name,
                            Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_diffusion( 0.0 )
{
  registerWrapper( viewKeyStruct::fieldNameString, &m_fieldName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the field variable" );

  registerWrapper( viewKeyStruct::diffusionString, &m_diffusion )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Diffusion coefficient" );
}

void DiffusionFEM::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    NodeManager * const nodes = mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getNodeManager();

    nodes->registerWrapper< real64_array >( m_fieldName )->
      setApplyDefaultValue( 0.0 )->
      setPlotLevel( PlotLevel::LEVEL_0 )->
      setDescription( "Primary field variable" );

    
    m_dFieldName = "delta_" + m_fieldName;

    nodes->registerWrapper< real64_array >( m_dFieldName )->
      setApplyDefaultValue( 0.0 )->
      setPlotLevel( PlotLevel::LEVEL_0 )->
      setDescription( "Change at the last time step of the primary field" );
  }
}

void DiffusionFEM::AssembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                   real64 const dt,
                                   DomainPartition & domain,
                                   DofManager const & dofManager,
                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                   arrayView1d< real64 > const & localRhs )
{
  MeshLevel * const mesh = domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  NodeManager & nodeManager = *(mesh->getNodeManager());

  arrayView1d< globalIndex const > const &
  dofIndex =  nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );


  finiteElement::
    regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                  constitutive::NullModel,
                                  CellElementSubRegion,
                                  DiffusionFEMKernel >( *mesh,
                                                        targetRegionNames(),
                                                        this->getDiscretizationName(),
                                                        arrayView1d< string const >(),
                                                        dofIndex,
                                                        dofManager.
                                                        rankOffset(),
                                                        localMatrix,
                                                        localRhs,
                                                        dt,
                                                        m_fieldName,
                                                        m_dFieldName,
                                                        m_diffusion );
}

void DiffusionFEM::ApplySystemSolution( DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localSolution,
                                        real64 const scalingFactor,
                                        DomainPartition & domain )
{
  // Apply the system solution to update the primary field.
  dofManager.addVectorToField( localSolution,
                               m_fieldName,
                               m_fieldName,
                               scalingFactor );


  // Apply the system solution to the primary field variation at the acctual time step.
  ResetStateToBeginningOfStep( domain );

  dofManager.addVectorToField( localSolution,
                               m_fieldName,
                               m_dFieldName,
                               scalingFactor );
  
  // Synchronize ghost nodes.
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( m_fieldName );
  fieldNames["node"].emplace_back( m_dFieldName );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors(),
                                         true );
}

void DiffusionFEM::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodeManager = *mesh.getNodeManager();
  arrayView1d< real64 > const & dPrimaryField = nodeManager.template getReference< array1d< real64 > >( m_dFieldName );
  forAll< parallelDevicePolicy< 32 > >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    dPrimaryField( a ) = 0.0;
  } );
}


/*
* Remainings are some generic stuffs
*/

DiffusionFEM::~DiffusionFEM(){}

real64 DiffusionFEM::SolverStep( real64 const & time_n,
                                 real64 const & dt,
                                 const int cycleNumber,
                                 DomainPartition & domain )
{
  real64 dtReturn = dt;
  dtReturn = this->LinearImplicitStep( time_n, dt, cycleNumber, domain );
  return dtReturn;
}


void DiffusionFEM::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                      real64 const & GEOSX_UNUSED_PARAM( dt ),
                                      DomainPartition & domain )
{
  // Computation of the sparsity pattern
  SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );
}

void DiffusionFEM::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                         real64 const & GEOSX_UNUSED_PARAM( dt ),
                                         DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void DiffusionFEM::SetupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                              DofManager & dofManager ) const
{
  dofManager.addField( m_fieldName,
                       DofManager::Location::Node );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );
}

void DiffusionFEM::SetupSystem( DomainPartition & domain,
                                DofManager & dofManager,
                                CRSMatrix< real64, globalIndex > & localMatrix,
                                array1d< real64 > & localRhs,
                                array1d< real64 > & localSolution,
                                bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  SolverBase::SetupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );

  MeshLevel * const mesh = domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  NodeManager const * const nodeManager = mesh->getNodeManager();
  arrayView1d< globalIndex const > const &
  dofIndex = nodeManager->getReference< globalIndex_array >( dofManager.getKey( m_fieldName ) );

  SparsityPattern< globalIndex > sparsityPattern( dofManager.numLocalDofs(),
                                                  dofManager.numGlobalDofs(),
                                                  8*8*3 );

  finiteElement::fillSparsity< CellElementSubRegion,
                               DiffusionFEMKernel >( *mesh,
                                                   targetRegionNames(),
                                                   this->getDiscretizationName(),
                                                   dofIndex,
                                                   dofManager.rankOffset(),
                                                   sparsityPattern );

  sparsityPattern.compress();
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( sparsityPattern ) );
}

void DiffusionFEM::ApplyBoundaryConditions( real64 const time_n,
                                            real64 const dt,
                                            DomainPartition & domain,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
{
  ApplyDirichletBC_implicit( time_n + dt, dofManager, domain, localMatrix, localRhs );
}

void DiffusionFEM::SolveSystem( DofManager const & dofManager,
                                ParallelMatrix & matrix,
                                ParallelVector & rhs,
                                ParallelVector & solution )
{
  rhs.scale( -1.0 );
  solution.zero();
  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void DiffusionFEM::ApplyDirichletBC_implicit( real64 const time,
                                              DofManager const & dofManager,
                                              DomainPartition & domain,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();

  fsManager.Apply( time,
                   &domain,
                   "nodeManager",
                   m_fieldName,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const targetGroup,
                        string const & GEOSX_UNUSED_PARAM( fieldName ) )
  {
    bc->ApplyBoundaryConditionToSystem< FieldSpecificationEqual, parallelDevicePolicy< 32 > >( targetSet,
                                                                                               time,
                                                                                               targetGroup,
                                                                                               m_fieldName,
                                                                                               dofManager.getKey( m_fieldName ),
                                                                                               dofManager.rankOffset(),
                                                                                               localMatrix,
                                                                                               localRhs );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, DiffusionFEM, std::string const &, Group * const )

} /* namespace geosx */
