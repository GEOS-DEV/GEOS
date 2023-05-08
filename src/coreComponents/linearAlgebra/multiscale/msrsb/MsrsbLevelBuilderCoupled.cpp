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
 * @file MsrsbLevelBuilderCoupled.cpp
 */

#include "MsrsbLevelBuilderCoupled.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/multiscale/msrsb/MsrsbUtils.hpp"
#include "linearAlgebra/solvers/PreconditionerNull.hpp"

namespace geos
{

namespace multiscale
{

template< typename LAI >
MsrsbLevelBuilderCoupled< LAI >::MsrsbLevelBuilderCoupled( string name, LinearSolverParameters params )
  : Base( std::move( name ), std::move( params ) )
{
  GEOS_ASSERT( !m_params.block.subParams.empty() );
  for( LinearSolverParameters const * const p : m_params.block.subParams )
  {
    GEOS_ASSERT_MSG( p != nullptr, "Sub-preconditioner parameters for each block must be set by the solver" );
    string levelName = GEOS_FMT( "{}_{}", m_name, p->multiscale.label );
    m_builders.emplace_back( std::make_unique< MsrsbLevelBuilder< LAI > >( std::move( levelName ), *p ) );
    m_fields.push_back( p->multiscale.fieldName );
  }
  m_prolongationBlocks.resize( m_builders.size() );
  m_selectors.resize( m_builders.size() );
}

template< typename LAI >
void MsrsbLevelBuilderCoupled< LAI >::createSmoothers()
{
  auto const makeSmoother = [&]( auto const getSmoother ) -> std::unique_ptr< PreconditionerBase< LAI > >
  {
    if( m_params.multiscale.coupled.useBlockSmoother )
    {
      integer const numBlocks = m_params.block.subParams.size();
      GEOS_ERROR_IF_NE_MSG( numBlocks, 2, "Only 2 blocks are supported" );
      GEOS_ERROR_IF_NE_MSG( m_params.block.order.size(), numBlocks, "The order of blocks must be prescribed" );
      auto smoother = std::make_unique< BlockPreconditioner< LAI > >( m_params.block );
      for( integer i = 0; i < numBlocks; ++i )
      {
        LinearSolverParameters params;
        params.preconditionerType = m_params.block.subParams[i]->multiscale.smoother.type;
        geos::DofManager::SubComponent comp{ m_params.block.subParams[i]->multiscale.fieldName, { m_builders[i]->numComp(), true } };
        smoother->setupBlock( m_params.block.order[i], { comp }, getSmoother( *m_builders[i] ) );
        smoother->setProlongation( m_params.block.order[i], m_selectors[i] );

      }
      return smoother;
    }
    else
    {
      LinearSolverParameters params;
      params.preconditionerType = m_params.multiscale.smoother.type;
      return LAI::createPreconditioner( params );
    }
  };

  using PreOrPost = LinearSolverParameters::AMG::PreOrPost;
  PreOrPost const & preOrPost = m_params.multiscale.smoother.preOrPost;

  m_preSmoother = preOrPost == PreOrPost::pre || preOrPost == PreOrPost::both
                  ? makeSmoother( []( auto & b ) { return b.presmoother(); } )
                  : std::make_unique< PreconditionerNull< LAI > >();
  m_postSmoother = preOrPost == PreOrPost::post || preOrPost == PreOrPost::both
                   ? makeSmoother( []( auto & b ) { return b.postsmoother(); } )
                   : std::make_unique< PreconditionerNull< LAI > >();
}

template< typename LAI >
void MsrsbLevelBuilderCoupled< LAI >::initializeCommon( DomainPartition & domain,
                                                        MPI_Comm const & comm )
{
  // Populate the coupled multiscale DofManager
  m_dofManager.setDomain( domain );
  for( size_t i = 0; i < m_builders.size(); ++i )
  {
    m_dofManager.addField( m_fields[i], m_builders[i]->numComp(), m_builders[i]->manager() );
  }
  m_dofManager.reorderByRank();

  // Construct subproblem selectors
  for( size_t i = 0; i < m_builders.size(); ++i )
  {
    m_dofManager.makeRestrictor( m_fields[i], comm, true, m_selectors[i] );
  }

  localIndex const numLocalRows =
    std::accumulate( m_builders.begin(), m_builders.end(), localIndex{},
                     []( auto const s, auto const & b ){ return s + b->matrix().numLocalRows(); } );

  // Create a "fake" fine matrix (no data, just correct sizes/comms for use at coarse level init)
  m_matrix.createWithLocalSize( numLocalRows, numLocalRows, 0, comm );
}

template< typename LAI >
void MsrsbLevelBuilderCoupled< LAI >::initializeFineLevel( DomainPartition & domain,
                                                           geos::DofManager const & dofManager,
                                                           MPI_Comm const & comm )
{
  for( size_t i = 0; i < m_builders.size(); ++i )
  {
    m_builders[i]->initializeFineLevel( domain, dofManager, comm );
  }

  initializeCommon( domain, comm );
  createSmoothers();
}

std::unordered_map< globalIndex, globalIndex >
makeGhostDofMap( MeshObjectManager const & manager,
                 string const & oldDofKey,
                 string const & newDofKey )
{
  std::unordered_map< globalIndex, globalIndex > ghostDofMap;
  arrayView1d< globalIndex const > const oldDofNumber = manager.getReference< array1d< globalIndex > >( oldDofKey );
  arrayView1d< globalIndex const > const newDofNumber = manager.getReference< array1d< globalIndex > >( newDofKey );
  for( localIndex i = manager.numOwnedObjects(); i < manager.size(); ++i )
  {
    ghostDofMap[oldDofNumber[i]] = newDofNumber[i];
  }
  return ghostDofMap;
}

template< typename LAI >
void MsrsbLevelBuilderCoupled< LAI >::buildProlongationStructure( DofManager const & fineDofManager )
{
  GEOS_ASSERT_EQ( m_prolongationBlocks.size(), m_dofManager.numFields() );

  array1d< localIndex > rowLength( fineDofManager.numLocalDofs() );
  for( std::size_t blockIdx = 0; blockIdx < m_prolongationBlocks.size(); ++blockIdx )
  {
    localIndex const rowOffset = fineDofManager.localOffset( m_fields[blockIdx] );
    forAll< parallelHostPolicy >( m_prolongationBlocks[blockIdx].numRows(),
                                  [block = m_prolongationBlocks[blockIdx].toViewConst(),
                                   rowLength = rowLength.toView(),
                                   rowOffset]( localIndex const localRow )
    {
      rowLength[localRow + rowOffset] = block.numNonZeros( localRow );
    } );
  }

  m_localProlongation.resizeFromRowCapacities< parallelHostPolicy >( rowLength.size(),
                                                                     m_dofManager.numGlobalDofs(),
                                                                     rowLength.data() );

  for( std::size_t blockIdx = 0; blockIdx < m_prolongationBlocks.size(); ++blockIdx )
  {
    DofManager const & dofManager = m_builders[blockIdx]->dofManager();
    string const & fieldName = m_fields[blockIdx];

    globalIndex const minLocalDof = dofManager.rankOffset();
    globalIndex const maxLocalDof = minLocalDof + dofManager.numLocalDofs();

    globalIndex const colOffset = m_dofManager.globalOffset( fieldName ) - dofManager.globalOffset( fieldName );
    localIndex const rowOffset = fineDofManager.localOffset( fieldName );

    integer const numComp = m_dofManager.numComponents( fieldName );

    std::unordered_map< globalIndex, globalIndex > const ghostDofMap =
      makeGhostDofMap( m_builders[blockIdx]->manager(), dofManager.key( fieldName ), m_dofManager.key( fieldName ) );

    auto const mapGhostCol = [numComp, &ghostDofMap]( globalIndex const col )
    {
      if( numComp == 1 )
      {
        return ghostDofMap.at( col );
      }
      auto const d = std::div( col, static_cast< globalIndex >( numComp ) );
      return ghostDofMap.at( d.quot * numComp ) + d.rem;
    };

    auto const mapColumn = [mapGhostCol, colOffset, minLocalDof, maxLocalDof]( globalIndex const col )
    {
      return ( minLocalDof <= col && col < maxLocalDof ) ? colOffset + col : mapGhostCol( col );
    };

    forAll< parallelHostPolicy >( m_prolongationBlocks[blockIdx].numRows(),
                                  [block = m_prolongationBlocks[blockIdx].toViewConst(),
                                   prolongation = m_localProlongation.toView(),
                                   rowOffset, mapColumn]( localIndex const localRow )
    {
      auto const cols = block.getColumns( localRow );
      auto const vals = block.getEntries( localRow );
      for( localIndex k = 0; k < cols.size(); ++k )
      {
        prolongation.insertNonZero( localRow + rowOffset, mapColumn( cols[k] ), vals[k] );
      }
    } );
  }
}

template< typename LAI >
void MsrsbLevelBuilderCoupled< LAI >::initializeCoarseLevel( LevelBuilderBase< LAI > & fineLevel,
                                                             Matrix const & fineMatrix )
{
  MsrsbLevelBuilderCoupled< LAI > & fine = dynamicCast< MsrsbLevelBuilderCoupled< LAI > & >( fineLevel );
  GEOS_ASSERT( fine.m_builders.size() == m_builders.size() );
  m_fineLevel = &fine;

  for( size_t blockIdx = 0; blockIdx < m_builders.size(); ++blockIdx )
  {
    Matrix fineBlock;
    fineMatrix.multiplyPtAP( fine.m_selectors[blockIdx], fineBlock );
    m_builders[blockIdx]->initializeCoarseLevel( *fine.m_builders[blockIdx], fineBlock );
    m_prolongationBlocks[blockIdx] = m_builders[blockIdx]->prolongation().extract();
  }

  initializeCommon( fine.dofManager().domain(), fineMatrix.comm() );
  createSmoothers();

  buildProlongationStructure( fine.dofManager() );
  m_restriction = msrsb::makeRestriction( m_params.multiscale, m_prolongation );
}

template< typename LAI >
bool MsrsbLevelBuilderCoupled< LAI >::updateProlongation( Matrix const & fineMatrix )
{
  GEOS_MARK_FUNCTION;

  // Extract diagonal blocks, update and extract sub-block prolongations
  bool update = false;
  localIndex rowOffset = 0;
  for( size_t blockIdx = 0; blockIdx < m_builders.size(); ++blockIdx )
  {
    Matrix fineBlock;
    {
      GEOS_MARK_SCOPE( extract blocks );
      auto const & fine = dynamicCast< MsrsbLevelBuilderCoupled< LAI > const & >( *m_fineLevel );
      fineMatrix.multiplyPtAP( fine.m_selectors[blockIdx], fineBlock );
    }

    bool const updateBlock = m_builders[blockIdx]->updateProlongation( fineBlock );
    CRSMatrixView< real64, globalIndex const > const block = m_prolongationBlocks[blockIdx].toViewConstSizes();

    if( updateBlock )
    {
      GEOS_MARK_SCOPE( merge blocks );
      m_builders[blockIdx]->prolongation().extract( block );
      forAll< parallelHostPolicy >( block.numRows(),
                                    [block = block.toViewConst(),
                                     prolongation = m_localProlongation.toViewConstSizes(),
                                     rowOffset]( localIndex const localRow )
      {
        auto const blockVals = block.getEntries( localRow );
        auto const coupledVals = prolongation.getEntries( rowOffset + localRow );
        GEOS_ASSERT_EQ( blockVals.size(), coupledVals.size() );
        std::copy( blockVals.dataIfContiguous(),
                   blockVals.dataIfContiguous() + blockVals.size(),
                   coupledVals.dataIfContiguous() );
      } );
    }
    update = update || updateBlock;
    rowOffset += block.numRows();
  }

  if( update )
  {
    GEOS_MARK_SCOPE( create LAI object );
    GEOS_LOG_RANK_0_IF( m_params.multiscale.debugLevel >= 2, GEOS_FMT( "[MsRSB] {}: updating coupled prolongation", m_name ) );
    m_prolongation.create( m_localProlongation.toViewConst(),
                           m_dofManager.numLocalDofs(),
                           fineMatrix.comm() );
  }

  if( m_params.multiscale.debugLevel >= 4 )
  {
    m_prolongation.write( GEOS_FMT( "{}_P_{}.mtx", m_name, "conv" ), LAIOutputFormat::MATRIX_MARKET );
  }

  return update;
}

template< typename LAI >
std::unique_ptr< PreconditionerBase< LAI > >
MsrsbLevelBuilderCoupled< LAI >::makeCoarseSolver() const
{
  LinearSolverParameters params = m_params;
  if( params.multiscale.coarseType == LinearSolverParameters::PreconditionerType::direct )
  {
    params.solverType = LinearSolverParameters::SolverType::direct;
    return LAI::createSolver( params );
  }
  else if( params.multiscale.coarseType == LinearSolverParameters::PreconditionerType::block )
  {
    GEOS_ERROR_IF_NE_MSG( m_builders.size(), 2, "More than 2 blocks are not supported yet" );
    auto solver = std::make_unique< BlockPreconditioner< LAI > >( params.block );
    for( integer i = 0; i < static_cast< integer >( m_builders.size() ); ++i )
    {
      geos::DofManager::SubComponent comp{ params.block.subParams[i]->multiscale.fieldName, { m_builders[i]->numComp(), true } };
      solver->setupBlock( i, { comp }, m_builders[i]->makeCoarseSolver() );
      solver->setProlongation( i, m_selectors[i] );
    }
    return solver;
  }
  else
  {
    params.preconditionerType = params.multiscale.coarseType;
    return LAI::createPreconditioner( params );
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MsrsbLevelBuilderCoupled< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MsrsbLevelBuilderCoupled< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MsrsbLevelBuilderCoupled< PetscInterface >;
#endif

} // namespace multiscale
} // namespace geos
