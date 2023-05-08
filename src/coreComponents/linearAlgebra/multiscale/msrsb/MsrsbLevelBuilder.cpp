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
 * @file MsrsbLevelBuilder.cpp
 */

#include "MsrsbLevelBuilder.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/multiscale/mesh/MeshData.hpp"
#include "linearAlgebra/multiscale/mesh/MeshUtils.hpp"
#include "linearAlgebra/multiscale/msrsb/MsrsbUtils.hpp"
#include "linearAlgebra/solvers/PreconditionerNull.hpp"
#include "linearAlgebra/utilities/TransposeOperator.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{
namespace multiscale
{

template< typename LAI >
MsrsbLevelBuilder< LAI >::MsrsbLevelBuilder( string name, LinearSolverParameters params )
  : Base( std::move( name ), std::move( params ) ),
  m_mesh( m_name ),
  m_updateLag( m_params.multiscale.msrsb.updateFrequency )
{
  //GEOS_ASSERT( m_params.subParams.empty() );
}

template< typename LAI >
void MsrsbLevelBuilder< LAI >::createSmoothers()
{
  LinearSolverParameters smoother_params = m_params;
  smoother_params.preconditionerType = m_params.multiscale.smoother.type;
  smoother_params.dofsPerNode = numComp();

  using PreOrPost = LinearSolverParameters::AMG::PreOrPost;
  PreOrPost const & preOrPost = m_params.multiscale.smoother.preOrPost;

  m_preSmoother = preOrPost == PreOrPost::pre || preOrPost == PreOrPost::both
                ? LAI::createPreconditioner( smoother_params )
                : std::make_unique< PreconditionerNull< LAI > >();

  m_postSmoother = preOrPost == PreOrPost::post || preOrPost == PreOrPost::both
                 ? LAI::createPreconditioner( smoother_params )
                 : std::make_unique< PreconditionerNull< LAI > >();

  // TODO: pre/post smoother could be the same object
}

template< typename LAI >
void MsrsbLevelBuilder< LAI >::initializeFineLevel( DomainPartition & domain,
                                                    geos::DofManager const & dofManager,
                                                    MPI_Comm const & comm )
{
  GEOS_MARK_FUNCTION;

  m_mesh.buildFineMesh( domain, dofManager.support( m_params.multiscale.fieldName ) );

  m_location = dofManager.location( m_params.multiscale.fieldName );
  m_dofManager.setDomain( domain );
  m_dofManager.addField( m_params.multiscale.fieldName, dofManager.numComponents( m_params.multiscale.fieldName ), manager() );
  m_dofManager.reorderByRank();

  // Create a "fake" fine matrix (no data, just correct sizes/comms for use at coarse level init)
  m_matrix.createWithLocalSize( m_dofManager.numLocalDofs(), m_dofManager.numLocalDofs(), 0, comm );

  createSmoothers();
}

struct SupportDescription
{
  ArrayOfSets< localIndex > supports;
  array1d< integer > supportBoundaryIndicator;
  array1d< localIndex > initialPartition;
};

/**
 * @brief Build the basic sparsity pattern for prolongation.
 *
 * Support of a coarse nodal basis function is defined as the set of fine-scale nodes
 * that are adjacent exclusively to subdomains (coarse cells or boundaries) that are
 * also adjacent to that coarse node.
 */
static SupportDescription
buildMatchingNodalSupport( multiscale::MeshLevel const & fine,
                           multiscale::MeshLevel const & coarse,
                           arrayView1d< string const > const & boundaryNodeSets )
{
  GEOS_MARK_FUNCTION;

  ArrayOfSets< localIndex > const nodalConn =
    msrsb::buildLocalConnectivity( fine.nodeManager(), true, fine.cellManager(), false );
  ArrayOfSets< localIndex > const fineNodeToCellSubdomain =
    meshUtils::buildFineObjectToSubdomainMap( fine.nodeManager(),
                                              fine.cellManager().getField< fields::multiscale::CoarseCellLocalIndex >(),
                                              boundaryNodeSets );
  ArrayOfSets< localIndex > const coarseNodeToCellSubdomain =
    meshUtils::addBoundarySubdomains( coarse.nodeManager(),
                                      coarse.nodeManager().toDualRelation(),
                                      boundaryNodeSets );

  SupportDescription result;

  result.supports =
    msrsb::buildSupports( fineNodeToCellSubdomain.toViewConst(),
                          coarse.cellManager().toDualRelation().toViewConst(),
                          coarseNodeToCellSubdomain.toViewConst() );

  result.supportBoundaryIndicator =
    msrsb::findGlobalSupportBoundary( fineNodeToCellSubdomain.toViewConst() );

  result.initialPartition =
    msrsb::makeSeededPartition( nodalConn.toViewConst(),
                                coarse.nodeManager().getField< fields::multiscale::FineNodeLocalIndex >().toViewConst(),
                                result.supports.toViewConst() );

  return result;
}

/**
 * @brief Build the basic sparsity pattern for prolongation.
 *
 * Support of a coarse nodal basis function is defined as the set of fine-scale nodes
 * that are adjacent exclusively to subdomains (coarse cells or boundaries) that are
 * also adjacent to that coarse node.
 */
static SupportDescription
buildMatchingCellSupport( multiscale::MeshLevel const & fine,
                          multiscale::MeshLevel const & coarse,
                          arrayView1d< string const > const & boundaryNodeSets )
{
  GEOS_MARK_FUNCTION;

  // Build the nodal partition to act as dual "volumes"
  array1d< localIndex > const nodalPartition = [&]
  {
    ArrayOfSets< localIndex > const nodalConn =
      msrsb::buildLocalConnectivity( fine.nodeManager(), true, fine.cellManager(), false );
    ArrayOfSets< localIndex > const fineNodeToCellSubdomain =
      meshUtils::buildFineObjectToSubdomainMap( fine.nodeManager(),
                                                fine.cellManager().getField< fields::multiscale::CoarseCellLocalIndex >(),
                                                boundaryNodeSets );
    ArrayOfSets< localIndex > const coarseNodeToCellSubdomain =
      meshUtils::addBoundarySubdomains( coarse.nodeManager(),
                                        coarse.nodeManager().toDualRelation(),
                                        boundaryNodeSets );

    // Unfortunately, we need nodal supports (in order to limit nodal partition growth).
    ArrayOfSets< localIndex > const nodalSupports =
      msrsb::buildSupports( fineNodeToCellSubdomain.toViewConst(),
                            coarse.cellManager().toDualRelation(),
                            coarseNodeToCellSubdomain.toViewConst() );

    return msrsb::makeSeededPartition( nodalConn.toViewConst(),
                                       coarse.nodeManager().getField< fields::multiscale::FineNodeLocalIndex >(),
                                       nodalSupports.toViewConst() );
  }();

  ArrayOfSets< localIndex > const fineCellToNodalSubdomain =
    meshUtils::buildFineObjectToSubdomainMap( fine.cellManager(), nodalPartition.toViewConst(), {} );

  SupportDescription result;

  result.supports =
    msrsb::buildSupports( fineCellToNodalSubdomain.toViewConst(),
                          coarse.nodeManager().toDualRelation().toViewConst(),
                          coarse.cellManager().toDualRelation().toViewConst() );

  result.supportBoundaryIndicator =
    msrsb::findGlobalSupportBoundary( fineCellToNodalSubdomain.toViewConst() );

  result.initialPartition.resize( fine.cellManager().size() );
  result.initialPartition.setValues< parallelHostPolicy >( fine.cellManager().getField< fields::multiscale::CoarseCellLocalIndex >() );

  return result;
}

static SupportDescription
buildLayeredCellSupport( multiscale::MeshLevel const & fine,
                         multiscale::MeshLevel const & coarse,
                         integer const numLayers )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( coarse );

  ArrayOfSets< localIndex > const cellConn =
    msrsb::buildLocalConnectivity( fine.cellManager(), true, fine.nodeManager(), true );

  auto const coarseCells = fine.cellManager().getField< fields::multiscale::CoarseCellLocalIndex >().toViewConst();

  SupportDescription result;
  result.initialPartition.resize( fine.cellManager().size() );
  result.initialPartition.setValues< parallelHostPolicy >( coarseCells );
  result.supports = msrsb::buildLayeredSupport( numLayers, cellConn.toViewConst(), result.initialPartition.toViewConst() );
  result.supportBoundaryIndicator = msrsb::findLayeredSupportBoundary( cellConn.toViewConst(), result.supports.toViewConst() );

  return result;
}

static SupportDescription
buildLayeredNodalSupport( multiscale::MeshLevel const & fine,
                          multiscale::MeshLevel const & coarse,
                          integer const numLayers )
{
  GEOS_MARK_FUNCTION;

  ArrayOfSets< localIndex > const nodalConn =
    msrsb::buildLocalConnectivity( fine.nodeManager(), true, fine.cellManager(), true );

  auto const coarseNodes = coarse.nodeManager().getField< fields::multiscale::FineNodeLocalIndex >().toViewConst();

  SupportDescription result;
  result.initialPartition = msrsb::makeSeededPartition( nodalConn.toViewConst(), coarseNodes, {} );
  result.supports = msrsb::buildLayeredSupport( numLayers, nodalConn.toViewConst(), result.initialPartition.toViewConst() );
  result.supportBoundaryIndicator = msrsb::findLayeredSupportBoundary( nodalConn.toViewConst(), result.supports.toViewConst() );
  return result;
}

static SupportDescription
buildSupports( FieldLocation const loc,
               multiscale::MeshLevel const & fine,
               multiscale::MeshLevel const & coarse,
               LinearSolverParameters::Multiscale const & params )
{
  switch( params.msrsb.support )
  {
    case LinearSolverParameters::Multiscale::MsRSB::SupportType::matching:
    {
      auto const build = ( loc == FieldLocation::Node ) ? buildMatchingNodalSupport : buildMatchingCellSupport;
      return build( fine, coarse, params.boundarySets );
    }
    case LinearSolverParameters::Multiscale::MsRSB::SupportType::layers:
    {
      auto const build = ( loc == FieldLocation::Node ) ? buildLayeredNodalSupport : buildLayeredCellSupport;
      return build( fine, coarse, params.msrsb.numLayers );
    }
    default:
    {
      GEOS_ERROR( "Unrecognized support region algorithm" );
      return {};
    }
  }
}

template< typename LAI >
void MsrsbLevelBuilder< LAI >::initializeCoarseLevel( LevelBuilderBase< LAI > & fineLevel,
                                                      Matrix const & fineMatrix )
{
  GEOS_MARK_FUNCTION;

  auto & fine = dynamicCast< MsrsbLevelBuilder< LAI > & >( fineLevel );
  m_fineLevel = &fine;

  // Provide a local system matrix as a source for graph weights
  {
    CRSMatrix< real64, localIndex > localMat;
    auto & coarseningParams = m_params.multiscale.coarsening;
    if( coarseningParams.partitionType == LinearSolverParameters::Multiscale::Coarsening::PartitionType::graph
        && coarseningParams.graph.matrixWeights > 0 )
    {
      localMat = fineMatrix.extractLocal();
      coarseningParams.graph.localMatrix = &localMat;
    }

    m_mesh.buildCoarseMesh( fine.mesh(), m_params.multiscale.coarsening, m_params.multiscale.boundarySets );
    coarseningParams.graph.localMatrix = nullptr; // avoid a dangling pointer, just in case
  }

  m_location = fine.m_location;
  m_dofManager.setDomain( fine.dofManager().domain() );
  m_dofManager.addField( m_params.multiscale.fieldName, fine.dofManager().numComponents( m_params.multiscale.fieldName ), manager() );
  m_dofManager.reorderByRank();

  // Write data back to GEOSX for visualization and debug
  if( m_params.multiscale.debugLevel >= 1 )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "[MsRSB] {}: generated coarse grid with {} global cells and {} global nodes",
                               m_name,
                               m_mesh.cellManager().maxGlobalIndex() + 1,
                               m_mesh.nodeManager().maxGlobalIndex() + 1 ) );
  }

  if( m_params.multiscale.debugLevel >= 2 )
  {
    GEOS_LOG_RANK( GEOS_FMT( "[MsRSB] {}: generated coarse grid with {} local cells and {} local nodes",
                             m_name,
                             m_mesh.cellManager().numOwnedObjects(),
                             m_mesh.nodeManager().numOwnedObjects() ) );
  }

  if( m_params.multiscale.debugLevel >= 3 )
  {
    m_mesh.writeCellData( { ObjectManagerBase::viewKeyStruct::ghostRankString() } );
    m_mesh.writeNodeData( { ObjectManagerBase::viewKeyStruct::ghostRankString() } );
    fine.mesh().writeCellData( { fields::multiscale::CoarseCellLocalIndex::key(),
                                 fields::multiscale::CoarseCellGlobalIndex::key() } );
    fine.mesh().writeNodeData( { fields::multiscale::CoarseNodeLocalIndex::key(),
                                 fields::multiscale::CoarseNodeGlobalIndex::key() } );
  }

  // Build support region definitions and tentative partition of unity
  // For now, we only handle two types of basis functions - nodal and cell-centered - with specific algorithms.
  // In future this should be refactored into an extensible hierarchy of basis constructors.
  SupportDescription const desc = buildSupports( m_location,
                                                 fine.mesh(),
                                                 m_mesh,
                                                 m_params.multiscale );

  // Construct global internal/boundary DoF sets
  msrsb::makeGlobalDofLists( fine.dofManager(),
                             m_params.multiscale.fieldName,
                             desc.supportBoundaryIndicator,
                             m_boundaryDof,
                             m_interiorDof );

  // Convert the partitioning into an actual DoF-based local matrix
  CRSMatrix< real64, globalIndex > const localProlongation =
    msrsb::buildTentativeProlongation( fine.dofManager(),
                                       m_dofManager,
                                       m_params.multiscale.fieldName,
                                       desc.supports.toViewConst(),
                                       desc.initialPartition );

  // Assemble local pieces into a global prolongation operator and make restriction operator
  m_prolongation.create( localProlongation.toViewConst(),
                         m_dofManager.numLocalDofs(),
                         fine.matrix().comm() );
  m_restriction = msrsb::makeRestriction( m_params.multiscale, m_prolongation );

  // Create a "fake" coarse matrix (no data, just correct sizes/comms), to be computed later
  m_matrix.createWithLocalSize( m_dofManager.numLocalDofs(), m_dofManager.numLocalDofs(), 0, fine.matrix().comm() );

  createSmoothers();
}

namespace
{

template< typename Matrix >
Matrix filterMatrix( Matrix const & fineMatrix,
                     integer const numComp )
{
  GEOS_MARK_SCOPE( filter );

  // 1. Apply SC approximation
  Matrix filteredMatrix;
  fineMatrix.separateComponentFilter( filteredMatrix, numComp );

  // 2. Filter out positive off-diagonal elements (assumed positive diagonals)
  filteredMatrix.clampEntries( -LvArray::NumericLimits< real64 >::infinity, 0.0, true );

  // 3. Enforce rowsum = 0
  // 3.1. Compute rowsums
  typename Matrix::Vector rowSums;
  rowSums.create( fineMatrix.numLocalRows(), fineMatrix.comm() );
  filteredMatrix.getRowSums( rowSums, RowSumType::SumValues );

  // 3.2. Preserve Dirichlet rows by setting the diagonal update to zero
  typename Matrix::Vector diag;
  diag.create( fineMatrix.numLocalRows(), fineMatrix.comm() );
  filteredMatrix.extractDiagonal( diag );
  forAll< parallelHostPolicy >( diag.localSize(), [diagData = diag.values(),
                                                   rowSumData = rowSums.open()]( localIndex const localRow )
  {
    if( isEqual( diagData[localRow], rowSumData[localRow] ) )
    {
      rowSumData[localRow] = 0.0;
    }
  } );
  rowSums.close();

  // 3.3. Subtract the nonzero rowsums from diagonal elements
  filteredMatrix.addDiagonal( rowSums, -1.0 );

  return filteredMatrix;
}

template< typename MATRIX >
auto makeJacobiMatrix( MATRIX && fineMatrix,
                       real64 const omega )
{
  GEOS_MARK_SCOPE( jacobi );
  using Matrix = std::remove_const_t< TYPEOFREF( fineMatrix ) >;

  // 0. Copy or move input matrix into a new object
  Matrix iterMatrix( std::forward< MATRIX >( fineMatrix ) );

  // 1. Compute -w * Dinv * A;
  typename Matrix::Vector diag;
  diag.create( iterMatrix.numLocalRows(), iterMatrix.comm() );
  iterMatrix.extractDiagonal( diag );
  diag.reciprocal();
  diag.scale( -omega );
  iterMatrix.leftScale( diag );

  // 2. Compute I - w * Dinv * A by adding identity diagonal
  diag.set( 1.0 );
  iterMatrix.addDiagonal( diag, 1.0 );
  return iterMatrix;
}

template< typename Matrix >
Matrix makeIterationMatrix( Matrix const & fineMatrix,
                            integer const numComp,
                            real64 const omega,
                            integer const debugLevel,
                            string const & debugPrefix )
{
  GEOS_MARK_SCOPE( build iteration matrix );

  Matrix filteredMatrix = filterMatrix( fineMatrix, numComp );
  if( debugLevel >= 4 )
  {
    filteredMatrix.write( debugPrefix + "_filtered.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  Matrix jacobiMatrix = makeJacobiMatrix( std::move( filteredMatrix ), omega );
  if( debugLevel >= 4 )
  {
    jacobiMatrix.write( debugPrefix + "_jacobi.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  return jacobiMatrix;
}

template< typename Matrix >
integer iterateBasis( Matrix const & jacobiMatrix,
                      arrayView1d< globalIndex const > const & boundaryDof,
                      arrayView1d< globalIndex const > const & interiorDof,
                      integer const maxIter,
                      real64 const tolerance,
                      integer const checkFreq,
                      integer const debugLevel,
                      string const & name,
                      Matrix & prolongation )
{
  GEOS_MARK_SCOPE( iterate basis );

  auto const saveForDebug = [&]( string const & suffix, integer const minDebugLevel )
  {
    if( debugLevel >= minDebugLevel )
    {
      GEOS_MARK_SCOPE( writeProlongationMatrix );
      prolongation.write( GEOS_FMT( "{}_P_{}.mtx", name, suffix ), LAIOutputFormat::MATRIX_MARKET );
    }
  };

  Matrix P( prolongation );
  integer iter = 0;
  real64 norm = LvArray::NumericLimits< real64 >::max;

  auto const computeAndLogConvergenceNorm = [&]()
  {
    GEOS_MARK_SCOPE( check );
    P.addEntries( prolongation, MatrixPatternOp::Equal, -1.0 );
    norm = P.normMax( interiorDof );
    GEOS_LOG_RANK_0_IF( debugLevel >= 2, GEOS_FMT( "[MsRSB] {}: iter = {}, conv = {:e}", name, iter, norm ) );
  };

  saveForDebug( "init", 4 );
  while( iter < maxIter && norm > tolerance )
  {
    // Keep 1-based iteration index for convenience
    ++iter;

    // Perform a step of Jacobi
    Matrix Ptemp;
    {
      GEOS_MARK_SCOPE( multiply );
      jacobiMatrix.multiply( prolongation, Ptemp );
    }

    // Restrict to the predefined prolongation pattern
    {
      GEOS_MARK_SCOPE( restrict );
      P.zero();
      P.addEntries( Ptemp, MatrixPatternOp::Restrict, 1.0 );
    }

    // Rescale to preserve partition of unity
    {
      GEOS_MARK_SCOPE( rescale );
      P.rescaleRows( boundaryDof, RowSumType::SumValues );
    }

    // Switch over to new prolongation operator
    std::swap( P, prolongation );
    saveForDebug( std::to_string( iter ), 6 );

    // Compute update norm, check convergence
    if( iter % checkFreq == 0 )
    {
      computeAndLogConvergenceNorm();
    }
  }

  // Compute update norm and check convergence one final time if needed (in case we ran out of iterations)
  if( iter % checkFreq != 0 )
  {
    computeAndLogConvergenceNorm();
  }

  GEOS_LOG_RANK_0_IF( debugLevel >= 1,
                      GEOS_FMT( "[MsRSB] {}: {} in {} iterations",
                                name, norm <= tolerance ? "converged" : "failed to converge", iter ) );

  saveForDebug( "conv", 4 );
  return iter;
}

} // namespace

template< typename LAI >
bool MsrsbLevelBuilder< LAI >::updateProlongation( Matrix const & fineMatrix )
{
  GEOS_MARK_FUNCTION;

  GEOS_LOG_RANK_0_IF( m_params.multiscale.debugLevel >= 2, GEOS_FMT( "[MsRSB] {}: building iteration matrix", m_name ) );
  Matrix const jacobiMatrix = makeIterationMatrix( fineMatrix,
                                                   numComp(),
                                                   m_params.multiscale.msrsb.relaxation,
                                                   m_params.multiscale.debugLevel,
                                                   m_name );

  GEOS_LOG_RANK_0_IF( m_params.multiscale.debugLevel >= 2, GEOS_FMT( "[MsRSB] {}: performing basis iteration", m_name ) );
  m_lastNumIter = iterateBasis( jacobiMatrix,
                                m_boundaryDof,
                                m_interiorDof,
                                m_params.multiscale.msrsb.maxIter,
                                m_params.multiscale.msrsb.tolerance,
                                m_lastNumIter <= m_params.multiscale.msrsb.checkFrequency ? 1 : m_params.multiscale.msrsb.checkFrequency,
                                m_params.multiscale.debugLevel,
                                m_name,
                                m_prolongation );

  if( m_params.multiscale.debugLevel >= 5 )
  {
    writeProlongationForDebug();
  }

  // Track the total number of iterations since last "update".
  // If exceeds limit, signal to the caller it's time to "update"
  // (which can mean recomputing RAP, smoothers, etc.)
  m_updateLag += m_lastNumIter;
  if( m_updateLag >= m_params.multiscale.msrsb.updateFrequency )
  {
    m_updateLag = 0;
    return true;
  }
  return false;
}

template< typename LAI >
void MsrsbLevelBuilder< LAI >::writeProlongationForDebug() const
{
  GEOS_MARK_FUNCTION;

  CRSMatrix< real64, globalIndex > const prolongation = m_prolongation.extract();
  MeshObjectManager & manager = m_location == FieldLocation::Node
                              ? m_mesh.fineMesh()->nodeManager()
                              : m_mesh.fineMesh()->cellManager();
  auto const writeFunc = [location = m_location]( multiscale::MeshLevel & mesh, std::vector< string > const & names )
  {
    return location == FieldLocation::Node
    ? mesh.writeNodeData( names )
    : mesh.writeCellData( names );
  };

  msrsb::writeProlongation( prolongation.toViewConst(),
                            m_fineLevel->dofManager(),
                            m_params.multiscale.fieldName,
                            m_name,
                            *m_mesh.fineMesh(),
                            manager,
                            writeFunc );
}

template< typename LAI >
std::unique_ptr< PreconditionerBase< LAI > >
MsrsbLevelBuilder< LAI >::makeCoarseSolver() const
{
  LinearSolverParameters params = m_params;
  if( params.multiscale.coarseType == LinearSolverParameters::PreconditionerType::direct )
  {
    params.solverType = LinearSolverParameters::SolverType::direct;
    return LAI::createSolver( params );
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
template class MsrsbLevelBuilder< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MsrsbLevelBuilder< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MsrsbLevelBuilder< PetscInterface >;
#endif

} // namespace multiscale
} // namespace geos
