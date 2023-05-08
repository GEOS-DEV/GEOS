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
 * @file MsrsbUtils.hpp
 *
 * Contains various utility functions used by MsRSB level builders.
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBUTILS_HPP_
#define GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBUTILS_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/multiscale/mesh/MeshObjectManager.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/TransposeOperator.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{
namespace multiscale
{

class DofManager;
class MeshObjectManager;
class MeshLevel;

namespace msrsb
{

/**
 * @brief Construct a local adjacency map between points (e.g. vertices) connected by duals (e.g. cells).
 * @param nodeManager object manager for points (nodes)
 * @param ghostNodes whether to include ghosted points in the map
 * @param dualManager object manager for duals (cells)
 * @param ghostDuals whether to include connections via ghosted duals
 * @return the constructed map
 */
ArrayOfSets< localIndex >
buildLocalConnectivity( MeshObjectManager const & nodeManager,
                        bool ghostNodes,
                        MeshObjectManager const & dualManager,
                        bool ghostDuals );

/**
 * @brief Construct a non-overlapping partitioning (clustering) of points using point connectivity and a set of partition seeds.
 * @param connectivity the point connectivity graph
 * @param seeds a list of seed point indices, corresponding to the number of desired partitions
 * @param supports support regions for each partition, i.e. a set of points that the partition is constrained to
 *                 (if empty, partition growth is unconstrained)
 * @return an array of partition indices for each point (of size connectivity.size())
 */
array1d< localIndex >
makeSeededPartition( ArrayOfSetsView< localIndex const > const & connectivity,
                     arrayView1d< localIndex const > const & seeds,
                     ArrayOfSetsView< localIndex const > const & supports );

/**
 * @brief Build a map of support regions through subdomain adjacency.
 * @param fineObjectToSubdomain map of fine-scale objects (points) to adjacent subdomains
 *        ("virtual" boundary subdomains, if included, must have negative indices)
 * @param subdomainToCoarseObject map of subdomains to adjacent coarse-scale objects (points).
 *        Does not need to include boundary subdomains.
 * @param coarseObjectToSubdomain map of coarse-scale objects (points) to adjacent subdomains
 *        (generally a transpose of @p subdomainToCoarseObject)
 * @return a map of fine objects (points) to a list of coarse objects (points) to the support of which they belong
 */
ArrayOfSets< localIndex >
buildSupports( ArrayOfSetsView< localIndex const > const & fineObjectToSubdomain,
               ArrayOfSetsView< localIndex const > const & subdomainToCoarseObject,
               ArrayOfSetsView< localIndex const > const & coarseObjectToSubdomain );

/**
 * @brief Compute an array of global support boundary indicators.
 * @param fineObjectToSubdomain an adjacency map of fine mesh support points to coarse subdomains
 * @return an array of size fineObjectToSubdomain.size(), with 1s indicating boundary points
 */
array1d< integer >
findGlobalSupportBoundary( ArrayOfSetsView< localIndex const > const & fineObjectToSubdomain );

/**
 * @brief Compute the sparsity pattern of the prolongation operator.
 * @param fineDofManager fine level DofManager
 * @param coarseDofManager coarse level DofManager
 * @param fieldName dof field name to extract info from DofManagers
 * @param supports support map (for each fine point, a list of coarse points that interpolate it)
 * @return the local sparsity pattern
 */
SparsityPattern< globalIndex >
buildProlongationSparsity( DofManager const & fineDofManager,
                           DofManager const & coarseDofManager,
                           string const & fieldName,
                           ArrayOfSetsView< localIndex const > const & supports );

/**
 * @brief Initialize a tentative prolongation operator according to a non-overlapping partitioning of points.
 * @param fineDofManager fine level DofManager
 * @param coarseDofManager coarse level DofManager
 * @param fieldName field name for the dof registered in DofManager
 * @param supports support map (for each fine point, a list of coarse points that interpolate it)
 * @param initPart non-overlapping partition of fine points that defines the initial interpolation
 * @return local part of the prolongation operator
 */
CRSMatrix< real64, globalIndex >
buildTentativeProlongation( DofManager const & fineDofManager,
                            DofManager const & coarseDofManager,
                            string const & fieldName,
                            ArrayOfSetsView< localIndex const > const & supports,
                            arrayView1d< localIndex const > const & initPart );

/**
 * @brief Build lists of boundary and interior dof indices to use in MsRSB basis smoothing iteration.
 * @param dofManager the current level DofManager
 * @param fieldName field name for the dof registered in DofManager
 * @param indicator indicator array with 1s indicating global support boundary points
 * @param boundaryDof an output array to be populated with a list of locally present boundary dof indices
 * @param interiorDof an output array to be populated with a list of locally present interior dof indices
 */
void makeGlobalDofLists( DofManager const & dofManager,
                         string const & fieldName,
                         arrayView1d< integer const > const & indicator,
                         array1d< globalIndex > & boundaryDof,
                         array1d< globalIndex > & interiorDof );

/**
 * @brief Build (possibly overlapping) supports consisting of a fixed number of layers added to initial (non-overlapping) partitions
 * @param numLayers number of layers to add
 * @param connectivity point connectivity map
 * @param initialPartition initial non-overlapping partition of the points
 * @return support map (for each fine point, a list of coarse points that interpolate it)
 */
ArrayOfSets< localIndex >
buildLayeredSupport( integer numLayers,
                     ArrayOfSetsView< localIndex const > const & connectivity,
                     arrayView1d< localIndex const > const & initialPartition );

/**
 * @brief Given a layered support, find the global support boundary
 * @param connectivity point connectivity map
 * @param support support map (for each fine point, a list of coarse points that interpolate it)
 * @return an array of size connectivity.size() with 1s indicating boundary points
 */
array1d< integer >
findLayeredSupportBoundary( ArrayOfSetsView< localIndex const > const & connectivity,
                            ArrayOfSetsView< localIndex const > const & support );

/**
 * @brief Make a restriction operator
 * @tparam MATRIX type of sparse matrix
 * @param params multiscale parameters
 * @param prolongation the prolongation operator
 * @return the linear operator representing restriction
 *
 * The operator returned may be an explicit matrix or an implicit operator wrapping @p prolongation, which must have appropriate lifetime.
 */
template< typename MATRIX >
std::unique_ptr< LinearOperator< typename MATRIX::Vector > >
makeRestriction( LinearSolverParameters::Multiscale const & params,
                 MATRIX const & prolongation )
{
  std::unique_ptr< LinearOperator< typename MATRIX::Vector > > restriction;
  if( params.galerkin )
  {
    // Make a transpose operator with a reference to P, which will be computed later
    restriction = std::make_unique< TransposeOperator< MATRIX > >( prolongation );
  }
  else
  {
    // Make an explicit transpose of tentative (initial) P
    MATRIX R;
    prolongation.transpose( R );
    restriction = std::make_unique< MATRIX >( std::move( R ) );
  }
  return restriction;
}

/**
 * @brief Drop entries from the matrix that are below @p relTol times max abs value in current row.
 * @param mat tne input matrix view
 * @param relTol relative tolerance parameter
 * @return the new matrix with remaining entries and compressed sparsity pattern
 */
CRSMatrix< real64, globalIndex >
dropEntries( CRSMatrixView< real64 const, globalIndex const > const & mat,
             real64 relTol );

/**
 * @brief Compute the triple product R*A*P with optional post-filtering
 * @tparam MATRIX type of sparse matrix
 * @param params multiscale parameters
 * @param fineMatrix the fine-scale system matrix
 * @param prolongation the prolongation operator
 * @param restrictionOp the restriction operator
 * @return the coarse-level matrix
 */
template< typename MATRIX >
MATRIX
computeRAP( LinearSolverParameters::Multiscale const & params,
            MATRIX const & fineMatrix,
            MATRIX const & prolongation,
            LinearOperator< typename MATRIX::Vector > const & restrictionOp )
{
  MATRIX coarseMatrix;
  if( params.galerkin )
  {
    fineMatrix.multiplyPtAP( prolongation, coarseMatrix );
  }
  else
  {
    MATRIX const & restriction = dynamicCast< MATRIX const & >( restrictionOp );
    fineMatrix.multiplyRAP( restriction, prolongation, coarseMatrix );
  }

  if( params.droptol > 0.0 )
  {
    CRSMatrix< real64, globalIndex > mat = coarseMatrix.extract();
    mat = dropEntries( mat.toViewConst(), params.droptol );
    coarseMatrix.create( mat.toViewConst(), coarseMatrix.numLocalCols(), coarseMatrix.comm() );
  }
  return coarseMatrix;
}

/**
 * @brief Write the prolongation operator as a collection of fine-level fields (one per coarse point)
 * @param prolongation the local prolongation operator
 * @param dofManager the fine-scale DofManager
 * @param fieldName the dof field name registered in DofManager
 * @param prefix string prefix to use for all field names
 * @param mesh the fine-scale mesh
 * @param fineManager the fine-scale object manager to write fields into
 * @param writeFunc the callback that will be triggered once all fields have been registered on @p mesh
 *                  and written to and before they are remove, and passed the list of field names (keys).
 *                  This can be used to propagate data further to finer level (or the original GEOS mesh).
 *
 * This function should only be used for debugging purposes (e.g. to dump basis functions to 3d plot outputs).
 * It is very slow as it has to write a massive amount of data.
 */
void writeProlongation( CRSMatrixView< real64 const, globalIndex const > const & prolongation,
                        multiscale::DofManager const & dofManager,
                        string const & fieldName,
                        string const & prefix,
                        multiscale::MeshLevel & mesh,
                        multiscale::MeshObjectManager & fineManager,
                        std::function< void ( multiscale::MeshLevel &, std::vector< string > const & ) > const & writeFunc );

} // namespace msrsb
} // namespace multiscale
} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBUTILS_HPP_
