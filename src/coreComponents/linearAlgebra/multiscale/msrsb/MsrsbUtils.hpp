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

namespace geosx
{
namespace multiscale
{

class DofManager;
class MeshObjectManager;
class MeshLevel;

namespace msrsb
{

ArrayOfSets< localIndex >
buildLocalConnectivity( MeshObjectManager const & nodeManager,
                        bool ghostNodes,
                        MeshObjectManager const & dualManager,
                        bool ghostDuals );

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
 * @param supportBoundaryIndicator auxiliary output array that marks points on global support boundary with a nonzero value (1).
 * @return a map of fine objects (points) to a list of coarse objects (points) to the support of which they belong
 */
ArrayOfSets< localIndex >
buildSupports( ArrayOfSetsView< localIndex const > const & fineObjectToSubdomain,
               ArrayOfSetsView< localIndex const > const & subdomainToCoarseObject,
               ArrayOfSetsView< localIndex const > const & coarseObjectToSubdomain );

array1d< integer >
findGlobalSupportBoundary( ArrayOfSetsView< localIndex const > const & fineObjectToSubdomain );

SparsityPattern< globalIndex >
buildProlongationSparsity( DofManager const & fineIdx,
                           DofManager const & coarseDofManager,
                           string const & fieldName,
                           ArrayOfSetsView< localIndex const > const & supports );

CRSMatrix< real64, globalIndex >
buildTentativeProlongation( DofManager const & fineIdx,
                            DofManager const & coarseDofManager,
                            string const & fieldName,
                            ArrayOfSetsView< localIndex const > const & supports,
                            arrayView1d< localIndex const > const & initPart );

void makeGlobalDofLists( DofManager const & objIdx,
                         string const & fieldName,
                         arrayView1d< integer const > const & indicator,
                         array1d< globalIndex > & boundaryDof,
                         array1d< globalIndex > & interiorDof );

ArrayOfSets< localIndex >
buildLayeredSupport( integer numLayers,
                     ArrayOfSetsView< localIndex const > const & connectivity,
                     arrayView1d< localIndex const > const & initialPartition );

array1d< integer >
findLayeredSupportBoundary( ArrayOfSetsView< localIndex const > const & connectivity,
                            ArrayOfSetsView< localIndex const > const & support );

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

CRSMatrix< real64, globalIndex >
dropEntries( CRSMatrixView< real64 const, globalIndex const > const & mat,
             real64 relTol );

template< typename MATRIX >
void computeRAP( LinearSolverParameters::Multiscale const & params,
                 MATRIX const & fineMatrix,
                 MATRIX const & prolongation,
                 LinearOperator< typename MATRIX::Vector > const & Rop,
                 MATRIX & coarseMatrix )
{
  if( params.galerkin )
  {
    fineMatrix.multiplyPtAP( prolongation, coarseMatrix );
  }
  else
  {
    MATRIX const & restriction = dynamicCast< MATRIX const & >( Rop );
    fineMatrix.multiplyRAP( restriction, prolongation, coarseMatrix );
  }

  if( params.droptol > 0.0 )
  {
    CRSMatrix< real64, globalIndex > mat = coarseMatrix.extract();
    mat = dropEntries( mat.toViewConst(), params.droptol );
    coarseMatrix.create( mat.toViewConst(), coarseMatrix.numLocalCols(), coarseMatrix.comm() );
  }
}

void writeProlongation( CRSMatrixView< real64 const, globalIndex const > const & prolongation,
                        multiscale::DofManager const & dofManager,
                        string const & fieldName,
                        string const & prefix,
                        multiscale::MeshLevel & mesh,
                        multiscale::MeshObjectManager & fineManager,
                        std::function< void ( multiscale::MeshLevel &, std::vector< string > const & ) > const & writeFunc );

} // namespace msrsb
} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MSRSBUTILS_HPP_
