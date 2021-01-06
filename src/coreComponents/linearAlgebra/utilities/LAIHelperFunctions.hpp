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

/*
 * LAIHelperFunctions.hpp
 */


#ifndef GEOSX_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geosx
{
namespace LAIHelperFunctions
{

/**
 * @brief Create a permutation matrix for a given nodal variable.
 * @param[in]  nodeManager       the node manager
 * @param[in]  nRows             number of local rows in the matrix
 * @param[in]  nCols             number of local columns in the matrix
 * @param[in]  nDofPerNode       number of degrees-of-freedom per node
 * @param[in]  dofKey            DofManager key used to access dof index array
 * @param[out] permutationMatrix the target matrix
 */
void CreatePermutationMatrix( NodeManager const * const nodeManager,
                              localIndex const nRows,
                              localIndex const nCols,
                              int const nDofPerNode,
                              string const dofKey,
                              ParallelMatrix & permutationMatrix );

/**
 * @brief Create a permutation matrix for a given nodal variable.
 * @param[in]  elemManager       the element region manager
 * @param[in]  nRows             number of local rows in the matrix
 * @param[in]  nCols             number of local columns in the matrix
 * @param[in]  nDofPerNode       number of degrees-of-freedom per node
 * @param[in]  dofKey            DofManager key used to access dof index array
 * @param[out] permutationMatrix the target matrix
 */
void CreatePermutationMatrix( ElementRegionManager const * const elemManager,
                              localIndex const nRows,
                              localIndex const nCols,
                              int const nDofPerNode,
                              string const DofKey,
                              ParallelMatrix & permutationMatrix );

/**
 * @brief Permute a vector.
 * @param[in] vector            the source vector
 * @param[in] permutationMatrix the permutation matrix
 * @return the permuted vector
 */
ParallelVector PermuteVector( ParallelVector const & vector,
                              ParallelMatrix const & permutationMatrix );

/**
 * @brief Permute rows and columns of a square matrix.
 * @param[in] matrix            the source matrix
 * @param[in] permutationMatrix permutation matrix
 * @return the permuted matrix
 */
ParallelMatrix PermuteMatrix( ParallelMatrix const & matrix,
                              ParallelMatrix const & permutationMatrix );

/**
 * Permute rows and columns of a rectangular matrix
 * @param[in] matrix                 the source matrix
 * @param[in] permutationMatrixLeft  left permutation matrix
 * @param[in] permutationMatrixRight right permutation matrix
 * @return the permuted matrix
 */
ParallelMatrix PermuteMatrix( ParallelMatrix const & matrix,
                              ParallelMatrix const & permutationMatrixLeft,
                              ParallelMatrix const & permutationMatrixRight );

void PrintPermutedVector( ParallelVector const & vector,
                          ParallelMatrix const & permuationMatrix,
                          std::ostream & os );


void PrintPermutedMatrix( ParallelMatrix const & matrix,
                          ParallelMatrix const & permutationMatrix,
                          std::ostream & os );

void PrintPermutedMatrix( ParallelMatrix const & matrix,
                          ParallelMatrix const & permutationMatrixLeft,
                          ParallelMatrix const & permutationMatrixRight,
                          std::ostream & os );

/**
 * @brief Apply a separate component approximation (filter) to a matrix.
 * @tparam MATRIX the type of matrices
 * @param src         the source matrix
 * @param dst         the target (filtered) matrix
 * @param dofsPerNode number of degrees-of-freedom per node
 */
template< typename MATRIX >
void SeparateComponentFilter( MATRIX const & src,
                              MATRIX & dst,
                              const localIndex dofsPerNode )
{
  GEOSX_ERROR_IF( dofsPerNode < 2, "Function requires dofsPerNode > 1" );

  const localIndex localRows  = src.numLocalRows();
  const localIndex maxEntries = src.maxRowLength();
  const localIndex maxDstEntries = maxEntries / dofsPerNode;

  dst.createWithLocalSize( localRows, maxEntries, MPI_COMM_WORLD );
  dst.open();

  array1d< real64 > srcValues;
  array1d< real64 > dstValues( maxDstEntries );

  array1d< globalIndex > srcIndices;
  array1d< globalIndex > dstIndices( maxDstEntries );

  for( globalIndex row=src.ilower(); row<src.iupper(); ++row )
  {
    const globalIndex rowComponent = row % dofsPerNode;
    const localIndex rowLength = src.globalRowLength( row );
    srcIndices.resize( rowLength );
    srcValues.resize( rowLength );

    src.getRowCopy( row, srcIndices, srcValues );

    localIndex k=0;
    for( localIndex col=0; col<rowLength; ++col )
    {
      const globalIndex colComponent = srcIndices[col] % dofsPerNode;
      if( rowComponent == colComponent )
      {
        dstValues[k] = srcValues[col];
        dstIndices[k] = srcIndices[col];
        k++;
      }
    }
    dst.insert( row, dstIndices.data(), dstValues.data(), k );
  }
  dst.close();
}

template< typename VECTOR >
void ComputeRigidBodyModes( MeshLevel const & mesh,
                            DofManager const & dofManager,
                            std::vector< std::string > const & selection,
                            array1d< VECTOR > & rigidBodyModes )
{
  NodeManager const & nodeManager = *mesh.getNodeManager();

  localIndex numComponents = 0;
  array1d< globalIndex > globalNodeList;
  for( localIndex k = 0; k < LvArray::integerConversion< localIndex >( selection.size() ); ++k )
  {
    if( dofManager.getLocation( selection[k] ) == DofManager::Location::Node )
    {
      string const & dispDofKey = dofManager.getKey( selection[k] );
      arrayView1d< globalIndex const > const & dofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
      localIndex const numComponentsField = dofManager.numComponents( selection[k] );
      numComponents = numComponents > 0 ? numComponents : numComponentsField;
      GEOSX_ERROR_IF( numComponents != numComponentsField, "Rigid body modes called with different number of components." );
      globalIndex const globalOffset = dofManager.globalOffset( selection[k] );
      globalIndex const numLocalDofs = LvArray::integerConversion< globalIndex >( dofManager.numLocalDofs( selection[k] ) );
      for( globalIndex i = 0; i < dofNumber.size(); ++i )
      {
        if( dofNumber[i] >= globalOffset && ( dofNumber[i] - globalOffset ) < numLocalDofs )
        {
          globalNodeList.emplace_back( ( dofNumber[i]-globalOffset )/numComponentsField );
        }
      }
    }
  }
  localIndex const numNodes = globalNodeList.size();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  localIndex const numRidigBodyModes = numComponents * ( numComponents + 1 ) / 2;
  rigidBodyModes.resize( numRidigBodyModes );
  for( localIndex k = 0; k < numComponents; ++k )
  {
    rigidBodyModes[k].createWithLocalSize( numNodes*numComponents, MPI_COMM_GEOSX );
    rigidBodyModes[k].open();
    for( localIndex i = 0; i < numNodes; ++i )
    {
      rigidBodyModes[k].set( rigidBodyModes[k].getGlobalRowID( numComponents*i+k ), 1.0 );
    }
    rigidBodyModes[k].close();
    rigidBodyModes[k].scale( 1.0/rigidBodyModes[k].norm2() );
  }
  switch( numComponents )
  {
    case 2:
    {
      localIndex const k = 2;
      rigidBodyModes[k].createWithLocalSize( numNodes*numComponents, MPI_COMM_GEOSX );
      rigidBodyModes[k].open();
      for( localIndex i = 0; i < numNodes; ++i )
      {
        rigidBodyModes[k].set( rigidBodyModes[k].getGlobalRowID( numComponents*i+0 ), -nodePosition[globalNodeList[i]][1] );
        rigidBodyModes[k].set( rigidBodyModes[k].getGlobalRowID( numComponents*i+1 ), +nodePosition[globalNodeList[i]][0] );
      }
      rigidBodyModes[k].close();
      for( localIndex j = 0; j < k; ++j )
      {
        rigidBodyModes[k].axpy( -rigidBodyModes[k].dot( rigidBodyModes[j] ), rigidBodyModes[j] );
      }
      rigidBodyModes[k].scale( 1.0/rigidBodyModes[k].norm2() );
      break;
    }
    case 3:
    {
      localIndex k = 3;
      rigidBodyModes[k].createWithLocalSize( numNodes*numComponents, MPI_COMM_GEOSX );
      rigidBodyModes[k].open();
      for( localIndex i = 0; i < numNodes; ++i )
      {
        rigidBodyModes[k].set( rigidBodyModes[k].getGlobalRowID( numComponents*i+0 ), +nodePosition[globalNodeList[i]][1] );
        rigidBodyModes[k].set( rigidBodyModes[k].getGlobalRowID( numComponents*i+1 ), -nodePosition[globalNodeList[i]][0] );
      }
      rigidBodyModes[k].close();
      for( localIndex j = 0; j < k; ++j )
      {
        rigidBodyModes[k].axpy( -rigidBodyModes[k].dot( rigidBodyModes[j] ), rigidBodyModes[j] );
      }
      rigidBodyModes[k].scale( 1.0/rigidBodyModes[k].norm2() );

      ++k;
      rigidBodyModes[k].createWithLocalSize( numNodes*numComponents, MPI_COMM_GEOSX );
      rigidBodyModes[k].open();
      for( localIndex i = 0; i < numNodes; ++i )
      {
        rigidBodyModes[k].set( rigidBodyModes[k].getGlobalRowID( numComponents*i+1 ), -nodePosition[globalNodeList[i]][2] );
        rigidBodyModes[k].set( rigidBodyModes[k].getGlobalRowID( numComponents*i+2 ), +nodePosition[globalNodeList[i]][1] );
      }
      rigidBodyModes[k].close();
      for( localIndex j = 0; j < k; ++j )
      {
        rigidBodyModes[k].axpy( -rigidBodyModes[k].dot( rigidBodyModes[j] ), rigidBodyModes[j] );
      }
      rigidBodyModes[k].scale( 1.0/rigidBodyModes[k].norm2() );

      ++k;
      rigidBodyModes[k].createWithLocalSize( numNodes*numComponents, MPI_COMM_GEOSX );
      rigidBodyModes[k].open();
      for( localIndex i = 0; i < numNodes; ++i )
      {
        rigidBodyModes[k].set( rigidBodyModes[k].getGlobalRowID( numComponents*i+0 ), +nodePosition[globalNodeList[i]][2] );
        rigidBodyModes[k].set( rigidBodyModes[k].getGlobalRowID( numComponents*i+2 ), -nodePosition[globalNodeList[i]][0] );
      }
      rigidBodyModes[k].close();
      for( localIndex j = 0; j < k; ++j )
      {
        rigidBodyModes[k].axpy( -rigidBodyModes[k].dot( rigidBodyModes[j] ), rigidBodyModes[j] );
      }
      rigidBodyModes[k].scale( 1.0/rigidBodyModes[k].norm2() );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Rigid body modes computation unsupported for " << numComponents << " components." );
    }
  }
}

} // LAIHelperFunctions namespace

} // geosx namespace

#endif /*GEOSX_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_*/
