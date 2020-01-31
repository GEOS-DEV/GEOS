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

/*
 * LAIHelperFunctions.cpp
 *
 */

#include "LAIHelperFunctions.hpp"

namespace geosx
{
namespace LAIHelperFunctions
{

void CreatePermutationMatrix(NodeManager const * const nodeManager,
                             localIndex const nRows,
                             localIndex const nCols,
                             int const nDofPerNode,
                             string const DofKey,
                             ParallelMatrix & permutationMatrix)
{
  /* Crates a permutation matrix for a given nodal variable specified by the DofKey. It consider that nDofPerNode
   * dofs are associated with each node (e.g., nDofPerNode = 3 for the displacement).
   *
   * The permutation matrix maps from the dofs ordering set by the DOFManager to the ordering based on the global
   * indexes of the nodes.
   */

  // Create permuation matrix based on size provided.
  permutationMatrix.createWithGlobalSize(nRows, nCols, 1, MPI_COMM_GEOSX);

  arrayView1d<globalIndex const> const &  DofNumber =  nodeManager->getReference<globalIndex_array>( DofKey );

  for( localIndex a=0 ; a<nodeManager->size() ; ++a )
      {
        if (DofNumber[a] >= 0)
        {
          for( int d=0 ; d<nDofPerNode ; ++d )
          {
            globalIndex const rowIndex    = nodeManager->m_localToGlobalMap[a]*nDofPerNode + d;
            globalIndex const columnIndex = DofNumber[a] + d;

            permutationMatrix.insert(rowIndex, columnIndex, 1.0);
          }
        }
      }
  permutationMatrix.close();
  permutationMatrix.set(1);
}

void CreatePermutationMatrix(ElementRegionManager const * const elemManager,
                             localIndex const nRows,
                             localIndex const nCols,
                             int const nDofPerCell,
                             string const DofKey,
                             ParallelMatrix & permutationMatrix)
{
  /* Crates a permutation matrix for a given cell centered variable specified by the DofKey. It consider that nDofPerNode
   * dofs are associated with each node (e.g., nDofPerNode = 3 for the displacement).
   *
   * The permutation matrix maps from the dofs ordering set by the DOFManager to the ordering based on the global
   * indexes of the cells.
   */

  // Create permuation matrix based on size provided.
  permutationMatrix.createWithGlobalSize(nRows, nCols, 1, MPI_COMM_GEOSX);

  elemManager->forElementSubRegions([&]( ElementSubRegionBase const * const elementSubRegion )
  {
    localIndex const numElems = elementSubRegion->size();
    arrayView1d<globalIndex const> const &
    DofNumber = elementSubRegion->getReference< array1d<globalIndex> >( DofKey );

    for( localIndex k=0 ; k<numElems ; ++k )
    {
      if (DofNumber[k] >= 0)
      {
        for( int d=0 ; d<nDofPerCell ; ++d )
        {
          globalIndex const rowIndex    = elementSubRegion->m_localToGlobalMap[k] * nDofPerCell + d;
          globalIndex const columnIndex = DofNumber[k]*nDofPerCell + d;

          permutationMatrix.insert(rowIndex, columnIndex, 1.0);
        }
      }
    }
  });
  permutationMatrix.close();
  permutationMatrix.set(1);
}

ParallelVector PermuteVector(ParallelVector const & vector,
                             ParallelMatrix const & permuationMatrix)
{
  ParallelVector permutedVector;

  permutedVector.createWithLocalSize(vector.localSize(), MPI_COMM_GEOSX);

  permuationMatrix.multiply(vector, permutedVector);

  permutedVector.close();

  return permutedVector;
}

/*
 *  permutedMatrix = permutationMatrix * matrix * permutationMatrix^T;
 */
ParallelMatrix PermuteMatrix(ParallelMatrix const & matrix,
                             ParallelMatrix const & permutationMatrix)
{
  ParallelMatrix temp;
  ParallelMatrix permutedMatrix;
  // The value 24 is hardcoded and should probably be changed (It s fine for displacement).
  temp.createWithLocalSize( matrix.localRows(),
                            matrix.localCols(),
                            24,
                            MPI_COMM_GEOSX );

  permutedMatrix.createWithLocalSize( matrix.localRows(),
                                      matrix.localCols(),
                                      24,
                                      MPI_COMM_GEOSX );

  permutationMatrix.multiply(matrix, temp, true);
  permutationMatrix.rightMultiplyTranspose(temp, permutedMatrix, true);

  return permutedMatrix;
}


/*
 *  permutedMatrix = permutationMatrixLeft * matrix * permutationMatrixRight^T;
 */
ParallelMatrix PermuteMatrix(ParallelMatrix const & matrix,
                             ParallelMatrix const & permutationMatrixLeft,
                             ParallelMatrix const & permutationMatrixRight)
{
  ParallelMatrix temp;
  ParallelMatrix permutedMatrix;

  temp.createWithLocalSize( matrix.localRows(),
                            matrix.localCols(),
                            24,
                            MPI_COMM_GEOSX );

  permutedMatrix.createWithLocalSize( matrix.localRows(),
                                      matrix.localCols(),
                                      24,
                                      MPI_COMM_GEOSX );

  permutationMatrixLeft.multiply( matrix, temp, true );
  permutationMatrixRight.rightMultiplyTranspose( temp, permutedMatrix, true );

  return permutedMatrix;
}

void PrintPermutedVector(ParallelVector const & vector,
                         ParallelMatrix const & permuationMatrix,
                         std::ostream & os )
{
  ParallelVector permutedVector = PermuteVector(vector, permuationMatrix);
  permutedVector.print(os);
}

void PrintPermutedMatrix(ParallelMatrix const & matrix,
                         ParallelMatrix const & permutationMatrix,
                         std::ostream & os)
{
  ParallelMatrix permutedMatrix = PermuteMatrix(matrix, permutationMatrix);
  permutedMatrix.print(os);
}

void PrintPermutedMatrix(ParallelMatrix const & matrix,
                         ParallelMatrix const & permutationMatrixLeft,
                         ParallelMatrix const & permutationMatrixRight,
                         std::ostream & os)

{
  ParallelMatrix permutedMatrix = PermuteMatrix(matrix, permutationMatrixLeft, permutationMatrixRight);
  permutedMatrix.print(os);
}


void SeparateComponentFilter(ParallelMatrix const & src,
                             ParallelMatrix & dst,
                             const localIndex dofsPerNode)
{
  GEOSX_ERROR_IF(dofsPerNode < 2,"Function requires dofsPerNode > 1");

  const localIndex  localRows  = src.localRows();
  const localIndex  maxEntries = src.maxRowLength();
  const localIndex  maxDstEntries = maxEntries / dofsPerNode;

  dst.createWithLocalSize(localRows,maxEntries,MPI_COMM_WORLD);
  dst.open();

  array1d<real64> srcValues;
  array1d<real64> dstValues( maxDstEntries );

  array1d<globalIndex> srcIndices;
  array1d<globalIndex> dstIndices( maxDstEntries );

  for(globalIndex row=src.ilower(); row<src.iupper(); ++row)
  {
     const globalIndex rowComponent = row % dofsPerNode;

     src.getRowCopy(row,srcIndices,srcValues);

     localIndex k=0;
     for(localIndex col=0; col<srcIndices.size(); ++col)
     {
        const globalIndex colComponent = srcIndices[col] % dofsPerNode;
        if( rowComponent == colComponent ) 
        {
          dstValues[k] = srcValues[col];
          dstIndices[k] = srcIndices[col];
          k++;
        } 
     }
     dst.insert(row,dstIndices.data(),dstValues.data(),k);
  }
  dst.close();
}

} // namespace LAIHelperFunctions

} // namespace geosx




