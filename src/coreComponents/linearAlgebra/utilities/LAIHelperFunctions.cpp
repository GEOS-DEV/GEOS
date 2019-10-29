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

void CreateNodalUnknownPermutationMatrix(NodeManager* const nodeManager,
                                         localIndex const nRows,
                                         localIndex const nCols,
                                         int const nDofPerNode,
                                         string const DofKey,
                                         ParallelMatrix & permutationMatrix)
{
  // Before outputting anything I build the permutation matrix
  permutationMatrix.createWithLocalSize(nRows, nCols, 1, MPI_COMM_GEOSX);

  arrayView1d<globalIndex> const &  DofNumber =  nodeManager->getReference<globalIndex_array>( DofKey );

  for( localIndex a=0 ; a<nodeManager->size() ; ++a )
      {
        for( int d=0 ; d<nDofPerNode ; ++d )
        {
          globalIndex const rowIndex = nodeManager->m_localToGlobalMap[a]*nDofPerNode + d;
          globalIndex const columnIndex = DofNumber[a] + d;

          permutationMatrix.insert(rowIndex, columnIndex, 1.0);
        }
      }
  permutationMatrix.close();
  permutationMatrix.set(1);
}

void CreateCellUnknownPermutationMatrix()
{
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
 *  permutedMatrix = permutationMatrixLeft * matrix * permutationMatrix^T;
 */
ParallelMatrix PermuteMatrix(ParallelMatrix const & matrix,
                             ParallelMatrix const & permutationMatrix)
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

  matrix.multiply( false, permutationMatrix, true, temp, false );
  temp.close();
  permutationMatrix.multiply( false, temp, false, permutedMatrix, false );
  permutedMatrix.close();
  // permutedMatrix.unwrappedPointer()->MakeDataContiguous();

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

  matrix.multiply( false, permutationMatrixRight, true, temp, false );
  temp.close();
  permutationMatrixLeft.multiply( false, temp, false, permutedMatrix, false );
  permutedMatrix.close();
  // permutedMatrix.unwrappedPointer()->MakeDataContiguous();

  return permutedMatrix;
}


} // namespace LAIHelperFunctions

} // namespace geosx




