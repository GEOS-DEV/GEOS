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
 * LAIHelperFunctions.hpp
 */


#ifndef SRC_CORECOMPONENTS_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geosx
{
namespace LAIHelperFunctions
{

/**
 * Create a permuation matrix for a given nodal variable.
 * @param[in]  nodeManager
 * @param[in]  nRows
 * @param[in]  nCols
 * @param[in]  nDofPerNode
 * @param[in]  DofKey
 * @param[out] permutationMatrix
 */
void CreatePermutationMatrix( NodeManager const * const nodeManager,
                              localIndex const nRows,
                              localIndex const nCols,
                              int const nDofPerNode,
                              string const DofKey,
                              ParallelMatrix & permutationMatrix);

/**
 * Create a permuation matrix for a given nodal variable.
 * @param[in]  elementRegionManager
 * @param[in]  nRows
 * @param[in]  nCols
 * @param[in]  nDofPerCell
 * @param[in]  DofKey
 * @param[out] permutationMatrix
 */
void CreatePermutationMatrix( ElementRegionManager const * const elemManager,
                              localIndex const nRows,
                              localIndex const nCols,
                              int const nDofPerNode,
                              string const DofKey,
                              ParallelMatrix & permutationMatrix);

/**
 * Create a permuation matrix for a given nodal variable.
 * @param[in]  nodeManager
 * @param[in]  nRows
 * @param[in]  nCols
 * @param[in]  nDofPerCell
 * @param[in]  DofKey
 * @param[out] permuationMatrix
 */
ParallelVector PermuteVector(ParallelVector const & vector,
                             ParallelMatrix const & permuationMatrix);

/**
 * Permute a square matrix
 * @param[in] matrix to be permuted
 * @param[in] permutation matrix
 * @param[out] permutedMatrix
 */
ParallelMatrix PermuteMatrix(ParallelMatrix const & matrix,
                             ParallelMatrix const & permutationMatrix);

/**
 * Permute a rectangular matrix
 * @param[in] matrix to be permuted
 * @param[in] left permutation matrix
 * @param[in] right permutation matrix
 * @param[out] permutedMatrix
 */
ParallelMatrix PermuteMatrix(ParallelMatrix const & matrix,
                             ParallelMatrix const & permuationMatrixLeft,
                             ParallelMatrix const & permutationMatrixRight);


void PrintPermutedVector(ParallelVector const & vector,
                         ParallelMatrix const & permuationMatrix,
                         std::ostream & os);


void PrintPermutedMatrix(ParallelMatrix const & matrix,
                         ParallelMatrix const & permutationMatrix,
                         std::ostream & os);

void PrintPermutedMatrix(ParallelMatrix const & matrix,
                         ParallelMatrix const & permutationMatrixLeft,
                         ParallelMatrix const & permutationMatrixRight,
                         std::ostream & os);


//void SeparateComponentFilter(ParallelMatrix const & src,
//                             ParallelMatrix & dst,
//                             const localIndex dofsPerNode);

template<typename LAI>
void SeparateComponentFilter(typename LAI::ParallelMatrix const & src,
                             typename LAI::ParallelMatrix & dst,
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

} // LAIHelperFunctions namespace

} // geosx namespace

#endif /* SRC_CORECOMPONENTS_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_ */
