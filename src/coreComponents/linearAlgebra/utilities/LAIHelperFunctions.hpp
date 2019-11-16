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


void SeparateComponentFilter(ParallelMatrix const & src,
                             ParallelMatrix & dst,
                             const localIndex dofsPerNode);

} // LAIHelperFunctions namespace

} // geosx namespace

#endif /* SRC_CORECOMPONENTS_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_ */
