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

/**
 * @file HypreSuiteSparse.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUITESPARSE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUITESPARSE_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/direct/SuiteSparse.hpp"

#include <HYPRE_utilities.h>
#include <seq_mv.h>

namespace geosx
{

/**
 * @brief Converts a matrix from Hypre to SuiteSparse format
 * @param[in] matrix the HypreMatrix object
 * @param[out] SSData the structure containing the matrix in SuiteSparse format
 */
void ConvertHypreToSuiteSparseMatrix( HypreMatrix const & matrix,
                                      SuiteSparse & SSData );

/**
 * @brief Solves a linear system with SuiteSparse (matrix has already been factorized)
 * @param[in,out] SSData the structure containing the matrix in SuiteSparse format
 * @param[in] b the right-hand side in Hypre format
 * @param[out] x the solution in Hypre format
 * @return info error code
 */
int SuiteSparseSolve( SuiteSparse & SSData,
                      HypreVector const & b,
                      HypreVector & x );

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUITESPARSE_HPP_*/
