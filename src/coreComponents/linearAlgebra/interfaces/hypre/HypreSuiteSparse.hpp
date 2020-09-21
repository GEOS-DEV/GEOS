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
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <HYPRE_utilities.h>
#include <seq_mv.h>
#include <umfpack.h>

namespace geosx
{

/**
 * SuiteSparse integer definition
 */
typedef SuiteSparse_long Int;

/**
 * SuiteSparse data
 */
struct SuiteSparseData
{
  hypre_CSRMatrix * CSRmatrix;        //!< global matrix in HypreMatrix format
  integer logLevel;                   //!< log level
  Int numRows;                        //!< number of rows
  Int numCols;                        //!< number of columns
  Int nonZeros;                       //!< number of entries
  Int * rowPtr;                       //!< row pointers
  Int * colIndices;                   //!< column indices
  real64 * data;                      //!< values
  real64 Info[UMFPACK_INFO];          //!< data structure to gather various info
  real64 Control[UMFPACK_CONTROL];    //!< SuiteSparse options
  void * Symbolic;                    //!< pointer to the symbolic factorization
  void * Numeric;                     //!< pointer to the numeric factorization
  MPI_Comm comm;                      //!< MPI communicator
  int workingRank;                    //!< MPI rank carring out the solution
};

/**
 * @brief Creates the SuiteSparse data structure
 * @param[in] matrix the HypreMatrix object
 * @param[in] params the linear solver parameters
 * @param[out] SSData the structure containing the matrix in SuiteSparse format
 */
void SuiteSparseCreate( HypreMatrix const & matrix,
                        LinearSolverParameters const & params,
                        SuiteSparseData & SSData );

/**
 * @brief Factorizes a linear system with SuiteSparse
 * @param[in,out] SSData the structure containing the matrix in SuiteSparse format
 * @param[out] time time spent in the factorization phase
 * @return info error code
 */
int SuiteSparseSetup( SuiteSparseData & SSData,
                      real64 & time );

/**
 * @brief Solves a linear system with SuiteSparse (matrix has already been factorized)
 * @param[in,out] SSData the structure containing the matrix in SuiteSparse format
 * @param[in] b the right-hand side in Hypre format
 * @param[out] x the solution in Hypre format
 * @param[out] time time spent in the solution phase
 * @return info error code
 */
int SuiteSparseSolve( SuiteSparseData & SSData,
                      HypreVector const & b,
                      HypreVector & x,
                      real64 & time );

/**
 * @brief Estimates the condition number of the matrix
 * @param[in] SSData the structure containing the matrix in SuiteSparse format
 * @return the estimated condition number
 */
real64 SuiteSparseCondEst( SuiteSparseData const & SSData );

/**
 * @brief Deallocates a SuiteSparse data structure
 * @param[in,out] SSData the structure containing the matrix in SuiteSparse format
 */
void SuiteSparseDestroy( SuiteSparseData & SSData );

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUITESPARSE_HPP_*/
