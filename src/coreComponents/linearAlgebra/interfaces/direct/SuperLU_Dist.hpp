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
 * @file SuperLU_Dist.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_SUPERLU_DIST_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_SUPERLU_DIST_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <superlu_ddefs.h>

namespace geosx
{

/**
 * @brief Convert GEOSX globalIndex value to SuperLU_Dist int_t
 * @param index the input value
 * @return the converted value
 */
inline int_t toSuperLU_intT( globalIndex const index )
{
  return LvArray::integerConversion< int_t >( index );
}

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to SuperLU_Dist int_t
 * @param[in] index the input array
 * @return the converted array
 */
inline int_t * toSuperLU_intT( globalIndex * const index )
{
  return reinterpret_cast< int_t * >( index );
}

/**
 * SuperLU_Dist data
 */
class SuperLU_Dist
{

public:

  /**
   * @brief Constructor
   */
  SuperLU_Dist();

  /**
   * @brief Constructor with parameters
   * @param[in] params the linear solver parameters
   */
  SuperLU_Dist( LinearSolverParameters const & params );

  /**
   * @brief Destructor
   */
  ~SuperLU_Dist();

  /**
   * @brief Creates the SuperLU_Dist data structure
   * @param[in] params the linear solver parameters
   */
  void create( LinearSolverParameters const & params );

  /**
   * @brief Factorizes a linear system with SuperLU_Dist
   * @return info error code
   */
  int setup();

  /**
   * @brief Solves a linear system with SuperLU_Dist (matrix has already been factorized)
   * @param[in] b the right-hand side
   * @param[out] x the solution vector
   * @return info error code
   */
  int solve( real64 const * b, real64 * x );

  /**
   * @brief Estimates the condition number of the matrix
   * @return the estimated condition number
   */
  real64 condEst();

  /**
   * @brief Estimates the relative tolerance for the matrix
   * @return the relative tolerance (condEst * eps)
   */
  real64 relativeTolerance();

  /**
   * @brief Deallocates a SuperLU_Dist data structure
   */
  void destroy();

  /**
   * @brief Sets the global number of rows
   * @param numGlobalRows the global number of rows
   */
  void setNumGlobalRows( int_t const numGlobalRows );

  /**
   * @brief Returns the global number of rows
   * @return the global number of rows
   */
  int_t numGlobalRows() const;

  /**
   * @brief Sets the local number of rows
   * @param numLocalRows the local number of rows
   */
  void setNumLocalRows( int_t const numLocalRows );

  /**
   * @brief Returns the local number of rows
   * @return the local number of rows
   */
  int_t numLocalRows() const;

  /**
   * @brief Sets the communicator
   * @param[in] comm the MPI communicator
   */
  void setComm( MPI_Comm const comm );

  /**
   * @brief Returns the communicator
   * @return the communicator
   */
  MPI_Comm getComm() const;

  /**
   * @brief Returns the matrix in SuperLU_Dist SuperMatrix format
   * @return the matrix in SuperLU_Dist SuperMatrix format
   */
  SuperMatrix & mat();

  /**
   * @brief Allocates the row pointers array
   * @param numRows the number of rows
   */
  void createRowPtr( localIndex const numRows );

  /**
   * @brief Sets the row pointers array
   * @param rowPtr the row pointers array
   */
  void setRowPtr( int_t * const rowPtr );

  /**
   * @brief Returns the array with the row pointers
   * @return the array with the row pointers
   */
  int_t * rowPtr();

  /**
   * @brief Allocates the column indices array
   * @param numNonzeros the number of entries
   */
  void createColIndices( localIndex const numNonzeros );

  /**
   * @brief Sets the column indices array
   * @param colIndices the column indices array
   */
  void setColIndices( int_t * const colIndices );

  /**
   * @brief Returns the array with the column indices
   * @return the array with the column indices
   */
  int_t * colIndices();

  /**
   * @brief Allocates the values array
   * @param numNonzeros the number of entries
   */
  void createValues( localIndex const numNonzeros );

  /**
   * @brief Sets the values array
   * @param values the values array
   */
  void setValues( real64 * const values );

  /**
   * @brief Returns the array with the values
   * @return the array with the values
   */
  real64 * values();

  /**
   * @brief Provides the setup time
   * @return the setup time
   */
  real64 setupTime() const;

  /**
   * @brief Provides the solve time
   * @return the solve time
   */
  real64 solveTime() const;

  /**
   * @brief Returns the parameters used to initialize this object
   * @return the parameters used to initialize this object
   */
  LinearSolverParameters getParameters() const;

  /**
   * @brief Returns the machine precision used in SuperLU_Dist class
   * @return the machine precision used in SuperLU_Dist class
   */
  real64 machinePrecision() const;

private:

  /// number of global rows
  int_t m_numGlobalRows;

  /// number of local rows
  int_t m_numLocalRows;

  /// row pointers
  int_t * m_rowPtr;

  /// column indices
  int_t * m_colIndices;

  /// values
  real64 * m_values;

  /// SuperLU_Dist matrix format
  SuperMatrix m_mat;

  /// data structure to scale and permute the matrix
  dScalePermstruct_t m_ScalePermstruct;

  /// data structure to store the LU factorization
  dLUstruct_t m_LUstruct;

  /// data structure to gather some statistics
  SuperLUStat_t m_stat;

  /// SuperLU_Dist MPI subdivision of load
  gridinfo_t m_grid;

  /// data structure to solve the matrix
  dSOLVEstruct_t m_SOLVEstruct;

  /// SuperLU_Dist options
  superlu_dist_options_t m_options;

  /// MPI communicator
  MPI_Comm m_comm;

  /// condition number estimation
  real64 m_condEst;

  /// setup time
  real64 m_setupTime;

  /// solve time
  real64 m_solveTime;

  /// linear solver parameters
  LinearSolverParameters m_params;

  /// Add two orders of magnitude to allow small error in condition number estimate
  real64 const m_machinePrecision = 100.0 * std::numeric_limits< real64 >::epsilon();
};

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_SUPERLU_DIST_HPP_*/
