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
 * @file SuiteSparse.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_SUITESPARSE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_SUITESPARSE_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <umfpack.h>

namespace geosx
{

/**
 * SuiteSparse integer definition
 */
typedef SuiteSparse_long Int;

/**
 * @brief Convert GEOSX globalIndex value to SuiteSparse Int
 * @param index the input value
 * @return the converted value
 */
inline Int toSuiteSparse_Int( globalIndex const index )
{
  return LvArray::integerConversion< Int >( index );
}

/**
 * @class SuiteSparse
 * @brief This class provides an interface for UMFPACK direct solver from
 *        SuiteSparse linear algebra package.
 */
class SuiteSparse
{

public:

  /**
   * @brief Constructor
   */
  SuiteSparse();

  /**
   * @brief Constructor with parameters
   * @param[in] params the linear solver parameters
   */
  SuiteSparse( LinearSolverParameters const & params );

  /**
   * @brief Destructor
   */
  ~SuiteSparse();

  /**
   * @brief Creates the SuiteSparse data structure
   * @param[in] params the linear solver parameters
   */
  void create( LinearSolverParameters const & params );

  /**
   * @brief Factorizes a linear system with SuiteSparse
   * @return info error code
   */
  int setup();

  /**
   * @brief Solves a linear system with SuiteSparse (matrix has already been factorized)
   * @param[in] b the right-hand side
   * @param[out] x the solution
   * @return info error code
   */
  int solveWorkingRank( real64 * b, real64 * x );

  /**
   * @brief Sycronizes times across ranks
   */
  void syncTimes();

  /**
   * @brief Estimates the condition number of the matrix
   * @return the estimated condition number
   */
  real64 condEst() const;

  /**
   * @brief Estimates the relative tolerance for the matrix
   * @return the relative tolerance (condEst * eps)
   */
  real64 relativeTolerance() const;

  /**
   * @brief Deallocates a SuiteSparse data structure
   */
  void destroy();

  /**
   * @brief Sets the working rank
   * @param[in] workingRank the working rank
   */
  void setWorkingRank( int const workingRank );

  /**
   * @brief Returns the working rank
   * @return the working rank
   */
  int workingRank() const;

  /**
   * @brief Sets the working rank in the sub-communicator
   * @param[in] subCommWorkingRank the working rank in the sub-communicator
   */
  void setSubCommWorkingRank( int const subCommWorkingRank );

  /**
   * @brief Returns the working rank in the sub-communicator
   * @return the working rank in the sub-communicator
   */
  int subCommWorkingRank() const;

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
   * @brief Sets the subcommunicator
   * @param[in] subComm the MPI subcommunicator
   */
  void setSubComm( MPI_Comm const subComm );

  /**
   * @brief Returns the subcommunicator
   * @return the subcommunicator
   */
  MPI_Comm getSubComm() const;

  /**
   * @brief Sets the number of rows
   * @param[in] numRows the number of rows
   */
  void setNumRows( Int const numRows );

  /**
   * @brief Returns the number of rows
   * @returns the number of rows
   */
  Int numRows() const;

  /**
   * @brief Sets the number of columns
   * @param[in] numCols the number of columns
   */
  void setNumCols( Int const numCols );

  /**
   * @brief Returns the number of columns
   * @return the number of columns
   */
  Int numCols() const;

  /**
   * @brief Sets the number of non zeros
   * @param[in] nonZeros the number of non zeros
   */
  void setNonZeros( Int const nonZeros );

  /**
   * @brief Returns the number of non zeros
   * @return the number of non zeros
   */
  Int nonZeros() const;

  /**
   * @brief Allocate the internal data storage arrays
   */
  void createInternalStorage();

  /**
   * @brief Returns the array with the row pointers
   * @return the array with the row pointers
   */
  array1d< Int > & rowPtr();

  /**
   * @brief Returns the array with the column indices
   * @return the array with the column indices
   */
  array1d< Int > & colIndices();

  /**
   * @brief Returns the array with the matrix values
   * @return the array with the matrix values
   */
  array1d< real64 > & values();

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

private:

  /// log level
  integer m_logLevel;

  /// number of rows
  Int m_numRows;

  /// number of columns
  Int m_numCols;

  /// number of entries
  Int m_nonZeros;

  /// row pointers
  array1d< Int > m_rowPtr;

  /// column indices
  array1d< Int > m_colIndices;

  /// values
  array1d< real64 > m_values;

  /// data structure to gather various info
  real64 m_Info[UMFPACK_INFO];

  /// SuiteSparse options
  real64 m_Control[UMFPACK_CONTROL];

  /// pointer to the symbolic factorization
  void * m_Symbolic;

  /// pointer to the numeric factorization
  void * m_Numeric;

  /// MPI communicator
  MPI_Comm m_comm;

  /// MPI sub-communicator for ranks that have parts of the matrix
  MPI_Comm m_subComm;

  /// flag to check if the sub-communicator is used
  bool m_usingSubComm;

  /// MPI rank carring out the solution
  int m_workingRank;

  /// MPI rank carring out the solution in the sub-communicator
  int m_subCommWorkingRank;

  /// condition number estimation
  real64 m_condEst;

  /// setup time
  real64 m_setupTime;

  /// solve time
  real64 m_solveTime;

};

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_SUITESPARSE_HPP_*/
