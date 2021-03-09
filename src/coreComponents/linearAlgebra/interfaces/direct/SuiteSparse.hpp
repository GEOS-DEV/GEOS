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

// Pre-define some suitesparse variables since they are not properly defined in the header for alternate index types.
#if GEOSX_GLOBALINDEX_TYPE_FLAG==0
#define SuiteSparse_long int
#define SuiteSparse_long_max 2147483647
#define SuiteSparse_long_idd "d"
#define SuiteSparse_long_id "%d"
#endif

#include <umfpack.h>

namespace geosx
{

/**
 * SuiteSparse integer definition
 */
using SSInt = SuiteSparse_long;

/**
 * @brief Convert GEOSX globalIndex value to SuiteSparse int
 * @param index the input value
 * @return the converted value
 */
inline SSInt toSuiteSparse_Int( globalIndex const index )
{
  return LvArray::integerConversion< SSInt >( index );
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
   * @param[in] transpose whether to solve for the original or the transpose matrix
   * @return info error code
   */
  int solveWorkingRank( real64 * b, real64 * x, bool transpose = false );

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
   * @brief Returns the number of rows
   * @returns the number of rows
   */
  SSInt numRows() const;

  /**
   * @brief Returns the number of columns
   * @return the number of columns
   */
  SSInt numCols() const;

  /**
   * @brief Returns the number of non zeros
   * @return the number of non zeros
   */
  SSInt nonZeros() const;

  /**
   * @brief Allocate the internal data storage arrays
   * @param[in] numRows the number of rows
   * @param[in] numCols the number of columns
   * @param[in] nonZeros the number of non zeros
   */
  void resize( SSInt const numRows, SSInt const numCols, SSInt const nonZeros );

  /**
   * @brief Returns the array with the row pointers
   * @return the array with the row pointers
   */
  arrayView1d< SSInt > rowPtr();

  /**
   * @brief Returns the array with the column indices
   * @return the array with the column indices
   */
  arrayView1d< SSInt > colIndices();

  /**
   * @brief Returns the array with the matrix values
   * @return the array with the matrix values
   */
  arrayView1d< real64 > values();

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
   * @brief Returns the precision tolarance used in SuiteSparse class
   * @return the precision tolerance used in SuiteSparse class
   */
  real64 precisionTolerance() const;

private:

  /// log level
  integer m_logLevel;

  /// number of rows
  SSInt m_numRows;

  /// number of columns
  SSInt m_numCols;

  /// number of entries
  SSInt m_nonZeros;

  /// row pointers
  array1d< SSInt > m_rowPtr;

  /// column indices
  array1d< SSInt > m_colIndices;

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

  /// Add two orders of magnitude to allow small error in condition number estimate
  real64 const m_precisionTolerance = 100.0 * std::numeric_limits< real64 >::epsilon();

};

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_SUITESPARSE_HPP_*/
