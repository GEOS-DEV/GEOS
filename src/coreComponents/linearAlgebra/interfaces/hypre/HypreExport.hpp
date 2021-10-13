/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreExport.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREEXPORT_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREEXPORT_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

/// Forward declaration
class HypreMatrix;

/// Forward declaration
class HypreVector;

/**
 * @brief Facilitates exporting Hypre matrix and associated vector objects (either in parallel or serial).
 */
class HypreExport
{
public:

  /**
   * @brief Parallel export constructor.
   *
   * Use this constructor when exporting local part of the matrix on each rank.
   */
  HypreExport();

  /**
   * @brief Single-rank export constructor
   * @param mat matrix to export, only used to infer parallel partitioning info (doesn't need to be filled yet)
   * @param targetRank the rank to gather matrix on
   *
   * Use this constructor when exporting the whole matrix on a single rank.
   */
  HypreExport( HypreMatrix const & mat, integer targetRank );

  /**
   * @brief Destructor.
   */
  ~HypreExport();

  /**
   * @brief Export the matrix into CRS arrays provided by the user.
   * @tparam OFFSET_TYPE row pointer offset type
   * @tparam COLUMN_TYPE column index type
   * @param mat the source matrix
   * @param rowOffsets pointer to user-managed array of row pointers
   * @param colIndices pointer to user-managed array of column indices
   * @param values pointer to user-managed array of matrix values
   *
   * This function must be called on all ranks in the matrix's communicator.
   * Only target rank needs to provide meaningful pointer parameters when doing single-rank export.
   */
  template< typename OFFSET_TYPE, typename COLUMN_TYPE >
  void exportCRS( HypreMatrix const & mat,
                  OFFSET_TYPE * rowOffsets,
                  COLUMN_TYPE * colIndices,
                  real64 * values ) const;

  /**
   * @brief Export the target vector into an array provided by the user.
   * @param vec the source vector, must be compatible with matrix row distribution
   * @param values pointer to user-managed array of vector values
   *
   * This method can be used to extract data from vectors associated with the original matrix.
   * It is mostly useful for single-rank export (gather) of vector values.
   * In parallel users can just call vec.extractLocalVector() and obtain values without copying.
   *
   * This function must be called on all ranks in the matrix/vector communicator.
   * Only target rank needs to provide meaningful pointer parameter when doing single-rank export.
   */
  void exportVector( HypreVector const & vec, real64 * values ) const;

  /**
   * @brief Import the target vector from an array provided by the user.
   * @param values pointer to user-managed array of vector values
   * @param vec the target vector, must be compatible with matrix row distribution
   *
   * This method can be used to populate data into vectors associated with the original matrix.
   * It is mostly useful for single-rank import (scatter) of vector values.
   * In parallel users can just call vec.extractLocalVector() and populate values without copying.
   *
   * This function must be called on all ranks in the matrix/vector communicator.
   * Only target rank needs to provide meaningful pointer parameter when doing single-rank import.
   */
  void importVector( real64 const * values, HypreVector & vec ) const;

private:

  /// Target rank for single-rank export/import
  integer m_targetRank = -1;

  /// Subcommunicator for ranks that have rows (needed by hypre's scatter op).
  MPI_Comm m_subComm = MPI_COMM_NULL;
};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_INTERFACES_HYPREEXPORT_HPP_
