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
 * @file EpetraExport.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_EPETRAEXPORT_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_EPETRAEXPORT_HPP_

#include "common/DataTypes.hpp"

#include <memory>

class Epetra_Map;
class Epetra_Import;

namespace geosx
{

/// Forward declaration
class EpetraMatrix;

/// Forward declaration
class EpetraVector;

/**
 * @brief Facilitates exporting Epetra matrix and associated vector objects (either in parallel or serial).
 */
class EpetraExport
{
public:

  /**
   * @brief Parallel export constructor.
   * @param mat matrix to export, only used to infer parallel partitioning info (doesn't need to be filled yet)
   *
   * Use this constructor when exporting local part of the matrix on each rank.
   */
  EpetraExport();

  /**
   * @brief Destructor.
   */
  ~EpetraExport();

  /**
   * @brief Single-rank export constructor
   * @param mat the matrix
   * @param targetRank the rank to gather matrix on
   *
   * Use this constructor when exporting the whole matrix on a single rank.
   */
  EpetraExport( EpetraMatrix const & mat, integer targetRank );

  /**
   * @brief Export the matrix into CRS arrays provided by the user.
   * @tparam OFFSET_TYPE row pointer offset type
   * @tparam COLUMN_TYPE column index type
   * @param mat the source matrix
   * @param rowOffsets array view to user-managed array of row pointers
   * @param colIndices array view to user-managed array of column indices
   * @param values array view to user-managed array of matrix values
   *
   * This function must be called on all ranks in the matrix's communicator.
   * Only target rank needs to provide meaningful pointer parameters when doing single-rank export.
   */
  template< typename OFFSET_TYPE, typename COLUMN_TYPE >
  void exportCRS( EpetraMatrix const & mat,
                  arrayView1d< OFFSET_TYPE > const & rowOffsets,
                  arrayView1d< COLUMN_TYPE > const & colIndices,
                  arrayView1d< real64 > const & values ) const;

  /**
   * @brief Export the target vector into an array provided by the user.
   * @param vec the source vector, must be compatible with matrix row distribution
   * @param values array view to user-managed array of vector values
   *
   * This method can be used to extract data from vectors associated with the original matrix.
   * It is mostly useful for single-rank export (gather) of vector values.
   * In parallel users can just call vec.extractLocalVector() and obtain values without copying.
   *
   * This function must be called on all ranks in the matrix/vector communicator.
   * Only target rank needs to provide meaningful pointer parameter when doing single-rank export.
   */
  void exportVector( EpetraVector const & vec,
                     arrayView1d< real64 > const & values ) const;

  /**
   * @brief Import the target vector from an array provided by the user.
   * @param values array view to user-managed array of vector values
   * @param vec the target vector, must be compatible with matrix row distribution
   *
   * This method can be used to populate data into vectors associated with the original matrix.
   * It is mostly useful for single-rank import (scatter) of vector values.
   * In parallel users can just call vec.extractLocalVector() and populate values without copying.
   *
   * This function must be called on all ranks in the matrix/vector communicator.
   * Only target rank needs to provide meaningful pointer parameter when doing single-rank import.
   */
  void importVector( arrayView1d< const real64 > const & values, EpetraVector & vec ) const;

private:

  /// Target rank for single-rank export/import
  integer m_targetRank = -1;

  /// Epetra Map for the local vector used in single-rank export/import
  std::unique_ptr< Epetra_Map > m_serialMap;

  /// Epetra Import object used in single-rank export/import
  std::unique_ptr< Epetra_Import > m_serialImport;
};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_INTERFACES_EPETRAEXPORT_HPP_
