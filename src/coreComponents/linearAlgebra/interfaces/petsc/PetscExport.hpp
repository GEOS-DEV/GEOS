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
 * @file PetscExport.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_PETSCEXPORT_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_PETSCEXPORT_HPP_

#include "common/DataTypes.hpp"

/// IS struct forward declaration
extern "C" struct _p_IS;

/// VecScatter struct forward declaration
extern "C" struct _p_VecScatter;

namespace geosx
{

/// Forward declaration
class PetscMatrix;

/// Forward declaration
class PetscVector;

/**
 * @brief Facilitates exporting Petsc matrix and associated vector objects (either in parallel or serial).
 */
class PetscExport
{
public:

  /**
   * @brief Parallel export constructor.
   *
   * Use this constructor when exporting local part of the matrix on each rank.
   */
  PetscExport();

  /**
   * @brief Single-rank export constructor
   * @param mat matrix to export, only used to infer parallel partitioning info (doesn't need to be filled yet)
   * @param targetRank the rank to gather matrix on
   *
   * Use this constructor when exporting the whole matrix on a single rank.
   */
  PetscExport( PetscMatrix const & mat, integer targetRank );

  /**
   * @brief Destructor.
   */
  ~PetscExport();

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
  void exportCRS( PetscMatrix const & mat,
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
  void exportVector( PetscVector const & vec, real64 * values ) const;

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
  void importVector( real64 const * values, PetscVector & vec ) const;

private:

  /// Alias for PETSc index set struct pointer
  using IS = struct _p_IS *;

  /// Alias for PETSc vector scatter struct pointer
  using VecScatter = struct _p_VecScatter *;

  /// Target rank for single-rank export/import
  integer m_targetRank = -1;

  /// Index set used in single-rank export/import
  IS m_indexSet{};

  /// PETSc vector scatter context used in single-rank export/import
  VecScatter m_scatter{};
};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_INTERFACES_PETSCEXPORT_HPP_
