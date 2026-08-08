/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DistributedMatrixStatistics.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_DISTRIBUTEDMATRIXSTATISTICS_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_DISTRIBUTEDMATRIXSTATISTICS_HPP_

#include "common/DataTypes.hpp"

#include <vector>

namespace geos
{

/**
 * @brief Structural and communication statistics for one rank of a row-distributed sparse matrix.
 */
struct MatrixRankStatistics
{
  /// First global row owned by this rank.
  globalIndex firstRow = 0;

  /// Number of rows owned by this rank.
  globalIndex numRows = 0;

  /// Number of stored entries in the owned rows.
  globalIndex numNonzeros = 0;

  /// Number of stored entries whose column is owned by another rank.
  globalIndex numOffRankNonzeros = 0;

  /// Number of unique vector entries received for one sparse matrix-vector product.
  globalIndex numHaloReceiveDofs = 0;

  /// Number of unique owned vector entries sent for one sparse matrix-vector product.
  globalIndex numHaloSendDofs = 0;

  /// Number of ranks from which vector entries are received.
  globalIndex numReceiveNeighbors = 0;

  /// Number of ranks to which vector entries are sent.
  globalIndex numSendNeighbors = 0;

  /// Maximum number of stored entries in an owned row.
  globalIndex maxRowNonzeros = 0;

  /// Number of owned rows with no stored entries.
  globalIndex numEmptyRows = 0;
};

/**
 * @brief Structural and communication statistics for every rank of a distributed sparse matrix.
 */
struct DistributedMatrixStatistics
{
  /// Statistics in MPI-rank order. The complete array is available on every rank.
  std::vector< MatrixRankStatistics > ranks;
};

/**
 * @brief Compute exact distribution and sparse matrix-vector communication statistics.
 * @param[in] localMatrix local rows with global column indices
 * @param[in] firstLocalRow first global row owned by this rank
 * @param[in] comm communicator over which rows and columns are distributed
 * @return statistics for every rank, replicated on every rank
 *
 * Each distinct off-rank column is one vector value that must be received for a sparse
 * matrix-vector product. The corresponding send counts are obtained by transposing the
 * rank-to-owner request counts. This measures the communication implied by the assembled
 * matrix pattern without depending on a particular linear algebra backend.
 */
DistributedMatrixStatistics
computeDistributedMatrixStatistics( CRSMatrixView< real64 const, globalIndex const > const & localMatrix,
                                    globalIndex firstLocalRow,
                                    MPI_Comm comm );

/**
 * @brief Write one MatrixMarket file per rank and a JSON ownership manifest.
 * @param[in] localMatrix local rows with global column indices
 * @param[in] firstLocalRow first global row owned by this rank
 * @param[in] numGlobalRows global matrix row count
 * @param[in] numGlobalColumns global matrix column count
 * @param[in] filenamePrefix common output prefix
 * @param[in] comm matrix communicator
 *
 * Each piece is a valid global-shape MatrixMarket coordinate matrix containing
 * only the rows owned by its rank. The manifest records the zero-based owned
 * row interval and local entry count for reconstruction and partition-aware
 * analysis.
 */
void writeDistributedMatrixMarket( CRSMatrixView< real64 const, globalIndex const > const & localMatrix,
                                   globalIndex firstLocalRow,
                                   globalIndex numGlobalRows,
                                   globalIndex numGlobalColumns,
                                   string const & filenamePrefix,
                                   MPI_Comm comm );

} // namespace geos

#endif // GEOS_LINEARALGEBRA_UTILITIES_DISTRIBUTEDMATRIXSTATISTICS_HPP_
