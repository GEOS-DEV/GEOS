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
 * @file SparsityPatternUtilities.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_SPARSITYPATTERNUTILITIES_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_SPARSITYPATTERNUTILITIES_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @brief Append all entries from one sparsity pattern to another.
 *
 * Coupled solvers build a diagonal pattern from the DofManager and then widen
 * it with their coupling entries. Keeping this small operation in one place
 * avoids subtly different row-copy loops in every such solver.
 *
 * Entries already present in @p target are discarded, so the result is the union
 * of the two patterns. @p target must have been sized for at least that union:
 * a SparsityPattern grows on demand, but growing it row by row is quadratic.
 *
 * @param target the pattern to append to
 * @param source the pattern whose entries are appended
 */
inline void appendSparsityPattern( SparsityPattern< globalIndex > & target,
                                   SparsityPattern< globalIndex > const & source )
{
  GEOS_ERROR_IF_NE( target.numRows(), source.numRows() );
  GEOS_ERROR_IF_NE( target.numColumns(), source.numColumns() );
  for( localIndex row = 0; row < source.numRows(); ++row )
  {
    globalIndex const * const columns = source.getColumns( row ).dataIfContiguous();
    target.insertNonZeros( row, columns, columns + source.numNonZeros( row ) );
  }
}

} // namespace geos

#endif // GEOS_LINEARALGEBRA_UTILITIES_SPARSITYPATTERNUTILITIES_HPP_
