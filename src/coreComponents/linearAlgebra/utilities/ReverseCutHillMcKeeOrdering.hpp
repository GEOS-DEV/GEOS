/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ReverseCutHillMcKeeOrdering.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_REVERSECUTHILLMCKEEORDERING_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_REVERSECUTHILLMCKEEORDERING_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

class ReverseCutHillMcKeeOrdering
{
public:

/*
 * @brief This function actually does the RCM ordering of a symmetric csr matrix (entire)
          The original implementation can be found in src/parcsr_ls/par_ilu.c in hypre
 * @param[in] offsets row offsets in the matrix to reorder
 * @param[in] columns column indices in the matrix to reorder
 * @param[in] rankOffset offset of this rank (assuming numComps = 1)
 * @param[out] perm the permutation array
 */
  static void
  computePermutation( localIndex const * const offsets,
                      globalIndex const * const columns,
                      localIndex const rankOffset,
                      arrayView1d< localIndex > const perm );

private:

/**
 * @brief This function finds the unvisited node with the minimum degree
 * @param[in] numRows number of rows in the matrix to reorder
 * @param[in] degree the degree array
 * @param[in] marker the marker array for unvisited node
 * @param[out] rootRef reference to the root (= node with the minimum degree)
 */
  static void
  findNodeWithMinDegree( localIndex const numRows,
                         arrayView1d< localIndex > const degree,
                         arrayView1d< localIndex > const marker,
                         localIndex & rootRef );

/*
 * @brief This function find a pseudo-peripheral node start from root
 * @param[in] offsets row offsets in the matrix to reorder
 * @param[in] columns column indices in the matrix to reorder
 * @param[in] numRows number of rows in the matrix to reorder
 * @param[in] rankOffset offset of this rank (assuming numComps = 1)
 * @param[out] rootRef on return will be a end of the pseudo-peripheral
 * @param[out] marker the marker array for unvisited node
 */
  static void
  findPseudoPeripheralNode( localIndex const * const offsets,
                            globalIndex const * const columns,
                            localIndex const numRows,
                            localIndex const rankOffset,
                            localIndex & rootRef,
                            arrayView1d< localIndex > const marker );


/*
 * @brief This function build level structure start from root
 * @param[in] offsets row offsets in the matrix to reorder
 * @param[in] columns column indices in the matrix to reorder
 * @param[in] numRows number of rows in the matrix to reorder
 * @param[in] rankOffset offset of this rank (assuming numComps = 1)
 * @param[in] root pointer to the root
 * @param[out] marker the marker array for unvisited node
 * @param[out] level_i points to the start/end of position on level_j, similar to CSR Matrix
 * @param[out] level_j store node number on each level
 * @param[out] numLevelsRef return the number of level on this level structure
 */
  static void
  buildLevel( localIndex const * const offsets,
              globalIndex const * const columns,
              localIndex const numRows,
              localIndex const rankOffset,
              localIndex const root,
              arrayView1d< localIndex > const marker,
              arrayView1d< localIndex > const level_i,
              arrayView1d< localIndex > const level_j,
              localIndex & numLevelsRef );


/*
 * @brief This qsort is very specialized, not worth to put into utilities
 * Sort a part of array perm based on degree value (ascend)
 * That is, if degree[perm[i]] < degree[perm[j]], we should have i < j
 * @param[in] perm the perm array
 * @param[in] start start in perm
 * @param[in] end end in perm
 * @param[out] degree degree array
 */
  static void
  qsort( arrayView1d< localIndex > const perm,
         localIndex const start,
         localIndex const end,
         arrayView1d< localIndex > const degree );

/*
 * @brief Last step in RCM, reverse it
 * @param[out] perm perm array
 * @param[in] start start position
 * @param[in] end end position
 */
  static void
  reversePermutation( arrayView1d< localIndex > const perm,
                      localIndex const start,
                      localIndex const end );

/*
 * @brief Utility function to swap two values in an array
 * @param[inout] v
 * @param[in] i
 * @param[in] j
 */
  static void
  swap( arrayView1d< localIndex > const v,
        localIndex const i,
        localIndex const j );

/*
 * @brief This function generate numbering for a connected component
 * @param[in] offsets row offsets in the matrix to reorder
 * @param[in] columns column indices in the matrix to reorder
 * @param[in] numRows number of rows in the matrix to reorder
 * @param[in] rankOffset offset of this rank (assuming numComps = 1)
 * @param[in] root pointer to the root
 * @param[out] marker the marker array for unvisited node
 * @param[out] perm permutation array
 * @param[out] currentNumRef number of nodes already have a perm value
 */
  static void
  computeNewOrdering( localIndex const * const offsets,
                      globalIndex const * const columns,
                      localIndex const numRows,
                      localIndex const rankOffset,
                      localIndex const root,
                      arrayView1d< localIndex > const marker,
                      arrayView1d< localIndex > const perm,
                      localIndex & currentNumRef );

};


} // namespace geos

#endif // GEOS_LINEARALGEBRA_UTILITIES_REVERSECUTHILLMCKEEORDERING_HPP_
