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

#include "ReverseCutHillMcKeeOrdering.hpp"

namespace geos
{

void
ReverseCutHillMcKeeOrdering::
  computePermutation( localIndex const * const offsets,
                      globalIndex const * const columns,
                      localIndex const rankOffset,
                      arrayView1d< localIndex > const perm )
{
  localIndex i = 0;
  localIndex root = 0;
  localIndex currentNum = 0;
  localIndex const numRows = perm.size();

  // get the degree for each node
  array1d< localIndex > degree( numRows );
  array1d< localIndex > marker( numRows );
  for( i = 0; i < numRows; i++ )
  {
    degree[i] = offsets[i + 1] - offsets[i];
    marker[i] = -1;
  }

  // start RCM loop
  currentNum = 0;
  while( currentNum < numRows )
  {
    // Find the unvisited node with the minimum degree
    findNodeWithMinDegree( numRows, degree, marker, root );

    // This is a new connected component
    findPseudoPeripheralNode( offsets, columns, numRows, rankOffset, root, marker );

    // Numbering of this component
    computeNewOrdering( offsets, columns, numRows, rankOffset, root, marker, perm, currentNum );
  }
}

void
ReverseCutHillMcKeeOrdering::
  findNodeWithMinDegree( localIndex const numRows,
                         arrayView1d< localIndex > const degree,
                         arrayView1d< localIndex > const marker,
                         localIndex & rootRef )
{
  localIndex i = 0;
  localIndex minDegree = numRows + 1;
  localIndex root = 0;

  for( i = 0; i < numRows; i++ )
  {
    if( marker[i] < 0 )
    {
      if( degree[i] < minDegree )
      {
        root = i;
        minDegree = degree[i];
      }
    }
  }
  rootRef = root;
}

void
ReverseCutHillMcKeeOrdering::
  findPseudoPeripheralNode( localIndex const * const offsets,
                            globalIndex const * const columns,
                            localIndex const numRows,
                            localIndex const rankOffset,
                            localIndex & rootRef,
                            arrayView1d< localIndex > const marker )
{
  localIndex i = 0;
  localIndex r1 = 0;
  localIndex r2 = 0;
  localIndex row = 0;
  localIndex minDegree = 0;
  localIndex levDegree = 0;
  localIndex numLevels = 0;
  localIndex newNumLevels = 0;
  localIndex root = rootRef;

  // at most numRows levels
  array1d< localIndex > level_i( numRows+1 );
  array1d< localIndex > level_j( numRows );

  // build initial level structure from root
  buildLevel( offsets, columns, numRows, rankOffset, root,
              marker, level_i, level_j, newNumLevels );

  numLevels  = newNumLevels - 1;
  while( numLevels < newNumLevels )
  {
    numLevels = newNumLevels;
    r1 =  level_i[numLevels - 1];
    r2 =  level_i[numLevels];
    minDegree = numRows;
    for( i = r1; i < r2; i++ )
    {
      // select the last level, pick min-degree node
      row = level_j[i];
      levDegree = offsets[row + 1] - offsets[row];
      if( minDegree > levDegree )
      {
        minDegree = levDegree;
        root = row;
      }
    }
    buildLevel( offsets, columns, numRows, rankOffset, root,
                marker, level_i, level_j, newNumLevels );
  }
  rootRef = root;
}

void
ReverseCutHillMcKeeOrdering::
  buildLevel( localIndex const * const offsets,
              globalIndex const * const columns,
              localIndex const numRows,
              localIndex const rankOffset,
              localIndex const root,
              arrayView1d< localIndex > const marker,
              arrayView1d< localIndex > const level_i,
              arrayView1d< localIndex > const level_j,
              localIndex & numLevelsRef )
{
  localIndex i = 0;
  localIndex j = 0;
  localIndex l1 = 0;
  localIndex l2 = 0;
  localIndex lCurrent = 0;
  localIndex r1 = 0;
  localIndex r2 = 0;
  localIndex row_i = 0;
  localIndex row_j = 0;
  localIndex numLevels = 0;

  // set first level first
  level_i[0] = 0;
  level_j[0] = root;
  marker[root] = 0;
  numLevels = 1;
  l1 = 0;
  l2 = 1;
  lCurrent = l2;

  // explore nbhds of all nodes in current level
  while( l2 > l1 )
  {
    level_i[numLevels++] = l2;

    // loop through last level
    for( i = l1; i < l2; i++ )
    {
      // the node to explore
      row_i = level_j[i];
      r1 = offsets[row_i];
      r2 = offsets[row_i + 1];
      for( j = r1; j < r2; j++ )
      {
        row_j = columns[j] - rankOffset;
        if( row_j >= 0 && row_j < numRows )
        {
          if( marker[row_j] < 0 )
          {
            // unmarked row
            marker[row_j] = 0;
            level_j[lCurrent++] = row_j;
          }
        }
      }
    }
    l1 = l2;
    l2 = lCurrent;
  }
  // after this we always have a "ghost" last level
  numLevels--;

  // reset marker
  for( i = 0; i < l2; i++ )
  {
    marker[level_j[i]] = -1;
  }

  numLevelsRef = numLevels;
}

void
ReverseCutHillMcKeeOrdering::
  qsort( arrayView1d< localIndex > const perm,
         localIndex const start,
         localIndex const end,
         arrayView1d< localIndex > const degree )
{
  localIndex i = 0;
  localIndex mid = 0;

  if( start >= end )
  {
    return; // TODO: throw an error here
  }

  swap( perm, start, (start + end) / 2 );
  mid = start;

  // loop to split
  for( i = start + 1; i <= end; i++ )
  {
    if( degree[perm[i]] < degree[perm[start]] )
    {
      swap( perm, ++mid, i );
    }
  }
  swap( perm, start, mid );
  qsort( perm, mid + 1, end, degree );
  qsort( perm, start, mid - 1, degree );
}

void
ReverseCutHillMcKeeOrdering::
  reversePermutation( arrayView1d< localIndex > const perm,
                      localIndex const start,
                      localIndex const end )
{
  localIndex i = 0;
  localIndex j = 0;
  localIndex const mid = (start + end + 1) / 2;

  for( i = start, j = end; i < mid; i++, j-- )
  {
    swap( perm, i, j );
  }
}

void
ReverseCutHillMcKeeOrdering::
  swap( arrayView1d< localIndex > const v,
        localIndex const i,
        localIndex const j )
{
  localIndex const temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}

void
ReverseCutHillMcKeeOrdering::
  computeNewOrdering( localIndex const * const offsets,
                      globalIndex const * const columns,
                      localIndex const numRows,
                      localIndex const rankOffset,
                      localIndex const root,
                      arrayView1d< localIndex > const marker,
                      arrayView1d< localIndex > const perm,
                      localIndex & currentNumRef )
{
  localIndex i = 0;
  localIndex j = 0;
  localIndex l1 = 0;
  localIndex l2 = 0;
  localIndex r1 = 0;
  localIndex r2 = 0;
  localIndex row_i = 0;
  localIndex row_j = 0;
  localIndex rowStart = 0;
  localIndex rowEnd = 0;
  localIndex currentNum = currentNumRef;

  marker[root]       = 0;
  l1                 = currentNum;
  perm[currentNum++] = root;
  l2                 = currentNum;

  while( l2 > l1 )
  {
    // loop through all nodes is current level
    for( i = l1; i < l2; i++ )
    {
      row_i = perm[i];
      r1 = offsets[row_i];
      r2 = offsets[row_i + 1];
      rowStart = currentNum;
      for( j = r1; j < r2; j++ )
      {
        row_j = columns[j] - rankOffset;
        if( row_j >= 0 && row_j < numRows )
        {
          if( marker[row_j] < 0 )
          {
            // save the degree in marker and add it to perm */
            marker[row_j] = offsets[row_j + 1] - offsets[row_j];
            perm[currentNum++] = row_j;
          }
        }
      }
      rowEnd = currentNum;
      qsort( perm, rowStart, rowEnd - 1, marker );
    }
    l1 = l2;
    l2 = currentNum;
  }

  // reverse
  reversePermutation( perm, currentNumRef, currentNum - 1 );
  currentNumRef = currentNum;
}


}
