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

#ifndef GEOS_NEWGLOBALNUMBERING_HPP
#define GEOS_NEWGLOBALNUMBERING_HPP

#include "Indices.hpp"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <set>

namespace geos::ghosting
{

using Edge = std::tuple< NodeGlbIdx, NodeGlbIdx >;
using Face = std::vector< NodeGlbIdx >;

struct Buckets
{
  std::map< std::set< MpiRank >, std::set< NodeGlbIdx > > nodes;
  std::map< std::set< MpiRank >, std::set< Edge > > edges;
  std::map< std::set< MpiRank >, std::set< Face > > faces;
};

inline void to_json( json & j,
                     const Buckets & v )
{
  j = json{ { "nodes", v.nodes },
            { "edges", v.edges },
            { "faces", v.faces } };
}

struct BucketOffsets
{
  std::map< std::set< MpiRank >, EdgeGlbIdx > edges;
  std::map< std::set< MpiRank >, FaceGlbIdx > faces;
};

inline void to_json( json & j,
                     const BucketOffsets & v )
{
  j = json{ { "edges", v.edges },
            { "faces", v.faces } };
}

inline void from_json( const json & j,
                       BucketOffsets & v )
{
  v.edges = j.at( "edges" ).get< std::map< std::set< MpiRank >, EdgeGlbIdx > >();
  v.faces = j.at( "faces" ).get< std::map< std::set< MpiRank >, FaceGlbIdx > >();
}

/**
 * @brief Order the nodes of the faces in a way that can be reproduced across the MPI ranks.
 * @param nodes The list of nodes as provided by the mesh.
 * @return A face with the nodes in the appropriate order
 * @details The nodes will be ordered in the following way.
 * First, we look for the lowest node index. It will become the first node.
 * Then we must pick the second node. We have to choices: just before or just after the first.
 * (we do not want to shuffle the nodes completely, we need to keep track of the order).
 * To do this, we select the nodes with the lowest index as well.
 * This also defines a direction in which we'll pick the other nodes.
 * For example, the face <tt>[2, 3, 1, 5, 9, 8]</tt> will become <tt>[1, 3, 2, 8, 9, 5]</tt>
 * because we'll start with @c 1 and then select the @c 3 over the @c 5.
 * Which defines the direction @c 2, @c 8, @c 9, @c 5.
 * @note This is the same pattern that we apply for edges.
 * Except that edges having only two nodes, it's not necessary to implement a dedicated function
 * and <tt>std::minmax</tt> is enough.
 */
Face reorderFaceNodes( std::vector< NodeGlbIdx > const & nodes );

std::tuple< Buckets, BucketOffsets > doTheNewGlobalNumbering( vtkSmartPointer< vtkDataSet > mesh,
                                                              std::set< MpiRank > const & neighbors );

} // end of namespace

#endif //GEOS_NEWGLOBALNUMBERING_HPP
