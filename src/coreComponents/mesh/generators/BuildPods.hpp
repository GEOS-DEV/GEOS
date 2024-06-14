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

#ifndef GEOS_BUILDPODS_HPP
#define GEOS_BUILDPODS_HPP

#include "Indices.hpp"

#include "include/MeshMappings.hpp"
#include "Pods.hpp"

#include <array>
#include <map>
#include <set>

namespace geos::ghosting
{

/**
 * @brief Complete set of information when a cell refers to a face.
 * @details When a cell refers to a face, w.r.t. the canonical ordering of the face,
 * the nodes of the face (in the cell) can be shifted,
 * and travelled in a different direction (clockwise or counter-clockwise).
 * The @c isFlipped parameter informs about the travel direction.
 * The @c start parameter informs about the shift.
 * @note This class does not refer to how multiple faces are ordered by a cell,
 * but to the precise information when refering to one given face.
 */
struct FaceInfo
{
  /// The global index of the face the cell is referring to.
  FaceGlbIdx index;
  /// Is the face travelled in the same direction as the canonical face.
  bool isFlipped;
  /// Where to start iterating over the nodes of the canonical faces to get back to the original face numbering in the cell.
  std::uint8_t start;
};

/**
 * @brief Complete set of information when a face refers to an edge.
 * @details When a face refers to an edge, w.r.t. the canonical ordering of the edge,
 * there's one node that comes before the other.
 * The @c start parameter (which will be equal to either 0 or 1) will hold this information.
 */
struct EdgeInfo
{
  /// The global index of the edge the face is referring to.
  EdgeGlbIdx index;
  /// Which nodes comes first (and therefore which comes second) in the original edge numbering in the cell.
  std::uint8_t start; // Where the initial index in the canonical face ordering was in the initial face.
};

struct MeshGraph
{
  std::map< CellGlbIdx, std::vector< FaceInfo > > c2f;
  std::map< FaceGlbIdx, std::vector< EdgeInfo > > f2e;
  std::map< EdgeGlbIdx, std::tuple< NodeGlbIdx, NodeGlbIdx > > e2n; // TODO use Edge here?
  std::map< NodeGlbIdx, std::array< double, 3 > > n2pos;
};

struct GhostSend
{
  std::map< NodeGlbIdx, std::set< MpiRank > > nodes;
  std::map< EdgeGlbIdx, std::set< MpiRank > > edges;
  std::map< FaceGlbIdx, std::set< MpiRank > > faces;
  std::map< CellGlbIdx, std::set< MpiRank > > cells;
};

struct GhostRecv
{
  std::map< NodeGlbIdx, MpiRank > nodes;
  std::map< EdgeGlbIdx, MpiRank > edges;
  std::map< FaceGlbIdx, MpiRank > faces;
  std::map< CellGlbIdx, MpiRank > cells;
};

void buildPods( MeshGraph const & owned,
                MeshGraph const & present,
                MeshGraph const & ghosts,
                GhostRecv const & recv,
                GhostSend const & send,
                MeshMappingImpl & meshMappings );

}

#endif //GEOS_BUILDPODS_HPP
