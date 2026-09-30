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
 * @file VTKRefinementTopology.hpp
 */

#ifndef GEOS_VTK_REFINEMENT_TOPOLOGY_HPP
#define GEOS_VTK_REFINEMENT_TOPOLOGY_HPP

#include <vtkType.h>

#include <array>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

namespace geos::vtk::refinement
{

using Coordinates = std::array< double, 3 >;
using Connectivity = std::vector< vtkIdType >;

enum class EntityKind : std::uint8_t
{
  vertex,
  edge,
  face,
  cell
};

/** Full integral identity; hashes are routing hints, never IDs or equality. */
struct EntityKey
{
  std::uint64_t meshNamespace{};
  EntityKind kind{};
  Connectivity corners;
  bool operator<( EntityKey const & other ) const;
  bool operator==( EntityKey const & other ) const;
};

/** Smallest rotation over both orientations. Rejects repeated corners. */
Connectivity canonicalCycle( Connectivity const & corners );
std::uint64_t stableHash( EntityKey const & key );

struct EntityKeyHash
{
  std::size_t operator()( EntityKey const & key ) const { return stableHash( key ); }
};

using Participants = std::vector< int >;
using Sharing = std::map< EntityKey, Participants >;

/** Build an exact vertex/edge/face key from global corner IDs. */
EntityKey entityKey( EntityKind kind, Connectivity corners, std::uint64_t meshNamespace = 0 );

/** Incidence on the previous level; cell supports are owned volumes, never replicas.
 */
struct EntitySupport
{
  EntityKey key;
  Connectivity localCorners;
  Participants participants;
};

struct ConnectivityHash
{
  std::size_t operator()( Connectivity const & corners ) const;
};

struct PointRecipe
{
  EntityKey key;
  Connectivity support; ///< Local indices of old support points, in canonical ID order.
  Coordinates position;
};

/** Rank-local registry shared by ALL volume cells and their surface traces. */
class PointRegistry
{
public:
  PointRegistry( std::vector< Coordinates > coordinates, Connectivity globalIds, std::uint64_t meshNamespace = 0 );

  vtkIdType edge( vtkIdType a, vtkIdType b );
  vtkIdType face( Connectivity const & corners );
  vtkIdType cell( vtkIdType globalCellId, Connectivity const & corners );
  /** Lookup an already planned point without extending the registry. */
  vtkIdType pointForKey( EntityKey const & key ) const;
  Coordinates const & position( vtkIdType point ) const { return m_points.at( point ).position; }
  std::vector< PointRecipe > const & points() const { return m_points; }
  vtkIdType originalSize() const { return m_globalIds.size(); }

private:
  vtkIdType insert( EntityKey key, Connectivity corners );
  std::uint64_t m_namespace;
  Connectivity m_globalIds;
  std::vector< PointRecipe > m_points;
  std::unordered_map< EntityKey, vtkIdType, EntityKeyHash > m_indices;
  std::unordered_map< Connectivity, Connectivity, ConnectivityHash > m_faceCycles;
};

/** Derive fine sharing from the smallest containing old topological entity.
 * Endpoint participant intersections are insufficient: ranks may meet only at
 * a vertex, without participating in an incident edge or face. This local
 * incidence lookup never runs a distributed directory on fine entities.
 */
class SharingInheritance
{
public:
  SharingInheritance( PointRegistry const & points, std::vector< EntitySupport > entities, int localRank );
  EntitySupport const & support( Connectivity const & fineCorners ) const;
  Participants const & participants( Connectivity const & fineCorners ) const { return support( fineCorners ).participants; }

private:
  std::vector< Connectivity > m_pointSupports;
  std::vector< EntitySupport > m_entities;
  std::vector< std::vector< std::size_t > > m_incident;
};

} // namespace geos::vtk::refinement
#endif
