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

/** @file VTKRefinementSharing.hpp */
#ifndef GEOS_VTK_REFINEMENT_SHARING_HPP
#define GEOS_VTK_REFINEMENT_SHARING_HPP

#include "VTKRefinementTemplates.hpp"

namespace geos::vtk::refinement
{

// Internal implementation of geos::vtk::refineUniformly.
/// @cond DO_NOT_DOCUMENT

struct CoarseBoundary
{
  std::vector< EntitySupport > entities; ///< Exposed volume faces and their edges/vertices.
  std::vector< EntityKey > volumeIds;    ///< ID-only records for duplicate volume ownership checks.
};

/** Local coarse boundary extraction for conforming manifold volumes.
 * Rejects duplicate IDs, crossed local face cycles, equal internal-face
 * orientations and more than two local incident volumes. Global manifold and
 * attachment validation remain the orchestrator's responsibility. ID-only
 * volume records are separate from boundary geometry and must be accounted for.
 */
CoarseBoundary coarseBoundary( std::vector< Cell > const & cells, Connectivity const & globalCellIds, Connectivity const & globalPointIds,
                               int localRank, std::uint64_t meshNamespace = 0 );

/** Sharing inheritance using only shared vertices/edges/faces of the old level.
 * Interior participants are implicit {localRank}. Templates supply valid fine
 * supports; this class does not validate arbitrary disconnected fine entities.
 * The point registry must outlive this object. Shared incidence must be closed:
 * every shared face's edges and every shared edge's vertices are present with
 * participant lists containing the parent list. No MPI or LvArray policy here.
 */
class InterfaceSharing
{
public:
  InterfaceSharing( PointRegistry & points, std::vector< EntitySupport > sharedEntities, int localRank );
  Participants const & participants( Connectivity const & fineCorners ) const;
  /** Build the next shared incidence by subdividing only old interfaces.
   * Templates must already have created all required points. This call cannot
   * extend the point registry after global IDs have been resolved.
   */
  std::vector< EntitySupport > fineSupports( Connectivity const & fineGlobalIds );

private:
  PointRegistry & m_points;
  Participants m_local;
  std::vector< EntitySupport > m_entities;
  std::unordered_map< EntityKey, std::size_t, EntityKeyHash > m_indices;
  std::unordered_map< vtkIdType, std::vector< std::size_t > > m_incident;
};

/// @endcond

} // namespace geos::vtk::refinement
#endif
