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

/** @file VTKRefinementAssociations.hpp */
#ifndef GEOS_VTK_REFINEMENT_ASSOCIATIONS_HPP
#define GEOS_VTK_REFINEMENT_ASSOCIATIONS_HPP

#include "VTKRefinementTopology.hpp"

namespace geos::vtk::refinement
{

// Internal implementation of geos::vtk::refineUniformly.
/// @cond DO_NOT_DOCUMENT
struct MainFace
{
  Connectivity globalCorners;
  Participants owners;
};

struct SurfaceSide
{
  EntityKey mainFace;
  Connectivity mainCornersBySurface; ///< Main GIDs in the surface cell's corner order.
  Participants owners;               ///< Actual volume-face owners, not endpoint allocators.
};

/** Match collocation buckets only against actual incident volume faces.
 * Input faces may include remote records obtained during coarse association
 * discovery. No coordinates or Cartesian products of buckets establish identity.
 */
class SurfaceAssociations
{
public:
  explicit SurfaceAssociations( std::vector< MainFace > faces, std::uint64_t mainNamespace = 0 );
  std::vector< SurfaceSide > match( Connectivity const & surfaceCorners, std::vector< Connectivity > const & collocationBuckets ) const;
  /** Project one actual surface recipe onto one matched side. */
  static EntityKey supportKey( EntityKind kind, Connectivity const & surfaceSupport, Connectivity const & surfaceFace,
                               SurfaceSide const & side );

private:
  std::uint64_t m_namespace;
  std::vector< MainFace > m_faces;
  std::unordered_map< vtkIdType, std::vector< std::size_t > > m_incident;
};
/// @endcond

} // namespace geos::vtk::refinement
#endif
