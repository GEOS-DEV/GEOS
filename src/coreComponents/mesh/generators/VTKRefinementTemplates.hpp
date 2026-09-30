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
 * @file VTKRefinementTemplates.hpp
 */

#ifndef GEOS_VTK_REFINEMENT_TEMPLATES_HPP
#define GEOS_VTK_REFINEMENT_TEMPLATES_HPP

#include "VTKRefinementTopology.hpp"

class vtkCell;

namespace geos::vtk::refinement
{

/** GEOS' VTK import ordering throughout; prismSides identifies a polygonal
 * prism. */
struct Cell
{
  int vtkType;
  Connectivity points;
  int prismSides{};
};

struct Subdivision
{
  std::vector< Cell > children;
  /// Parent face order, each containing its oriented child-face trace.
  std::vector< std::vector< Connectivity > > faceChildren;
};

/** Recognize a closed supported incidence graph, without using GEOS ordering.
 */
Cell normalizeCell( vtkCell & cell );
/** Choose the initial reference frame using original IDs and orientation.
 * Call once on coarse input; descendants inherit their template frames.
 * Reapplying this to fine cells would introduce allocated-ID tie breaking.
 */
Cell normalizeCoarseCell( vtkCell & cell, PointRegistry const & points );
std::vector< Connectivity > cellFaces( Cell const & cell );
/** Complete local incidence; discovery may select exposed boundary entities. */
std::vector< EntitySupport > cellEntitySupports( Cell const & cell, vtkIdType globalCellId, Connectivity const & globalPointIds,
                                                 int localRank, std::uint64_t meshNamespace = 0 );
/** Deduplicate rank-local incidence with full-key hash equality, in expected O(N).
 * Cell IDs must be locally unique. Output order follows first local incidence;
 * it is not an ID-allocation order. No global mesh data or MPI is involved.
 */
std::vector< EntitySupport > meshEntitySupports( std::vector< Cell > const & cells, Connectivity const & globalCellIds,
                                                 Connectivity const & globalPointIds, int localRank, std::uint64_t meshNamespace = 0 );
std::vector< Connectivity > subdivideFace( Connectivity const & corners, PointRegistry & points );
Subdivision subdivideCell( Cell const & cell, vtkIdType globalCellId, PointRegistry & points );

/** Signed measure of first-order maps, with a common bilinear quad fan on
 * polygonal caps. */
double signedMeasure( Cell const & cell, PointRegistry const & points );
void validateGeometry( Cell const & cell, PointRegistry const & points );

struct CellCounts
{
  std::uint64_t hexahedra{}, tetrahedra{}, wedges{}, pyramids{};
  std::array< std::uint64_t, 12 > prisms{};
  std::uint64_t total() const;
  CellCounts next() const;
};

} // namespace geos::vtk::refinement
#endif
