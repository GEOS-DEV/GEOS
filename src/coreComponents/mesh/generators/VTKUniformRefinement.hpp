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

/** @file VTKUniformRefinement.hpp */
#ifndef GEOS_VTK_UNIFORM_REFINEMENT_HPP
#define GEOS_VTK_UNIFORM_REFINEMENT_HPP

#include "VTKRefinementCommunication.hpp"
#include "VTKRefinementFields.hpp"
#include "mesh/ElementType.hpp"
#include <optional>

class vtkDataSet;

namespace geos::vtk
{
class AllMeshes;

struct RefinementBlockDescriptor
{
  std::string sourceName;
  std::string name;
  ElementType sourceType;
  ElementType type;
  int attribute;
  refinement::Connectivity cells;
};

struct UniformRefinementOptions
{
  std::string regionAttribute = "attribute";
  refinement::TransferPolicies fields;
  std::map< std::string, refinement::TransferPolicies > faceBlockFields;
  /// Production import selects required arrays; unset preserves all fields for component callers.
  std::optional< std::set< std::string > > requiredPointArrays;
  std::optional< std::set< std::string > > requiredCellArrays;
  std::optional< std::set< std::string > > requiredFaceBlockCellArrays;
  std::uint64_t chunkBytes = UINT64_C( 1 ) << 30;
  bool reportStatistics = false;
};

/** Local forecasts. Point/field bounds allow no inter-cell reuse; byte models
 * exclude allocator/runtime overhead and are not a physical-memory budget.
 */
struct RefinementResourceEstimate
{
  std::uint64_t volumeCells{}, surfaceCellCopies{}, pointCopiesUpperBound{}, connectivityEntries{};
  std::uint64_t fieldBytesUpperBound{}, vtkBytesUpperBound{}, geosOwnedConnectivityBytes{};
  double modeledRefinerPeakBytes{}, exchangeBytesUpperBound{};
  /// Available with reporting; models whole datasets from all other ranks.
  std::optional< double > geosGhostConnectivityBytesModel;
};

struct RefinementLevelStatistics
{
  std::uint64_t ownedVolumeCells{}, mainPointCopies{}, ownedMainPoints{}, sharedMainPointCopies{};
  refinement::CommunicationStatistics communication;
};

struct UniformRefinementResult
{
  std::vector< RefinementBlockDescriptor > blocks;
  refinement::Participants neighbors;
  refinement::CommunicationStatistics communication;
  refinement::CommunicationStatistics coarseCommunication;
  std::vector< RefinementResourceEstimate > resources;
  std::vector< RefinementLevelStatistics > levels;
};

/** Refine the coupled, already partitioned datasets transactionally.
 * A zero-level call returns immediately, with no communicator or mesh allocation.
 * All new communication uses comm; coarse volume owners never change.
 */
UniformRefinementResult refineUniformly( AllMeshes & meshes, int levels, UniformRefinementOptions const & options, MPI_Comm comm );

/** Check physical coordinates without applying the import transform twice. */
void validateRefinedTransform( vtkDataSet & mesh, refinement::Coordinates const & translation,
                               refinement::Coordinates const & scale, MPI_Comm comm );
} // namespace geos::vtk
#endif
