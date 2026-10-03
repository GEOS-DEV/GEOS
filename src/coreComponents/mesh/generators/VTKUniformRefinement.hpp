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

/// Final cell block produced by uniform refinement from one coarse (type, attribute) block.
struct RefinementBlockDescriptor
{
  std::string sourceName;     ///< Name of the coarse cell block the cells descend from.
  std::string name;           ///< Name of the refined cell block.
  ElementType sourceType;     ///< Element type of the coarse cell block.
  ElementType type;           ///< Element type of the refined cells (templates may split one type into several).
  int attribute;              ///< Region attribute value, or -1 when the mesh has no region attribute.
  refinement::Connectivity cells; ///< Local indices of the refined cells in the final mesh.
};

/// Inputs controlling uniform refinement and field transfer.
struct UniformRefinementOptions
{
  std::string regionAttribute = "attribute"; ///< Cell array holding the region attribute.
  refinement::TransferPolicies fields;       ///< Field transfer policies for the main mesh.
  std::map< std::string, refinement::TransferPolicies > faceBlockFields; ///< Field transfer policies per face block.
  /// Production import selects required arrays; unset preserves all fields for component callers.
  std::optional< std::set< std::string > > requiredPointArrays;
  /// Main mesh cell arrays to keep; unset keeps all of them.
  std::optional< std::set< std::string > > requiredCellArrays;
  /// Face block cell arrays to keep; unset keeps all of them.
  std::optional< std::set< std::string > > requiredFaceBlockCellArrays;
  std::uint64_t chunkBytes = UINT64_C( 1 ) << 30; ///< Largest MPI message chunk, in bytes.
  bool reportStatistics = false;                  ///< Log forecasts and per-level statistics.
  /// Component callers can inspect full lineage; production import needs only the root ID.
  bool diagnosticLineage = true;
};

/** Local forecasts. Point/field bounds allow no inter-cell reuse; byte models
 * exclude allocator/runtime overhead and are not a physical-memory budget.
 */
struct RefinementResourceEstimate
{
  std::uint64_t volumeCells{};                ///< Volume cells on this rank.
  std::uint64_t surfaceCellCopies{};          ///< Surface cell copies on this rank.
  std::uint64_t pointCopiesUpperBound{};      ///< Upper bound on point copies on this rank.
  std::uint64_t connectivityEntries{};        ///< Cell connectivity entries on this rank.
  std::uint64_t fieldBytesUpperBound{};       ///< Upper bound on point and cell field bytes.
  std::uint64_t vtkBytesUpperBound{};         ///< Upper bound on the VTK dataset bytes.
  std::uint64_t geosOwnedConnectivityBytes{}; ///< Connectivity bytes GEOS will own after import.
  double modeledRefinerPeakBytes{};           ///< Modeled peak bytes held by the refiner.
  double exchangeBytesUpperBound{};           ///< Upper bound on bytes exchanged with neighbors.
  /// Available with reporting; models whole datasets from all other ranks.
  std::optional< double > geosGhostConnectivityBytesModel;
};

/// Statistics of one refinement level on this rank.
struct RefinementLevelStatistics
{
  std::uint64_t ownedVolumeCells{};      ///< Volume cells owned by this rank.
  std::uint64_t mainPointCopies{};       ///< Main mesh points held by this rank, shared ones included.
  std::uint64_t ownedMainPoints{};       ///< Main mesh points owned by this rank.
  std::uint64_t sharedMainPointCopies{}; ///< Main mesh points this rank shares with other ranks.
  refinement::CommunicationStatistics communication; ///< Communication spent on this level.
};

/// Output of uniform refinement.
struct UniformRefinementResult
{
  std::vector< RefinementBlockDescriptor > blocks;  ///< Refined cell blocks.
  refinement::Participants neighbors;               ///< Ranks sharing mesh entities with this rank.
  refinement::CommunicationStatistics communication; ///< Total communication.
  refinement::CommunicationStatistics coarseCommunication; ///< Communication spent on the coarse mesh.
  std::vector< RefinementResourceEstimate > resources; ///< Forecast per level.
  std::vector< RefinementLevelStatistics > levels;  ///< Statistics per level.
};

/**
 * @brief Refine the coupled, already partitioned datasets transactionally.
 * @details A zero-level call returns immediately, with no communicator or mesh allocation.
 * All new communication uses comm; coarse volume owners never change.
 * @param meshes The main mesh and face blocks; replaced by their refined versions.
 * @param levels Number of uniform refinement levels.
 * @param options Refinement and field transfer options.
 * @param comm The MPI communicator.
 * @return The refined blocks, neighbors and statistics.
 */
UniformRefinementResult refineUniformly( AllMeshes & meshes, int levels, UniformRefinementOptions const & options, MPI_Comm comm );

/**
 * @brief Check the import transform of a refined mesh without applying it.
 * @details The transform must be finite, nonsingular and orientation preserving,
 * and every transformed coordinate must be finite. Check the transformed
 * geometry because rounding and underflow can collapse distinct vertices.
 * @param mesh The refined mesh.
 * @param translation The import translation.
 * @param scale The import scaling.
 * @param comm The MPI communicator.
 */
void validateRefinedTransform( vtkDataSet & mesh, refinement::Coordinates const & translation,
                               refinement::Coordinates const & scale, MPI_Comm comm );
} // namespace geos::vtk
#endif
