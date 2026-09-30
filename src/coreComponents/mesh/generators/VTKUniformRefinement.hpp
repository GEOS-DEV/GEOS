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
  std::uint64_t chunkBytes = UINT64_C( 1 ) << 30;
};

struct UniformRefinementResult
{
  std::vector< RefinementBlockDescriptor > blocks;
  refinement::Participants neighbors;
  refinement::CommunicationStatistics communication;
};

/** Refine the coupled, already partitioned datasets transactionally.
 * A zero-level call returns immediately, with no communicator or mesh allocation.
 * All new communication uses comm; coarse volume owners never change.
 */
UniformRefinementResult refineUniformly( AllMeshes & meshes, int levels, UniformRefinementOptions const & options, MPI_Comm comm );
} // namespace geos::vtk
#endif
