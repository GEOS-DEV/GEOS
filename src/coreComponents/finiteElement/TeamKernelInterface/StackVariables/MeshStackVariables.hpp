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
 * @file MeshStackVariables.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_MESH_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_MESH_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/BasisStackVariables.hpp"

namespace geosx
{

template < localIndex num_dofs_mesh_1d, localIndex num_quads_1d, localIndex dim, localIndex batch_size >
struct MeshStackVariables
{
  BasisStackVariables< num_dofs_mesh_1d, num_quads_1d > basis;

  GEOSX_HOST_DEVICE
  MeshStackVariables( LaunchContext & ctx ) : basis( ctx )
  {
    localIndex const batch_index = GEOSX_THREAD_ID(z);
    // Mesh nodes
    GEOSX_STATIC_SHARED real64 s_mesh_nodes[batch_size][num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim];
    mesh_nodes = &s_mesh_nodes[batch_index];

    // Mesh jacobians
    GEOSX_STATIC_SHARED real64 s_jacobians[batch_size][num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim][dim];
    jacobians = &s_jacobians[batch_index];
  }

  // Mesh nodes
  real64 ( * mesh_nodes )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim]; // Could be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getNodes() const )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim]
  {
    return *mesh_nodes;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getNodes() )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim]
  {
    return *mesh_nodes;
  }

  // Mesh jacobians
  real64 ( * jacobians )[num_quads_1d][num_quads_1d][num_quads_1d][dim][dim]; // Can be in registers

  GEOSX_HOST_DEVICE
  real64 ( & getJacobians() )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim][dim]
  {
    return *jacobians;
  }

  GEOSX_HOST_DEVICE
  real64 const ( & getJacobians() const )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim][dim]
  {
    return *jacobians;
  }
};

} // namespace geosx

#endif // GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_MESH_HPP_
