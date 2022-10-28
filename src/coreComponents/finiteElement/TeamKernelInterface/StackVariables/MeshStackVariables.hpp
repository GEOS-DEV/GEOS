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
#include "tensor/tensor_types.hpp"

namespace geosx
{

namespace stackVariables
{

// Distributed on threads x & y
template < localIndex num_dofs_mesh_1d, localIndex num_quads_1d, localIndex dim, localIndex batch_size >
struct StackMesh
{
  StackBasis< num_dofs_mesh_1d, num_quads_1d > basis;

  template < localIndex... Sizes >
  using Tensor = tensor::StaticDTensor< Sizes... >; // TODO generalize

  GEOSX_HOST_DEVICE
  StackMesh( LaunchContext & ctx ) : basis( ctx )
  { }

  // Mesh nodes
  Tensor< num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, dim > mesh_nodes;

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, dim > const & getNodes() const
  {
    return mesh_nodes;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, dim > & getNodes()
  {
    return mesh_nodes;
  }

  // Mesh jacobians
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, dim, dim > jacobians;

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, dim, dim > & getJacobians()
  {
    return jacobians;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, dim, dim > const & getJacobians() const
  {
    return jacobians;
  }
};

// Distributed on threads x & y
template < localIndex num_dofs_mesh_1d, localIndex num_quads_1d, localIndex dim, localIndex batch_size >
struct Distributed2DMesh
{
  StackBasis< num_dofs_mesh_1d, num_quads_1d > basis;

  template < localIndex... Sizes >
  using Tensor = tensor::Static2dThreadDTensor< Sizes... >; // TODO generalize

  GEOSX_HOST_DEVICE
  Distributed2DMesh( LaunchContext & ctx ) : basis( ctx )
  { }

  // Mesh nodes
  Tensor< num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, dim > mesh_nodes;

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, dim > const & getNodes() const
  {
    return mesh_nodes;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_mesh_1d, num_dofs_mesh_1d, num_dofs_mesh_1d, dim > & getNodes()
  {
    return mesh_nodes;
  }

  // Mesh jacobians
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, dim, dim > jacobians;

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, dim, dim > & getJacobians()
  {
    return jacobians;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, dim, dim > const & getJacobians() const
  {
    return jacobians;
  }
};

// Distributed on 3D threads
template < localIndex num_dofs_mesh_1d, localIndex num_quads_1d, localIndex dim, localIndex batch_size >
struct Distributed3DMesh
{
  BasisStackVariables< num_dofs_mesh_1d, num_quads_1d > basis;

  GEOSX_HOST_DEVICE
};

template < localIndex num_dofs_mesh_1d, localIndex num_quads_1d, localIndex dim, localIndex batch_size >
struct SharedMesh
{
  SharedBasis< num_dofs_mesh_1d, num_quads_1d > basis;

  GEOSX_HOST_DEVICE
  SharedMesh( LaunchContext & ctx ) : basis( ctx )
  {
    // Mesh nodes
    GEOSX_STATIC_SHARED real64 s_mesh_nodes[batch_size][num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim];
    // Mesh jacobians
    GEOSX_STATIC_SHARED real64 s_jacobians[batch_size][num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim][dim];

    loop<thread_z> (ctx, RAJA::RangeSegment(0, batch_size), [&] (localIndex batch_index) {
      mesh_nodes = &s_mesh_nodes[batch_index];
      jacobians = &s_jacobians[batch_index];
    } );
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
  real64 ( & getJacobians() )[num_quads_1d][num_quads_1d][num_quads_1d][dim][dim]
  {
    return *jacobians;
  }

  GEOSX_HOST_DEVICE
  real64 const ( & getJacobians() const )[num_quads_1d][num_quads_1d][num_quads_1d][dim][dim]
  {
    return *jacobians;
  }
};

// Generic type alias
template < Location location,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex batch_size >
struct Mesh_t;

template < localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex batch_size >
struct Mesh_t< Location::Stack, num_dofs_1d, num_quads_1d, dim, batch_size >
{
  using type = StackMesh< num_dofs_1d, num_quads_1d, dim, batch_size >;
};

template < localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex batch_size >
struct Mesh_t< Location::Shared, num_dofs_1d, num_quads_1d, dim, batch_size >
{
  using type = SharedMesh< num_dofs_1d, num_quads_1d, dim, batch_size >;
};

template < localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex batch_size >
struct Mesh_t< Location::Distributed2D, num_dofs_1d, num_quads_1d, dim, batch_size >
{
  using type = Distributed2DMesh< num_dofs_1d, num_quads_1d, dim, batch_size >;
};


template < Location location,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex batch_size >
using Mesh = typename Mesh_t< location, num_dofs_1d, num_quads_1d, dim, batch_size >::type;

} // namespace stackVariables

} // namespace geosx

#endif // GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_MESH_HPP_
