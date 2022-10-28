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
 * @file VectorElementStackVariables.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_VECTORELEMENT_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_VECTORELEMENT_HPP_

#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/BasisStackVariables.hpp"
#include "tensor/tensor_types.hpp"

namespace geosx
{

namespace stackVariables
{

template < localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex num_comp,
           localIndex batch_size >
struct StackVectorElement
{
  StackBasis< num_dofs_1d, num_quads_1d > basis;

  template < localIndex... Sizes >
  using Tensor = tensor::StaticDTensor< Sizes... >; // TODO generalize

  GEOSX_HOST_DEVICE
  StackVectorElement( LaunchContext & ctx ) : basis( ctx )
  {
  }

  // Element input dofs of the primary field
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > dofs_in;

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & getDofsIn() const
  {
    return dofs_in;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > & getDofsIn()
  {
    return dofs_in;
  }

  // Element primary field gradients at quadrature points
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > q_gradient_values;

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > const & getGradientValues() const
  {
    return q_gradient_values;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > & getGradientValues()
  {
    return q_gradient_values;
  }

  // Element "geometric factors"
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > Du;

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > const & getQuadValues() const
  {
    return Du;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > & getQuadValues()
  {
    return Du;
  }

  // Element contribution to the residual
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > dofs_out;

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & getDofsOut() const
  {
    return dofs_out;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > & getDofsOut()
  {
    return dofs_out;
  }
};

template < localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex num_comp,
           localIndex batch_size >
struct Distributed2DVectorElement
{
  StackBasis< num_dofs_1d, num_quads_1d > basis;

  template < localIndex... Sizes >
  using Tensor = tensor::Static2dThreadDTensor< Sizes... >; // TODO generalize

  GEOSX_HOST_DEVICE
  Distributed2DVectorElement( LaunchContext & ctx ) : basis( ctx )
  {
  }

  // Element input dofs of the primary field
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > dofs_in;

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & getDofsIn() const
  {
    return dofs_in;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > & getDofsIn()
  {
    return dofs_in;
  }

  // Element primary field gradients at quadrature points
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > q_gradient_values;

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > const & getGradientValues() const
  {
    return q_gradient_values;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > & getGradientValues()
  {
    return q_gradient_values;
  }

  // Element "geometric factors"
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > Du;

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > const & getQuadValues() const
  {
    return Du;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, dim > & getQuadValues()
  {
    return Du;
  }

  // Element contribution to the residual
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > dofs_out;

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & getDofsOut() const
  {
    return dofs_out;
  }

  GEOSX_HOST_DEVICE
  Tensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > & getDofsOut()
  {
    return dofs_out;
  }
};


template < localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex num_comp,
           localIndex batch_size >
struct SharedVectorElement
{
  SharedBasis< num_dofs_1d, num_quads_1d > basis;

  GEOSX_HOST_DEVICE
  SharedVectorElement( LaunchContext & ctx ) : basis( ctx )
  {
    // Element input dofs of the primary field
    GEOSX_STATIC_SHARED real64 s_dofs_in[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d][dim];
    // Element primary field gradients at quadrature points
    GEOSX_STATIC_SHARED real64 s_q_gradient_values[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d][dim][dim];
    // Element "geometric factors"
    GEOSX_STATIC_SHARED real64 s_Du[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d][dim][dim];
    // Element contribution to the residual
    GEOSX_STATIC_SHARED real64 s_dofs_out[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d][dim];

    loop<thread_z> (ctx, RAJA::RangeSegment(0, batch_size), [&] (localIndex batch_index) {
      dofs_in = &s_dofs_in[batch_index];
      q_gradient_values = &s_q_gradient_values[batch_index];
      Du = &s_Du[batch_index];
      dofs_out = &s_dofs_out[batch_index];
    } );
  }

  // Element input dofs of the primary field
  real64 ( * dofs_in )[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp]; // Could be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getDofsIn() const )[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp]
  {
    return *dofs_in;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getDofsIn() )[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp]
  {
    return *dofs_in;
  }

  // Element primary field gradients at quadrature points
  real64 ( * q_gradient_values )[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][dim]; // Can be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getGradientValues() const )[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][dim]
  {
    return *q_gradient_values;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getGradientValues() )[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][dim]
  {
    return *q_gradient_values;
  }

  // Element "geometric factors"
  real64 ( * Du )[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][dim]; // Could be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getQuadValues() const )[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][dim]
  {
    return *Du;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getQuadValues() )[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][dim]
  {
    return *Du;
  }

  // Element contribution to the residual
  real64 ( * dofs_out )[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp]; // Can be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getDofsOut() const )[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp]
  {
    return *dofs_out;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getDofsOut() )[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp]
  {
    return *dofs_out;
  }
};

// Generic type alias
template < Location location,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex num_comp,
           localIndex batch_size >
struct VectorElement_t;

template < localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex num_comp,
           localIndex batch_size >
struct VectorElement_t< Location::Stack, num_dofs_1d, num_quads_1d, dim, num_comp, batch_size >
{
  using type = StackVectorElement< num_dofs_1d, num_quads_1d, dim, num_comp, batch_size >;
};

template < localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex num_comp,
           localIndex batch_size >
struct VectorElement_t< Location::Shared, num_dofs_1d, num_quads_1d, dim, num_comp, batch_size >
{
  using type = SharedVectorElement< num_dofs_1d, num_quads_1d, dim, num_comp, batch_size >;
};

template < localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex num_comp,
           localIndex batch_size >
struct VectorElement_t< Location::Distributed2D, num_dofs_1d, num_quads_1d, dim, num_comp, batch_size >
{
  using type = Distributed2DVectorElement< num_dofs_1d, num_quads_1d, dim, num_comp, batch_size >;
};


template < Location location,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex dim,
           localIndex num_comp,
           localIndex batch_size >
using VectorElement = typename VectorElement_t< location, num_dofs_1d, num_quads_1d, dim, num_comp, batch_size >::type;

} // namespace stackVariables

} // namespace geosx

#endif // GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_VECTORELEMENT_HPP_
