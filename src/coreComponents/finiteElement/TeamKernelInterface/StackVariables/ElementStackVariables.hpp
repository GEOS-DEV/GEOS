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
 * @file ElementStackVariables.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_ELEMENT_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_ELEMENT_HPP_

#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/BasisStackVariables.hpp"

namespace geosx
{

template < localIndex num_dofs_1d, localIndex num_quads_1d, localIndex dim, localIndex batch_size >
struct ElementStackVariables
{
  BasisStackVariables< num_dofs_1d, num_quads_1d > basis;

  GEOSX_HOST_DEVICE
  ElementStackVariables( LaunchContext & ctx ) : basis( ctx )
  {
    localIndex const batch_index = GEOSX_THREAD_ID(z);
    // Element input dofs of the primary field
    GEOSX_STATIC_SHARED real64 s_dofs_in[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d];
    dofs_in = &s_dofs_in[batch_index];

    // Element primary field gradients at quadrature points
    GEOSX_STATIC_SHARED real64 s_q_gradient_values[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d][dim];
    q_gradient_values = &s_q_gradient_values[batch_index];

    // Element "geometric factors"
    GEOSX_STATIC_SHARED real64 s_Du[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d][dim];
    Du = &s_Du[batch_index];

    // Element contribution to the residual
    GEOSX_STATIC_SHARED real64 s_dofs_out[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d];
    dofs_out = &s_dofs_out[batch_index];
  }

  // Element input dofs of the primary field
  real64 ( * dofs_in )[num_dofs_1d][num_dofs_1d][num_dofs_1d]; // Could be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getDofsIn() const )[num_dofs_1d][num_dofs_1d][num_dofs_1d]
  {
    return *dofs_in;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getDofsIn() )[num_dofs_1d][num_dofs_1d][num_dofs_1d]
  {
    return *dofs_in;
  }

  // Element primary field gradients at quadrature points
  real64 ( * q_gradient_values )[num_quads_1d][num_quads_1d][num_quads_1d][dim]; // Can be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getGradientValues() const )[num_quads_1d][num_quads_1d][num_quads_1d][dim]
  {
    return *q_gradient_values;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getGradientValues() )[num_quads_1d][num_quads_1d][num_quads_1d][dim]
  {
    return *q_gradient_values;
  }

  // Element "geometric factors"
  real64 ( * Du )[num_quads_1d][num_quads_1d][num_quads_1d][dim]; // Could be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getQuadValues() const )[num_quads_1d][num_quads_1d][num_quads_1d][dim]
  {
    return *Du;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getQuadValues() )[num_quads_1d][num_quads_1d][num_quads_1d][dim]
  {
    return *Du;
  }

  // Element contribution to the residual
  real64 ( * dofs_out )[num_dofs_1d][num_dofs_1d][num_dofs_1d]; // Can be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getDofsOut() const )[num_dofs_1d][num_dofs_1d][num_dofs_1d]
  {
    return *dofs_out;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getDofsOut() )[num_dofs_1d][num_dofs_1d][num_dofs_1d]
  {
    return *dofs_out;
  }
};

} // namespace geosx

#endif // GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_ELEMENT_HPP_
