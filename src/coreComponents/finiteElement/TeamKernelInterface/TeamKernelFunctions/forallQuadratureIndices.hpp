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
 * @file forallQuadratureIndices.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_FORALL_QUAD_PT_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_FORALL_QUAD_PT_HPP_

#include "finiteElement/TeamKernelInterface/common.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

template < typename Stack, typename Lambda >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void forallQuadratureIndices( Stack & stack, Lambda && lambda )
{
  constexpr localIndex num_quads_1d = Stack::num_quads_1d;
  loop3D( stack, num_quads_1d, num_quads_1d, num_quads_1d,
          [&]( localIndex ind_x, localIndex ind_y, localIndex ind_z )
  {
    TensorIndex quad_index { ind_x, ind_y, ind_z };
    lambda( quad_index );
  } );
}

} // namespace finiteElement
} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_FORALL_QUAD_PT_HPP_ */
