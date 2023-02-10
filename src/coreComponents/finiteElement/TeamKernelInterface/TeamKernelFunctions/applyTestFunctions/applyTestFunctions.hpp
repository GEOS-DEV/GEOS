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
 * @file applyTestFunctions.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_HPP_

// Legacy implementation of 2D distributed applyTestFunctions
#include "legacy.hpp"

// Implementation of non-distributed applyTestFunctions
#include "stack.hpp"

// Implementation of 2D distributed applyTestFunctions
#include "distributed_2d.hpp"

// Implementation of 3D distributed applyTestFunctions
#include "distributed_3d.hpp"

#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

template < typename StackVariables,
           typename Basis,
           typename QValues,
           typename Dofs >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyTestFunctions( StackVariables & stack,
                         Basis const & basis,
                         QValues const & q_values,
                         Dofs& dofs )
{
  impl::applyTestFunctions( stack,
                            basis.getValuesAtQuadPts(),
                            basis.getGradientValuesAtQuadPts(),
                            q_values,
                            dofs );
}

} // namespace finiteElement

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_HPP_ */
