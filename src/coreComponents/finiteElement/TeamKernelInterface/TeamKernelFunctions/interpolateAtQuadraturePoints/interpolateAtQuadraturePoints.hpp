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
 * @file interpolateAtQuadraturePoints.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_HPP_

// Legacy implementation of 2D distributed interpolateAtQuadraturePoints
#include "legacy.hpp"

// Implementation of non-distributed interpolateAtQuadraturePoints
#include "stack.hpp"

// Implementation of 2D distributed interpolateAtQuadraturePoints
#include "distributed_2d.hpp"

// Implementation of 3D distributed interpolateAtQuadraturePoints
#include "distributed_3d.hpp"

#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

} // namespace finiteElement

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_HPP_ */
