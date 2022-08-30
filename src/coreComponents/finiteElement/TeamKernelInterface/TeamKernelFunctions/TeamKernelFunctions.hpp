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
 * @file TeamKernelFunction.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_HPP_

/// Functions to load a field from global memory to local memory.
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/readField.hpp"
/// Functions to write a field in local memory to global memory.
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/writeField.hpp"
/// Functions to write-add using atomics a field in local memory to global memory.
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/writeAddField.hpp"


/// Functions to interpolate a field at quadrature points.
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/interpolateAtQuadraturePoints.hpp"
/// Functions to interpolate gradient of a field at quadrature points.
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/interpolateGradientAtQuadraturePoints.hpp"
/// Functions to apply test functions.
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/applyTestFunctions.hpp"
/// Functions to apply gradient of test functions.
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/applyGradientTestFunctions.hpp"
/// Functions to compute the diagonal of a "grad-grad" operator.
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/gradGradDiagonal.hpp"

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_HPP_ */
