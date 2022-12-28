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
 * @file TeamKernelInterface.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNEINTERFACE_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNEINTERFACE_HPP_

// The base class for Team kernels.
#include "TeamKernelBase.hpp"

// Main functions to interpolate fields and apply test functions.
#include "TeamKernelFunctions/TeamKernelFunctions.hpp"

// Stack variables meant to be used in Team kernels.
#include "StackVariables/StackVariables.hpp"

// Helper functions for computation at quadrature points.
#include "QuadraturePointKernelFunctions/QuadratureFunctionsHelper.hpp"
#include "QuadraturePointKernelFunctions/qLocalLoad.hpp"
#include "QuadraturePointKernelFunctions/qLocalWrite.hpp"

#endif /* GEOSX_FINITEELEMENT_TEAMKERNEINTERFACE_HPP_ */
