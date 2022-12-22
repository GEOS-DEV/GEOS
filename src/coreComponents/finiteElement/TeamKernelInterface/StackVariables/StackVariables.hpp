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
 * @file StackVariables.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_HPP_

// A simple structure representing the shape functions and their derivatives evaluated at quadrature points.
#include "finiteElement/TeamKernelInterface/StackVariables/Basis.hpp"
// A simple structure containing quadrature weights.
#include "finiteElement/TeamKernelInterface/StackVariables/QuadratureWeights.hpp"

// An object containing shared mewmory.
#include "finiteElement/TeamKernelInterface/StackVariables/SharedMem.hpp"
// A set of shared memory buffers.
#include "finiteElement/TeamKernelInterface/StackVariables/SharedMemBuffers.hpp"

#endif
