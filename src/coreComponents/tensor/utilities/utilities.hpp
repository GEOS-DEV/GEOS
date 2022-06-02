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
 * @file utilities.hpp
 */

#ifndef GEOSX_TENSOR_UTILITIES
#define GEOSX_TENSOR_UTILITIES

/**
 * Utilities regroup important tensor abstractions to implement algorithms on
 * tensors.
 * */

/// Helper constants
#include "helper_constants.hpp"
/// Factories for kernel configurations
#include "config.hpp"
/// Function to iterate on tensor dimensions
#include "foreach.hpp"
/// GEOSX_FORALL_CONFIG selecting threading strategy according to the KernelConfig
#include "forall.hpp"
/// TMP patterns to manipulate integer lists at compilation
#include "int_list.hpp"
/// Utility functions, and TMP patterns
#include "util.hpp"

#endif // GEOSX_TENSOR_UTILITIES
