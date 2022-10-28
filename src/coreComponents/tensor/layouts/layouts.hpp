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
 * @file layouts.hpp
 */

#ifndef GEOSX_LAYOUTS_HPP_
#define GEOSX_LAYOUTS_HPP_

/**
 * Layouts are simple data structures (often empty) representing the mapping
 * from a rank N index to a linear index. The linear index representing the
 * tensor's associated index in the container for the value.
 * Their main purpose is to differentiate statically/dynamically known sizes
 * for the tensors, and to abstract threading models when utilizing tensors.
 * Layouts also defines the thread block sizes to run a kernel, and handle it
 * automatically for the user.
 * */

/// A statically sized layout
#include "static_layout.hpp"

/// Distributed layouts statically sized
#include "distributed_2d_layout.hpp"

#endif // GEOSX_LAYOUTS_HPP_
