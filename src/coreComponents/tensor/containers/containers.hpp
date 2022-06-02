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
 * @file containers.hpp
 */

#ifndef GEOSX_TENSOR_CONTAINERS_HPP_
#define GEOSX_TENSOR_CONTAINERS_HPP_

/**
 *  The container classes represent different types of linear memories
 *  storing the values of a tensor.
 * */

/// A read/write access pointer container
#include "device_container.hpp"
/// A read only access pointer container
#include "read_container.hpp"
/// A statically sized container allocated on the stack
#include "stack_container.hpp"
/// A view container (reference) to another container
#include "view_container.hpp"

#endif // GEOSX_TENSOR_CONTAINERS_HPP_
