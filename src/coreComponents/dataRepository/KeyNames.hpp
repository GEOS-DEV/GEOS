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
 * @file KeyNames.hpp
 */

#ifndef GEOS_DATAREPOSITORY__KEYNAMES_HPP_
#define GEOS_DATAREPOSITORY__KEYNAMES_HPP_

#include <string>

namespace geos
{
namespace dataRepository
{
namespace keys
{

/// @cond DO_NOT_DOCUMENT

static constexpr auto ProblemManager = "Problem";
static constexpr auto cellManager = "cellManager";
static constexpr auto particleManager = "particleManager";

/// @endcond

}
}
}
#endif /* GEOS_DATAREPOSITORY__KEYNAMES_HPP_ */
