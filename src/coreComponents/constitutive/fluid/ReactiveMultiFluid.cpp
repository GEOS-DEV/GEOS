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
 * @file ReactiveMultiFluid.cpp
 */
#include "ReactiveMultiFluid.hpp"
#include "ReactiveMultiFluidExtrinsicData.hpp"


namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


REGISTER_CATALOG_ENTRY( ConstitutiveBase, , string const &, Group * const )
} //namespace constitutive

} //namespace geosx
