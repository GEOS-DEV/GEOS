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
 * @file NegativeTwoPhaseFlashModelImpl.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODELIMPL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODELIMPL_HPP_

#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

using NegativeTwoPhaseFlashPRPR = NegativeTwoPhaseFlashModel< PengRobinsonEOS, PengRobinsonEOS >;
using NegativeTwoPhaseFlashSRKSRK = NegativeTwoPhaseFlashModel< SoaveRedlichKwongEOS, SoaveRedlichKwongEOS >;

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODELIMPL_HPP_
