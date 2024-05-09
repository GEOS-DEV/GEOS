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
 * @file DataTypes.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_DATATYPES_HPP
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_DATATYPES_HPP

#include "Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{
namespace constitutive
{
namespace multifluid
{
template< localIndex MAXSIZE >
using CompositionStackArray = StackArray< real64, 1, MAXSIZE >;

using CompositionSlice = arraySlice1d< real64 >;
using CompositionSliceConst = arraySlice1d< real64 const >;

using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;
using FluidProp = MultiFluidVar< real64, 2, multifluid::LAYOUT_FLUID, multifluid::LAYOUT_FLUID_DC >;


} // namespace multifluid
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_DATATYPES_HPP
