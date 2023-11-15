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

/*
 * PVTDriverRunTestCompositionalMultiphaseFluid.cpp
 */

#include "constitutive/fluid/multifluid/PVTDriverRunTest.hpp"
#ifdef GEOSX_USE_PVTPackage
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluidPVTPackage.hpp"
#endif
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"

namespace geos
{
#ifdef GEOSX_USE_PVTPackage
template void PVTDriver::runTest< constitutive::CompositionalMultiphaseFluidPVTPackage >( constitutive::CompositionalMultiphaseFluidPVTPackage &, arrayView2d< real64 > const & );
#endif
template void PVTDriver::runTest< constitutive::CompositionalTwoPhasePengRobinsonConstantViscosity >(
  constitutive::CompositionalTwoPhasePengRobinsonConstantViscosity &, arrayView2d< real64 > const & );
template void PVTDriver::runTest< constitutive::CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity >(
  constitutive::CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity &, arrayView2d< real64 > const & );
}
