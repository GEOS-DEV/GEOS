/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SurfaceAreaFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SURFACEAREA_SURFACEAREAFIELDS_HPP_
#define GEOS_CONSTITUTIVE_SURFACEAREA_SURFACEAREAFIELDS_HPP_

#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace surfacearea
{

using array3dLayoutSpecies = array3d< real64, constitutive::reactivefluid::LAYOUT_SPECIES >;

// NOPLOT: the number of minerals is zero for a purely equilibrium chemical system, and the VTK
// writer rejects a field with no component.
DECLARE_FIELD( surfaceArea,
               "surfaceArea",
               array3dLayoutSpecies,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Reactive surface area of each mineral" );

DECLARE_FIELD( initialSurfaceArea,
               "initialSurfaceArea",
               array3dLayoutSpecies,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Initial reactive surface area of each mineral" );

}

}

}

#endif // GEOS_CONSTITUTIVE_SURFACEAREA_SURFACEAREAFIELDS_HPP_
