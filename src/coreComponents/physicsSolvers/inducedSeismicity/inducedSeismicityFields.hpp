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
 * @file inducedSeismicityFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace inducedSeismicity
{

DECLARE_FIELD( t_a,
               "t_a",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Dieterich constitutive relaxation time of seismicity rate" );

DECLARE_FIELD( aSigma,
               "aSigma",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Product of direct effect parameter a and initial meanNormal stress" );

DECLARE_FIELD( pressureRate,
               "pressureRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Pore pressure rate");

DECLARE_FIELD( initialmeanNormalStress,
               "initialmeanNormalStress",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Initial meanNormal stress acting on the fault");

DECLARE_FIELD( meanNormalStress,
               "meanNormalStress",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "meanNormal stress acting on the fault");

DECLARE_FIELD( meanNormalStress_n,
               "meanNormalStress_n",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "meanNormal stress acting on the fault at the previous converged time step");

DECLARE_FIELD( meanNormalStressRate,
               "meanNormalStressRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "meanNormal stress rate acting on the fault");
               
DECLARE_FIELD( initialmeanShearStress,
               "initialmeanShearStress",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Initial meanShear stress acting on the fault");

DECLARE_FIELD( meanShearStress,
               "meanShearStress",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "meanShear stress acting on the fault");

DECLARE_FIELD( meanShearStress_n,
               "meanShearStress_n",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "meanShear stress acting on the fault at the previous converged time step");

DECLARE_FIELD( meanShearStressRate,
               "meanShearStressRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "meanShear stress rate acting on the fault");

DECLARE_FIELD( seismicityRate,
               "seismicityRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Seismicity rate");

DECLARE_FIELD( h,
               "h",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Log of the seismicity rate");

DECLARE_FIELD( h_n,
               "h_n",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Log of the seismicity rate at the previous converged time step");

DECLARE_FIELD( logDenom,
               "logDenom",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Log of the denominator of the integral form of the seismicity rate");

DECLARE_FIELD( logDenom_n,
               "logDenom_n",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Log of the denominator of the integral form of the seismicity rate at the previous converged time step");

}

}

}

#endif // GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_
