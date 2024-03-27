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
 * @file ElasticVTIFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICVTIFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICVTIFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace elasticvtifields
{
  DECLARE_FIELD( Delta,
               "delta",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Delta thomsen anisotropy parameter" );

  DECLARE_FIELD( Epsilon,
               "epsilon",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Epsilon thomsen anisotropy parameter" );

  DECLARE_FIELD( Gamma,
               "gamma",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Gamma thomsen anisotropy parameter" );

}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_ELASTICVTIFIELDS */
