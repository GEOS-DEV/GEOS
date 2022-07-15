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
 * @file QuadratureWeightsStackVariables.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_WEIGHTS_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_WEIGHTS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

template < localIndex num_quads_1d >
struct QuadratureWeightsStackVariables
{
  real64 ( * weights )[num_quads_1d];

  QuadratureWeightsStackVariables( LaunchContext & ctx )
  {
    // Initialize quadrature weights
    GEOSX_STATIC_SHARED real64 s_weights[num_quads_1d];
    weights = &s_weights;
    // TODO generalize/use threads
    s_weights[0] = 1.0;
    s_weights[1] = 1.0;
  }

  real64 operator()( localIndex quad_x, localIndex quad_y, localIndex quad_z )
  {
    return (*weights)[ quad_x ] * (*weights)[ quad_y ] * (*weights)[ quad_z ];
  }
};

} // namespace geosx

#endif // GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_WEIGHTS_HPP_
