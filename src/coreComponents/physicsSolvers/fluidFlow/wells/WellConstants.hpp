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
 * @file WellConstants.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTANTS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTANTS_HPP

namespace geos
{

/**
 * @struct WellConstants
 * @brief A container for constants used in the well solver
 */
struct WellConstants
{
  /// The default BHP for a rate controlled producer when the BHP is not specified
  static constexpr real64 defaultProducerBHP = 1.01325e5;

  /// The default BHP for a rate controlled injector when the BHP is not specified
  static constexpr real64 defaultInjectorBHP = 1.01325e8;
};

} //namespace geos

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTANTS_HPP
