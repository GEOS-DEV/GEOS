/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file singlePhaseSelector.hpp
 */
#ifndef GEOSX_CONSTITUTIVE_FLUID_SINGLEPHASESELECTOR_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_SINGLEPHASESELECTOR_HPP_

#include "constitutive/fluid/CompressibleSinglePhaseFluid.hpp"

namespace geosx
{

namespace constitutive
{

template< typename LAMBDA >
bool constitutiveUpdatePassThru( SingleFluidBase const & fluid,
                                 LAMBDA && lambda )
{
  bool rval = true;
  if( dynamicCast< CompressibleSinglePhaseFluid const * >( &fluid ) )
  {
    lambda( static_cast< CompressibleSinglePhaseFluid const & >( fluid ) );
  }
  else
  {
    rval = false;
  }

  return rval;
}

template< typename LAMBDA >
bool constitutiveUpdatePassThru( SingleFluidBase & fluid,
                                 LAMBDA && lambda )
{
  bool rval = true;
  if( dynamicCast< CompressibleSinglePhaseFluid * >( &fluid ) )
  {
    lambda( static_cast< CompressibleSinglePhaseFluid & >( fluid ) );
  }
  else
  {
    rval = false;
  }

  return rval;
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_SINGLEPHASESELECTOR_HPP_
