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
 * @file capillaryPressureSelector.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSURESELECTOR_HPP
#define GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSURESELECTOR_HPP

#include "constitutive/capillaryPressure/BrooksCoreyCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/VanGenuchtenCapillaryPressure.hpp"

namespace geosx
{

namespace constitutive
{

#ifndef PASSTHROUGH_HANDLE_CASE
#define PASSTHROUGH_HANDLE_CASE( MODEL ) \
  if( dynamicCast< MODEL * >( &capPres ) ) \
  { \
    lambda( static_cast< MODEL & >( capPres ) ); \
    return true; \
  }

template< typename LAMBDA >
bool constitutiveUpdatePassThru( CapillaryPressureBase const & capPres,
                                 LAMBDA && lambda )
{
  PASSTHROUGH_HANDLE_CASE( BrooksCoreyCapillaryPressure const )
  PASSTHROUGH_HANDLE_CASE( VanGenuchtenCapillaryPressure const )
  return false;
}

template< typename LAMBDA >
bool constitutiveUpdatePassThru( CapillaryPressureBase & capPres,
                                 LAMBDA && lambda )
{
  PASSTHROUGH_HANDLE_CASE( BrooksCoreyCapillaryPressure )
  PASSTHROUGH_HANDLE_CASE( VanGenuchtenCapillaryPressure )
  return false;
}

#undef PASSTHROUGH_HANDLE_CASE
#endif

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSURESELECTOR_HPP
