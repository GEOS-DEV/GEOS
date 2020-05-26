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
 * @file multiFluidSelector.hpp
 */
#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_

#include "constitutive/fluid/MultiFluidPVTPackageWrapper.hpp"
#include "constitutive/fluid/MultiPhaseMultiComponentFluid.hpp"

namespace geosx
{

namespace constitutive
{

#ifndef PASSTHROUGH_HANDLE_CASE
#define PASSTHROUGH_HANDLE_CASE( MODEL ) \
  if( dynamicCast< MODEL * >( &fluid ) ) \
  { \
    lambda( static_cast< MODEL & >( fluid ) ); \
    return true; \
  }

template< typename LAMBDA >
bool constitutiveUpdatePassThru( MultiFluidBase const & fluid,
                                 LAMBDA && lambda )
{
  PASSTHROUGH_HANDLE_CASE( MultiFluidPVTPackageWrapper const )
  PASSTHROUGH_HANDLE_CASE( MultiPhaseMultiComponentFluid const )
  return false;
}

template< typename LAMBDA >
bool constitutiveUpdatePassThru( MultiFluidBase & fluid,
                                 LAMBDA && lambda )
{
  PASSTHROUGH_HANDLE_CASE( MultiFluidPVTPackageWrapper )
  PASSTHROUGH_HANDLE_CASE( MultiPhaseMultiComponentFluid )
  return false;
}

#undef PASSTHROUGH_HANDLE_CASE
#endif

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
