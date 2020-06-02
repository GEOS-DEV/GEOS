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
 * @file relativePermeabilitySelector.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYSELECTOR_HPP
#define GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYSELECTOR_HPP

#include "constitutive/relativePermeability/BrooksCoreyRelativePermeability.hpp"
#include "constitutive/relativePermeability/BrooksCoreyBakerRelativePermeability.hpp"
#include "constitutive/relativePermeability/VanGenuchtenBakerRelativePermeability.hpp"

namespace geosx
{

namespace constitutive
{

#ifndef PASSTHROUGH_HANDLE_CASE
#define PASSTHROUGH_HANDLE_CASE( MODEL ) \
  if( dynamicCast< MODEL * >( &relPerm ) ) \
  { \
    lambda( static_cast< MODEL & >( relPerm ) ); \
    return true; \
  }

template< typename LAMBDA >
bool constitutiveUpdatePassThru( RelativePermeabilityBase const & relPerm,
                                 LAMBDA && lambda )
{
  PASSTHROUGH_HANDLE_CASE( BrooksCoreyRelativePermeability const )
  PASSTHROUGH_HANDLE_CASE( BrooksCoreyBakerRelativePermeability const )
  PASSTHROUGH_HANDLE_CASE( VanGenuchtenBakerRelativePermeability const )
  return false;
}

template< typename LAMBDA >
bool constitutiveUpdatePassThru( RelativePermeabilityBase & relPerm,
                                 LAMBDA && lambda )
{
  PASSTHROUGH_HANDLE_CASE( BrooksCoreyRelativePermeability )
  PASSTHROUGH_HANDLE_CASE( BrooksCoreyBakerRelativePermeability )
  PASSTHROUGH_HANDLE_CASE( VanGenuchtenBakerRelativePermeability )
  return false;
}

#undef PASSTHROUGH_HANDLE_CASE
#endif

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYSELECTOR_HPP
