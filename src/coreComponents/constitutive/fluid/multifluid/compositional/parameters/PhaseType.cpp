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
 * @file PhaseType.cpp
 */

#include "PhaseType.hpp"
#include "common/format/StringUtilities.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

static stdUnorderedMap< PhaseType, std::string > const phase_aliases{
  {PhaseType::LIQUID, "liquid,liq,oil"},
  {PhaseType::VAPOUR, "gas,vap,vapor,vapour"},
  {PhaseType::AQUEOUS, "wat,water,aqueous"}
};

PhaseType getPhaseTypeFromName( string const & name )
{
  std::string const lowerCaseName = stringutilities::toLower( name );
  for( auto const & it : phase_aliases )
  {
    auto const nameContainer = stringutilities::tokenize( it.second, ",", true, false );
    if( std::find( nameContainer.begin(), nameContainer.end(), lowerCaseName ) != nameContainer.end())
    {
      return it.first;
    }
  }
  // If the name has not been found, throw
  GEOS_THROW( GEOS_FMT( "Phase name {} is not valid for compositional fluid model. "
                        "Allowed values are liquid, vapor and water.", name ), InputError );
  return PhaseType::LIQUID;
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
