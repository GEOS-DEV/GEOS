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
 * @file SubstrateCoverageSurfaceArea.cpp
 */

#include "constitutive/surfaceArea/SubstrateCoverageSurfaceArea.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SubstrateCoverageSurfaceArea::SubstrateCoverageSurfaceArea( string const & name, Group * const parent ):
  SurfaceAreaBase( name, parent )
{
  registerWrapper( viewKeyStruct::characteristicVolumeFractionString(), &m_characteristicVolumeFraction ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Characteristic precipitate volume fraction needed to significantly coat or block "
                    "the substrate, one value per mineral" );

  registerWrapper( viewKeyStruct::criticalPorosityString(), &m_criticalPorosity ).
    setApplyDefaultValue( 1.0e-4 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Porosity below which the rock does not react anymore" );

  registerWrapper( viewKeyStruct::substrateSpecificSurfaceAreaString(), &m_substrateSpecificSurfaceArea ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Specific surface area of the substrate available to precipitation" );

  registerWrapper( viewKeyStruct::mineralExponentString(), &m_mineralExponent ).
    setApplyDefaultValue( 2.0 / 3.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Exponent applied to the ratio of the current to the initial mineral volume fraction" );
}

void SubstrateCoverageSurfaceArea::postInputInitialization()
{
  SurfaceAreaBase::postInputInitialization();

  GEOS_THROW_IF_NE_MSG( m_characteristicVolumeFraction.size(), m_defaultInitialSurfaceArea.size(),
                        GEOS_FMT( "{}: attributes '{}' and '{}' must have the same size",
                                  getFullName(),
                                  viewKeyStruct::characteristicVolumeFractionString(),
                                  viewKeyStruct::defaultInitialSurfaceAreaString() ),
                        InputError );

  for( localIndex r = 0; r < m_characteristicVolumeFraction.size(); ++r )
  {
    GEOS_THROW_IF_LE_MSG( m_characteristicVolumeFraction[r], 0.0,
                          GEOS_FMT( "{}: invalid value of attribute '{}'",
                                    getFullName(), viewKeyStruct::characteristicVolumeFractionString() ),
                          InputError );
  }

  GEOS_THROW_IF_LT_MSG( m_criticalPorosity, 0.0,
                        GEOS_FMT( "{}: invalid value of attribute '{}'",
                                  getFullName(), viewKeyStruct::criticalPorosityString() ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( m_substrateSpecificSurfaceArea, 0.0,
                        GEOS_FMT( "{}: invalid value of attribute '{}'",
                                  getFullName(), viewKeyStruct::substrateSpecificSurfaceAreaString() ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( m_mineralExponent, 0.0,
                        GEOS_FMT( "{}: invalid value of attribute '{}'",
                                  getFullName(), viewKeyStruct::mineralExponentString() ),
                        InputError );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SubstrateCoverageSurfaceArea, string const &, Group * const )

}
} /* namespace geos */
