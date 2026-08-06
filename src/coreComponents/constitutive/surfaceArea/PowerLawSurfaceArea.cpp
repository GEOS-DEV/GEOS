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
 * @file PowerLawSurfaceArea.cpp
 */

#include "constitutive/surfaceArea/PowerLawSurfaceArea.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

PowerLawSurfaceArea::PowerLawSurfaceArea( string const & name, Group * const parent ):
  SurfaceAreaBase( name, parent )
{
  registerWrapper( viewKeyStruct::mineralExponentString(), &m_mineralExponent ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Exponent applied to the ratio of the current to the initial mineral volume fraction" );

  registerWrapper( viewKeyStruct::porosityExponentString(), &m_porosityExponent ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Exponent applied to the ratio of the current to the initial porosity" );
}

void PowerLawSurfaceArea::postInputInitialization()
{
  SurfaceAreaBase::postInputInitialization();

  GEOS_THROW_IF_LT_MSG( m_mineralExponent, 0.0,
                        GEOS_FMT( "{}: invalid value of attribute '{}'",
                                  getFullName(), viewKeyStruct::mineralExponentString() ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( m_porosityExponent, 0.0,
                        GEOS_FMT( "{}: invalid value of attribute '{}'",
                                  getFullName(), viewKeyStruct::porosityExponentString() ),
                        InputError );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PowerLawSurfaceArea, string const &, Group * const )

}
} /* namespace geos */
