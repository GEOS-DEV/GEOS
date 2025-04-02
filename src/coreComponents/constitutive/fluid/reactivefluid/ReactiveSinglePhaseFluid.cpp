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
 * @file ReactiveSinglePhaseFluid.cpp
 */
#include "ReactiveSinglePhaseFluid.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

namespace reactivefluid
{

template< typename BASE >
ReactiveSinglePhaseFluid< BASE >::
  ReactiveSinglePhaseFluid( string const & name, Group * const parent ) :
  BASE( name, parent )
{

  registerWrapper( viewKeyStruct::chemicalSystemNameString(), &m_chemicalSystemType ).
  setInputFlag( InputFlags::REQUIRED ).
  setDescription( "Chemical System type. Available options are: "
                  "``" + EnumStrings< ChemicalSystemType >::concat( "|" ) + "``" );
  // For now this is being hardcoded. We will see where this should come from.
  m_numPrimarySpecies = 7;
  m_numSecondarySpecies = 11;

  registerField( fields::reactivefluid::primarySpeciesConcentration{}, &m_primarySpeciesConcentration );
  registerField( fields::reactivefluid::secondarySpeciesConcentration{}, &m_secondarySpeciesConcentration );
  registerField( fields::reactivefluid::primarySpeciesTotalConcentration{}, &m_primarySpeciesTotalConcentration );
}

template< typename BASE >
std::unique_ptr< ConstitutiveBase > ReactiveSinglePhaseFluid< BASE >::
  deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = BASE::deliverClone( name, parent );

  ReactiveSinglePhaseFluid & newConstitutiveRelation = dynamicCast< ReactiveSinglePhaseFluid & >( *clone );

  return clone;
}

template< typename BASE >
void ReactiveSinglePhaseFluid< BASE >::postInputInitialization()
{}

template< typename BASE >
void ReactiveSinglePhaseFluid< BASE >::resizeFields( localIndex const size, localIndex const numPts )
{
  BASE::resizeFields( size, numPts );

  integer const numPrimarySpecies = this->numPrimarySpecies();
  integer const numSecondarySpecies = this->numSecondarySpecies();
  
  m_primarySpeciesConcentration.resize( size, numPrimarySpecies );
  m_secondarySpeciesConcentration.resize( size, numSecondarySpecies );
  m_primarySpeciesTotalConcentration.resize( size, numPrimarySpecies );
}

} // namespace reactivefluid

} // namespace constitutive

} // namespace geos
