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
#include "ReactiveFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

namespace reactivefluid
{

using namespace hpcReact::bulkGeneric;

template< typename BASE >
ReactiveSinglePhaseFluid< BASE >::
ReactiveSinglePhaseFluid( string const & name, Group * const parent ):
  BASE( name, parent )
{

  this->registerWrapper( viewKeyStruct::chemicalSystemNameString(), &m_chemicalSystemType ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Chemical System type. Available options are: "
                    "``" + EnumStrings< ChemicalSystemType >::concat( "|" ) + "``" );
  // For now this is being hardcoded. We will see where this should come from.
  m_numPrimarySpecies = 5;
  m_numSecondarySpecies = 11;

  this->registerField( fields::reactivefluid::primarySpeciesConcentration{}, &m_primarySpeciesConcentration );
  this->registerField( fields::reactivefluid::secondarySpeciesConcentration{}, &m_secondarySpeciesConcentration );
  this->registerField( fields::reactivefluid::primarySpeciesAggregateConcentration{}, &m_primarySpeciesAggregateConcentration );
  this->registerField( fields::reactivefluid::primarySpeciesAggregateConcentration_n{}, &m_primarySpeciesAggregateConcentration_n );
  this->registerField( fields::reactivefluid::dPrimarySpeciesAggregateConcentration_dLogPrimaryConc{}, &m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc );
}

template< typename BASE >
std::unique_ptr< ConstitutiveBase > ReactiveSinglePhaseFluid< BASE >::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = BASE::deliverClone( name, parent );

  ReactiveSinglePhaseFluid & newConstitutiveRelation = dynamicCast< ReactiveSinglePhaseFluid & >( *clone );

  GEOS_UNUSED_VAR( newConstitutiveRelation );

  return clone;
}

template< typename BASE >
void ReactiveSinglePhaseFluid< BASE >::postInputInitialization()
{
  BASE::postInputInitialization();
}

template< typename BASE >
void ReactiveSinglePhaseFluid< BASE >::saveConvergedState() const
{
  BASE::saveConvergedState();

  m_primarySpeciesAggregateConcentration_n.setValues< parallelDevicePolicy<> >( m_primarySpeciesAggregateConcentration.toViewConst() );
}

template< typename BASE >
void ReactiveSinglePhaseFluid< BASE >::resizeFields( localIndex const size, localIndex const numPts )
{
  GEOS_UNUSED_VAR( numPts );
  integer const numPrimarySpecies = this->numPrimarySpecies();
  integer const numSecondarySpecies = this->numSecondarySpecies();

  m_primarySpeciesConcentration.resize( size, numPrimarySpecies );
  m_secondarySpeciesConcentration.resize( size, numSecondarySpecies );
  m_primarySpeciesAggregateConcentration.resize( size, numPrimarySpecies );
  m_primarySpeciesAggregateConcentration_n.resize( size, numPrimarySpecies );
  m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc.resize( size, numPrimarySpecies, numPrimarySpecies );
}

template< typename BASE >
void ReactiveSinglePhaseFluid< BASE >::allocateConstitutiveData( dataRepository::Group & parent,
                                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  BASE::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );
}

template class ReactiveSinglePhaseFluid< CompressibleSinglePhaseFluid >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveCompressibleSinglePhaseFluid, string const &, Group * const )

template class ReactiveSinglePhaseFluid< ThermalCompressibleSinglePhaseFluid >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveThermalCompressibleSinglePhaseFluid, string const &, Group * const )

} // namespace reactivefluid

} // namespace constitutive

} // namespace geos
