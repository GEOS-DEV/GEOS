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
using namespace hpcReact::geochemistry;
using namespace hpcReact::ChainGeneric;
using namespace hpcReact::MoMasBenchmark;

template< typename BASE >
ReactiveSinglePhaseFluid< BASE >::
ReactiveSinglePhaseFluid( string const & name, Group * const parent ):
  BASE( name, parent )
{

  this->registerWrapper( viewKeyStruct::chemicalSystemNameString(), &m_chemicalSystemType ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Chemical System type. Available options are: "
                    "``" + EnumStrings< ChemicalSystemType >::concat( "|" ) + "``" );

  this->registerWrapper( viewKeyStruct::activityModelNameString(), &m_activityModelType ).
    setApplyDefaultValue( ActivityModelType::identity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Activity model applied to the chemical system. Available options are: "
                    "``" + EnumStrings< ActivityModelType >::concat( "|" ) + "``. "
                    "``bdot`` requires ion size and b-dot parameters, which only the geochemical "
                    "systems carry." );

  this->registerWrapper( viewKeyStruct::solventMassPerSolutionVolumeString(), &m_solventMassPerSolutionVolume ).
    setApplyDefaultValue( 1000.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Mass of solvent per unit volume of solution [kg/m^3], used to convert species "
                    "molality [mol/kg solvent] to molarity [mol/m^3 solution]. The default of 1000 "
                    "approximates an aqueous solution by the density of pure water." );

  this->template registerField< fields::reactivefluid::initialPrimarySpeciesConcentration >( &m_initialPrimarySpeciesConcentration );
  this->template registerField< fields::reactivefluid::secondarySpeciesConcentration >( &m_secondarySpeciesConcentration );
  this->template registerField< fields::reactivefluid::primarySpeciesAggregateConcentration >( &m_primarySpeciesAggregateConcentration );
  this->template registerField< fields::reactivefluid::primarySpeciesAggregateConcentration_n >( &m_primarySpeciesAggregateConcentration_n );
  this->template registerField< fields::reactivefluid::primarySpeciesMobileAggregateConcentration >( &m_primarySpeciesMobileAggregateConcentration );
  this->template registerField< fields::reactivefluid::dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations >(
    &m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations );
  this->template registerField< fields::reactivefluid::dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations >(
    &m_dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations );
  this->template registerField< fields::reactivefluid::kineticReactionRates >( &m_kineticReactionRates );
  this->template registerField< fields::reactivefluid::aggregateSpeciesRates >( &m_aggregateSpeciesRates );
  this->template registerField< fields::reactivefluid::dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations >( &m_dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations );
}

template< typename BASE >
std::unique_ptr< ConstitutiveBase > ReactiveSinglePhaseFluid< BASE >::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = BASE::deliverClone( name, parent );

  ReactiveSinglePhaseFluid & newConstitutiveRelation = dynamicCast< ReactiveSinglePhaseFluid & >( *clone );

  newConstitutiveRelation.m_chemicalSystemType = m_chemicalSystemType;
  newConstitutiveRelation.m_activityModelType = m_activityModelType;
  newConstitutiveRelation.m_numPrimarySpecies = m_numPrimarySpecies;
  newConstitutiveRelation.m_numSecondarySpecies = m_numSecondarySpecies;
  newConstitutiveRelation.m_numKineticReactions = m_numKineticReactions;
  newConstitutiveRelation.m_solventMassPerSolutionVolume = m_solventMassPerSolutionVolume;

  return clone;
}

template< typename BASE >
void ReactiveSinglePhaseFluid< BASE >::postInputInitialization()
{
  BASE::postInputInitialization();

  GEOS_THROW_IF( !isSupportedReactionSystem( m_chemicalSystemType, m_activityModelType ),
                 GEOS_FMT( "{}: chemical system '{}' has no '{}' activity model.",
                           this->getDataContext(),
                           EnumStrings< ChemicalSystemType >::toString( m_chemicalSystemType ),
                           EnumStrings< ActivityModelType >::toString( m_activityModelType ) ),
                 InputError );

  // taken from the reaction system itself
  forEachReactionSystem( [&]( auto system )
  {
    using System = decltype( system );
    if( System::chemicalSystem == m_chemicalSystemType && System::activityModel == m_activityModelType )
    {
      using ReactionParams = typename System::ReactionParamsType;
      m_numPrimarySpecies = ReactionParams::numPrimarySpecies();
      m_numSecondarySpecies = ReactionParams::numSecondarySpecies();
      m_numKineticReactions = ReactionParams::numKineticReactions();
    }
  } );

  GEOS_THROW_IF_LE_MSG( m_solventMassPerSolutionVolume, 0.0,
                        GEOS_FMT( "invalid value of attribute '{}'",
                                  viewKeyStruct::solventMassPerSolutionVolumeString() ),
                        InputError, this->getDataContext() );
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
  integer const numPrimarySpecies = this->numPrimarySpecies();
  integer const numSecondarySpecies = this->numSecondarySpecies();
  integer const numKineticReactions = this->numKineticReactions();

  m_initialPrimarySpeciesConcentration.resize( size, numPts, numPrimarySpecies );
  m_secondarySpeciesConcentration.resize( size, numPts, numSecondarySpecies );
  m_primarySpeciesAggregateConcentration.resize( size, numPts, numPrimarySpecies );
  m_primarySpeciesAggregateConcentration_n.resize( size, numPts, numPrimarySpecies );
  m_primarySpeciesMobileAggregateConcentration.resize( size, numPts, numPrimarySpecies );
  m_dPrimarySpeciesAggregateConcentration_dLogPrimarySpeciesConcentrations.resize( size, numPts, numPrimarySpecies, numPrimarySpecies );
  m_dPrimarySpeciesMobileAggregateConcentration_dLogPrimarySpeciesConcentrations.resize( size, numPts, numPrimarySpecies, numPrimarySpecies );
  m_kineticReactionRates.resize( size, numPts, numKineticReactions );
  m_aggregateSpeciesRates.resize( size, numPts, numPrimarySpecies );
  m_dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations.resize( size, numPts, numPrimarySpecies, numPrimarySpecies );
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
