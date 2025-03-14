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
 * @file ReactiveSingleFluid.cpp
 */
#include "ReactiveSingleFluid.hpp"

#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/multifluid/reactive/ReactiveMultiFluidFields.hpp"


namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ReactiveSingleFluid::
  ReactiveSingleFluid( string const & name, Group * const parent ):
  SingleFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::primarySpeciesNamesString(), &m_primarySpeciesNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of primary species names" );
    
  // For now this is being hardcoded. We will see where this should come from.
  m_numPrimarySpecies = 3;
  m_numSecondarySpecies = 11;
  m_numKineticReactions = 2;

  registerField( fields::reactivefluid::primarySpeciesConcentration{}, &m_primarySpeciesConcentration );
  registerField( fields::reactivefluid::secondarySpeciesConcentration{}, &m_secondarySpeciesConcentration );
  registerField( fields::reactivefluid::primarySpeciesAggregateConcentration{}, &m_primarySpeciesAggregateConcentration );
  registerField( fields::reactivefluid::primarySpeciesAggregateConcentration_n{}, &m_primarySpeciesAggregateConcentration_n );
  registerField( fields::reactivefluid::dPrimarySpeciesAggregateConcentration_dLogPrimaryConc{}, &m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc );
  registerField( fields::reactivefluid::kineticReactionRates{}, &m_kineticReactionRates );
}

void ReactiveSingleFluid::postInputInitialization()
{
  SingleFluidBase::postInputInitialization();

  // // call to correctly set member array tertiary sizes on the 'main' material object
  // resizeFields( 0, 0 );

  // createChemicalReactions();
}

void ReactiveSingleFluid::saveConvergedState() const
{
  SingleFluidBase::saveConvergedState();

  m_primarySpeciesAggregateConcentration_n.setValues< parallelDevicePolicy<> >( m_primarySpeciesAggregateConcentration.toViewConst() );
}

void ReactiveSingleFluid::resizeFields( localIndex const size, localIndex const numPts )
{
  GEOS_UNUSED_VAR( numPts );
  
  integer const numPrimarySpecies = this->numPrimarySpecies();
  integer const numSecondarySpecies = this->numSecondarySpecies();
  integer const numKineticReactions = this->numKineticReactions();

  m_primarySpeciesConcentration.resize( size, numPrimarySpecies );
  m_secondarySpeciesConcentration.resize( size, numSecondarySpecies );
  m_primarySpeciesAggregateConcentration.resize( size, numPrimarySpecies );
  m_primarySpeciesAggregateConcentration_n.resize( size, numPrimarySpecies );
  m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc.resize( size, numPrimarySpecies, numPrimarySpecies );
  m_kineticReactionRates.resize( size, numKineticReactions );
}

void ReactiveSingleFluid::allocateConstitutiveData( dataRepository::Group & parent,
                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  SingleFluidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );
}

// void ReactiveSingleFluid::createChemicalReactions()
// {
//   // instantiate reactions objects
//   m_equilibriumReactions = std::make_unique< chemicalReactions::EquilibriumReactions >( getName() + "_equilibriumReactions", m_numPrimarySpecies, m_numSecondarySpecies );
//   m_kineticReactions = std::make_unique< chemicalReactions::KineticReactions >( getName() + "_kineticReactions", m_numPrimarySpecies, m_numSecondarySpecies, m_numKineticReactions );
// }

} //namespace constitutive

} //namespace geos
