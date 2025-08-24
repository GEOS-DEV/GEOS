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
 * @file ReactiveMultiFluid.cpp
 */
#include "ReactiveMultiFluid.hpp"
#include "ReactiveMultiFluidFields.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ReactiveMultiFluid::
  ReactiveMultiFluid( string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{
  // For now this is being hardcoded. We will see where this should come from.
  m_numPrimarySpecies = 7;
  m_numSecondarySpecies = 11;
  m_numKineticReactions = 2;
}

bool ReactiveMultiFluid::isThermal() const
{
  return true;
}

std::unique_ptr< ConstitutiveBase > ReactiveMultiFluid::
  deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  ReactiveMultiFluid & newConstitutiveRelation = dynamicCast< ReactiveMultiFluid & >( *clone );

  newConstitutiveRelation.createChemicalReactions();

  return clone;
}

void ReactiveMultiFluid::postInputInitialization()
{
  MultiFluidBase::postInputInitialization();

  GEOS_THROW_IF_NE_MSG( numFluidPhases(), 1,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );

  createChemicalReactions();
}

void ReactiveMultiFluid::allocateConstitutiveData( dataRepository::Group & parent,
                                                   localIndex const numPts )
{
  auto subregion = dynamic_cast< ElementSubRegionBase * >(&parent); // TODO remove

  integer const numPrimarySpecies = this->numPrimarySpecies();
  integer const numSecondarySpecies = this->numSecondarySpecies();
  integer const numKineticReactions = this->numKineticReactions();

  subregion->registerField< fields::reactivefluid::primarySpeciesConcentration >( getName(), &m_primarySpeciesConcentration ).
    reference().resizeDimension< 1 >( numPrimarySpecies );
  subregion->registerField< fields::reactivefluid::secondarySpeciesConcentration >( getName(), &m_secondarySpeciesConcentration ).
    reference().resizeDimension< 1 >( numSecondarySpecies );
  subregion->registerField< fields::reactivefluid::primarySpeciesTotalConcentration >( getName(), &m_primarySpeciesTotalConcentration ).
    reference().resizeDimension< 1 >( numPrimarySpecies );
  subregion->registerField< fields::reactivefluid::kineticReactionRates >( getName(), &m_kineticReactionRates ).
    reference().resizeDimension< 1 >( numKineticReactions );

  MultiFluidBase::allocateConstitutiveData( parent, numPts );
}

void ReactiveMultiFluid::createChemicalReactions()
{
  // instantiate reactions objects
  m_equilibriumReactions = std::make_unique< chemicalReactions::EquilibriumReactions >( getName() + "_equilibriumReactions", m_numPrimarySpecies, m_numSecondarySpecies );
  m_kineticReactions = std::make_unique< chemicalReactions::KineticReactions >( getName() + "_kineticReactions", m_numPrimarySpecies, m_numSecondarySpecies, m_numKineticReactions );
}

} //namespace constitutive

} //namespace geos
