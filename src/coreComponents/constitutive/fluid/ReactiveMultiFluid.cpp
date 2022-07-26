/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ReactiveMultiFluid.cpp
 */
#include "ReactiveMultiFluid.hpp"
#include "ReactiveMultiFluidExtrinsicData.hpp"


namespace geosx
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
  
  registerExtrinsicData( extrinsicMeshData::reactivefluid::primarySpeciesConcentration{}, &m_primarySpeciesConcentration);
  registerExtrinsicData( extrinsicMeshData::reactivefluid::secondarySpeciesConcentration{}, &m_secondarySpeciesConcentration );
  registerExtrinsicData( extrinsicMeshData::reactivefluid::primarySpeciesTotalConcentration{}, &m_primarySpeciesTotalConcentration );
  registerExtrinsicData( extrinsicMeshData::reactivefluid::kineticReactionRates{}, &m_kineticReactionRates );
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

void ReactiveMultiFluid::postProcessInput()
{
  MultiFluidBase::postProcessInput();

  GEOSX_THROW_IF_NE_MSG( numFluidPhases(), 1,
                         GEOSX_FMT( "{}: invalid number of phases", getFullName() ),
                         InputError );

  createChemicalReactions();
}

void ReactiveMultiFluid::resizeFields( localIndex const size, localIndex const numPts )
{ 
  MultiFluidBase::resizeFields(size, numPts);

  integer const numPrimarySpecies = this->numPrimarySpecies();
  integer const numSecondarySpecies = this->numSecondarySpecies();
  integer const numKineticReactions = this->numKineticReactions();

  m_primarySpeciesConcentration.resize( size, numPrimarySpecies );
  m_secondarySpeciesConcentration.resize( size, numSecondarySpecies );
  m_primarySpeciesTotalConcentration.resize( size, numPrimarySpecies );
  m_kineticReactionRates.resize( size, numKineticReactions );
}

void ReactiveMultiFluid::createChemicalReactions()
{
  // instantiate reactions objects
  m_equilibriumReactions = std::make_unique< chemicalReactions::EquilibriumReactions >( getName() + "_equilibriumReactions", m_numPrimarySpecies, m_numSecondarySpecies );
  m_kineticReactions = std::make_unique< chemicalReactions::KineticReactions >( getName() + "_kineticReactions", m_numPrimarySpecies, m_numSecondarySpecies, m_numKineticReactions );
}

ReactiveMultiFluid::KernelWrapper
ReactiveMultiFluid::createKernelWrapper()
{
  return KernelWrapper( m_componentMolarWeight.toViewConst(),
                        m_useMass,
                        isThermal(),
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView(),
                        m_numPrimarySpecies,
                        *m_equilibriumReactions,
                        *m_kineticReactions,
                        m_primarySpeciesConcentration.toView(),
                        m_secondarySpeciesConcentration.toView(),
                        m_primarySpeciesTotalConcentration.toView(),
                        m_kineticReactionRates.toView() );
}

ReactiveMultiFluid::KernelWrapper::
  KernelWrapper( arrayView1d< geosx::real64 const > componentMolarWeight,
                 bool const useMass,
                 bool const isThermal,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseProp::ViewType phaseEnthalpy,
                 PhaseProp::ViewType phaseInternalEnergy,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity,
                 integer const numPrimarySpecies,
                 chemicalReactions::EquilibriumReactions const & equilibriumReactions,
                 chemicalReactions::KineticReactions const & kineticReactions,
                 arrayView2d< real64 > const & primarySpeciesConcentration,
                 arrayView2d< real64 > const & secondarySpeciesConcentration,
                 arrayView2d< real64 > const & primarySpeciesTotalConcentration,
                 arrayView2d< real64 > const & kineticReactionRates )
  : MultiFluidBase::KernelWrapper( std::move( componentMolarWeight ),
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseEnthalpy ),
                                   std::move( phaseInternalEnergy ),
                                   std::move( phaseCompFraction ),
                                   std::move( totalDensity ) ),
  m_numPrimarySpecies( numPrimarySpecies ),
  m_equilibriumReactions( equilibriumReactions.createKernelWrapper() ),
  m_kineticReactions( kineticReactions.createKernelWrapper() ),
  m_primarySpeciesConcentration( primarySpeciesConcentration ),
  m_secondarySpeciesConcentration( secondarySpeciesConcentration ),
  m_primarySpeciesTotalConcentration( primarySpeciesTotalConcentration ),
  m_kineticReactionRates( kineticReactionRates )
{}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ReactiveMultiFluid, string const &, Group * const )
} //namespace constitutive

} //namespace geosx
