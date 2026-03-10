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
 * @file CompositionalMultiphaseFluid.cpp
 */

#include "CompositionalMultiphaseFluid.hpp"

#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PhaseType.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/format/StringUtilities.hpp"

namespace geos
{

namespace constitutive
{

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::
CompositionalMultiphaseFluid( string const & name, Group * const parent )
  : MultiFluidBase( name, parent ),
  m_componentProperties( std::make_unique< compositional::ComponentProperties >( m_componentNames, m_componentMolarWeight ) ),
  m_parameters( createModelParameters() )
{
  using InputFlags = dataRepository::InputFlags;
  using RestartFlags = dataRepository::RestartFlags;

  getWrapperBase( viewKeyStruct::componentNamesString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::componentMolarWeightString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString() ).setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::componentCriticalPressureString(), &m_componentProperties->m_componentCriticalPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical pressures" );

  registerWrapper( viewKeyStruct::componentCriticalTemperatureString(), &m_componentProperties->m_componentCriticalTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical temperatures" );

  registerWrapper( viewKeyStruct::componentAcentricFactorString(), &m_componentProperties->m_componentAcentricFactor ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component acentric factors" );

  registerWrapper( viewKeyStruct::componentVolumeShiftString(), &m_componentProperties->m_componentVolumeShift ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component volume shifts" );

  registerWrapper( viewKeyStruct::componentBinaryCoeffString(), &m_componentProperties->m_componentBinaryCoeff ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Table of binary interaction coefficients" );

  registerField< fields::multifluid::kValues >( &m_kValues );

  // Link parameters specific to each model
  m_parameters->registerParameters( this );

  // Register extra wrappers to enable auto-cloning
  registerWrapper( "componentType", &m_componentProperties->m_componentType )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
integer CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::getWaterPhaseIndex() const
{
  auto const phaseTypes = getPhaseTypes();
  integer const aqueous = static_cast< integer >(compositional::PhaseType::AQUEOUS);
  for( integer ip = 0; ip < numFluidPhases(); ++ip )
  {
    if( phaseTypes[ip] == aqueous )
    {
      return ip;
    }
  }
  return -1;
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
string CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::catalogName()
{
  return GEOS_FMT( "Compositional{}Fluid{}",
                   FLASH::catalogName(),
                   PHASE1::Viscosity::catalogName() );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::initializeState() const
{
  // Zero k-Values to force re-initialisation
  m_kValues.zero();

  MultiFluidBase::initializeState();
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  m_kValues.resize( 0, numPts, numFluidPhases()-1, numFluidComponents() );

  MultiFluidBase::allocateConstitutiveData( parent, numPts );

  // Zero k-Values to force re-initialisation
  m_kValues.zero();
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::postInputInitialization()
{
  MultiFluidBase::postInputInitialization();

  integer const NC = numFluidComponents();
  integer const NP = numFluidPhases();

  GEOS_THROW_IF_NE_MSG( NP, NUM_PHASES,
                        GEOS_FMT( "{}: invalid number of phases in '{}'. There should be {} phases",
                                  getFullName(), viewKeyStruct::phaseNamesString(), NUM_PHASES ),
                        InputError );

  // Phase types should not be repeated
  auto const phaseTypes = getPhaseTypes();
  std::set< integer > uniquePhases;
  for( integer ip = 0; ip < NP; ++ip )
  {
    string const type_name = EnumStrings< compositional::PhaseType >::toString( static_cast< compositional::PhaseType >(phaseTypes[ip]));
    GEOS_THROW_IF ( uniquePhases.find( phaseTypes[ip] ) != uniquePhases.end(),
                    GEOS_FMT( "{}: phase with name {} is of type {} which is repeated. "
                              "Phase types should be unique.", getFullName(), m_phaseNames[ip],
                              type_name ), InputError );
    uniquePhases.insert( phaseTypes[ip] );
  }

  auto const checkInputSize = [&]( auto const & array, integer const expected, string const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( array.size(), expected,
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                          InputError );

  };
  checkInputSize( m_componentProperties->m_componentCriticalPressure, NC, viewKeyStruct::componentCriticalPressureString() );
  checkInputSize( m_componentProperties->m_componentCriticalTemperature, NC, viewKeyStruct::componentCriticalTemperatureString() );
  checkInputSize( m_componentProperties->m_componentAcentricFactor, NC, viewKeyStruct::componentAcentricFactorString() );

  m_componentProperties->classifyComponents();

  if( m_componentProperties->m_componentVolumeShift.empty() )
  {
    m_componentProperties->m_componentVolumeShift.resize( NC );
    m_componentProperties->m_componentVolumeShift.zero();
  }
  checkInputSize( m_componentProperties->m_componentVolumeShift, NC, viewKeyStruct::componentVolumeShiftString() );

  array2d< real64 > & componentBinaryCoeff = m_componentProperties->m_componentBinaryCoeff;
  if( componentBinaryCoeff.empty() )
  {
    componentBinaryCoeff.resize( NC, NC );
    componentBinaryCoeff.zero();
  }
  checkInputSize( componentBinaryCoeff, NC * NC, viewKeyStruct::componentBinaryCoeffString() );

  // Binary interaction coefficients should be symmetric and have zero diagonal
  GEOS_THROW_IF_NE_MSG( componentBinaryCoeff.size( 0 ), NC,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::componentBinaryCoeffString() ),
                        InputError );
  GEOS_THROW_IF_NE_MSG( componentBinaryCoeff.size( 1 ), NC,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::componentBinaryCoeffString() ),
                        InputError );
  for( integer ic = 0; ic < NC; ++ic )
  {
    GEOS_THROW_IF_GT_MSG( LvArray::math::abs( componentBinaryCoeff( ic, ic )), MultiFluidConstants::epsilon,
                          GEOS_FMT( "{}: {} entry at ({},{}) is {}: should be zero", getFullName(), viewKeyStruct::componentBinaryCoeffString(),
                                    ic, ic, componentBinaryCoeff( ic, ic ) ),
                          InputError );
    for( integer jc = ic + 1; jc < NC; ++jc )
    {
      real64 const difference = LvArray::math::abs( componentBinaryCoeff( ic, jc )-componentBinaryCoeff( jc, ic ));
      GEOS_THROW_IF_GT_MSG( difference, MultiFluidConstants::epsilon,
                            GEOS_FMT( "{}: {} entry at ({},{}) is {} and is different from entry at ({},{}) which is {}",
                                      getFullName(), viewKeyStruct::componentBinaryCoeffString(),
                                      ic, jc, componentBinaryCoeff( ic, jc ), jc, ic, componentBinaryCoeff( jc, ic ) ),
                            InputError );
    }
  }

  m_parameters->postInputInitialization( this, *m_componentProperties );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();

  // Create the fluid models
  createModels();
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
std::unique_ptr< ConstitutiveBase >
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::deliverClone( string const & name,
                                                                             Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );
  CompositionalMultiphaseFluid & newFluid = dynamicCast< CompositionalMultiphaseFluid & >( *clone );
  newFluid.createModels();
  return clone;
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
typename CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::KernelWrapper
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::createKernelWrapper()
{
  return KernelWrapper( *m_componentProperties,
                        *m_flash,
                        *m_phase1,
                        *m_phase2,
                        *m_phase3,
                        m_phaseOrder.toViewConst(),
                        m_componentMolarWeight,
                        m_useMass,
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView(),
                        m_kValues.toView() );
}

// Create the fluid models
template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::createModels()
{
  m_phaseType = getPhaseTypes();

  // Determine the phase ordering
  m_phaseOrder.resize( m_phaseType.size() );
  FlashModel::calculatePhaseOrdering( m_phaseType.toViewConst(), m_phaseOrder );

  m_flash = std::make_unique< FLASH >( getName() + '_' + FLASH::catalogName(),
                                       *m_componentProperties,
                                       *m_parameters,
                                       m_phaseType.toViewConst() );

  m_phase1 = std::make_unique< PHASE1 >( GEOS_FMT( "{}_PhaseModel1", getName() ),
                                         *m_componentProperties,
                                         0,
                                         *m_parameters );

  m_phase2 = std::make_unique< PHASE2 >( GEOS_FMT( "{}_PhaseModel2", getName() ),
                                         *m_componentProperties,
                                         1,
                                         *m_parameters );

  m_phase3 = std::make_unique< PHASE3 >( GEOS_FMT( "{}_PhaseModel3", getName() ),
                                         *m_componentProperties,
                                         2,
                                         *m_parameters );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
array1d< integer > CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::getPhaseTypes() const
{
  integer const numPhases = numFluidPhases();
  array1d< integer > phaseTypes( numPhases );
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    phaseTypes[ip] = static_cast< integer >(compositional::getPhaseTypeFromName( m_phaseNames[ip] ));
  }
  return phaseTypes;
}

// Create the fluid models
template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
std::unique_ptr< compositional::ModelParameters >
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::createModelParameters()
{
  std::unique_ptr< compositional::ModelParameters > parameters;
  parameters = FLASH::createParameters( std::move( parameters ));
  parameters = PHASE1::createParameters( std::move( parameters ));
  parameters = PHASE2::createParameters( std::move( parameters ));
  parameters = PHASE3::createParameters( std::move( parameters ));
  return parameters;
}

// Explicit instantiation of the model template.
template class CompositionalMultiphaseFluid<
    compositional::NegativeTwoPhaseFlashModel,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::ConstantViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::ConstantViscosity, compositional::NullModel > >;
template class CompositionalMultiphaseFluid<
    compositional::NegativeTwoPhaseFlashModel,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel > >;
template class CompositionalMultiphaseFluid<
    compositional::NegativeTwoPhaseFlashModel,
    compositional::PhaseModel< compositional::PhillipsBrineDensity, compositional::PhillipsBrineViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel > >;
template class CompositionalMultiphaseFluid<
    compositional::ImmiscibleWaterFlashModel,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::ImmiscibleWaterDensity, compositional::ImmiscibleWaterViscosity, compositional::NullModel > >;
template class CompositionalMultiphaseFluid<
    compositional::KValueFlashModel< 2 >,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel > >;
template class CompositionalMultiphaseFluid<
    compositional::KValueFlashModel< 2 >,
    compositional::PhaseModel< compositional::PhillipsBrineDensity, compositional::PhillipsBrineViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel > >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalTwoPhaseConstantViscosity,
                        string const &,
                        dataRepository::Group * const )

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalTwoPhaseLohrenzBrayClarkViscosity,
                        string const &,
                        dataRepository::Group * const )

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalTwoPhasePhillipsBrine,
                        string const &,
                        dataRepository::Group * const )

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalThreePhaseLohrenzBrayClarkViscosity,
                        string const &,
                        dataRepository::Group * const )

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalKValueLohrenzBrayClarkViscosity,
                        string const &,
                        dataRepository::Group * const )

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalKValuePhillipsBrine,
                        string const &,
                        dataRepository::Group * const )

} // namespace constitutive

} // namespace geos
