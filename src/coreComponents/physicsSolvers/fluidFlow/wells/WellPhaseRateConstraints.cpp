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

/*
 * @file WellConstraint.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellPhaseRateConstraints.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;


template< typename WellConstraintType >
PhaseConstraint< WellConstraintType >::PhaseConstraint( string const & name, Group * const parent )
  : WellConstraintType( name, parent )
{
  this->setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  this->registerWrapper( viewKeyStruct::phaseRateString(), &this->m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Phase rate,  (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s]) " );

  this->registerWrapper( viewKeyStruct::phaseNameString(), &this->m_phaseName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setDefaultValue( "" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Name of the target phase" );
}

template< typename WellConstraintType >
PhaseConstraint< WellConstraintType >::~PhaseConstraint()
{}

template< typename WellConstraintType >
void PhaseConstraint< WellConstraintType >::postInputInitialization()
{
  // Validate value and table options
  WellConstraintType::postInputInitialization();
}

template< typename WellConstraintType >
bool PhaseConstraint< WellConstraintType >::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return this->isViolated( currentConstraint.phaseVolumeRates()[this->m_phaseIndex], this->getConstraintValue( currentTime ));
}

namespace
{

typedef PhaseConstraint< InjectionConstraint > PhaseInjectionConstraint;
REGISTER_CATALOG_ENTRY( WellConstraintBase, PhaseInjectionConstraint, string const &, Group * const )
typedef PhaseConstraint< ProductionConstraint > PhaseProductionConstraint;
REGISTER_CATALOG_ENTRY( WellConstraintBase, PhaseProductionConstraint, string const &, Group * const )

}
// Explicit template instantiations to ensure constructors are emitted for registration
template class PhaseConstraint< InjectionConstraint >;
template class PhaseConstraint< ProductionConstraint >;

} //namespace geos
