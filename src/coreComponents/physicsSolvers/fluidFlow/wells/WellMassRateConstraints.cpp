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
 * @file WellMassRateConstraints.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellMassRateConstraints.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

template< typename WellConstraintType >
MassConstraint< WellConstraintType >::MassConstraint( string const & name, Group * const parent )
  : WellConstraintType( name, parent )
{
  this->setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  this->registerWrapper( constraintViewStruct::constraintValueKey::constraintValueString(), &this->m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Maximum mass rate (kg/s)" );
}
template< typename WellConstraintType >
MassConstraint< WellConstraintType >::~MassConstraint()
{}

template< typename WellConstraintType >
void MassConstraint< WellConstraintType >::postInputInitialization()
{
  // Validate value and table options
  WellConstraintBase::postInputInitialization();

}

template< typename WellConstraintType >
bool MassConstraint< WellConstraintType >::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime )const
{
  // isViolated is defined as a static method on the specific WellConstraintType (Injection/Production)
  return WellConstraintType::isViolated( currentConstraint.massRate(), this->getConstraintValue( currentTime ));
}

typedef MassConstraint< InjectionConstraint > MassInjectionConstraint;
REGISTER_CATALOG_ENTRY( WellConstraintBase, MassInjectionConstraint, string const &, Group * const )
typedef MassConstraint< ProductionConstraint > MassProductionConstraint;
REGISTER_CATALOG_ENTRY( WellConstraintBase, MassProductionConstraint, string const &, Group * const )

// Explicit template instantiations to ensure symbols are emitted for the concrete types
template class MassConstraint< InjectionConstraint >;
template class MassConstraint< ProductionConstraint >;

} //namespace geos
