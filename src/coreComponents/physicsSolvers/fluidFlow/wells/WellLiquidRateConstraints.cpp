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
#include "WellLiquidRateConstraints.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

template< typename WellConstraintType >
LiquidConstraint< WellConstraintType >::LiquidConstraint( string const & name, Group * const parent )
  : WellConstraintType( name, parent )
{
  this->registerWrapper( viewKeyStruct::liquidRateString(), &this->m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Phase rate,  (if useSurfaceCondSitions: [surface m^3/s]; else [reservoir m^3/s]) " );

  this->registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setDefaultValue( "" ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Name of the target phase" );
}
template< typename WellConstraintType >
LiquidConstraint< WellConstraintType >::~LiquidConstraint()
{}
template< typename WellConstraintType >
void LiquidConstraint< WellConstraintType >::postInputInitialization()
{
  // Validate value and table options
  WellConstraintBase::postInputInitialization();
}

template< typename WellConstraintType >
bool LiquidConstraint< WellConstraintType >::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return WellConstraintType::isViolated( currentConstraint.liquidRate(), this->getConstraintValue( currentTime ));
}

namespace
{

typedef LiquidConstraint< InjectionConstraint > LiquidInjectionConstraint;
REGISTER_CATALOG_ENTRY( WellConstraintBase, LiquidInjectionConstraint, string const &, Group * const )
typedef LiquidConstraint< ProductionConstraint > LiquidProductionConstraint;
REGISTER_CATALOG_ENTRY( WellConstraintBase, LiquidProductionConstraint, string const &, Group * const )

}

// Explicit template instantiations to ensure constructors are emitted for registration
template class LiquidConstraint< InjectionConstraint >;
template class LiquidConstraint< ProductionConstraint >;

} //namespace geos
