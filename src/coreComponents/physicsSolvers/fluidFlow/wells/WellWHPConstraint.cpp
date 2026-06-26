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
#include "WellWHPConstraint.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

WHPConstraint::WHPConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent ),
  m_refElevation( 0.0 ),
  m_refGravCoef( 0.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::targetWHPString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Minimun bottom-hole production pressure [Pa]" );

  registerWrapper( viewKeyStruct::refElevString(), &m_refElevation ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference elevation where WHP control is enforced [m]" );

  registerWrapper( viewKeyStruct::flowTableNameString(), &m_flowTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the well WHP TABLE. \n" );


}


WHPConstraint::~WHPConstraint()
{}

void WHPConstraint::postInputInitialization()
{

  WellConstraintBase::postInputInitialization();

}

MinimumWHPConstraint::MinimumWHPConstraint( string const & name, Group * const parent )
  : WHPConstraint( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::targetWHPString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Minimun well head pressure [Pa]" );
}


MinimumWHPConstraint::~MinimumWHPConstraint()
{}

void MinimumWHPConstraint::postInputInitialization()
{
  WHPConstraint::postInputInitialization();
}

bool MinimumWHPConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return currentConstraint.bottomHolePressure() < getConstraintValue( currentTime );
}

MaximumWHPConstraint::MaximumWHPConstraint( string const & name, Group * const parent )
  : WHPConstraint( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}


MaximumWHPConstraint::~MaximumWHPConstraint()
{}

void MaximumWHPConstraint::postInputInitialization()
{
  // Validate value and table options
  WHPConstraint::postInputInitialization();
}
bool MaximumWHPConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return currentConstraint.bottomHolePressure() > getConstraintValue( currentTime );
}

REGISTER_CATALOG_ENTRY( WellConstraintBase, MinimumWHPConstraint, string const &, Group * const )
REGISTER_CATALOG_ENTRY( WellConstraintBase, MaximumWHPConstraint, string const &, Group * const )

} //namespace geos
