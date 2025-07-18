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
#include "WellPressureConstraints.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

BHPConstraint::BHPConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent ),
  m_refElevation( 0.0 ),
  m_refGravCoef( 0.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::targetBHPString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Minimun bottom-hole production pressure [Pa]" );

  registerWrapper( viewKeyStruct::refElevString(), &m_refElevation ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference elevation where BHP control is enforced [m]" );

}


BHPConstraint::~BHPConstraint()
{}

void BHPConstraint::postInputInitialization()
{

  WellConstraintBase::postInputInitialization();

}

MinimumBHPConstraint::MinimumBHPConstraint( string const & name, Group * const parent )
  : BHPConstraint( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::targetBHPString(), &m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Minimun bottom-hole production pressure [Pa]" );
}


MinimumBHPConstraint::~MinimumBHPConstraint()
{}

void MinimumBHPConstraint::postInputInitialization()
{

  BHPConstraint::postInputInitialization();

}

bool MinimumBHPConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return currentConstraint.bottomHolePressure()  > getConstraintValue( currentTime );
}

MaximumBHPConstraint::MaximumBHPConstraint( string const & name, Group * const parent )
  : BHPConstraint( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

}


MaximumBHPConstraint::~MaximumBHPConstraint()
{}

void MaximumBHPConstraint::postInputInitialization()
{
  // Validate value and table options
  BHPConstraint::postInputInitialization();

}
bool MaximumBHPConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return currentConstraint.bottomHolePressure() > getConstraintValue( currentTime );
}

} //namespace geos
