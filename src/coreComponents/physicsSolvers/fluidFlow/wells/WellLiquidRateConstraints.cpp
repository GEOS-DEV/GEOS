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

// *** Liquid Constraint for Production Well  ***************************************************************
LiquidConstraint::LiquidConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

}

LiquidConstraint::~LiquidConstraint()
{}

void LiquidConstraint::postInputInitialization()
{
  // Validate value and table options
  WellConstraintBase::postInputInitialization();
}


// *** Liquid Constraint for Production Well  ***************************************************************
LiquidProductionConstraint::LiquidProductionConstraint( string const & name, Group * const parent )
  : LiquidConstraint( name, parent )
{
  m_rateSign=-1.0;
}

LiquidProductionConstraint::~LiquidProductionConstraint()
{}

void LiquidProductionConstraint::postInputInitialization()
{
  // Validate value and table options
  LiquidConstraint::postInputInitialization();

}

bool LiquidProductionConstraint::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  return currentConstraint.liquidRate() < getConstraintValue( currentTime );
}



} //namespace geos
