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
 * @file WellBHPConstraints.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellBHPConstraints.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

template< BHPConstraintTypeId T >
BHPConstraint< T >::BHPConstraint( string const & name, Group * const parent )
  : WellConstraintBase( name, parent ),
  m_refElevation( 0.0 ),
  m_refGravCoef( 0.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::refElevString(), &m_refElevation ).
    setDefaultValue( -1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference elevation where BHP control is enforced [m]" );
  if constexpr (T == BHPConstraintTypeId::MAX)
  {
    // override the description for minimum BHP constraint
    registerWrapper( viewKeyStruct::targetBHPString(), &m_constraintValue ).
      setDefaultValue( 0.0 ).
      setInputFlag( InputFlags::OPTIONAL ).
      setRestartFlags( RestartFlags::WRITE_AND_READ ).
      setDescription( "Maximum bottom-hole production pressure [Pa]" );
  }
  else
  {
    registerWrapper( viewKeyStruct::targetBHPString(), &m_constraintValue ).
      setDefaultValue( 0.0 ).
      setInputFlag( InputFlags::OPTIONAL ).
      setRestartFlags( RestartFlags::WRITE_AND_READ ).
      setDescription( "Minimum bottom-hole production pressure [Pa]" );
  }
}

template< BHPConstraintTypeId T >
BHPConstraint< T >::~BHPConstraint()
{}

template< BHPConstraintTypeId T >
void BHPConstraint< T >::postInputInitialization()
{
  WellConstraintBase::postInputInitialization();
}


template< BHPConstraintTypeId T >
bool BHPConstraint< T >::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime ) const
{
  if constexpr (T == BHPConstraintTypeId::MAX)
  {
    return currentConstraint.bottomHolePressure() > getConstraintValue( currentTime );
  }
  else
  {
    return currentConstraint.bottomHolePressure() < getConstraintValue( currentTime );
  }
}

template class BHPConstraint< BHPConstraintTypeId::MIN >;
template class BHPConstraint< BHPConstraintTypeId::MAX >;

namespace
{
REGISTER_CATALOG_ENTRY( WellConstraintBase, MinimumBHPConstraint, string const &, Group * const )
REGISTER_CATALOG_ENTRY( WellConstraintBase, MaximumBHPConstraint, string const &, Group * const )
}

} //namespace geos
