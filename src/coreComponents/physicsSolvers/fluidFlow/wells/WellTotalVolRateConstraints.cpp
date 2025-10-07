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
 * @file WellTotalVolRateConstraints.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellTotalVolRateConstraints.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

template< typename WellConstraintType >
TotalVolConstraint< WellConstraintType >::TotalVolConstraint( string const & name, Group * const parent )
  : WellConstraintType( name, parent )
{
  this->setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  this->registerWrapper( viewKeyStruct::volumeRateString(), &this->m_constraintValue ).
    setDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Volumetric rate (if useSurfaceConditions: [surface m^3/s]; else [reservoir m^3/s])" );

}

template< typename WellConstraintType >
TotalVolConstraint< WellConstraintType >::~TotalVolConstraint()
{}

template< typename WellConstraintType >
void TotalVolConstraint< WellConstraintType >::postInputInitialization()
{
  WellConstraintBase::postInputInitialization();
}

template< typename WellConstraintType >
bool TotalVolConstraint< WellConstraintType >::checkViolation( WellConstraintBase const & currentConstraint, real64 const & currentTime )const
{
  return WellConstraintType::isViolated( currentConstraint.totalVolumeRate(), this->getConstraintValue( currentTime ));
}

typedef TotalVolConstraint< InjectionConstraint > TotalVolInjectionConstraint;
REGISTER_CATALOG_ENTRY( WellConstraintBase, TotalVolInjectionConstraint, string const &, Group * const )
typedef TotalVolConstraint< ProductionConstraint > TotalVolProductionConstraint;
REGISTER_CATALOG_ENTRY( WellConstraintBase, TotalVolProductionConstraint, string const &, Group * const )


// Explicit template instantiations to ensure constructors are emitted for registration
template class TotalVolConstraint< InjectionConstraint >;
template class TotalVolConstraint< ProductionConstraint >;

} //namespace geos
