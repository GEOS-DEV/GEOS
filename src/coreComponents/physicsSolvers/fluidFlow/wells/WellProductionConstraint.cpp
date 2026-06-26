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
 * @file WellProductionConstraint.cpp
 */

#include "LogLevelsInfo.hpp"
#include "WellProductionConstraint.hpp"
#include "WellConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "functions/FunctionManager.hpp"

#include "WellLiquidRateConstraint.hpp"
#include "WellMassRateConstraint.hpp"
#include "WellPhaseVolumeRateConstraint.hpp"
#include "WellVolumeRateConstraint.hpp"

namespace geos
{

template< typename ConstraintRateType >
ProductionConstraint< ConstraintRateType >::ProductionConstraint( string const & name, Group * const parent )
  : ConstraintRateType( name, parent )
{
  // set rate sign for producers (base class member)
  this->m_rateSign = -1.0;
}
template< typename ConstraintRateType >
ProductionConstraint< ConstraintRateType >::~ProductionConstraint()
{}

template< typename ConstraintRateType >
void ProductionConstraint< ConstraintRateType >::postInputInitialization()
{
  // Validate value and table options
  ConstraintRateType::postInputInitialization();

}
// Register concrete wrapper constraint types and instantiate templates.

template class ProductionConstraint< LiquidRateConstraint >;
using ProductionLiquidRateConstraint = ProductionConstraint< LiquidRateConstraint >;
REGISTER_CATALOG_ENTRY( WellConstraintBase, ProductionLiquidRateConstraint, string const &, Group * const )

template class ProductionConstraint< MassRateConstraint >;
using ProductionMassRateConstraint = ProductionConstraint< MassRateConstraint >;
REGISTER_CATALOG_ENTRY( WellConstraintBase, ProductionMassRateConstraint, string const &, Group * const )

template class ProductionConstraint< PhaseVolumeRateConstraint >;
using ProductionPhaseVolumeRateConstraint = ProductionConstraint< PhaseVolumeRateConstraint >;
REGISTER_CATALOG_ENTRY( WellConstraintBase, ProductionPhaseVolumeRateConstraint, string const &, Group * const )

template class ProductionConstraint< VolumeRateConstraint >;
using ProductionVolumeRateConstraint = ProductionConstraint< VolumeRateConstraint >;
REGISTER_CATALOG_ENTRY( WellConstraintBase, ProductionVolumeRateConstraint, string const &, Group * const )

} //namespace geos
