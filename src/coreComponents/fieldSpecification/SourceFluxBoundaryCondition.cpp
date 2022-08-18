/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * SourceFluxBoundaryCondition.cpp
 *
 */

#include "SourceFluxBoundaryCondition.hpp"

namespace geosx
{
using namespace dataRepository;

SourceFluxBoundaryCondition::SourceFluxBoundaryCondition( string const & name, Group * const parent ):
  FieldSpecificationBase( name, parent )
{}

SourceFluxBoundaryCondition::~SourceFluxBoundaryCondition()
{
  // TODO Auto-generated destructor stub
}

void SourceFluxBoundaryCondition::initializePreSubGroups()
{
  this->setFieldName( viewKeyStruct::fluxBoundaryConditionString() );

  FieldSpecificationBase::initializePreSubGroups();
}


REGISTER_CATALOG_ENTRY( FieldSpecificationBase, SourceFluxBoundaryCondition, string const &, Group * const )

} /* namespace geosx */
