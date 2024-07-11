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
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"

namespace geos
{
using namespace dataRepository;

SourceFluxBoundaryCondition::SourceFluxBoundaryCondition( string const & name, Group * const parent ):
  FieldSpecificationBase( name, parent )
{
  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  addLogLevel( "logLevel >= 1 and first newton iteration", "Print the log message issued by the solver if the boundary condition is called" );
  addLogLevel( "logLevel >= 1 and regions with no dof", "Warnings about non-simulated region intersecting, that can cancel sourceFlux effects" );
  addLogLevel( "logLevel >= 1 and first nonlinear iteration", "Information abonout the Dirichlet pressure, temperature boundary conditions" );
  addLogLevel( "logLevel >= 1 and first nonlinear iteration and this is a thermal simulation", "Information on single phase thermal simulation" );
}

REGISTER_CATALOG_ENTRY( FieldSpecificationBase, SourceFluxBoundaryCondition, string const &, Group * const )

} /* namespace geos */
