/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file DirichletBoundaryCondition.cpp
 */

#include "DirichletBoundaryCondition.hpp"

namespace geosx
{
using namespace dataRepository;

DirichletBoundaryCondition::DirichletBoundaryCondition( string const & name, Group *const parent ):
  FieldSpecificationBase( name, parent )
{
  // TODO Auto-generated constructor stub

}

DirichletBoundaryCondition::~DirichletBoundaryCondition()
{
  // TODO Auto-generated destructor stub
}



REGISTER_CATALOG_ENTRY( FieldSpecificationBase, DirichletBoundaryCondition, string const &, Group * const )

} /* namespace geosx */
