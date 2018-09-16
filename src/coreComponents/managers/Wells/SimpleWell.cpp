/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * @file SimpleWell.cpp
 *
 */

#include "SimpleWell.hpp"

namespace geosx
{

using namespace dataRepository;

SimpleWell::SimpleWell(string const & name, dataRepository::ManagedGroup * const parent)
  : WellBase(name, parent)
{

}

SimpleWell::~SimpleWell()
{

}

void SimpleWell::FillDocumentationNode()
{
  WellBase::FillDocumentationNode();
}

const string SimpleWell::getCatalogName() const
{
  return CatalogName();
}

REGISTER_CATALOG_ENTRY( WellBase, SimpleWell, string const &, ManagedGroup * const )

} //namespace geosx