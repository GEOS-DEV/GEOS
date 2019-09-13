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

#include "LASWellGenerator.hpp"

namespace geosx
{
using namespace dataRepository;

LASWellGenerator::LASWellGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent )
{
}

LASWellGenerator::~LASWellGenerator()
{
}

void LASWellGenerator::PostProcessInput()
{
}

Group * LASWellGenerator::CreateChild( string const & childKey, string const & childName )
{
  return nullptr;
}

void LASWellGenerator::GenerateMesh( DomainPartition * const domain )
{
}

} // namespace
