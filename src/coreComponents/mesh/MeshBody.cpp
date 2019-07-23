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
 * @file MeshBody.cpp
 */

#include "MeshBody.hpp"
#include "MeshLevel.hpp"

namespace geosx
{
using namespace dataRepository;

MeshBody::MeshBody( string const & name,
                    ManagedGroup * const parent ):
  ManagedGroup(name,parent),
  m_globalLengthScale(0)
{
  RegisterViewWrapper<integer>( viewKeys.meshLevels );
}

MeshBody::~MeshBody()
{
  // TODO Auto-generated destructor stub
}



MeshLevel * MeshBody::CreateMeshLevel( integer const newLevel )
{
  return this->RegisterGroup<MeshLevel>( "Level0" );
}

void MeshBody::setGlobalLengthScale( real64 scale )
{
  m_globalLengthScale = scale;
}

} /* namespace geosx */
