// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * MeshBody.cpp
 *
 *  Created on: Sep 13, 2017
 *      Author: settgast
 */

#include "MeshBody.hpp"
#include "MeshLevel.hpp"

namespace geosx
{
using namespace dataRepository;

MeshBody::MeshBody( string const & name,
                    ManagedGroup * const parent ):
  ManagedGroup(name,parent)
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

} /* namespace geosx */
