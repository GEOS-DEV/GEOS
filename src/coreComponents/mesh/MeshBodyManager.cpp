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
 * MeshBodyManager.cpp
 *
 *  Created on: Oct 4, 2018
 *      Author: Antoine Mazuyer
 */

#include "MeshBodyManager.hpp"
#include "MeshBody.hpp"

namespace geosx
{

using namespace dataRepository;

MeshBodyManager::MeshBodyManager( const std::string& name,
                                        ManagedGroup * const parent ):
  ManagedGroup( name, parent )
{}

void MeshBodyManager::CreateChild( string const & childKey, string const & childName  ) {
    std::cout << "Using a MeshBody provided by " << childKey << " named " << childName << std::endl;
    std::unique_ptr<MeshBody> meshBody =
        MeshBody::CatalogInterface::Factory( childKey, childName, this );
    this->RegisterGroup<MeshBody>( childKey, std::move(meshBody) );
}

} //namespace
