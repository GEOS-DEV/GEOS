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

/*
 * GeometricObjectManager.cpp
 *
 *  Created on: Oct 18, 2017
 *      Author: sherman
 */

#include "GeometricObjectManager.hpp"

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

GeometricObjectManager::GeometricObjectManager( std::string const & name,
                                                ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{}

GeometricObjectManager::~GeometricObjectManager()
{}

ManagedGroup * GeometricObjectManager::CreateChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0("Adding Geometric Object: " << childKey << ", " << childName);
  std::unique_ptr<SimpleGeometricObjectBase> geometriObject = SimpleGeometricObjectBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup<SimpleGeometricObjectBase>( childName, std::move(geometriObject) );
}



} /* namespace geosx */
