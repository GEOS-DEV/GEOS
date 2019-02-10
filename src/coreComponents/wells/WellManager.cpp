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
 * @file WellManager.cpp
 *
 */

#include "WellManager.hpp"
#include "Well.hpp"

namespace geosx
{

using namespace dataRepository;


WellManager::WellManager(string const & name,
			 dataRepository::ManagedGroup * const parent)
  : ManagedGroup(name, parent)
{


}

WellManager::~WellManager()
{

}

ManagedGroup * WellManager::CreateChild(string const & childKey, string const & childName)
{
  std::cout << "Adding Well: " << childKey << ", " << childName << std::endl;
  std::unique_ptr<Well> well = Well::CatalogInterface::Factory( childKey, childName, this );
  return RegisterGroup<Well>( childName, std::move( well ) );
}

Well * WellManager::getWell(string const & name)
{
  return GetGroup<Well>(name);
}

void WellManager::setGravityVector(R1Tensor const & gravity, bool gravityFlag)
{
  m_gravityVector = gravity;
  m_gravityFlag = gravityFlag;
}

} //namespace geosx
