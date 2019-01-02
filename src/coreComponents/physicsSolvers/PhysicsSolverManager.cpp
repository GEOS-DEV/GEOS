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
 * PhysicsSolverManager.cpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#include "PhysicsSolverManager.hpp"

#include "SolverBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

PhysicsSolverManager::PhysicsSolverManager( std::string const & name,
                                            ManagedGroup * const parent ):
  ManagedGroup( name, parent),
  m_gravityVector( R1Tensor(0.0) ),
  m_blockSystemRepository()
{
  this->RegisterViewWrapper( viewKeyStruct::gravityVectorString, &m_gravityVector, 0 )->
    setApplyDefaultValue({0,0,0})->
    setInputFlag(InputFlags::OPTIONAL);

  this->RegisterViewWrapper( viewKeyStruct::blockSystemRepositoryString, &m_blockSystemRepository, 0 )->setRestartFlags( RestartFlags::NO_WRITE );
}

PhysicsSolverManager::~PhysicsSolverManager()
{}


ManagedGroup * PhysicsSolverManager::CreateChild( string const & childKey, string const & childName )
{
  ManagedGroup * rval = nullptr;
  if( SolverBase::CatalogInterface::hasKeyName(childKey) )
  {
    GEOS_LOG_RANK_0("Adding Solver of type " << childKey << ", named " << childName);
    rval = RegisterGroup( childName,
                          SolverBase::CatalogInterface::Factory( childKey, childName, this ) );
  }
  return rval;
}



} /* namespace geosx */
