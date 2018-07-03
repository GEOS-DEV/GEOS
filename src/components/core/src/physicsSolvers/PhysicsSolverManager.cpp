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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
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

#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "SolverBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

PhysicsSolverManager::PhysicsSolverManager( std::string const & name,
                                            ManagedGroup * const parent ):
  ManagedGroup( name, parent),
  m_gravityVector( R1Tensor(0.0) )
{
  this->RegisterViewWrapper( viewKeyStruct::gravityVectorString, &m_gravityVector, 0 );
}

PhysicsSolverManager::~PhysicsSolverManager()
{}


void PhysicsSolverManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("Solvers");
  docNode->setSchemaType("UniqueNode");
  docNode->setShortDescription("Solver manager");

  docNode->AllocateChildNode( viewKeyStruct::gravityVectorString,
                              viewKeyStruct::gravityVectorString,
                              -1,
                              "R1Tensor",
                              "R1Tensor",
                              "Number of Nodes Per Element",
                              "Number of Nodes Per Element",
                              "",
                              "",
                              0,
                              1,
                              0);
}


void PhysicsSolverManager::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Adding Solver: " << childKey << ", " << childName << std::endl;
  std::unique_ptr<SolverBase> solver = SolverBase::CatalogInterface::Factory( childKey, childName, this );
  this->RegisterGroup<SolverBase>( childName, std::move(solver) );
}



} /* namespace geosx */
