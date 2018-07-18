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
 * RestartOutput.cpp
 *
 *  Created on: Jan 31, 2018
 *      Author: sherman
 */

#include "RestartOutput.hpp"
#include "DocumentationNode.hpp"
#include "fileIO/silo/SiloFile.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "physicsSolvers/BoundaryConditions/BoundaryConditionManager.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

RestartOutput::RestartOutput( std::string const & name,
                              ManagedGroup * const parent ):
  OutputBase( name, parent)
{
}

RestartOutput::~RestartOutput()
{}

void RestartOutput::FillDocumentationNode()
{
  OutputBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("Restart");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Outputs Restart format files");

}

void RestartOutput::Execute(real64 const& time_n,
                            real64 const& dt,
                            const int cycleNumber,
                            ManagedGroup * domain)
{
  DomainPartition* domainPartition = ManagedGroup::group_cast<DomainPartition*>(domain);
  
  // #ifdef USE_ATK
  //   char fileName[200] = {0};
  //   sprintf(fileName, "%s_%09d", "restart", cycle);

  //   NewFunctionManager * functionManager = NewFunctionManager::Instance();
  //   BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
    
  //   domainPartition->getParent()->prepareToWrite();
  //   functionManager->prepareToWrite();
  //   bcManager->prepareToWrite();
  //   SidreWrapper::writeTree( 1, fileName, "sidre_hdf5", MPI_COMM_WORLD );
  //   domainPartition->getParent()->finishWriting();
  //   functionManager->finishWriting();
  //   bcManager->finishWriting();
  // #endif

  #ifdef USE_ATK
    ViewWrapper<std::string>::rtype problemName = this->getName();
    char fileName[200] = {0};
    sprintf(fileName, "%s_%s_%09d", problemName.data(), "restart", cycleNumber);

    domainPartition->getParent()->prepareToWrite();
    m_functionManager->prepareToWrite();
    BoundaryConditionManager::get()->prepareToWrite();
    SidreWrapper::writeTree( 1, fileName, "sidre_hdf5", MPI_COMM_WORLD );
    domainPartition->getParent()->finishWriting();
    m_functionManager->finishWriting();
    BoundaryConditionManager::get()->finishWriting();
  #endif
}


REGISTER_CATALOG_ENTRY( OutputBase, RestartOutput, std::string const &, ManagedGroup * const )
} /* namespace geosx */
