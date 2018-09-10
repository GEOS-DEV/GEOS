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

/**
 * @file RestartOutput.cpp
 */

#include "RestartOutput.hpp"
#include "DocumentationNode.hpp"
#include "fileIO/silo/SiloFile.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"


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
                            const integer cycleNumber,
                            real64 const & eventPosition,
                            ManagedGroup * domain)
{
  DomainPartition* domainPartition = ManagedGroup::group_cast<DomainPartition*>(domain);
  
  // #ifdef GEOSX_USE_ATK
  //   char fileName[200] = {0};
  //   sprintf(fileName, "%s_%09d", "restart", cycle);

  //   NewFunctionManager * functionManager = NewFunctionManager::Instance();
  //   BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
    
  //   domainPartition->getParent()->prepareToWrite();
  //   functionManager->prepareToWrite();
  //   bcManager->prepareToWrite();
  //   SidreWrapper::writeTree( 1, fileName, "sidre_hdf5", MPI_COMM_GEOSX );
  //   domainPartition->getParent()->finishWriting();
  //   functionManager->finishWriting();
  //   bcManager->finishWriting();
  // #endif

#ifdef GEOSX_USE_ATK
  integer eventPositionPercent = static_cast<integer>(eventPosition);

  string const problemName = this->getName();
  char fileName[200] = {0};
  sprintf(fileName, "%s_%s__%03d_%09d", problemName.data(), "restart", eventPositionPercent, cycleNumber);

  domainPartition->getParent()->prepareToWrite();
  NewFunctionManager::Instance()->prepareToWrite();
  BoundaryConditionManager::get()->prepareToWrite();
  SidreWrapper::writeTree( 1, fileName, "sidre_hdf5", MPI_COMM_GEOSX );
  domainPartition->getParent()->finishWriting();
  NewFunctionManager::Instance()->finishWriting();
  BoundaryConditionManager::get()->finishWriting();
#endif
}


REGISTER_CATALOG_ENTRY( OutputBase, RestartOutput, std::string const &, ManagedGroup * const )
} /* namespace geosx */
