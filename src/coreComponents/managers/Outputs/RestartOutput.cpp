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
 * @file RestartOutput.cpp
 */

#include "RestartOutput.hpp"
#include "fileIO/silo/SiloFile.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"


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

void RestartOutput::Execute(real64 const time_n,
                            real64 const dt,
                            integer const cycleNumber,
                            integer const eventCounter,
                            real64 const eventProgress,
                            ManagedGroup * domain)
{
#ifdef GEOSX_USE_ATK
  GEOSX_MARK_FUNCTION;

  DomainPartition* domainPartition = ManagedGroup::group_cast<DomainPartition*>(domain);
  ProblemManager* problemManager = ManagedGroup::group_cast<ProblemManager*>(domainPartition->getParent());

  // Ignoring the eventProgress indicator for now to be compliant with the integrated test repo
  // integer const eventProgressPercent = static_cast<integer const>(eventProgress * 100.0);
  char fileName[200] = {0};
  sprintf(fileName, "%s_%s_%09d", problemManager->getProblemName().c_str(), "restart", cycleNumber);

  problemManager->prepareToWrite();
  NewFunctionManager::Instance()->prepareToWrite();
  FieldSpecificationManager::get()->prepareToWrite();
  int numFiles;
  MPI_Comm_size( MPI_COMM_GEOSX, &numFiles );
  SidreWrapper::writeTree( numFiles, fileName, "sidre_hdf5", MPI_COMM_GEOSX );
  problemManager->finishWriting();
  NewFunctionManager::Instance()->finishWriting();
  FieldSpecificationManager::get()->finishWriting();
#endif
}


REGISTER_CATALOG_ENTRY( OutputBase, RestartOutput, std::string const &, ManagedGroup * const )
} /* namespace geosx */
