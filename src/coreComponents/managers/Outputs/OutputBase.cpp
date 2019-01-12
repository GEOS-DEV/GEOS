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
 * @file OutputBase.cpp
 */

#include "OutputBase.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

OutputBase::OutputBase( std::string const & name,
                        ManagedGroup * const parent ):
  ExecutableGroup( name, parent),
  m_slaveDirectory(),
  m_parallelThreads(1)
{
  RegisterViewWrapper(viewKeysStruct::slaveDirectoryString, &m_slaveDirectory, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("slave directory path");

  RegisterViewWrapper(viewKeysStruct::parallelThreadsString, &m_parallelThreads, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Number of plot files.");

}

OutputBase::~OutputBase()
{}

OutputBase::CatalogInterface::CatalogType& OutputBase::GetCatalog()
{
  static OutputBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


void OutputBase::InitializePreSubGroups( ManagedGroup * const group )
{
  // This command doesn't seem to work anymore
  // SetupDirectoryStructure();
}


void OutputBase::SetupDirectoryStructure()
{
  string slaveDirectory = m_slaveDirectory;

  int rank;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
  if (rank  == 0)
  {
    if (!slaveDirectory.empty())
    {
      string cmd = "mkdir -p " + slaveDirectory;
      int ret = std::system(cmd.c_str());
      GEOS_ERROR_IF(ret != 0, "Command '" << cmd << "' exited with code " << std::to_string(ret));
    }
  }
}


} /* namespace geosx */
