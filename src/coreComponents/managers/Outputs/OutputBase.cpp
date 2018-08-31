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
 * @file OutputBase.cpp
 */

#include "OutputBase.hpp"
#include "DocumentationNode.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

OutputBase::OutputBase( std::string const & name,
                        ManagedGroup * const parent ):
  ExecutableGroup( name, parent)
{
}

OutputBase::~OutputBase()
{}

OutputBase::CatalogInterface::CatalogType& OutputBase::GetCatalog()
{
  static OutputBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void OutputBase::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("OutputBase");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Outputs Base Class");

  docNode->AllocateChildNode( outputBaseViewKeys.slaveDirectory.Key(),
                              outputBaseViewKeys.slaveDirectory.Key(),
                              -1,
                              "string",
                              "string",
                              "slave directory path",
                              "slave directory path",
                              "",
                              "",
                              0,
                              1,
                              0 );

    docNode->AllocateChildNode( outputBaseViewKeys.parallelThreads.Key(),
                                outputBaseViewKeys.parallelThreads.Key(),
                              -1,
                              "integer",
                              "integer",
                              "plot file threads",
                              "plot file threads",
                              "1",
                              "",
                              0,
                              1,
                              0 );
}

void OutputBase::Initialize( ManagedGroup * const group )
{
  // This command doesn't seem to work anymore
  // SetupDirectoryStructure();
}


void OutputBase::SetupDirectoryStructure()
{
  string slaveDirectory = this->getReference<string>(outputBaseViewKeys.slaveDirectory);

  int rank;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
  if (rank  == 0)
  {
    if (!slaveDirectory.empty())
    {
      string cmd = "mkdir -p " + slaveDirectory;
      std::system(cmd.c_str());
    }
  }
}


} /* namespace geosx */
