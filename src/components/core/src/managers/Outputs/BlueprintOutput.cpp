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
 * BlueprintOutput.cpp
 *
 *  Created on: Jan 31, 2018
 *      Author: sherman
 */

#include "BlueprintOutput.hpp"
#include "DocumentationNode.hpp"
#include "mesh/MeshLevel.hpp"
#include "managers/DomainPartition.hpp"
#include "fileIO/blueprint/Blueprint.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

BlueprintOutput::BlueprintOutput( std::string const & name,
                        ManagedGroup * const parent ):
  OutputBase( name, parent)
{
}

BlueprintOutput::~BlueprintOutput()
{}

void BlueprintOutput::FillDocumentationNode()
{
  OutputBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("Blueprint");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Outputs Blueprint format files");
}

void BlueprintOutput::Execute(real64 const& time_n,
                              real64 const& dt,
                              const int cycleNumber,
                              ManagedGroup * domain)
{
  DomainPartition* domainPartition = ManagedGroup::group_cast<DomainPartition*>(domain);
  
  const MeshLevel * meshLevel = domainPartition->getMeshBody(0)->getMeshLevel(0);
  Blueprint bpWriter(*meshLevel->getNodeManager(),
                     *meshLevel->getElemManager(),
                     "bp_plot", MPI_COMM_WORLD);
  
  bpWriter.write(cycleNumber);
}


REGISTER_CATALOG_ENTRY( OutputBase, BlueprintOutput, std::string const &, ManagedGroup * const )
} /* namespace geosx */
