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
 * @file BlueprintOutput.cpp
 */
#include "BlueprintOutput.hpp"
#include "mesh/MeshLevel.hpp"
#include "managers/DomainPartition.hpp"
#include "fileIO/blueprint/Blueprint.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

BlueprintOutput::BlueprintOutput( std::string const & name,
                        Group * const parent ):
  OutputBase( name, parent)
{
}

BlueprintOutput::~BlueprintOutput()
{}


void BlueprintOutput::Execute(real64 const GEOSX_UNUSED_ARG( time_n ),
                              real64 const GEOSX_UNUSED_ARG( dt ),
                              integer const cycleNumber,
                              integer const eventCounter,
                              real64 const GEOSX_UNUSED_ARG( eventProgress ),
                              Group * domain)
{
  DomainPartition* domainPartition = Group::group_cast<DomainPartition*>(domain);
  
  const MeshLevel * meshLevel = domainPartition->getMeshBody(0)->getMeshLevel(0);
  Blueprint bpWriter(*meshLevel->getNodeManager(),
                     *meshLevel->getElemManager(),
                     "bp_plot", MPI_COMM_GEOSX);
  
  bpWriter.write(cycleNumber, eventCounter );
}


REGISTER_CATALOG_ENTRY( OutputBase, BlueprintOutput, std::string const &, Group * const )
} /* namespace geosx */
