/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
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


void BlueprintOutput::Execute(real64 const time_n,
                              real64 const dt,
                              integer const cycleNumber,
                              integer const eventCounter,
                              real64 const eventProgress,
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
