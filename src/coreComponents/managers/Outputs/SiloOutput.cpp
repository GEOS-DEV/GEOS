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
 * @file SiloOutput.cpp
 */

#include "SiloOutput.hpp"
#include "DocumentationNode.hpp"
#include "fileIO/silo/SiloFile.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

SiloOutput::SiloOutput( std::string const & name,
                        ManagedGroup * const parent ):
  OutputBase( name, parent)
{
}

SiloOutput::~SiloOutput()
{}

void SiloOutput::FillDocumentationNode()
{
  OutputBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("Silo");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Outputs SILO format files");



  docNode->AllocateChildNode( siloOutputViewKeys.plotLevel.Key(),
                              siloOutputViewKeys.plotLevel.Key(),
                              -1,
                              "integer",
                              "integer",
                              "output level",
                              "output level",
                              "1",
                              "",
                              0,
                              1,
                              0 );

//  docNode->AllocateChildNode( viewKeys.writeFEMFaces.Key(),
//                              viewKeys.writeFEMFaces.Key(),
//                              -1,
//                              "integer",
//                              "integer",
//                              "flag to write FEM faces",
//                              "flag to write FEM faces",
//                              "0",
//                              "",
//                              0,
//                              1,
//                              0 );
}

void SiloOutput::Execute(real64 const& time_n,
                         real64 const& dt,
                         integer const cycleNumber,
                         integer const eventCounter,
                         real64 const & eventProgress,
                         ManagedGroup * domain)
{
  DomainPartition* domainPartition = ManagedGroup::group_cast<DomainPartition*>(domain);
  SiloFile silo;

  integer rank;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
  MPI_Barrier( MPI_COMM_GEOSX );

  integer numFiles = this->getReference<integer>( siloOutputViewKeys.parallelThreads);

  silo.setPlotLevel( getReference<integer>( siloOutputViewKeys.plotLevel ) );

  silo.Initialize(PMPIO_WRITE , numFiles );
  silo.WaitForBatonWrite(rank, cycleNumber, eventCounter, false );
  silo.WriteDomainPartition( *domainPartition, cycleNumber,  time_n + dt * eventProgress, 0);
  silo.HandOffBaton();
  silo.ClearEmptiesFromMultiObjects(cycleNumber);
  silo.Finish();

}


REGISTER_CATALOG_ENTRY( OutputBase, SiloOutput, std::string const &, ManagedGroup * const )
} /* namespace geosx */
