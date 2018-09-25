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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


#include "DummySolver.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include <thread>
#include <chrono>

namespace geosx
{


using namespace dataRepository;


DummySolver::DummySolver( const std::string& name,
                                                  ManagedGroup * const parent ):
  SolverBase( name, parent )
{}



DummySolver::~DummySolver()
{
  // TODO Auto-generated destructor stub
}


void DummySolver::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  SolverBase::FillDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Dummy solver for testing time-stepping behavior");

  docNode->AllocateChildNode( dummyViewKeys.rand_scale.Key(),
                              dummyViewKeys.rand_scale.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Scale for modifying requested dt",
                              "Scale for modifying requested dt",
                              "1e-9",
                              "",
                              1,
                              1,
                              0 );

}


void DummySolver::Initialize( ManagedGroup * const problemManager )
{
  integer rank = 0;
  #ifdef GEOSX_USE_MPI
    MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
  #endif
  std::srand(rank * 12345);
}


real64 DummySolver::SolverStep( real64 const& time_n,
                                        real64 const& dt,
                                        const int cycleNumber,
                                        DomainPartition * domain )
{
  std::this_thread::sleep_for(std::chrono::seconds(1));
  return dt;
}


real64 DummySolver::GetTimestepRequest(real64 const time)
{
  integer rank = 0;
  #ifdef GEOSX_USE_MPI
    MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
  #endif

  real64 const rand_scale = this->getReference<real64>(dummyViewKeys.rand_scale);
  real64 dt_request = std::rand() * rand_scale;

  std::cout << "time=" << time << ", solver=" << this->getName() << ", rank=" << rank << ", dt_r=" << dt_request << std::endl;

  return dt_request;
}


REGISTER_CATALOG_ENTRY( SolverBase, DummySolver, std::string const &, ManagedGroup * const )
} /* namespace ANST */
