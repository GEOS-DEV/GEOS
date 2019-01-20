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


#include "DummySolver.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include <thread>
#include <chrono>

namespace geosx
{


using namespace dataRepository;


DummySolver::DummySolver( const std::string& name,
                                                  ManagedGroup * const parent ):
  SolverBase( name, parent )
{

  RegisterViewWrapper<real64>( dummyViewKeys.rand_scale.Key() )->
    setApplyDefaultValue(1e-9)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Scale for modifying requested dt");
}



DummySolver::~DummySolver()
{
  // TODO Auto-generated destructor stub
}


void DummySolver::InitializePreSubGroups( ManagedGroup * const problemManager )
{
  integer const rank = CommunicationTools::MPI_Rank(MPI_COMM_GEOSX );
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
  integer const rank = CommunicationTools::MPI_Rank(MPI_COMM_GEOSX );

  real64 const rand_scale = this->getReference<real64>(dummyViewKeys.rand_scale);
  real64 dt_request = std::rand() * rand_scale;

  std::cout << "time=" << time << ", solver=" << this->getName() << ", rank=" << rank << ", dt_r=" << dt_request << std::endl;

  return dt_request;
}


REGISTER_CATALOG_ENTRY( SolverBase, DummySolver, std::string const &, ManagedGroup * const )
} /* namespace ANST */
