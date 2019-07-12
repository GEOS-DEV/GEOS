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
  SolverBase( name, parent ),
  m_randScale(0.0),
  m_randSeed(0)
{
  RegisterViewWrapper(viewKeyStruct::randScaleString, &m_randScale, false )->
    setApplyDefaultValue(1e-9)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Scale for modifying requested dt");

  RegisterViewWrapper(viewKeyStruct::randSeedString, &m_randSeed, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Scale for modifying requested dt");
}



DummySolver::~DummySolver()
{
  // TODO Auto-generated destructor stub
}


void DummySolver::InitializePreSubGroups( ManagedGroup * const problemManager )
{
  if (m_randSeed > 0)
  {
    integer const rank = CommunicationTools::MPI_Rank( MPI_COMM_GEOSX );
    std::srand((1 + rank) * m_randSeed);
  }
}


real64 DummySolver::SolverStep( real64 const& time_n,
                                        real64 const& dt,
                                        const int cycleNumber,
                                        DomainPartition * domain )
{
  return dt;
}


real64 DummySolver::GetTimestepRequest(real64 const time)
{
  real64 dt_request = 1.0 + std::rand() * m_randScale;
  return dt_request;
}


REGISTER_CATALOG_ENTRY( SolverBase, DummySolver, std::string const &, ManagedGroup * const )
} /* namespace ANST */
