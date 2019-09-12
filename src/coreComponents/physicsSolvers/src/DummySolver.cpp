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


#include "DummySolver.hpp"
#include <thread>
#include <chrono>

#include "mpiCommunications/CommunicationTools.hpp"
#include "dataRepository/Group.hpp"

namespace geosx
{


using namespace dataRepository;


DummySolver::DummySolver( const std::string& name,
                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_randScale(0.0),
  m_randSeed(0)
{
  registerWrapper(viewKeyStruct::randScaleString, &m_randScale, false )->
    setApplyDefaultValue(1e-9)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Scale for modifying requested dt");

  registerWrapper(viewKeyStruct::randSeedString, &m_randSeed, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Scale for modifying requested dt");
}



DummySolver::~DummySolver()
{
  // TODO Auto-generated destructor stub
}


void DummySolver::InitializePreSubGroups( Group * const GEOSX_UNUSED_ARG( problemManager ) )
{
  if (m_randSeed > 0)
  {
    integer const rank = CommunicationTools::MPI_Rank( MPI_COMM_GEOSX );
    std::srand((1 + rank) * m_randSeed);
  }
}


real64 DummySolver::SolverStep( real64 const& GEOSX_UNUSED_ARG( time_n ),
                                real64 const& dt,
                                const int GEOSX_UNUSED_ARG( cycleNumber ),
                                DomainPartition * GEOSX_UNUSED_ARG( domain ) )
{
  return dt;
}


real64 DummySolver::GetTimestepRequest( real64 const GEOSX_UNUSED_ARG( time ) )
{
  real64 dt_request = 1.0 + std::rand() * m_randScale;
  return dt_request;
}


REGISTER_CATALOG_ENTRY( SolverBase, DummySolver, std::string const &, Group * const )
} /* namespace ANST */
