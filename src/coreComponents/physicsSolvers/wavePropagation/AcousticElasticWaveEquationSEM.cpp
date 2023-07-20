/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file AcousticElasticWaveEquation.cpp
 */

#include "AcousticElasticWaveEquationSEM.hpp"
#include "dataRepository/Group.hpp"

namespace geos
{
using namespace dataRepository;

void AcousticElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  auto acousNodesSet = acousticSolver()->getSolverNodesSet();
  auto elasNodesSet = elasticSolver()->getSolverNodesSet();

  for( auto val : acousNodesSet )
  {
    if( elasNodesSet.contains( val ))
      m_interfaceNodesSet.insert( val );
  }

  std::cout << "\t[AcousticElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups] "
            << "m_interfaceNodesSet.size()=" << m_interfaceNodesSet.size() << std::endl;
  GEOS_THROW_IF( m_interfaceNodesSet.size() == 0, "Failed to compute interface: check xml input (solver order)", std::runtime_error );
}

real64 AcousticElasticWaveEquationSEM::solverStep( real64 const & time_n,
                                                   real64 const & dt,
                                                   int const cycleNumber,
                                                   DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  std::cout << "\t[AcousticElasticWaveEquationSEM::solverStep]" << std::endl;

  auto elasSolver = elasticSolver();
  auto acousSolver = acousticSolver();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    elasSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    // TODO: coupling - get and write data into elastic fields

    elasSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    acousSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    // TODO: coupling - get and write data into acoustic fields

    acousSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );
  } );

  return dt;
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticElasticWaveEquationSEM, string const &, Group * const )

} /* namespace geos */
