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

  SortedArrayView< localIndex const > const & interfaceNodesSet = m_interfaceNodesSet.toViewConst();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 > const p_n = nodeManager.getField< fields::Pressure_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_np1 >();

    arrayView1d< real32 > const ux_nm1 = nodeManager.getField< fields::Displacementx_nm1 >();
    arrayView1d< real32 > const uy_nm1 = nodeManager.getField< fields::Displacementy_nm1 >();
    arrayView1d< real32 > const uz_nm1 = nodeManager.getField< fields::Displacementz_nm1 >();
    arrayView1d< real32 > const ux_n   = nodeManager.getField< fields::Displacementx_n >();
    arrayView1d< real32 > const uy_n   = nodeManager.getField< fields::Displacementy_n >();
    arrayView1d< real32 > const uz_n   = nodeManager.getField< fields::Displacementz_n >();
    arrayView1d< real32 > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();

    elasSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    real64 rhof = 1020;  // hardcoded until github.com/GEOS-DEV/GEOS/pull/2548 is merged

    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = interfaceNodesSet[n];
      // TODO: obtain basis functions
      ux_np1[a] += 1e-3 * (p_n[a] / rhof);
      uy_np1[a] += 1e-3 * (p_n[a] / rhof);
      uz_np1[a] += 1e-3 * (p_n[a] / rhof);
      /*
         ux_np1[a] += Afsx * (p_n[a] / ρf);
         uy_np1[a] += Afsy * (p_n[a] / ρf);
         uz_np1[a] += Afsz * (p_n[a] / ρf);
       */
    } );

    elasSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    acousSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    real64 const dt2 = pow( dt, 2 );

    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = interfaceNodesSet[n];
      // TODO: obtain basis functions
      p_np1[a] += rhof * (
        1e-3 * (( ux_np1[a] - 2.0 * ux_n[a] + ux_nm1[a] )) +
        1e-3 * (( uy_np1[a] - 2.0 * uy_n[a] + uy_nm1[a] )) +
        1e-3 * (( uz_np1[a] - 2.0 * uz_n[a] + uz_nm1[a] ))
        ) / dt2; // debug
      /*
         p_np1[a] += (
           Asfx * (ρf * ( ux_np1[a] - 2.0 * ux_n[a] + ux_nm1[a] ) / dt2 ) +
           Asfy * (ρf * ( uy_np1[a] - 2.0 * uy_n[a] + uy_nm1[a] ) / dt2 ) +
           Asfz * (ρf * ( uz_np1[a] - 2.0 * uz_n[a] + uz_nm1[a] ) / dt2 )
         );
       */
    } );

    acousSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

  } );

  return dt;
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticElasticWaveEquationSEM, string const &, Group * const )

} /* namespace geos */
