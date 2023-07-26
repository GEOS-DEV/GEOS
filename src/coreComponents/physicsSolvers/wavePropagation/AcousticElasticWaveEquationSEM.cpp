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
 * @file AcousticElasticWaveEquationSEM.cpp
 */

#include "AcousticElasticWaveEquationSEM.hpp"
#include "AcousticElasticWaveEquationSEMKernel.hpp"
#include "dataRepository/Group.hpp"

namespace geos
{
using namespace dataRepository;

void AcousticElasticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    nodeManager.registerField< fields::CouplingVector >( this->getName() );
  } );
}

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

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X32 = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    ArrayOfArraysView< localIndex const > const faceToNode  = faceManager.nodeList().toViewConst();
    arrayView2d< localIndex const > const & faceToRegion    = faceManager.elementRegionList();
    arrayView2d< localIndex const > const & faceToSubRegion = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const & faceToElement   = faceManager.elementList();

    arrayView1d< real32 > const coupling = nodeManager.getField< fields::CouplingVector >();
    coupling.zero();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticElasticWaveEquationSEMKernels::CouplingKernel< FE_TYPE > kernelC( finiteElement );

        kernelC.template launch< EXEC_POLICY, ATOMIC_POLICY >( faceManager.size(),
                                                               X32,
                                                               faceToRegion,
                                                               faceToSubRegion,
                                                               faceToElement,
                                                               faceToNode,
                                                               coupling );
      } );
    } );
  } );
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
  // auto acous2elasCoupling = acousticElasticWaveEquationSEMKernels_old::AcousticElasticSEMFactory( interfaceNodesSet, dt );
  // auto elas2acousCoupling = acousticElasticWaveEquationSEMKernels_old::ElasticAcousticSEMFactory( interfaceNodesSet, dt );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const p_n      = nodeManager.getField< fields::Pressure_n >();
    arrayView1d< real32 const > const ux_nm1   = nodeManager.getField< fields::Displacementx_nm1 >();
    arrayView1d< real32 const > const uy_nm1   = nodeManager.getField< fields::Displacementy_nm1 >();
    arrayView1d< real32 const > const uz_nm1   = nodeManager.getField< fields::Displacementz_nm1 >();
    arrayView1d< real32 const > const ux_n     = nodeManager.getField< fields::Displacementx_n >();
    arrayView1d< real32 const > const uy_n     = nodeManager.getField< fields::Displacementy_n >();
    arrayView1d< real32 const > const uz_n     = nodeManager.getField< fields::Displacementz_n >();
    arrayView1d< real32 const > const coupling = nodeManager.getField< fields::CouplingVector >();

    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_np1 >();
    arrayView1d< real32 > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();

    elasSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    real64 rhof = 1020;  // hardcoded until github.com/GEOS-DEV/GEOS/pull/2548 is merged

    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = interfaceNodesSet[n];
      ux_np1[a] += coupling[a] * (p_n[a] / rhof);
      uy_np1[a] += coupling[a] * (p_n[a] / rhof);
      uz_np1[a] += coupling[a] * (p_n[a] / rhof);
    } );

    elasSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    acousSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    real64 const dt2 = pow( dt, 2 );

    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = interfaceNodesSet[n];
      p_np1[a] += coupling[a] * rhof * (
        ( ux_np1[a] - 2.0 * ux_n[a] + ux_nm1[a] ) +
        ( uy_np1[a] - 2.0 * uy_n[a] + uy_nm1[a] ) +
        ( uz_np1[a] - 2.0 * uz_n[a] + uz_nm1[a] )
        ) / dt2;
    } );

    acousSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    /*
    elasSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    finiteElement::
      regionBasedKernelApplication< EXEC_POLICY,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            getDiscretizationName(),
                                                            "",
                                                            acous2elasCoupling );

    acousSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    finiteElement::
      regionBasedKernelApplication< EXEC_POLICY,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            getDiscretizationName(),
                                                            "",
                                                            elas2acousCoupling );

    elasSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );
    acousSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );
    */
  } );

  return dt;
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticElasticWaveEquationSEM, string const &, Group * const )

} /* namespace geos */
