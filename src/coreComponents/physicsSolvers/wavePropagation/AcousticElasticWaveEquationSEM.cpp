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
    nodeManager.registerField< fields::CouplingVectorx >( this->getName() );
    nodeManager.registerField< fields::CouplingVectory >( this->getName() );
    nodeManager.registerField< fields::CouplingVectorz >( this->getName() );
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

    ArrayOfArraysView< localIndex const > const faceToNode = faceManager.nodeList().toViewConst();
    arrayView2d< localIndex const > const faceToRegion     = faceManager.elementRegionList();
    arrayView2d< localIndex const > const faceToSubRegion  = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const faceToElement    = faceManager.elementList();
    arrayView2d< real64 const > const faceNormals          = faceManager.faceNormal().toViewConst();

    arrayView1d< real32 > const atoex = nodeManager.getField< fields::CouplingVectorx >();
    atoex.zero();

    arrayView1d< real32 > const atoey = nodeManager.getField< fields::CouplingVectory >();
    atoey.zero();

    arrayView1d< real32 > const atoez = nodeManager.getField< fields::CouplingVectorz >();
    atoez.zero();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const targetIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticElasticWaveEquationSEMKernels::CouplingKernel< FE_TYPE > kernelC;

        // arrayView1d< real32 const > const density = elementSubRegion.getField< fields::MediumDensity >();

        kernelC.template launch< EXEC_POLICY, ATOMIC_POLICY >( faceManager.size(),
                                                               targetIndex,
                                                               X32,
                                                               // density,
                                                               faceToRegion,
                                                               faceToSubRegion,
                                                               faceToElement,
                                                               faceToNode,
                                                               faceNormals,
                                                               atoex,
                                                               atoey,
                                                               atoez );
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

  auto acousSolver = acousticSolver();
  auto elasSolver = elasticSolver();

  SortedArrayView< localIndex const > const & interfaceNodesSet = m_interfaceNodesSet.toViewConst();

#if 1
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const p_n    = nodeManager.getField< fields::Pressure_n >();
    arrayView1d< real32 const > const ux_nm1 = nodeManager.getField< fields::Displacementx_nm1 >();
    arrayView1d< real32 const > const uy_nm1 = nodeManager.getField< fields::Displacementy_nm1 >();
    arrayView1d< real32 const > const uz_nm1 = nodeManager.getField< fields::Displacementz_nm1 >();
    arrayView1d< real32 const > const ux_n   = nodeManager.getField< fields::Displacementx_n >();
    arrayView1d< real32 const > const uy_n   = nodeManager.getField< fields::Displacementy_n >();
    arrayView1d< real32 const > const uz_n   = nodeManager.getField< fields::Displacementz_n >();
    arrayView1d< real32 const > const atoex  = nodeManager.getField< fields::CouplingVectorx >();
    arrayView1d< real32 const > const atoey  = nodeManager.getField< fields::CouplingVectory >();
    arrayView1d< real32 const > const atoez  = nodeManager.getField< fields::CouplingVectorz >();

    arrayView1d< real32 > const p_np1  = nodeManager.getField< fields::Pressure_np1 >();
    arrayView1d< real32 > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();

    elasSolver->computeUnknowns( time_n, dt, cycleNumber, domain, mesh, regionNames );

    real64 const rho0 = 1020.0;  // hardcoded until github.com/GEOS-DEV/GEOS/pull/2548 is merged
    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = interfaceNodesSet[n];
      
      real32 const localIncrementx = atoex[a] * (p_n[a] / rho0);
      real32 const localIncrementy = atoey[a] * (p_n[a] / rho0);
      real32 const localIncrementz = atoez[a] * (p_n[a] / rho0);

      // printf(
      //   "\t[AcousticElasticWaveEquationSEM::solverStep] atoex=%g atoey=%g atoez=%g incx=%g incy=%g incz=%g\n",
      //   atoex[a], atoey[a], atoez[a],
      //   localIncrementx, localIncrementy, localIncrementz
      // );

      RAJA::atomicAdd< ATOMIC_POLICY >( &ux_np1[a], localIncrementx );
      RAJA::atomicAdd< ATOMIC_POLICY >( &uy_np1[a], localIncrementy );
      RAJA::atomicAdd< ATOMIC_POLICY >( &uz_np1[a], localIncrementz );
    } );

    acousSolver->computeUnknowns( time_n, dt, cycleNumber, domain, mesh, regionNames );

    real64 const dt2 = pow( dt, 2 );
    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = interfaceNodesSet[n];

      real32 const localIncrement = rho0 * (
        atoex[a] * ( ux_np1[a] - 2.0 * ux_n[a] + ux_nm1[a] ) +
        atoey[a] * ( uy_np1[a] - 2.0 * uy_n[a] + uy_nm1[a] ) +
        atoez[a] * ( uz_np1[a] - 2.0 * uz_n[a] + uz_nm1[a] )
      ) / dt2;

      // printf(
      //   "\t[AcousticElasticWaveEquationSEM::solverStep] atoex=%g atoey=%g atoez=%g inc=%g\n",
      //   atoex[a], atoey[a], atoez[a],
      //   localIncrement
      // );

      RAJA::atomicAdd< ATOMIC_POLICY >( &p_np1[a], localIncrement );
    } );

    acousSolver->synchronizeUnknowns( time_n, dt, cycleNumber, domain, mesh, regionNames );
    acousSolver->prepareNextTimestep( mesh );
    elasSolver->synchronizeUnknowns( time_n, dt, cycleNumber, domain, mesh, regionNames );
    elasSolver->prepareNextTimestep( mesh );
  } );

#else

  auto acous2elasCoupling2 = acousticElasticWaveEquationSEMKernels2::AcousticElasticSEMFactory( interfaceNodesSet, dt );
  auto elas2acousCoupling2 = acousticElasticWaveEquationSEMKernels2::ElasticAcousticSEMFactory( interfaceNodesSet, dt );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elementRegionManager = mesh.getElemManager();

    elasSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                         CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        auto kernel = acousticElasticWaveEquationSEMKernels3::AcousticToElasticSEM<
          TYPEOFREF( elementSubRegion ),
          TYPEOFREF( finiteElement )
        >( nodeManager, elementSubRegion, finiteElement, interfaceNodesSet, dt );
        kernel.template kernelLaunch< EXEC_POLICY >( elementSubRegion.size() );
      } );
    } );
    // finiteElement::
    //   regionBasedKernelApplication< EXEC_POLICY,
    //                                 constitutive::NullModel,
    //                                 CellElementSubRegion >( mesh,
    //                                                         regionNames,
    //                                                         getDiscretizationName(),
    //                                                         "",
    //                                                         acous2elasCoupling2 );

    acousSolver->unknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );

    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                         CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        auto kernel = acousticElasticWaveEquationSEMKernels3::ElasticToAcousticSEM<
          TYPEOFREF( elementSubRegion ),
          TYPEOFREF( finiteElement )
        >( nodeManager, elementSubRegion, finiteElement, interfaceNodesSet, dt );
        kernel.template kernelLaunch< EXEC_POLICY >( elementSubRegion.size() );
      } );
    } );

    // finiteElement::
    //   regionBasedKernelApplication< EXEC_POLICY,
    //                                 constitutive::NullModel,
    //                                 CellElementSubRegion >( mesh,
    //                                                         regionNames,
    //                                                         getDiscretizationName(),
    //                                                         "",
    //                                                         elas2acousCoupling2 );

    elasSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );
    acousSolver->postUnknownsUpdate( time_n, dt, cycleNumber, domain, mesh, regionNames );
  } );

#endif

  return dt;
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticElasticWaveEquationSEM, string const &, Group * const )

} /* namespace geos */
