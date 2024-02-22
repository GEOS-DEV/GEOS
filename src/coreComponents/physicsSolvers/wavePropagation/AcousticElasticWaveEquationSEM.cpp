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
#include <typeinfo>
#include <limits>

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
    nodeManager.registerField< fields::CouplingVectorx >( getName() );
    nodeManager.registerField< fields::CouplingVectory >( getName() );
    nodeManager.registerField< fields::CouplingVectorz >( getName() );
  } );
}

void AcousticElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  auto acousSolver = acousticSolver();
  auto elasSolver = elasticSolver();

  auto acousNodesSet = acousSolver->getSolverNodesSet();
  auto elasNodesSet = elasSolver->getSolverNodesSet();

  for( auto val : acousNodesSet )
  {
    if( elasNodesSet.contains( val ) )
      m_interfaceNodesSet.insert( val );
  }
  localIndex const numInterfaceNodes = MpiWrapper::sum( m_interfaceNodesSet.size() );
  GEOS_THROW_IF( numInterfaceNodes == 0, "Failed to compute interface: check xml input (solver order)", std::runtime_error );

  m_acousRegions = acousSolver->getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );
  m_elasRegions = elasSolver->getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & GEOS_UNUSED_PARAM( regionNames ) )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    arrayView2d< real64 const > const faceNormals          = faceManager.faceNormal().toViewConst();
    ArrayOfArraysView< localIndex const > const faceToNode = faceManager.nodeList().toViewConst();
    arrayView2d< localIndex const > const faceToSubRegion  = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const faceToRegion     = faceManager.elementRegionList();
    arrayView2d< localIndex const > const faceToElement    = faceManager.elementList();

    arrayView1d< real32 > const couplingVectorx = nodeManager.getField< fields::CouplingVectorx >();
    couplingVectorx.zero();

    arrayView1d< real32 > const couplingVectory = nodeManager.getField< fields::CouplingVectory >();
    couplingVectory.zero();

    arrayView1d< real32 > const couplingVectorz = nodeManager.getField< fields::CouplingVectorz >();
    couplingVectorz.zero();

    elemManager.forElementRegions( m_acousRegions, [&] ( localIndex const regionIndex, ElementRegionBase const & elemRegion )
    {
      elemRegion.forElementSubRegionsIndex( [&]( localIndex const subRegionIndex, ElementSubRegionBase const & elementSubRegion )
      {
        finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

        finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
        {
          using FE_TYPE = TYPEOFREF( finiteElement );

          acousticElasticWaveEquationSEMKernels::CouplingKernel< FE_TYPE > kernelC;
          kernelC.template launch< EXEC_POLICY, ATOMIC_POLICY >( faceManager.size(),
                                                                 nodeCoords,
                                                                 regionIndex,
                                                                 subRegionIndex,
                                                                 faceToSubRegion,
                                                                 faceToRegion,
                                                                 faceToElement,
                                                                 faceToNode,
                                                                 faceNormals,
                                                                 couplingVectorx,
                                                                 couplingVectory,
                                                                 couplingVectorz );
        } );
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

  auto acousSolver = acousticSolver();
  auto elasSolver = elasticSolver();

  SortedArrayView< localIndex const > const interfaceNodesSet = m_interfaceNodesSet.toViewConst();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & GEOS_UNUSED_PARAM( regionNames ) )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const acousticMass = nodeManager.getField< fields::AcousticMassVector >();
    arrayView1d< real32 const > const elasticMass = nodeManager.getField< fields::ElasticMassVector >();
    arrayView1d< localIndex > const acousticFSNodeIndicator = nodeManager.getField< fields::AcousticFreeSurfaceNodeIndicator >();
    arrayView1d< localIndex > const elasticFSNodeIndicator = nodeManager.getField< fields::ElasticFreeSurfaceNodeIndicator >();

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

    real32 const dt2 = pow( dt, 2 );

    elasSolver->computeUnknowns( time_n, dt, cycleNumber, domain, mesh, m_elasRegions );

    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = interfaceNodesSet[n];
      if( elasticFSNodeIndicator[a] == 1 )
        return;

      real32 const aux = -p_n[a] / elasticMass[a];
      real32 const localIncrementx = dt2 * atoex[a] * aux;
      real32 const localIncrementy = dt2 * atoey[a] * aux;
      real32 const localIncrementz = dt2 * atoez[a] * aux;

      RAJA::atomicAdd< ATOMIC_POLICY >( &ux_np1[a], localIncrementx );
      RAJA::atomicAdd< ATOMIC_POLICY >( &uy_np1[a], localIncrementy );
      RAJA::atomicAdd< ATOMIC_POLICY >( &uz_np1[a], localIncrementz );
    } );

    elasSolver->synchronizeUnknowns( time_n, dt, cycleNumber, domain, mesh, m_elasRegions );

    acousSolver->computeUnknowns( time_n, dt, cycleNumber, domain, mesh, m_acousRegions );

    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = interfaceNodesSet[n];
      if( acousticFSNodeIndicator[a] == 1 )
        return;

      real32 const localIncrement = (
        atoex[a] * ( ux_np1[a] - 2.0 * ux_n[a] + ux_nm1[a] ) +
        atoey[a] * ( uy_np1[a] - 2.0 * uy_n[a] + uy_nm1[a] ) +
        atoez[a] * ( uz_np1[a] - 2.0 * uz_n[a] + uz_nm1[a] )
        ) / acousticMass[a];

      RAJA::atomicAdd< ATOMIC_POLICY >( &p_np1[a], localIncrement );
    } );

    acousSolver->synchronizeUnknowns( time_n, dt, cycleNumber, domain, mesh, m_acousRegions );

    acousSolver->prepareNextTimestep( mesh );
    elasSolver->prepareNextTimestep( mesh );
  } );

  return dt;
}

void AcousticElasticWaveEquationSEM::cleanup( real64 const time_n,
                                              integer const cycleNumber,
                                              integer const eventCounter,
                                              real64 const eventProgress,
                                              DomainPartition & domain )
{
  elasticSolver()->cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );
  acousticSolver()->cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticElasticWaveEquationSEM, string const &, Group * const )

} /* namespace geos */
