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
    nodeManager.registerField< fields::CouplingDensity>( this->getName() );
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
    arrayView2d< localIndex const > const faceToElement    = faceManager.elementList();
    arrayView2d< real64 const > const faceNormals          = faceManager.faceNormal().toViewConst();

    arrayView1d< real32 > const couplingVectorx = nodeManager.getField< fields::CouplingVectorx >();
    couplingVectorx.zero();

    arrayView1d< real32 > const couplingVectory = nodeManager.getField< fields::CouplingVectory >();
    couplingVectory.zero();

    arrayView1d< real32 > const couplingVectorz = nodeManager.getField< fields::CouplingVectorz >();
    couplingVectorz.zero();

    arrayView1d< real32 > const couplingDensity = nodeManager.getField< fields::CouplingDensity >();
    couplingDensity.zero();

    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    auto regions = elementRegionManager.getRegions();

#if 0
    auto fluid = regions["Fluid"];
    std::cout << "region=" << fluid << " type=" << typeid(fluid).name() << std::endl;  // geos::dataRepository::Group*
#endif

    // TODO: make this generic by looping through `getGroupByPath`

    localIndex const fluid_index = regions.getIndex( "Fluid" );
    localIndex const solid_index = regions.getIndex( "Solid" );
    printf("\t[AcousticElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups] fluid_index=%i solid_index=%i\n", fluid_index, solid_index);

    m_acousRegions.resize(1);
    m_acousRegions[0] = "Fluid";
    m_elasRegions.resize(1);
    m_elasRegions[0] = "Solid";

    elementRegionManager.forElementSubRegions< CellElementSubRegion >( m_acousRegions, [&]( localIndex const targetIndex,  // GEOS_UNUSED_PARAM(targetIndex)
                                                                                          CellElementSubRegion & elementSubRegion )
    {
#if 0
      auto & parent_name = elementSubRegion.getParent().getParent().getName();
      printf(
        "\t[AcousticElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups] targetIndex=%i elementSubRegion.getName()=%s parent=%s\n",
        targetIndex, elementSubRegion.getName().c_str(), parent_name.c_str()
      );

      if (parent_name != "Fluid") return;  // only loop over the fluid region
#endif
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      arrayView1d< real32 const > const fluid_density = elementSubRegion.getField< fields::MediumDensityA >();

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticElasticWaveEquationSEMKernels::CouplingKernel< FE_TYPE > kernelC;

        kernelC.template launch< EXEC_POLICY, ATOMIC_POLICY >( faceManager.size(),
                                                               X32,
                                                               fluid_index,
                                                               fluid_density,
                                                               faceToRegion,
                                                               faceToElement,
                                                               faceToNode,
                                                               faceNormals,
                                                               couplingVectorx,
                                                               couplingVectory,
                                                               couplingVectorz,
                                                               couplingDensity );
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
  // arrayView1d< string const > const & acousRegions = m_acousRegions.toViewConst();
  // arrayView1d< string const > const & elasRegions = m_elasRegions.toViewConst();

#if 1
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const massA = nodeManager.getField< fields::MassVectorA >();
    arrayView1d< real32 const > const massE = nodeManager.getField< fields::MassVectorE >();
    arrayView1d< localIndex > const freeSurfaceNodeIndicatorA = nodeManager.getField< fields::FreeSurfaceNodeIndicatorA >();
    arrayView1d< localIndex > const freeSurfaceNodeIndicatorE = nodeManager.getField< fields::FreeSurfaceNodeIndicatorE >();

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
    arrayView1d< real32 const > const rho0   = nodeManager.getField< fields::CouplingDensity >();

    arrayView1d< real32 > const p_np1  = nodeManager.getField< fields::Pressure_np1 >();
    arrayView1d< real32 > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();

    // source ok when using (elas=false, acous=true, coupling=false) or (elas=true, acous=true, coupling=false)
    bool const elas = !helpers::benv("NO_ELAS");
    bool const acous = !helpers::benv("NO_ACOUS");
    bool const coupling = !helpers::benv("NO_COUPLING");
    bool const dump = helpers::ienv("DUMP") > 0;
    // std::cout << "elas=" << (elas ? 'T' : 'F') << " acous=" << (acous ? 'T' : 'F') << " coupling=" << (coupling ? 'T' : 'F')  << " dump=" << (dump ? 'T' : 'F') << std::endl;

    for (auto nm : regionNames)
      std::cout << "\t[AcousticElasticWaveEquationSEM::solverStep] regionName=" << nm << std::endl;

    for (auto nm : m_acousRegions)
      std::cout << "\t[AcousticElasticWaveEquationSEM::solverStep] acousRegion=" << nm << std::endl;

    for (auto nm : m_elasRegions)
      std::cout << "\t[AcousticElasticWaveEquationSEM::solverStep] elasRegion=" << nm << std::endl;

    if (elas) elasSolver->computeUnknowns( time_n, dt, cycleNumber, domain, mesh, m_elasRegions );

    if (coupling)
    {
      forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
      {
        localIndex const a = interfaceNodesSet[n];
        if (freeSurfaceNodeIndicatorE[a] == 1) return;

        real32 const aux = (p_n[a] / rho0[a]) / massE[a];
        real32 const localIncrementx = atoex[a] * aux;
        real32 const localIncrementy = atoey[a] * aux;
        real32 const localIncrementz = atoez[a] * aux;

        if (dump) printf(
          "\t[AcousticElasticWaveEquationSEM::solverStep] atoex=%g atoey=%g atoez=%g p_n=%g rho0=%g massE=%g aux=%g incx=%g incy=%g incz=%g\n",
          atoex[a], atoey[a], atoez[a], p_n[a], rho0[a], massE[a], aux,
          localIncrementx, localIncrementy, localIncrementz
        );

        RAJA::atomicAdd< ATOMIC_POLICY >( &ux_np1[a], localIncrementx );
        RAJA::atomicAdd< ATOMIC_POLICY >( &uy_np1[a], localIncrementy );
        RAJA::atomicAdd< ATOMIC_POLICY >( &uz_np1[a], localIncrementz );
      } );
    }

    if (elas) elasSolver->synchronizeUnknowns( time_n, dt, cycleNumber, domain, mesh, m_elasRegions );

    if (acous) acousSolver->computeUnknowns( time_n, dt, cycleNumber, domain, mesh, m_acousRegions );

    if (coupling)
    {
      forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
      {
        localIndex const a = interfaceNodesSet[n];
        if (freeSurfaceNodeIndicatorA[a] == 1) return;

        real32 const localIncrement = rho0[a] * (
          -atoex[a] * ( ux_np1[a] - 2.0 * ux_n[a] + ux_nm1[a] ) +
          -atoey[a] * ( uy_np1[a] - 2.0 * uy_n[a] + uy_nm1[a] ) +
          -atoez[a] * ( uz_np1[a] - 2.0 * uz_n[a] + uz_nm1[a] )
        ) / massA[a];

        if (dump) printf(
          "\t[AcousticElasticWaveEquationSEM::solverStep] atoex=%g atoey=%g atoez=%g rho0=%g massA=%g inc=%g\n",
          atoex[a], atoey[a], atoez[a], rho0[a], massA[a], localIncrement
        );

        RAJA::atomicAdd< ATOMIC_POLICY >( &p_np1[a], localIncrement );
      } );
    }

    if (acous) acousSolver->synchronizeUnknowns( time_n, dt, cycleNumber, domain, mesh, m_acousRegions );

    if (acous) acousSolver->prepareNextTimestep( mesh );
    if (elas) elasSolver->prepareNextTimestep( mesh );
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
