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
 * @file SinglePhasePoromechanicsSolver.cpp
 */

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "SinglePhasePoromechanicsSolver.hpp"

#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"

namespace geosx
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SinglePhasePoromechanicsSolver::SinglePhasePoromechanicsSolver( const string & name,
                                                                Group * const parent )
  : Base( name, parent )
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void SinglePhasePoromechanicsSolver::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< string >( viewKeyStruct::porousMaterialNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );
    } );
  } );
}

void SinglePhasePoromechanicsSolver::setupCoupling( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                    DofManager & dofManager ) const
{
  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void SinglePhasePoromechanicsSolver::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                       [&]( localIndex const,
                                                                            ElementSubRegionBase & subRegion )
    {
      string & porousName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      porousName = getConstitutiveName< CoupledSolidBase >( subRegion );
      GEOSX_ERROR_IF( porousName.empty(), GEOSX_FMT( "Solid model not found on subregion {}", subRegion.getName() ) );
    } );
  } );
}

void SinglePhasePoromechanicsSolver::setupSystem( DomainPartition & domain,
                                                  DofManager & dofManager,
                                                  CRSMatrix< real64, globalIndex > & localMatrix,
                                                  ParallelVector & rhs,
                                                  ParallelVector & solution,
                                                  bool const setSparsity )
{
  if( m_precond )
  {
    m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner();
  }
}

void SinglePhasePoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{
  if( flowSolver()->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics;
  }
}

real64 SinglePhasePoromechanicsSolver::solverStep( real64 const & time_n,
                                                   real64 const & dt,
                                                   int const cycleNumber,
                                                   DomainPartition & domain )
{
  real64 dt_return = dt;

  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_rhs,
               m_solution );

  implicitStepSetup( time_n, dt, domain );

  dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void SinglePhasePoromechanicsSolver::assembleSystem( real64 const time_n,
                                                     real64 const dt,
                                                     DomainPartition & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{

  GEOSX_MARK_FUNCTION;

  real64 poromechanicsMaxForce = 0.0;
  real64 mechanicsMaxForce = 0.0;


  // step 1: apply the full poromechanics coupling on the target regions on the poromechanics solver

  set< string > poromechanicsRegionNames;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    poromechanicsRegionNames.insert( regionNames.begin(), regionNames.end() );

    string const flowDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    poromechanicsMaxForce =
      assemblyLaunch< constitutive::PorousSolidBase,
                      poromechanicsKernels::SinglePhaseKernelFactory >( mesh,
                                                                        dofManager,
                                                                        regionNames,
                                                                        viewKeyStruct::porousMaterialNamesString(),
                                                                        localMatrix,
                                                                        localRhs,
                                                                        flowDofKey,
                                                                        FlowSolverBase::viewKeyStruct::fluidNamesString() );



  } );


  // step 2: apply mechanics solver on its target regions not included in the poromechanics solver target regions

  solidMechanicsSolver()->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                        MeshLevel & mesh,
                                                                                        arrayView1d< string const > const & regionNames )
  {

    // collect the target region of the mechanics solver not included in the poromechanics target regions
    array1d< string > filteredRegionNames;
    filteredRegionNames.reserve( regionNames.size() );
    for( string const & regionName : regionNames )
    {
      // if the mechanics target region is not included in the poromechanics target region, save the string
      if( poromechanicsRegionNames.count( regionName ) == 0 )
      {
        filteredRegionNames.emplace_back( regionName );
      }
    }

    // if the array is empty, the mechanics and poromechanics solver target regions overlap perfectly, there is nothing to do
    if( filteredRegionNames.empty() )
    {
      return;
    }

    mechanicsMaxForce =
      assemblyLaunch< constitutive::SolidBase,
                      solidMechanicsLagrangianFEMKernels::QuasiStaticFactory >( mesh,
                                                                                dofManager,
                                                                                filteredRegionNames.toViewConst(),
                                                                                SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                                                localMatrix,
                                                                                localRhs );
  } );

  solidMechanicsSolver()->getMaxForce() = LvArray::math::max( mechanicsMaxForce, poromechanicsMaxForce );


  // step 3: compute the fluxes (face-based contributions)

  flowSolver()->assemblePoroelasticFluxTerms( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              " " );
}

void SinglePhasePoromechanicsSolver::createPreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( solidMechanicsSolver()->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { solidMechanics::totalDisplacement::key(), { 3, true } } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( flowSolver()->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { flow::pressure::key(), { 1, true } } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void SinglePhasePoromechanicsSolver::updateState( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion & subRegion )
    {
      flowSolver()->updateFluidState( subRegion );

    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
