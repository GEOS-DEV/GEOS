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
 *
 */


#include "SinglePhasePoromechanicsSolver.hpp"

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "SinglePhasePoromechanicsKernel.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhasePoromechanicsSolver::SinglePhasePoromechanicsSolver( const string & name,
                                                                Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName()

{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver to use in the poromechanics solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fluid mechanics solver to use in the poromechanics solver" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void SinglePhasePoromechanicsSolver::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forMeshTargets( meshBodies, [&] ( string const &,
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

void SinglePhasePoromechanicsSolver::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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

void SinglePhasePoromechanicsSolver::setupDofs( DomainPartition const & domain,
                                                DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->setupDofs( domain, dofManager );
  m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          extrinsicMeshData::flow::pressure::key(),
                          DofManager::Connector::Elem );
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

void SinglePhasePoromechanicsSolver::implicitStepSetup( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_solidSolver->implicitStepSetup( time_n, dt, domain );
}

void SinglePhasePoromechanicsSolver::implicitStepComplete( real64 const & time_n,
                                                           real64 const & dt,
                                                           DomainPartition & domain )
{
  m_solidSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );
}

void SinglePhasePoromechanicsSolver::postProcessInput()
{
  SolverBase::postProcessInput();

  m_flowSolver = &this->getParent().getGroup< SinglePhaseBase >( m_flowSolverName );
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
}

void SinglePhasePoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{
  if( m_flowSolver->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics;
  }
}

SinglePhasePoromechanicsSolver::~SinglePhasePoromechanicsSolver()
{
  // TODO Auto-generated destructor stub
}

void SinglePhasePoromechanicsSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_solidSolver->resetStateToBeginningOfStep( domain );
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
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();

    string const dofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
    arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

    string const pDofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );

//  m_solidSolver->resetStressToBeginningOfStep( domain );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    PoromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                                  pDofKey,
                                                                  dofManager.rankOffset(),
                                                                  localMatrix,
                                                                  localRhs,
                                                                  gravityVectorData,
                                                                  FlowSolverBase::viewKeyStruct::fluidNamesString() );

    // Cell-based contributions
    m_solidSolver->getMaxForce() =
      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                      constitutive::PorousSolidBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              this->getDiscretizationName(),
                                                              viewKeyStruct::porousMaterialNamesString(),
                                                              kernelFactory );

  } );

  m_flowSolver->assemblePoroelasticFluxTerms( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              " " );
}

void SinglePhasePoromechanicsSolver::applyBoundaryConditions( real64 const time_n,
                                                              real64 const dt,
                                                              DomainPartition & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  m_solidSolver->applyBoundaryConditions( time_n, dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  m_flowSolver->applyBoundaryConditions( time_n, dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

real64 SinglePhasePoromechanicsSolver::calculateResidualNorm( DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm = m_solidSolver->calculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = m_flowSolver->calculateResidualNorm( domain, dofManager, localRhs );

  GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "    ( Rsolid, Rfluid ) = ( {:4.2e}, {:4.2e} )", momementumResidualNorm, massResidualNorm ) );

  return sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm );
}

void SinglePhasePoromechanicsSolver::createPreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( m_solidSolver->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { keys::TotalDisplacement, { 3, true } } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( m_flowSolver->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { extrinsicMeshData::flow::pressure::key(), { 1, true } } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void SinglePhasePoromechanicsSolver::solveSystem( DofManager const & dofManager,
                                                  ParallelMatrix & matrix,
                                                  ParallelVector & rhs,
                                                  ParallelVector & solution )
{
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void SinglePhasePoromechanicsSolver::applySystemSolution( DofManager const & dofManager,
                                                          arrayView1d< real64 const > const & localSolution,
                                                          real64 const scalingFactor,
                                                          DomainPartition & domain )
{
  // update displacement field
  m_solidSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

void SinglePhasePoromechanicsSolver::updateState( DomainPartition & domain )
{
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion & subRegion )
    {
      m_flowSolver->updateFluidState( subRegion );

    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
