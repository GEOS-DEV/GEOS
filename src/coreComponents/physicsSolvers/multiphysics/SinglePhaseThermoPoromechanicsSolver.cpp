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
 * @file SinglePhaseThermoPoromechanicsSolver.cpp
 *
 */

#include "SinglePhaseThermoPoromechanicsSolver.hpp"
#include "SinglePhaseThermoPoromechanicsKernel.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhaseThermoPoromechanicsSolver::SinglePhaseThermoPoromechanicsSolver( string const & name,
                                                                            Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName()
{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver to use in the ThermoPoromechanics solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fluid mechanics solver to use in the ThermoPoromechanics solver" );

  registerWrapper( viewKeyStruct::porousMaterialNamesString(), &m_porousMaterialNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the material that should be used in the constitutive updates" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhaseThermoPoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}


void SinglePhaseThermoPoromechanicsSolver::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    NodeManager & nodeManager = meshBody.getMeshLevel( 0 ).getNodeManager();

    nodeManager.registerWrapper< array1d< real64 > >( viewKeyStruct::temperatureString() ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setDescription( "An array that holds the nodal temperature" );

    nodeManager.registerWrapper< array1d< real64 > >( viewKeyStruct::newDeltaTemperatureString() ).
      setPlotLevel( PlotLevel::LEVEL_3 ).
      setDescription( "An array that holds the nodal temperature increment" );

    nodeManager.registerWrapper< array1d< real64 > >( viewKeyStruct::oldDeltaTemperatureString() ).
      setPlotLevel( PlotLevel::LEVEL_3 ).
      setDescription( "An array that holds the nodal temperature increment of the previous time step" );

  } );

}


void SinglePhaseThermoPoromechanicsSolver::setupDofs( DomainPartition const & domain,
                                                      DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->setupDofs( domain, dofManager );
  m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString(),
                          DofManager::Connector::Elem );

  // setup dof for the temperature field

  dofManager.addField( viewKeyStruct::temperatureString(),
                       DofManager::Location::Node,
                       1,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::temperatureString(),
                          viewKeyStruct::temperatureString(),
                          DofManager::Connector::Elem );
}

void SinglePhaseThermoPoromechanicsSolver::setupSystem( DomainPartition & domain,
                                                        DofManager & dofManager,
                                                        CRSMatrix< real64, globalIndex > & localMatrix,
                                                        array1d< real64 > & localRhs,
                                                        array1d< real64 > & localSolution,
                                                        bool const setSparsity )
{
  if( m_precond )
  {
    m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner();
  }
}

void SinglePhaseThermoPoromechanicsSolver::implicitStepSetup( real64 const & time_n,
                                                              real64 const & dt,
                                                              DomainPartition & domain )
{
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_solidSolver->implicitStepSetup( time_n, dt, domain );
}

void SinglePhaseThermoPoromechanicsSolver::implicitStepComplete( real64 const & time_n,
                                                                 real64 const & dt,
                                                                 DomainPartition & domain )
{
  m_solidSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );
}

void SinglePhaseThermoPoromechanicsSolver::postProcessInput()
{
  SolverBase::postProcessInput();

  m_flowSolver = &this->getParent().getGroup< SinglePhaseBase >( m_flowSolverName );
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
}

void SinglePhaseThermoPoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{
/**
   if( m_flowSolver->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
   {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhaseThermoPoromechanics;
   }
 */
}

SinglePhaseThermoPoromechanicsSolver::~SinglePhaseThermoPoromechanicsSolver()
{}

void SinglePhaseThermoPoromechanicsSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_solidSolver->resetStateToBeginningOfStep( domain );

  // Reset thermal state
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real64 > const &
  temperature = nodeManager.template getReference< array1d< real64 > >( viewKeyStruct::temperatureString() );

  arrayView1d< real64 > const &
  newDeltaTemperature = nodeManager.template getReference< array1d< real64 > >( viewKeyStruct::newDeltaTemperatureString() );

  arrayView1d< real64 > const &
  oldDeltaTemperature = nodeManager.template getReference< array1d< real64 > >( viewKeyStruct::oldDeltaTemperatureString() );

  forAll< parallelDevicePolicy< 32 > >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    temperature( a ) -= newDeltaTemperature( a );
    oldDeltaTemperature( a ) = newDeltaTemperature( a );
    newDeltaTemperature( a ) = 0.0;
  } );
}

real64 SinglePhaseThermoPoromechanicsSolver::solverStep( real64 const & time_n,
                                                         real64 const & dt,
                                                         int const cycleNumber,
                                                         DomainPartition & domain )
{
  real64 dt_return = dt;

  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_localRhs,
               m_localSolution );

  implicitStepSetup( time_n, dt, domain );

  dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void SinglePhaseThermoPoromechanicsSolver::assembleSystem( real64 const time_n,
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBodies().getGroup< MeshBody >( 0 ).getMeshLevel( 0 );

  NodeManager const & nodeManager = mesh.getNodeManager();

  string const dofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

  string const pressureDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  string const temperatureDofKey = dofManager.getKey( viewKeyStruct::temperatureString() );
  arrayView1d< globalIndex const > const & temperatureDofNumber = nodeManager.getReference< globalIndex_array >( temperatureDofKey );

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  ThermoPoromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                                      pressureDofKey,
                                                                      temperatureDofNumber,
                                                                      dofManager.rankOffset(),
                                                                      localMatrix,
                                                                      localRhs,
                                                                      gravityVectorData,
                                                                      m_flowSolver->fluidModelNames(),
                                                                      dt );

  // Cell-based contributions
  m_solidSolver->getMaxForce() =
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                    constitutive::PorousSolidBase,
                                    CellElementSubRegion >( mesh,
                                                            targetRegionNames(),
                                                            this->getDiscretizationName(),
                                                            porousMaterialNames(),
                                                            kernelFactory );

  // Contributions of the fluid flux term
  m_flowSolver->assemblePoroelasticFluxTerms( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              " " );
}

void SinglePhaseThermoPoromechanicsSolver::applyBoundaryConditions( real64 const time_n,
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

  // Apply boundary temperature
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply( time_n + dt,
                   domain,
                   "nodeManager",
                   viewKeyStruct::temperatureString(),
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group & targetGroup,
                        string const & GEOSX_UNUSED_PARAM( fieldName ) )
  {
    bc.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                       parallelDevicePolicy< 32 > >( targetSet,
                                                                     time_n + dt,
                                                                     targetGroup,
                                                                     viewKeyStruct::temperatureString(),
                                                                     dofManager.getKey( viewKeyStruct::temperatureString() ),
                                                                     dofManager.rankOffset(),
                                                                     localMatrix,
                                                                     localRhs );
  } );
}

real64 SinglePhaseThermoPoromechanicsSolver::calculateResidualNorm( DomainPartition const & domain,
                                                                    DofManager const & dofManager,
                                                                    arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm = m_solidSolver->calculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = m_flowSolver->calculateResidualNorm( domain, dofManager, localRhs );

  GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "( Rsolid, Rfluid ) = ( {:4.2e}, {:4.2e} )", momementumResidualNorm, massResidualNorm ) );

  return sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm );
}

void SinglePhaseThermoPoromechanicsSolver::createPreconditioner()
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
                         { { SinglePhaseBase::viewKeyStruct::pressureString(), { 1, true } } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void SinglePhaseThermoPoromechanicsSolver::solveSystem( DofManager const & dofManager,
                                                        ParallelMatrix & matrix,
                                                        ParallelVector & rhs,
                                                        ParallelVector & solution )
{
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void SinglePhaseThermoPoromechanicsSolver::applySystemSolution( DofManager const & dofManager,
                                                                arrayView1d< real64 const > const & localSolution,
                                                                real64 const scalingFactor,
                                                                DomainPartition & domain )
{
  // update displacement field
  m_solidSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );

  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );


  // update temperature field
  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::temperatureString(),
                               viewKeyStruct::temperatureString(),
                               -scalingFactor );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::temperatureString(),
                               viewKeyStruct::newDeltaTemperatureString(),
                               -scalingFactor );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( string( viewKeyStruct::temperatureString() ) );
  fieldNames["node"].emplace_back( string( viewKeyStruct::newDeltaTemperatureString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
}

void SinglePhaseThermoPoromechanicsSolver::updateState( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  this->template forTargetSubRegions< CellElementSubRegion >( mesh, [&] ( localIndex const targetIndex,
                                                                          auto & subRegion )
  {
    m_flowSolver->updateFluidState( subRegion, targetIndex );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseThermoPoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
