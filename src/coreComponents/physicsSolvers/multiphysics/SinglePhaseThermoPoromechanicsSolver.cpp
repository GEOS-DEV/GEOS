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
 */

#include "SinglePhaseThermoPoromechanicsSolver.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhaseThermoPoromechanicsKernel.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhaseThermoPoromechanicsSolver::SinglePhaseThermoPoromechanicsSolver( const string & name,
                                                                Group * const parent )
  : Base( name, parent )
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void SinglePhaseThermoPoromechanicsSolver::registerDataOnMesh( Group & meshBodies )
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

	// Register thermal fields
	NodeManager & nodeManager = mesh.getNodeManager();

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

void SinglePhaseThermoPoromechanicsSolver::setupCoupling( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                    DofManager & dofManager ) const
{
  dofManager.addCoupling( keys::TotalDisplacement,
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void SinglePhaseThermoPoromechanicsSolver::setupDofs( DomainPartition const & domain,
                                                      DofManager & dofManager ) const
{
  Base::setupDofs(domain, dofManager);

  // setup dof for the temperature field
  dofManager.addField( viewKeyStruct::temperatureString(),
                       FieldLocation::Node,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::temperatureString(),
                          viewKeyStruct::temperatureString(),
                          DofManager::Connector::Elem );
}

void SinglePhaseThermoPoromechanicsSolver::initializePreSubGroups()
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

void SinglePhaseThermoPoromechanicsSolver::setupSystem( DomainPartition & domain,
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

void SinglePhaseThermoPoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{
  if( flowSolver()->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics;
  }
}

void SinglePhaseThermoPoromechanicsSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  Base::resetStateToBeginningOfStep( domain );

  // Reset thermal state
  //MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
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
               m_rhs,
               m_solution );

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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )

  {
    NodeManager const & nodeManager = mesh.getNodeManager();

    string const dofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
    arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

    string const pDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

	string const temperatureDofKey = dofManager.getKey( viewKeyStruct::temperatureString() );
    arrayView1d< globalIndex const > const & temperatureDofNumber = nodeManager.getReference< globalIndex_array >( temperatureDofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    thermoPoromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                                  pDofKey,
																  temperatureDofNumber,
                                                                  dofManager.rankOffset(),
                                                                  localMatrix,
                                                                  localRhs,
                                                                  gravityVectorData,
                                                                  FlowSolverBase::viewKeyStruct::fluidNamesString(),
																  dt );

    // Cell-based contributions
    solidMechanicsSolver()->getMaxForce() =
      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                      constitutive::PorousSolidBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              this->getDiscretizationName(),
                                                              viewKeyStruct::porousMaterialNamesString(),
                                                              kernelFactory );

  } );

  flowSolver()->assemblePoroelasticFluxTerms( time_n, dt,
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
  Base::applyBoundaryConditions( time_n, dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  // Apply boundary temperature
/**
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
*/
FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & )
  {

    fsManager.apply< NodeManager >( time_n + dt,
                                    mesh,
                                    viewKeyStruct::temperatureString(),
                                    [&]( FieldSpecificationBase const & bc,
                                         string const &,
                                         SortedArrayView< localIndex const > const & targetSet,
                                         NodeManager & targetGroup,
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
  } );
}

void SinglePhaseThermoPoromechanicsSolver::createPreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( solidMechanicsSolver()->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { keys::TotalDisplacement, { 3, true } } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( flowSolver()->getLinearSolverParameters() );
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

void SinglePhaseThermoPoromechanicsSolver::applySystemSolution( DofManager const & dofManager,
                                                                arrayView1d< real64 const > const & localSolution,
                                                                real64 const scalingFactor,
                                                                DomainPartition & domain )
{
  // update poromechanical fields
  Base::applySystemSolution( dofManager, localSolution, scalingFactor, domain );

  // update temperature field
  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::temperatureString(),
                               viewKeyStruct::temperatureString(),
                               -scalingFactor );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::temperatureString(),
                               viewKeyStruct::newDeltaTemperatureString(),
                               -scalingFactor );
/**
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( string( viewKeyStruct::temperatureString() ) );
  fieldNames["node"].emplace_back( string( viewKeyStruct::newDeltaTemperatureString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
*/

forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )

  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { viewKeyStruct::temperatureString() } );
	fieldsToBeSync.addFields( FieldLocation::Node, { viewKeyStruct::newDeltaTemperatureString() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );

/**
dofManager.addVectorToField( localSolution,
                               m_fieldName,
                               m_fieldName,
                               scalingFactor );


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )

  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { m_fieldName } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
*/
}

void SinglePhaseThermoPoromechanicsSolver::updateState( DomainPartition & domain )
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

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseThermoPoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
