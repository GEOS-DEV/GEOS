/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanicsLagrangianContactSolver.cpp
 *
 */


#include "SinglePhasePoromechanicsLagrangianContactSolver.hpp"

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
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
//#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/multiphysics/LagrangianContactSolver.hpp"
#include "physicsSolvers/multiphysics/LagrangianContactFlowSolver.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "SinglePhasePoromechanicsKernel.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhasePoromechanicsLagrangianContactSolver::SinglePhasePoromechanicsLagrangianContactSolver( const string & name,
                                                                Group * const parent ):
  SolverBase( name, parent ),
  m_contactSolverName(),
  m_flowSolverName()

{
  registerWrapper( viewKeyStruct::contactSolverNameString(), &m_contactSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the contact mechanics solver to use in the poromechanics solver" );

  registerWrapper( viewKeyStruct::contactFlowSolverNameString(), &m_contactFlowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the contact mechanics and flow solver to use in the poromechanics solver" );

  registerWrapper( viewKeyStruct::flowSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fluid mechanics solver to use in the poromechanics solver" );

  registerWrapper( viewKeyStruct::porousMaterialNamesString(), &m_porousMaterialNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the material that should be used in the constitutive updates" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void SinglePhasePoromechanicsLagrangianContactSolver::setupDofs( DomainPartition const & domain,
                                                DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  //m_contactFlowSolver->setupDofs( domain, dofManager );
  //m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addField( keys::TotalDisplacement,
                       DofManager::Location::Node,
                       3,
                       targetRegionNames() );

  dofManager.addCoupling( keys::TotalDisplacement,
                          keys::TotalDisplacement,
                          DofManager::Connector::Elem );

  // restrict coupling to fracture regions only
  ElementRegionManager const & elemManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getElemManager();
  string_array fractureRegions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elementRegion )
  {
    fractureRegions.emplace_back( elementRegion.getName() );
  } );

  dofManager.addField( LagrangianContactSolver::viewKeyStruct::tractionString(),
                       DofManager::Location::Elem,
                       3,
                       fractureRegions );
  dofManager.addCoupling( LagrangianContactSolver::viewKeyStruct::tractionString(),
                          LagrangianContactSolver::viewKeyStruct::tractionString(),
                          DofManager::Connector::Face,
                          fractureRegions );
  dofManager.addCoupling( keys::TotalDisplacement,
                          LagrangianContactSolver::viewKeyStruct::tractionString(),
                          DofManager::Connector::Elem,
                          fractureRegions );

  m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString(),
                          DofManager::Connector::Elem );
  dofManager.addCoupling( FlowSolverBase::viewKeyStruct::pressureString(),
                          LagrangianContactSolver::viewKeyStruct::tractionString(),
                          DofManager::Connector::None );
//  dofManager.addCoupling( LagrangianContactSolver::viewKeyStruct::tractionString(),
//                          FlowSolverBase::viewKeyStruct::pressureString(),
//                          DofManager::Connector::None );

//  dofManager.addCoupling( FlowSolverBase::viewKeyStruct::pressureString(),
//                          LagrangianContactSolver::viewKeyStruct::tractionString(),
//                          DofManager::Connector::None );
//  dofManager.addCoupling( LagrangianContactSolver::viewKeyStruct::tractionString(),
//                          FlowSolverBase::viewKeyStruct::pressureString(),
//                          DofManager::Connector::None );

  dofManager.printFieldInfo();
}

void SinglePhasePoromechanicsLagrangianContactSolver::setupSystem( DomainPartition & domain,
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

  m_contactFlowSolver->setupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );

  //ParallelMatrix MM;
  //MM.create( localMatrix.toViewConst(), MPI_COMM_GEOSX );
  //MM.print();

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner();
  }
}

void SinglePhasePoromechanicsLagrangianContactSolver::implicitStepSetup( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_contactFlowSolver->implicitStepSetup( time_n, dt, domain );
}

void SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete( real64 const & time_n,
                                                           real64 const & dt,
                                                           DomainPartition & domain )
{
  m_contactFlowSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, porousMaterialNames()[targetIndex] );
    porousMaterial.saveConvergedState();
  } );
}

void SinglePhasePoromechanicsLagrangianContactSolver::postProcessInput()
{
  SolverBase::postProcessInput();

  m_flowSolver = &this->getParent().getGroup< SinglePhaseBase >( m_flowSolverName );
  //m_contactFlowSolver = &this->getParent().getGroup< LagrangianContactSolver >( m_contactSolverName );
  m_contactSolver = &this->getParent().getGroup< LagrangianContactSolver >( m_contactSolverName );
  m_contactFlowSolver = &this->getParent().getGroup< LagrangianContactFlowSolver >( m_contactFlowSolverName );
}

void SinglePhasePoromechanicsLagrangianContactSolver::initializePostInitialConditionsPreSubGroups()
{
  if( m_flowSolver->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics;
  }
}

SinglePhasePoromechanicsLagrangianContactSolver::~SinglePhasePoromechanicsLagrangianContactSolver()
{
  // TODO Auto-generated destructor stub
}

void SinglePhasePoromechanicsLagrangianContactSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_contactFlowSolver->resetStateToBeginningOfStep( domain );
}

real64 SinglePhasePoromechanicsLagrangianContactSolver::solverStep( real64 const & time_n,
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

void SinglePhasePoromechanicsLagrangianContactSolver::assembleSystem( real64 const time_n,
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

  string const pDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  // Assemble K twice!!!
  //m_contactFlowSolver->assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );

  m_contactSolver->assembleForceResidualDerivativeWrtTraction( domain, dofManager, localMatrix, localRhs );
  m_contactSolver->assembleTractionResidualDerivativeWrtDisplacementAndTraction( domain, dofManager, localMatrix, localRhs );
  m_contactSolver->assembleStabilization( domain, dofManager, localMatrix, localRhs );



  m_contactFlowSolver->updateOpeningForFlow( domain );

  //m_flowSolver->assembleSystem( time,
  //                              dt,
  //                              domain,
  //                              dofManager,
  //                              localMatrix,
  //                              localRhs );

  ElementRegionManager const & elemManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getElemManager();
  string_array fractureRegions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elementRegion )
  {
    fractureRegions.emplace_back( elementRegion.getName() );
  } );

  m_flowSolver->assembleAccumulationTerms< parallelDevicePolicy<> >( domain,
                                                                     dofManager,
                                                                     localMatrix,
                                                                     localRhs );
  m_flowSolver->assembleHydrofracFluxTerms( time_n,
                                            dt,
                                            domain,
                                            dofManager,
                                            localMatrix,
                                            localRhs,
                                            m_contactFlowSolver->getDerivativeFluxResidual_dAperture(),
                                            fractureRegions );

  m_contactFlowSolver->assembleForceResidualDerivativeWrtPressure( domain, dofManager, localMatrix, localRhs );
  m_contactFlowSolver->assembleFluidMassResidualDerivativeWrtDisplacement( domain, dofManager, localMatrix, localRhs );
  m_contactFlowSolver->assembleStabilization( domain, dofManager, localMatrix, localRhs );

  m_contactFlowSolver->getRefDerivativeFluxResidual_dAperture()->zero();


//  m_contactFlowSolver->resetStressToBeginningOfStep( domain );

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  PoromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                                pDofKey,
                                                                dofManager.rankOffset(),
                                                                localMatrix,
                                                                localRhs,
                                                                gravityVectorData,
                                                                m_flowSolver->fluidModelNames() );

  // Cell-based contributions
  //m_contactFlowSolver->getMaxForce() =
  real64 maxForce = 
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                    constitutive::PorousSolidBase,
                                    CellElementSubRegion >( mesh,
                                                            targetRegionNames(),
                                                            this->getDiscretizationName(),
                                                            porousMaterialNames(),
                                                            kernelFactory );

  GEOSX_UNUSED_VAR( maxForce );

  m_flowSolver->assemblePoroelasticFluxTerms( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              dofManager.getKey( SinglePhaseBase::viewKeyStruct::pressureString() ) );
                                              //" " );
                                              //dofManager.getKey( LagrangianContactSolver::viewKeyStruct::tractionString() ) );
}

void SinglePhasePoromechanicsLagrangianContactSolver::applyBoundaryConditions( real64 const time_n,
                                                              real64 const dt,
                                                              DomainPartition & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  m_contactSolver->applyBoundaryConditions( time_n, dt,
                                            domain,
                                            dofManager,
                                            localMatrix,
                                            localRhs );

  //ParallelMatrix MM;
  //MM.create( localMatrix.toViewConst(), MPI_COMM_GEOSX );
  //std::cout << MM << std::endl;


  m_flowSolver->applyBoundaryConditions( time_n, dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

real64 SinglePhasePoromechanicsLagrangianContactSolver::calculateResidualNorm( DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm = m_contactSolver->calculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = m_flowSolver->calculateResidualNorm( domain, dofManager, localRhs );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output, "    ( Rsolid, Rfluid ) = ( %4.2e, %4.2e )", momementumResidualNorm, massResidualNorm );
    std::cout << output << std::endl;
  }

  return sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm );
}

void SinglePhasePoromechanicsLagrangianContactSolver::createPreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( m_contactFlowSolver->getLinearSolverParameters() );
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

void SinglePhasePoromechanicsLagrangianContactSolver::solveSystem( DofManager const & dofManager,
                                                  ParallelMatrix & matrix,
                                                  ParallelVector & rhs,
                                                  ParallelVector & solution )
{
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );

  int rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  if( rank == 0 )
  {
    string str;
    std::getline( std::cin, str );
    if( str.length() > 0 )
    {
      debugOutputSolution( 99, 99, 99, solution );
      GEOSX_ERROR( "STOP" );
    }
  }
  MpiWrapper::barrier( MPI_COMM_GEOSX );
}

void SinglePhasePoromechanicsLagrangianContactSolver::applySystemSolution( DofManager const & dofManager,
                                                          arrayView1d< real64 const > const & localSolution,
                                                          real64 const scalingFactor,
                                                          DomainPartition & domain )
{
  // update displacement field
  m_contactSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

void SinglePhasePoromechanicsLagrangianContactSolver::updateState( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  this->template forTargetSubRegions< CellElementSubRegion >( mesh, [&] ( localIndex const targetIndex,
                                                                          auto & subRegion )
  {
    m_flowSolver->updateFluidState( subRegion, targetIndex );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsLagrangianContactSolver, string const &, Group * const )

} /* namespace geosx */
