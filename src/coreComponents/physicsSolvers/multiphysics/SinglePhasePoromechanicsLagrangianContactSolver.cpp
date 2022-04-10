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
//#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
//#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
//#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/contact/LagrangianContactSolver.hpp"
#include "physicsSolvers/multiphysics/LagrangianContactFlowSolver.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "SinglePhasePoromechanicsKernel.hpp"
#include "math/interpolation/Interpolation.hpp"

#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/solvers/PreconditionerJacobi.hpp"
#include "linearAlgebra/solvers/PreconditionerBlockJacobi.hpp"
#include "linearAlgebra/solvers/BlockPreconditionerGeneral.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace interpolation;

SinglePhasePoromechanicsLagrangianContactSolver::SinglePhasePoromechanicsLagrangianContactSolver( const string & name,
                                                                                                  Group * const parent ):
  SinglePhasePoromechanicsSolver( name, parent ),
  m_contactSolverName()
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

  // TODO: stampare m_meshTargets
  // question: from here...
  dofManager.addField( keys::TotalDisplacement,
                       DofManager::Location::Node,
                       3,
                       m_meshTargets );

  dofManager.addCoupling( keys::TotalDisplacement,
                          keys::TotalDisplacement,
                          DofManager::Connector::Elem );
  // ... to here. Can we replace simply by m_solidSolver->setupDofs( domain, dofManager );   ??
  // OR, considering also the coupling u-t can we simply call m_contactSolver->setupDofs( domain, dofManager ); ??

  // restrict coupling to fracture regions only
  map< string, array1d< string > > meshTargets;
  forMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                MeshLevel const & meshLevel,
                                                arrayView1d< string const > const & regionNames )
  {
    array1d< string > regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< SurfaceElementRegion >( regionNames,
                                                                    [&]( localIndex const,
                                                                         SurfaceElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    meshTargets[meshBodyName] = std::move( regions );
  } );

  dofManager.addField( extrinsicMeshData::contact::traction::key(),
                       DofManager::Location::Elem,
                       3,
                       meshTargets );
//                       fractureRegions );
  dofManager.addCoupling( extrinsicMeshData::contact::traction::key(),
                          extrinsicMeshData::contact::traction::key(),
                          DofManager::Connector::Face,
                          meshTargets );
//                          fractureRegions );
  dofManager.addCoupling( keys::TotalDisplacement,
                          extrinsicMeshData::contact::traction::key(),
                          DofManager::Connector::Elem,
                          meshTargets );
//                          fractureRegions );

  m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          extrinsicMeshData::flow::pressure::key(),
                          DofManager::Connector::Elem );
  dofManager.addCoupling( extrinsicMeshData::flow::pressure::key(),
                          extrinsicMeshData::contact::traction::key(),
                          DofManager::Connector::None );
//  dofManager.addCoupling( extrinsicMeshData::contact::traction::key(),
//                          FlowSolverBase::viewKeyStruct::pressureString(),
//                          DofManager::Connector::None );

  dofManager.printFieldInfo();
}

bool SinglePhasePoromechanicsLagrangianContactSolver::updateConfiguration( DomainPartition & domain )
{
  return m_contactSolver->updateConfiguration( domain );
}

bool SinglePhasePoromechanicsLagrangianContactSolver::resetConfigurationToDefault( DomainPartition & domain ) const
{
  return m_contactSolver->resetConfigurationToDefault( domain );
}

void SinglePhasePoromechanicsLagrangianContactSolver::setupSystem( DomainPartition & domain,
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
  SinglePhasePoromechanicsSolver::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  m_contactFlowSolver->setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  //ParallelMatrix MM;
  //MM.create( localMatrix.toViewConst(), MPI_COMM_GEOSX );
  //MM.print();

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner( domain );
  }
}

void SinglePhasePoromechanicsLagrangianContactSolver::implicitStepSetup( real64 const & time_n,
                                                                         real64 const & dt,
                                                                         DomainPartition & domain )
{
  m_contactSolver->implicitStepSetup( time_n, dt, domain );
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  //m_contactFlowSolver->implicitStepSetup( time_n, dt, domain );
}

void SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete( real64 const & time_n,
                                                                            real64 const & dt,
                                                                            DomainPartition & domain )
{
  m_contactSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      string const porousMaterialName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      ConstitutiveBase const & porousMaterial = subRegion.getConstitutiveModel< ConstitutiveBase >( porousMaterialName );
      porousMaterial.saveConvergedState();
    } );

    // Laura print max pressure and displacement
    // pressure
    elemManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion & region )
    {
      region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
      {
        if( subRegion.hasWrapper( extrinsicMeshData::flow::pressure::key() ) )
        {
          arrayView1d< real64 > pres = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
          double * max_pres = std::max_element( pres.begin(), pres.end());
          GEOSX_LOG_RANK_0( GEOSX_FMT( "SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete -- max pres {:15.6e}", *max_pres ) );
          double * min_pres = std::min_element( pres.begin(), pres.end());
          GEOSX_LOG_RANK_0( GEOSX_FMT( "SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete -- min pres {:15.6e}", *min_pres ) );
        }
      } );
    } );
    // displacement
    NodeManager & nodeManager = mesh.getNodeManager();
    arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const disp = nodeManager.totalDisplacement();
    double * min_disp = std::min_element( disp.begin(), disp.end());
    GEOSX_LOG_RANK_0( GEOSX_FMT( "SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete -- min disp {:15.6e}", *min_disp ) );

    real64 const totalFlux = m_flowSolver->computeFluxFaceDirichlet( time_n, dt, m_dofManager, domain );
    GEOSX_LOG_RANK_0( GEOSX_FMT( "SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete -- total flux through Dirichlet faces {:15.6e}", totalFlux ) );
    // end Laura
  } );
}

void SinglePhasePoromechanicsLagrangianContactSolver::postProcessInput()
{
  SinglePhasePoromechanicsSolver::postProcessInput();

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
// Laura - this is call is already in the m_contactFlowSolver->resetStateToBeginningOfStep( domain );
//  m_flowSolver->resetStateToBeginningOfStep( domain );
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
               m_rhs,
               m_solution );

  implicitStepSetup( time_n, dt, domain );

  dt_return = SolverBase::nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

bool SinglePhasePoromechanicsLagrangianContactSolver::lineSearch( real64 const & time_n,
                                                                  real64 const & dt,
                                                                  integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                                                  DomainPartition & domain,
                                                                  DofManager const & dofManager,
                                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                  ParallelVector & rhs,
                                                                  ParallelVector & solution,
                                                                  real64 const scaleFactor,
                                                                  real64 & lastResidual )
{
  bool lineSearchSuccess = true;

  integer const maxNumberLineSearchCuts = m_nonlinearSolverParameters.m_lineSearchMaxCuts;

  real64 const sigma1 = 0.5;
  real64 const alpha = 1.e-4;

  real64 localScaleFactor = scaleFactor;
  real64 lamm = scaleFactor;
  real64 lamc = localScaleFactor;
  integer lineSearchIteration = 0;

  // get residual norm
  real64 residualNorm0 = lastResidual;

  applySystemSolution( dofManager, solution.values(), scaleFactor, domain );

  // re-assemble system
  localMatrix.zero();
  rhs.zero();

  {
    arrayView1d< real64 > const localRhs = rhs.open();
    assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
    rhs.close();
  }

  // get residual norm
  real64 residualNormT = calculateResidualNorm( domain, dofManager, rhs.values() );

  real64 ff0 = residualNorm0*residualNorm0;
  real64 ffT = residualNormT*residualNormT;
  real64 ffm = ffT;
  real64 cumulativeScale = scaleFactor;

  while( residualNormT >= (1.0 - alpha*localScaleFactor)*residualNorm0 )
  {
    real64 const previousLocalScaleFactor = localScaleFactor;
    // Apply the three point parabolic model
    if( lineSearchIteration == 0 )
    {
      localScaleFactor *= sigma1;
    }
    else
    {
      localScaleFactor = parabolicInterpolationThreePoints( lamc, lamm, ff0, ffT, ffm );
    }

    // Update x; keep the books on lambda
    real64 const deltaLocalScaleFactor = ( localScaleFactor - previousLocalScaleFactor );
    cumulativeScale += deltaLocalScaleFactor;

    if( !checkSystemSolution( domain, dofManager, solution.values(), deltaLocalScaleFactor ) )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search " << lineSearchIteration << ", solution check failed" );
      continue;
    }

    applySystemSolution( dofManager, solution.values(), deltaLocalScaleFactor, domain );
    lamm = lamc;
    lamc = localScaleFactor;

    // Keep the books on the function norms
    // re-assemble system
    // TODO: add a flag to avoid a completely useless Jacobian computation: rhs is enough
    localMatrix.zero();
    rhs.zero();

    {
      arrayView1d< real64 > const localRhs = rhs.open();
      assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
      applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
      rhs.close();
    }

    if( getLogLevel() >= 1 && logger::internal::rank==0 )
    {
      char output[100];
      sprintf( output, "        Line search @ %0.3f:      ", cumulativeScale );
      std::cout<<output;
    }

    // get residual norm
    residualNormT = calculateResidualNorm( domain, dofManager, rhs.values() );
    ffm = ffT;
    ffT = residualNormT*residualNormT;
    lineSearchIteration += 1;

    if( lineSearchIteration > maxNumberLineSearchCuts )
    {
      lineSearchSuccess = false;
      break;
    }
  }

  lastResidual = residualNormT;

  return lineSearchSuccess;
}

void SinglePhasePoromechanicsLagrangianContactSolver::assembleSystem( real64 const time_n,
                                                                      real64 const dt,
                                                                      DomainPartition & domain,
                                                                      DofManager const & dofManager,
                                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                      arrayView1d< real64 > const & localRhs )
{

  GEOSX_MARK_FUNCTION;

  // Need to synchronize the two iteration counters
  m_contactSolver->getNonlinearSolverParameters().m_numNewtonIterations = m_nonlinearSolverParameters.m_numNewtonIterations;

  // TODO: synchronizeFractureState ? ripristinare seguente istruzione
  //m_contactSolver->synchronizeFractureState( domain );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();

    string const dofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
    arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

    string const pDofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );

    // TODO TMP to be used if no coupling terms are assembled
    //GEOSX_UNUSED_VAR( regionNames );
    //m_contactSolver->assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );

    m_contactSolver->assembleForceResidualDerivativeWrtTraction( domain, dofManager, localMatrix, localRhs );
    m_contactSolver->assembleTractionResidualDerivativeWrtDisplacementAndTraction( domain, dofManager, localMatrix, localRhs );
    m_contactSolver->assembleStabilization( domain, dofManager, localMatrix, localRhs );

    m_contactFlowSolver->getRefDerivativeFluxResidual_dAperture()->zero();
    m_contactFlowSolver->updateOpeningForFlow( domain );

    m_contactFlowSolver->assembleForceResidualDerivativeWrtPressure( domain, dofManager, localMatrix, localRhs );
    m_contactFlowSolver->assembleFluidMassResidualDerivativeWrtDisplacement( domain, dofManager, localMatrix, localRhs );
    // No need for stabilization
    //m_contactFlowSolver->assembleStabilization( domain, dofManager, localMatrix, localRhs );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    poromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                                  pDofKey,
                                                                  dofManager.rankOffset(),
                                                                  localMatrix,
                                                                  localRhs,
                                                                  gravityVectorData,
                                                                  FlowSolverBase::viewKeyStruct::fluidNamesString() );

    // Cell-based contributions: coupling terms (u/p)
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                    constitutive::PorousSolidBase,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            this->getDiscretizationName(),
                                                            viewKeyStruct::porousMaterialNamesString(),
                                                            kernelFactory );

  } );

  // Transmissibility 3D/2D
  // FIXME: changes in assembleHydrofracFluxTerms to be checked
  m_flowSolver->assembleHydrofracFluxTerms( time_n,
                                            dt,
                                            domain,
                                            dofManager,
                                            localMatrix,
                                            localRhs,
                                            m_contactFlowSolver->getDerivativeFluxResidual_dAperture() );

  // Assemble transmissibility entries
  m_flowSolver->assemblePoroelasticFluxTerms( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              dofManager.getKey( extrinsicMeshData::flow::pressure::key() ) );

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
    sprintf( output, "    ( Rsolid, Rfluid, Rtotal ) = ( %4.2e, %4.2e, %4.2e )", momementumResidualNorm, massResidualNorm,
             sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm ) );
    std::cout << output << std::endl;
  }

  return sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm );
}

void SinglePhasePoromechanicsLagrangianContactSolver::createPreconditioner( DomainPartition const & domain )
{
  //m_contactFlowSolver->createPreconditioner( domain );
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    // TODO: move among inputs (xml)
    string const leadingBlockApproximation = "blockJacobi";
    //string const leadingBlockApproximation = "jacobi";

    LinearSolverParameters mechParams = m_contactSolver->getSolidSolver()->getLinearSolverParameters();
    // Because of boundary conditions
    mechParams.isSymmetric = false;

    std::vector< SchurComplementOption > schurOptions( 2 );
    std::unique_ptr< BlockPreconditionerGeneral< LAInterface > > precond;
    std::unique_ptr< PreconditionerBase< LAInterface > > tracPrecond;
    std::unique_ptr< PreconditionerBase< LAInterface > > flowPrecond;

    if( leadingBlockApproximation == "jacobi" )
    {
      schurOptions[0] = SchurComplementOption::Diagonal;

      // Using LAI implementation of Jacobi preconditioner
      LinearSolverParameters tracParams;
      tracParams.preconditionerType = LinearSolverParameters::PreconditionerType::jacobi;
      tracPrecond = LAInterface::createPreconditioner( tracParams );
    }
    else if( leadingBlockApproximation == "blockJacobi" )
    {
      schurOptions[0] = SchurComplementOption::UserDefined;

      tracPrecond = std::make_unique< PreconditionerBlockJacobi< LAInterface > >( mechParams.dofsPerNode );
    }
    else
    {
      GEOSX_ERROR( "SinglePhasePoromechanicsLagrangianContactSolver::CreatePreconditioner leadingBlockApproximation option " << leadingBlockApproximation << " not supported" );
    }

    // Flow + Jacobi: using LAI implementation of Jacobi preconditioner
    schurOptions[1] = SchurComplementOption::Diagonal;

    LinearSolverParameters flowParams = m_flowSolver->getLinearSolverParameters();
    flowPrecond = LAInterface::createPreconditioner( flowParams );

    precond = std::make_unique< BlockPreconditionerGeneral< LAInterface > >( 3,
                                                                             BlockShapeOption::LowerUpperTriangular,
                                                                             schurOptions );

    precond->setupBlock( 0,
                         { { extrinsicMeshData::contact::traction::key(), { 3, 0, 3 } } },
                         std::move( tracPrecond ) );

    precond->setupBlock( 1,
                         { { extrinsicMeshData::flow::pressure::key(), { 1, 0, 1 } } },
                         std::move( flowPrecond ) );

    if( mechParams.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
    {
      if( m_contactSolver->getSolidSolver()->getRigidBodyModes().empty() )
      {
        MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
        LAIHelperFunctions::computeRigidBodyModes( mesh,
                                                   m_dofManager,
                                                   { keys::TotalDisplacement },
                                                   m_contactSolver->getSolidSolver()->getRigidBodyModes() );
      }
    }

    // Preconditioner for the Schur complement: mechPrecond
    std::unique_ptr< PreconditionerBase< LAInterface > > mechPrecond = LAInterface::createPreconditioner( mechParams, m_contactSolver->getSolidSolver()->getRigidBodyModes() );
    precond->setupBlock( 2,
                         { { keys::TotalDisplacement, { 3, 0, 3 } } },
                         std::move( mechPrecond ) );

    m_precond = std::move( precond );
  }
  else if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block_fs )
  {
    // TODO: move among inputs (xml)
    string const leadingBlockApproximation = "blockJacobi";
    //string const leadingBlockApproximation = "jacobi";

    LinearSolverParameters mechParams = m_contactSolver->getSolidSolver()->getLinearSolverParameters();
    // Because of boundary conditions
    mechParams.isSymmetric = false;

    std::vector< SchurComplementOption > schurOptions( 2 );
    std::unique_ptr< BlockPreconditionerGeneral< LAInterface > > precond;
    std::unique_ptr< PreconditionerBase< LAInterface > > tracPrecond;
    std::unique_ptr< PreconditionerBase< LAInterface > > flowPrecond;

    if( leadingBlockApproximation == "jacobi" )
    {
      schurOptions[0] = SchurComplementOption::Diagonal;

      // Using LAI implementation of Jacobi preconditioner
      LinearSolverParameters tracParams;
      tracParams.preconditionerType = LinearSolverParameters::PreconditionerType::jacobi;
      tracPrecond = LAInterface::createPreconditioner( tracParams );
    }
    else if( leadingBlockApproximation == "blockJacobi" )
    {
      schurOptions[0] = SchurComplementOption::UserDefined;

      tracPrecond = std::make_unique< PreconditionerBlockJacobi< LAInterface > >( mechParams.dofsPerNode );
    }
    else
    {
      GEOSX_ERROR( "SinglePhasePoromechanicsLagrangianContactSolver::CreatePreconditioner leadingBlockApproximation option " << leadingBlockApproximation << " not supported" );
    }

    // Mechanics + Jacobi: using LAI implementation of Jacobi preconditioner
    schurOptions[1] = SchurComplementOption::Diagonal;

    LinearSolverParameters flowParams = m_flowSolver->getLinearSolverParameters();
    flowPrecond = LAInterface::createPreconditioner( flowParams );

    precond = std::make_unique< BlockPreconditionerGeneral< LAInterface > >( 3,
                                                                             BlockShapeOption::LowerUpperTriangular,
                                                                             schurOptions );

    precond->setupBlock( 0,
                         { { extrinsicMeshData::contact::traction::key(), { 3, 0, 3 } } },
                         std::move( tracPrecond ) );

    if( mechParams.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
    {
      if( m_contactSolver->getSolidSolver()->getRigidBodyModes().empty() )
      {
        MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
        LAIHelperFunctions::computeRigidBodyModes( mesh,
                                                   m_dofManager,
                                                   { keys::TotalDisplacement },
                                                   m_contactSolver->getSolidSolver()->getRigidBodyModes() );
      }
    }

    // Preconditioner for the Schur complement: mechPrecond
    std::unique_ptr< PreconditionerBase< LAInterface > > mechPrecond = LAInterface::createPreconditioner( mechParams, m_contactSolver->getSolidSolver()->getRigidBodyModes() );
    precond->setupBlock( 1,
                         { { keys::TotalDisplacement, { 3, 0, 3 } } },
                         std::move( mechPrecond ) );

    precond->setupBlock( 2,
                         { { extrinsicMeshData::flow::pressure::key(), { 1, 0, 1 } } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
  /*
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
                         { { extrinsicMeshData::flow::pressure::key(), { 1, true } } },
                         std::move( flowPrecond ) );

     m_precond = std::move( precond );
     }
     else
     {
     //TODO: Revisit this part such that is coherent across physics solver
     //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
     }
   */
}

void SinglePhasePoromechanicsLagrangianContactSolver::solveLinearSystem( DofManager const & dofManager,
                                                                         ParallelMatrix & matrix,
                                                                         ParallelVector & rhs,
                                                                         ParallelVector & solution )
{
  SinglePhasePoromechanicsSolver::solveLinearSystem( dofManager, matrix, rhs, solution );

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
  /*
     MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
     ElementRegionManager & elemManager = mesh.getElemManager();
     elemManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion & region )
     {
     region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
     {
      if( subRegion.hasWrapper( SinglePhaseBase::viewKeyStruct::pressureString() ) )
      {
        arrayView1d< real64 > pres = subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::pressureString() );
        arrayView1d< real64 > dpres = subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaPressureString() );
        std::cout << "pres0 " << pres[0] << " " << dpres[0] << std::endl;
      }
     } );
     } );
   */

  // update displacement field
  m_contactSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

void SinglePhasePoromechanicsLagrangianContactSolver::updateState( DomainPartition & domain )
{
//  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
//  this->template forTargetSubRegions< CellElementSubRegion >( mesh, [&] ( localIndex const targetIndex,
//                                                                          auto & subRegion )
//  {
//    m_flowSolver->updateFluidState( subRegion, targetIndex );
//  } );

  m_contactSolver->updateState( domain );

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

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsLagrangianContactSolver, string const &, Group * const )

} /* namespace geosx */
