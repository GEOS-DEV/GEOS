/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsModalAnalysis.cpp
 *
 * Modal analysis (vibration modes) of SolidMechanicsLagrangianFEM:
 * the generalized eigenproblem K phi = lambda M phi, solved with a spectral transformation that
 * requires solving linear systems with K - sigma M. For a negative shift sigma = -alpha, this is
 * the shifted elasticity operator K + alpha M, which is symmetric positive definite even for
 * structures with rigid-body modes.
 */

#include "SolidMechanicsLagrangianFEM.hpp"

#include "common/Stopwatch.hpp"
#include "common/TimingMacros.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "linearAlgebra/solvers/EigenSolverBase.hpp"
#include "linearAlgebra/solvers/KrylovSolver.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/utilities/DiagonalOperator.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

#include <cmath>
#include <functional>

namespace geos
{

using namespace dataRepository;
using namespace fields;

namespace
{

/**
 * @brief Operator applying the inverse of a matrix through a solver that was set up once.
 */
class SolverInverseOperator : public LinearOperator< ParallelVector >
{
public:

  /// Signature of the solve function: solve( rhs, solution )
  using SolveFunction = std::function< void ( ParallelVector const &, ParallelVector & ) >;

  SolverInverseOperator( ParallelMatrix const & matrix, SolveFunction solve ):
    m_matrix( matrix ),
    m_solve( std::move( solve ) )
  {}

  virtual void apply( ParallelVector const & src, ParallelVector & dst ) const override
  {
    dst.zero();
    m_solve( src, dst );
  }

  virtual globalIndex numGlobalRows() const override { return m_matrix.numGlobalRows(); }

  virtual globalIndex numGlobalCols() const override { return m_matrix.numGlobalCols(); }

  virtual localIndex numLocalRows() const override { return m_matrix.numLocalRows(); }

  virtual localIndex numLocalCols() const override { return m_matrix.numLocalCols(); }

  virtual MPI_Comm comm() const override { return m_matrix.comm(); }

private:

  ParallelMatrix const & m_matrix;
  SolveFunction m_solve;
};

} // namespace

globalIndex SolidMechanicsLagrangianFEM::computeModalDiagonals( real64 const time,
                                                                DomainPartition & domain,
                                                                ParallelVector & freeMask,
                                                                ParallelVector & massDiag )
{
  GEOS_MARK_FUNCTION;

  string const dofKey = m_dofManager.getKey( solidMechanics::totalDisplacement::key() );
  globalIndex const rankOffset = m_dofManager.rankOffset();

  freeMask.create( m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
  massDiag.create( m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
  freeMask.set( 1.0 );
  massDiag.zero();

  arrayView1d< real64 > const maskView = freeMask.open();
  arrayView1d< real64 > const massView = massDiag.open();

  // Lumped (row-sum) nodal mass, repeated for the three components
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    NodeManager const & nodes = mesh.getNodeManager();
    arrayView1d< globalIndex const > const dofNumber = nodes.getReference< globalIndex_array >( dofKey );
    arrayView1d< integer const > const ghostRank = nodes.ghostRank();
    arrayView1d< real64 const > const mass = nodes.getField< solidMechanics::mass >();

    forAll< parallelDevicePolicy<> >( nodes.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      if( ghostRank[a] < 0 && dofNumber[a] >= 0 )
      {
        localIndex const row = LvArray::integerConversion< localIndex >( dofNumber[a] - rankOffset );
        for( integer c = 0; c < 3; ++c )
        {
          massView[row + c] = mass[a];
        }
      }
    } );
  } );

  // Constrained degrees of freedom: homogeneous displacement boundary conditions
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    fsManager.apply< NodeManager >( time,
                                    mesh,
                                    solidMechanics::totalDisplacement::key(),
                                    [&]( FieldSpecification const & bc,
                                         string const &,
                                         SortedArrayView< localIndex const > const & targetSet,
                                         NodeManager & targetGroup,
                                         string const GEOS_UNUSED_PARAM( fieldName ) )
    {
      integer const component = bc.getComponent();
      GEOS_ERROR_IF_LT_MSG( component, 0, "Component index required for displacement BC ",
                            getDataContext(), bc.getDataContext() );
      arrayView1d< globalIndex const > const dofNumber = targetGroup.getReference< globalIndex_array >( dofKey );

      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        // Nodes outside the target regions have no degree of freedom (number -1): adding the component to this
        // sentinel would select the degree of freedom of another node
        globalIndex const nodeDof = dofNumber[ targetSet[i] ];
        globalIndex const row = nodeDof + component - rankOffset;
        if( nodeDof >= 0 && row >= 0 && row < maskView.size() )
        {
          maskView[row] = 0.0;
          massView[row] = 0.0;
        }
      } );
    } );
  } );

  freeMask.close();
  massDiag.close();

  return freeMask.globalSize() - LvArray::integerConversion< globalIndex >( std::llround( freeMask.norm1() ) );
}

real64 SolidMechanicsLagrangianFEM::modalAnalysisStep( real64 const & time_n,
                                                       real64 const & dt,
                                                       integer const cycleNumber,
                                                       DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  Stopwatch totalWatch;

  // ---- Linear system: sparsity, tangent stiffness at the current state ----
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName );
  Timestamp const meshTimestamp = LvArray::math::max( getMeshModificationTimestamp( domain ),
                                                      meshLevel.getModificationTimestamp() );
  m_dofManager.clear();
  setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution, true );
  setSystemSetupTimestamp( meshTimestamp );

  implicitStepSetup( time_n, dt, domain );

  m_localMatrix.zero();
  m_rhs.zero();
  {
    arrayView1d< real64 > const localRhs = m_rhs.open();
    m_isModalAssembly = true;
    assembleSystem( time_n, dt, domain, m_dofManager, m_localMatrix.toViewConstSizes(), localRhs );
    m_isModalAssembly = false;
    m_rhs.close();
  }
  m_matrix.create( m_localMatrix.toViewConst(), m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
  m_matrix.setDofManager( &m_dofManager );
  // GEOS assembles the Jacobian of the residual, which is minus the stiffness
  m_matrix.scale( -1.0 );

  // ---- Lumped mass and constrained degrees of freedom ----
  ParallelVector freeMask;
  ParallelVector massDiag;
  globalIndex const numConstrained = computeModalDiagonals( time_n + dt, domain, freeMask, massDiag );

  // Constrained rows and columns are removed symmetrically: K <- D K D, with D = diag( freeMask ).
  // They get a decoupled diagonal entry (the original stiffness diagonal) in the shifted matrix, and a
  // zero mass, so that the eigenvectors vanish there.
  ParallelVector constrainedDiagonal;
  if( numConstrained > 0 )
  {
    ParallelVector stiffnessDiagonal;
    stiffnessDiagonal.create( m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
    m_matrix.extractDiagonal( stiffnessDiagonal );
    m_matrix.leftRightScale( freeMask, freeMask );

    constrainedDiagonal.create( m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
    constrainedDiagonal.copy( freeMask );
    constrainedDiagonal.scale( -1.0 );
    arrayView1d< real64 > const values = constrainedDiagonal.open();
    arrayView1d< real64 const > const kDiagonal = stiffnessDiagonal.values();
    forAll< parallelDevicePolicy<> >( values.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      // (1 - mask) * K_ii
      values[i] = ( 1.0 + values[i] ) * kDiagonal[i];
    } );
    constrainedDiagonal.close();
  }

  // ---- Shifted operator K - sigma M ----
  real64 const shiftFrequency = m_modalShiftFrequency;
  real64 const shift = ( shiftFrequency < 0.0 ? -1.0 : 1.0 ) * ( 2.0 * M_PI * shiftFrequency ) * ( 2.0 * M_PI * shiftFrequency );

  ParallelMatrix shiftedMatrix( m_matrix );
  shiftedMatrix.setDofManager( &m_dofManager );
  shiftedMatrix.addDiagonal( massDiag, -shift );
  if( numConstrained > 0 )
  {
    shiftedMatrix.addDiagonal( constrainedDiagonal, 1.0 );
  }

  // ---- Linear solver for the shifted operator, set up once and applied for every Krylov vector ----
  LinearSolverParameters const & linParams = m_linearSolverParameters.get();
  bool const isDirectSolver = ( linParams.solverType == LinearSolverParameters::SolverType::direct );

  // Arnoldi needs the shift-and-invert operator, i.e. accurate solves with the shifted matrix.
  // LOBPCG only needs a preconditioner: a standalone one (e.g. one AMG cycle) unless the linear solver is direct.
  bool const needsInverse = ( m_modalSolverType == EigenSolverParameters::SolverType::arnoldi ) || isDirectSolver;

  std::unique_ptr< KrylovSolver< ParallelVector > > krylovSolver;
  std::unique_ptr< PreconditionerBase< LAInterface > > standalonePreconditioner;
  integer numLinearSolves = 0;
  integer numLinearIterations = 0;
  integer numLinearFailures = 0;
  real64 linearSolveTime = 0.0;
  {
    Stopwatch setupWatch;
    if( !needsInverse )
    {
      if( m_precond )
      {
        m_precond->setup( shiftedMatrix );
      }
      else
      {
        standalonePreconditioner = LAInterface::createPreconditioner( linParams, getLinearSolverNearNullKernel() );
        standalonePreconditioner->setup( shiftedMatrix );
      }
    }
    else if( isDirectSolver || !m_precond )
    {
      m_linearSolver = LAInterface::createSolver( linParams );

      LinearSolverExecutionContext executionContext;
      executionContext.solverName = getName();
      executionContext.cycleNumber = cycleNumber;
      executionContext.timeStepAttempt = m_nonlinearSolverParameters.m_numTimeStepAttempts;
      executionContext.configurationAttempt = m_nonlinearSolverParameters.m_numConfigurationAttempts;
      executionContext.nonlinearIteration = 0;
      executionContext.systemSetupTimestamp = getSystemSetupTimestamp();
      m_linearSolver->setExecutionContext( executionContext );
      m_linearSolver->setNearNullKernel( getLinearSolverNearNullKernel() );
      m_linearSolver->setup( shiftedMatrix );
    }
    else
    {
      m_precond->setup( shiftedMatrix );
      krylovSolver = KrylovSolver< ParallelVector >::create( linParams, shiftedMatrix, *m_precond );
    }
    GEOS_LOG_RANK_0_IF( getLogLevel() >= 1,
                        GEOS_FMT( "{}: modal analysis linear solver setup: {:.3f} s", getName(), setupWatch.elapsedTime() ) );
  }

  SolverInverseOperator shiftedInverse( shiftedMatrix, [&]( ParallelVector const & rhs, ParallelVector & sol )
  {
    Stopwatch solveWatch;
    LinearSolverResult result;
    if( krylovSolver )
    {
      krylovSolver->solve( rhs, sol );
      result = krylovSolver->result();
    }
    else
    {
      m_linearSolver->solve( rhs, sol );
      result = m_linearSolver->result();
    }
    linearSolveTime += solveWatch.elapsedTime();
    ++numLinearSolves;
    numLinearIterations += result.numIterations;
    if( !result.success() )
    {
      ++numLinearFailures;
    }
  } );

  // ---- Eigensolve ----
  EigenSolverParameters eigenParams;
  eigenParams.solverType = m_modalSolverType;
  eigenParams.numEigenvalues = m_modalNumModes;
  eigenParams.shift = shift;
  eigenParams.tolerance = m_modalTolerance;
  eigenParams.maxIterations = m_modalMaxIterations;
  eigenParams.subspaceSize = m_modalSubspaceSize;
  eigenParams.blockSize = m_modalBlockSize;
  eigenParams.completenessCheck = m_modalCompletenessCheck;
  eigenParams.seed = m_modalSeed;
  eigenParams.logLevel = getLogLevel() >= 2 ? 2 : ( getLogLevel() >= 1 ? 1 : 0 );

  DiagonalOperator< ParallelVector > massOperator( massDiag );
  // Preconditioned methods use one application of a set-up preconditioner (e.g. one multigrid cycle)
  LinearOperator< ParallelVector > const * preconditioner = nullptr;
  if( !needsInverse )
  {
    preconditioner = standalonePreconditioner ? standalonePreconditioner.get() : m_precond.get();
  }
  GeneralizedEigenProblem< ParallelVector > problem{ m_matrix, massOperator, needsInverse ? &shiftedInverse : nullptr, preconditioner,
                                                     m_solution.globalSize() - numConstrained, &freeMask };

  // Rigid-body modes of a free structure, deflated from the eigensolve
  array1d< ParallelVector > rigidBodyModes;
  if( m_modalDeflateRigidBodyModes != 0 )
  {
    GEOS_ERROR_IF( numConstrained > 0,
                   "Rigid-body modes cannot be deflated when displacement boundary conditions are applied",
                   getWrapperDataContext( viewKeyStruct::modalDeflateRigidBodyModesString() ) );
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  string_array const & )
    {
      if( rigidBodyModes.empty() )
      {
        NodeManager const & nodes = mesh.getNodeManager();
        arrayView1d< globalIndex const > const dofNumber =
          nodes.getReference< globalIndex_array >( m_dofManager.getKey( solidMechanics::totalDisplacement::key() ) );
        rigidBodyModes = LAIHelperFunctions::computeRigidBodyModes< ParallelVector >( nodes.referencePosition(),
                                                                                      dofNumber,
                                                                                      m_dofManager.rankOffset(),
                                                                                      m_dofManager.numLocalDofs() );
      }
    } );
    for( ParallelVector const & mode : rigidBodyModes )
    {
      problem.constraints.push_back( &mode );
    }
  }

  std::unique_ptr< GeneralizedEigenSolver< ParallelVector > > eigenSolver =
    GeneralizedEigenSolver< ParallelVector >::create( eigenParams );

  std::vector< ParallelVector > modes;
  EigenSolverResult const eigenResult = eigenSolver->solve( problem, m_solution, modes );

  GEOS_WARNING_IF( !eigenResult.converged,
                   GEOS_FMT( "Modal analysis: only {} of {} eigenpairs converged to the tolerance {:.1e} after {} restarts",
                             eigenResult.numConverged, m_modalNumModes, m_modalTolerance, eigenResult.numIterations ),
                   getDataContext() );
  // The convergence test of the eigensolver trusts the linear solves. If none of them converged, the modes are
  // meaningless, whatever the Ritz estimates say. This happens for example with cg when the shift is above the
  // lowest eigenvalue: K - sigma M is then indefinite.
  GEOS_ERROR_IF( numLinearSolves > 0 && numLinearFailures == numLinearSolves,
                 GEOS_FMT( "Modal analysis: none of the {} linear solves converged, so the modes are not valid. "
                           "With the cg solver, use a shift frequency below the lowest mode (a negative shift frequency) "
                           "so that K - sigma M is positive definite, or use gmres or a direct solver. "
                           "Otherwise tighten krylovTol or strengthen the preconditioner.",
                           numLinearSolves ),
                 getDataContext() );
  GEOS_WARNING_IF( numLinearFailures > 0,
                   GEOS_FMT( "Modal analysis: {} of {} linear solves did not converge, so the eigenpairs may be inaccurate; "
                             "tighten the linear solver tolerance (krylovTol) or strengthen the preconditioner",
                             numLinearFailures, numLinearSolves ),
                   getDataContext() );

  // ---- Store results ----
  integer const numModes = LvArray::integerConversion< integer >( modes.size() );
  m_modalEigenvalues.resize( numModes );
  m_modalFrequencies.resize( numModes );
  m_modalResiduals.resize( numModes );
  m_modalParticipationFactors.resize( numModes, 3 );
  m_modalParticipationFactors.zero();

  stdVector< string > shapeFieldNames;
  for( integer k = 0; k < numModes; ++k )
  {
    real64 const lambda = eigenResult.eigenvalues[k];
    m_modalEigenvalues[k] = lambda;
    m_modalFrequencies[k] = ( lambda < 0.0 ? -1.0 : 1.0 ) * std::sqrt( std::fabs( lambda ) ) / ( 2.0 * M_PI );
    m_modalResiduals[k] = eigenResult.residuals[k];

    // Participation factors: Gamma_d = sum_i m_i phi_i e_d
    arrayView1d< real64 const > const phi = modes[k].values();
    arrayView1d< real64 const > const mass = massDiag.values();
    RAJA::ReduceSum< parallelDeviceReduce, real64 > sumX( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > sumY( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > sumZ( 0.0 );
    forAll< parallelDevicePolicy<> >( phi.size() / 3, [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      sumX += mass[3 * n] * phi[3 * n];
      sumY += mass[3 * n + 1] * phi[3 * n + 1];
      sumZ += mass[3 * n + 2] * phi[3 * n + 2];
    } );
    m_modalParticipationFactors( k, 0 ) = MpiWrapper::sum( sumX.get(), MPI_COMM_GEOS );
    m_modalParticipationFactors( k, 1 ) = MpiWrapper::sum( sumY.get(), MPI_COMM_GEOS );
    m_modalParticipationFactors( k, 2 ) = MpiWrapper::sum( sumZ.get(), MPI_COMM_GEOS );

    // Mode shape to the nodal field
    string const fieldName = modeShapeFieldName( k + 1 );
    m_dofManager.copyVectorToField( modes[k].values(),
                                    solidMechanics::totalDisplacement::key(),
                                    fieldName,
                                    1.0 );
    shapeFieldNames.emplace_back( fieldName );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, shapeFieldNames );
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );

  // ---- Report ----
  if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
  {
    GEOS_LOG( GEOS_FMT( "\n{}: modal analysis, {} eigenpairs (shift {:.6e} Hz, sigma = {:.6e}), {} solver",
                        getName(), numModes, shiftFrequency, shift,
                        EnumStrings< EigenSolverParameters::SolverType >::toString( m_modalSolverType ) ) );
    if( numConstrained > 0 )
    {
      GEOS_LOG( GEOS_FMT( "  {} constrained degrees of freedom removed", numConstrained ) );
    }
    GEOS_LOG( "  mode   frequency [Hz]      eigenvalue [1/s^2]    residual     MPF-x         MPF-y         MPF-z" );
    for( integer k = 0; k < numModes; ++k )
    {
      GEOS_LOG( GEOS_FMT( "  {:4}   {:16.8e}   {:16.8e}   {:9.2e}   {:12.5e}  {:12.5e}  {:12.5e}",
                          k + 1, m_modalFrequencies[k], m_modalEigenvalues[k], m_modalResiduals[k],
                          m_modalParticipationFactors( k, 0 ),
                          m_modalParticipationFactors( k, 1 ),
                          m_modalParticipationFactors( k, 2 ) ) );
    }
    GEOS_LOG( GEOS_FMT( "  {} restarts/iterations, {} operator applications, {} linear solves ({} linear iterations, "
                        "{:.3f} s in linear solves), {:.3f} s eigensolve, {:.3f} s total",
                        eigenResult.numIterations, eigenResult.numOperatorApplications, numLinearSolves,
                        numLinearIterations, linearSolveTime, eigenResult.solveTime, totalWatch.elapsedTime() ) );
  }

  // The preconditioner and the linear solver refer to the shifted matrix, which is local to this function, while
  // the solver members m_precond and m_linearSolver outlive it. Release the matrix before it goes out of scope:
  // some backends (e.g. Trilinos/ML) crash when the matrix is deleted before the preconditioner.
  if( standalonePreconditioner )
  {
    standalonePreconditioner->clear();
  }
  if( m_precond )
  {
    m_precond->clear();
  }
  if( m_linearSolver )
  {
    m_linearSolver->clear();
  }

  return dt;
}

} // namespace geos
