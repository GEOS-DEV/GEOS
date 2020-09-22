/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TrilinosTpetraSolver.cpp
 */

#include "TrilinosTpetraSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/trilinos/TpetraMatrix.hpp"
#include "linearAlgebra/interfaces/trilinos/TpetraVector.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosTpetraPreconditioner.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <Amesos2.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

namespace geosx
{

TrilinosTpetraSolver::TrilinosTpetraSolver( LinearSolverParameters parameters ):
  m_parameters( std::move( parameters ) )
{}

TrilinosTpetraSolver::~TrilinosTpetraSolver() = default;

void TrilinosTpetraSolver::solve( TpetraMatrix & mat,
                                  TpetraVector & sol,
                                  TpetraVector & rhs,
                                  DofManager const * const GEOSX_UNUSED_PARAM( dofManager ) )
{
  GEOSX_LAI_ASSERT( mat.ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );

  if( m_parameters.solverType == LinearSolverParameters::SolverType::direct )
  {
    solve_direct( mat, sol, rhs );
  }
  else
  {
    solve_krylov( mat, sol, rhs );
  }
}

namespace
{

string getAmesos2SolverName( bool const parallel )
{
  // Define a list of direct solvers to try ordered by priority
  static std::array< string, 1 > parallelSolvers =
  {
    "superlu_dist",
  };

  static std::array< string, 2 > serialSolvers =
  {
    "umfpack",
    "Klu2"
  };

  // First try parallel solvers, if requested (which is normally the default)
  if( parallel )
  {
    for( string const & solver : parallelSolvers )
    {
      if( Amesos2::query( solver ) )
      {
        return solver;
      }
    }
  }

  // Fall back to serial if parallel not available or user doesn't want it
  for( string const & solver : serialSolvers )
  {
    if( Amesos2::query( solver ) )
    {
      return solver;
    }
  }
  GEOSX_ERROR( "No direct solver available in Trilinos/Amesos2" );
  return "";
}

string getAmesos2SuperLUDistColPerm( LinearSolverParameters::Direct::ColPerm const value )
{
  static std::map< LinearSolverParameters::Direct::ColPerm, string > const optionMap =
  {
    { LinearSolverParameters::Direct::ColPerm::none, "NATURAL" },
    { LinearSolverParameters::Direct::ColPerm::parmetis, "PARMETIS" },
  };

  // Since not all valid options can be accepted by Amesos2, return default value instead of error.
  string result;
  if( optionMap.count( value ) > 0 )
  {
    result = optionMap.at( value );
  }
  else
  {
    LinearSolverParameters::Direct::ColPerm const defaultValue = LinearSolverParameters::Direct::ColPerm::none;
    GEOSX_WARNING( "Unsupported Trilinos/Amesos2/SuperLU_Dist column permutation option: " << value << ". "
                                                                                                       "Default value '" << defaultValue << "' will be used." );
    result = optionMap.at( defaultValue );
  }
  return result;
}

} // namespace

void TrilinosTpetraSolver::solve_direct( TpetraMatrix & mat,
                                         TpetraVector & sol,
                                         TpetraVector & rhs )
{
  // Time setup and solve
  Stopwatch watch;

  // Create Epetra linear problem and instantiate solver.
  using Tpetra_Matrix = TpetraMatrix::Tpetra_CrsMatrix;
  using Tpetra_MultiVector = TpetraVector::Tpetra_MultiVector;
  using Amesos2_Solver = Amesos2::Solver< Tpetra_Matrix, Tpetra_MultiVector >;

  // Instantiate the solver
  // Note: we explicitly provide template arguments to create() because
  // Tpetra::Vector is not actually recognized by some Amesos2 internals
  string const solverName = getAmesos2SolverName( m_parameters.direct.parallel );
  Teuchos::RCP< Amesos2_Solver > solver =
    Amesos2::create< Tpetra_Matrix, Tpetra_MultiVector >( solverName.c_str(),
                                                          Teuchos::rcp( &mat.unwrapped(), false ),
                                                          Teuchos::rcp( &sol.unwrapped(), false ),
                                                          Teuchos::rcp( &rhs.unwrapped(), false ) );

  // Set additional solver options
  Teuchos::ParameterList params( "Amesos2" );
  if( solverName == "superlu_dist" )
  {
    // Only these two options are supported by Amesos2_Superludist
    params.set( "ColPerm", getAmesos2SuperLUDistColPerm( m_parameters.direct.colPerm ) );
    params.set( "ReplaceTinyPivot", bool( m_parameters.direct.replaceTinyPivot ) );
  }
  solver->setParameters( Teuchos::rcp( &params, false ) );

  // Factorize the matrix
  solver->symbolicFactorization();
  solver->numericFactorization();
  m_result.setupTime = watch.elapsedTime();

  // Solve the system
  watch.zero();
  solver->solve();
  m_result.solveTime = watch.elapsedTime();

  // Basic output
  if( m_parameters.logLevel > 0 )
  {
    solver->printTiming( *Teuchos::getFancyOStream( Teuchos::rcp( &std::cout, false ) ) );
  }

  m_result.status = LinearSolverResult::Status::Success;
  m_result.numIterations = 1;
  m_result.residualReduction = NumericTraits< real64 >::eps;
}

namespace
{

string getBelosSolverName( LinearSolverParameters::SolverType const & type )
{
  static std::map< LinearSolverParameters::SolverType, string > typeMap =
  {
    { LinearSolverParameters::SolverType::cg, "CG" },
    { LinearSolverParameters::SolverType::bicgstab, "BICGSTAB" },
    { LinearSolverParameters::SolverType::gmres, "GMRES" },
    { LinearSolverParameters::SolverType::fgmres, "FLEXIBLE GMRES" }
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Trilinos/Belos solver option: " << type );
  return typeMap.at( type );
}

Belos::MsgType getBelosVerbosity( integer const value )
{
  static std::array< Belos::MsgType, 5 > const optionMap =
  {
    Belos::MsgType( Belos::Errors ),
    Belos::MsgType( Belos::Errors | Belos::Warnings | Belos::FinalSummary ),
    Belos::MsgType( Belos::Errors | Belos::Warnings | Belos::FinalSummary | Belos::IterationDetails ),
    Belos::MsgType( Belos::Errors | Belos::Warnings | Belos::FinalSummary | Belos::IterationDetails | Belos::TimingDetails ),
    Belos::MsgType( Belos::Errors | Belos::Warnings | Belos::FinalSummary | Belos::IterationDetails | Belos::TimingDetails | Belos::Debug )
  };

  return optionMap[ std::max( std::min( value, integer( optionMap.size() - 1 ) ), 0 ) ];
}

integer getBelosOutputFrequency( integer const value )
{
  static std::array< integer, 5 > const optionMap =
  {
    -1,
    1,
    1,
    1,
    1
  };

  return optionMap[ std::max( std::min( value, integer( optionMap.size() - 1 ) ), 0 ) ];
}

} // namespace

void TrilinosTpetraSolver::solve_krylov( TpetraMatrix & mat,
                                         TpetraVector & sol,
                                         TpetraVector & rhs )
{
  // Time setup and solve
  Stopwatch watch;

  using Tpetra_MultiVector = TpetraVector::Tpetra_MultiVector;
  using Tpetra_Operator = TrilinosTpetraPreconditioner::Tpetra_Operator;
  using Belos_LinearProblem = Belos::LinearProblem< real64, Tpetra_MultiVector, Tpetra_Operator >;
  using Belos_Factory = Belos::SolverFactory< real64, Tpetra_MultiVector, Tpetra_Operator >;
  using Belos_Solver = Belos::SolverManager< real64, Tpetra_MultiVector, Tpetra_Operator >;

  // Deal with separate component approximation
  TpetraMatrix separateComponentMatrix;
  if( m_parameters.amg.separateComponents )
  {
    LAIHelperFunctions::SeparateComponentFilter( mat, separateComponentMatrix, m_parameters.dofsPerNode );
  }
  TpetraMatrix & precondMat = m_parameters.amg.separateComponents ? separateComponentMatrix : mat;

  // Create and compute preconditioner
  TrilinosTpetraPreconditioner precond( m_parameters );
  precond.compute( precondMat );

  // Create a linear problem.
  Teuchos::RCP< Belos_LinearProblem > problem =
    Teuchos::rcp( new Belos_LinearProblem( Teuchos::rcp( &mat.unwrapped(), false ),
                                           Teuchos::rcp( &sol.unwrapped(), false ),
                                           Teuchos::rcp( &rhs.unwrapped(), false ) ) );

  // Set as right preconditioner.
  // TODO: do we need an option for left/right?
  problem->setRightPrec( Teuchos::rcp( &precond.unwrapped(), false ) );
  problem->setProblem();

  Teuchos::ParameterList list;
  list.set( "Implicit Residual Scaling", "Norm of RHS" );
  list.set( "Explicit Residual Scaling", "Norm of RHS" );
  list.set( "Maximum Iterations", m_parameters.krylov.maxIterations );
  list.set( "Num Blocks", m_parameters.krylov.maxRestart );
  list.set( "Maximum Restarts", m_parameters.krylov.maxIterations / m_parameters.krylov.maxRestart - 1 );
  list.set( "Orthogonalization", "DGKS" ); // this one seems more robust
  list.set( "Convergence Tolerance", m_parameters.krylov.relTolerance );
  list.set( "Verbosity", getBelosVerbosity( m_parameters.logLevel ) );
  list.set( "Output Frequency", getBelosOutputFrequency( m_parameters.logLevel ) );

  Belos_Factory factory;
  Teuchos::RCP< Belos_Solver > solver = factory.create( getBelosSolverName( m_parameters.solverType ),
                                                        Teuchos::rcp( &list, false ) );
  solver->setProblem( problem );

  m_result.setupTime = watch.elapsedTime();

  // Do solve
  watch.zero();
  Belos::ReturnType const result = solver->solve();
  m_result.solveTime = watch.elapsedTime();

  m_result.status = ( result == Belos::Converged ) ? LinearSolverResult::Status::Success : LinearSolverResult::Status::NotConverged;
  m_result.numIterations = solver->getNumIters();
  m_result.residualReduction = solver->achievedTol();
}

} // end geosx namespace
