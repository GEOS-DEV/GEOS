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
 * @file TrilinosSolver.cpp
 */

#include "TrilinosSolver.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraMatrix.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraVector.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosPreconditioner.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraSuperLU_Dist.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraSuiteSparse.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <Epetra_Map.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Epetra_Import.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <ml_MultiLevelPreconditioner.h>

#include <memory>

namespace geosx
{

TrilinosSolver::TrilinosSolver( LinearSolverParameters parameters ):
  m_parameters( std::move( parameters ) )
{}

TrilinosSolver::~TrilinosSolver() = default;

void TrilinosSolver::solve( EpetraMatrix & mat,
                            EpetraVector & sol,
                            EpetraVector & rhs,
                            DofManager const * const dofManager )
{
  GEOSX_LAI_ASSERT( mat.ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );

  GEOSX_UNUSED_VAR( dofManager );

  if( m_parameters.scaling.useRowScaling )
  {
    Epetra_FECrsMatrix & mat_raw = mat.unwrapped();
    Epetra_MultiVector & rhs_raw = rhs.unwrapped();

    Epetra_Vector scaling( mat_raw.RowMap() );
    mat_raw.InvRowSums( scaling );
    mat_raw.LeftScale( scaling );

    Epetra_MultiVector tmp( rhs_raw );
    rhs_raw.Multiply( 1.0, scaling, tmp, 0.0 );
  }

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

void solve_parallelDirect( LinearSolverParameters const & parameters,
                           EpetraMatrix & mat,
                           EpetraVector & sol,
                           EpetraVector & rhs,
                           LinearSolverResult & result )
{
  // To be able to use SuperLU_Dist solver we need to disable floating point exceptions
  LvArray::system::FloatingPointExceptionGuard guard;

  SuperLU_Dist SLUDData( parameters );
  EpetraConvertToSuperMatrix( mat, SLUDData );

  GEOSX_LAI_CHECK_ERROR( SLUDData.setup() );
  sol.copy( rhs );
  GEOSX_LAI_CHECK_ERROR( SLUDData.solve( sol.localSize(), sol.extractLocalVector() ) );

  // Save setup and solution times
  result.setupTime = SLUDData.setupTime();
  result.solveTime = SLUDData.solveTime();

  EpetraVector res( rhs );
  mat.gemv( -1.0, sol, 1.0, res );
  result.residualReduction = res.norm2() / rhs.norm2();

  result.status = LinearSolverResult::Status::Success;
  if( result.residualReduction < SLUDData.relativeTolerance() )
  {
    result.status = LinearSolverResult::Status::Success;
    result.numIterations = 1;
  }
  else
  {
    result.status = LinearSolverResult::Status::Breakdown;
  }
}

#ifdef GEOSX_USE_SUITESPARSE
void solve_serialDirect( LinearSolverParameters const & parameters,
                         EpetraMatrix & mat,
                         EpetraVector & sol,
                         EpetraVector & rhs,
                         LinearSolverResult & result )
{
  // To be able to use UMFPACK direct solver we need to disable floating point exceptions
  LvArray::system::FloatingPointExceptionGuard guard;

  SuiteSparse SSData( parameters );

  Epetra_Map * serialMap = NULL;
  Epetra_Import * importToSerial = NULL;
  ConvertEpetraToSuiteSparseMatrix( mat, SSData, serialMap, importToSerial );

  int info = 0;
  info = SSData.setup();

  info += SuiteSparseSolve( SSData, serialMap, importToSerial, rhs, sol );

  delete serialMap;
  delete importToSerial;

  // Save setup and solution times
  result.setupTime = SSData.setupTime();
  result.solveTime = SSData.solveTime();

  if( info == 0 )
  {
    EpetraVector res( rhs );
    mat.gemv( -1.0, sol, 1.0, res );
    result.residualReduction = res.norm2() / rhs.norm2();
  }

  if( info == 0 && result.residualReduction < SSData.relativeTolerance() )
  {
    result.status = LinearSolverResult::Status::Success;
    result.numIterations = 1;
  }
  else
  {
    result.status = LinearSolverResult::Status::Breakdown;
  }
}
#endif

void CreateTrilinosKrylovSolver( LinearSolverParameters const & params, AztecOO & solver )
{
  switch( params.solverType )
  {
    case LinearSolverParameters::SolverType::gmres:
    {
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_solver, AZ_gmres ) );
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_kspace, params.krylov.maxRestart ) );
      break;
    }
    case LinearSolverParameters::SolverType::bicgstab:
    {
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_solver, AZ_bicgstab ) );
      break;
    }
    case LinearSolverParameters::SolverType::cg:
    {
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_solver, AZ_cg ) );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Solver type not supported in Trilinos interface: " << params.solverType );
    }
  }

  // Ask for a convergence normalized by the right hand side
  GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_conv, AZ_rhs ) );

  // Control output
  switch( params.logLevel )
  {
    case 1:
    {
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_output, AZ_summary ) );
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_diagnostics, AZ_all ) );
      break;
    }
    case 2:
    {
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_output, AZ_all ) );
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_diagnostics, AZ_all ) );
      break;
    }
    default:
    {
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_output, AZ_none ) );
    }
  }
}

} // namespace

void TrilinosSolver::solve_direct( EpetraMatrix & mat,
                                   EpetraVector & sol,
                                   EpetraVector & rhs )
{
  if( m_parameters.direct.parallel )
  {
    solve_parallelDirect( m_parameters, mat, sol, rhs, m_result );
  }
  else
  {
#ifdef GEOSX_USE_SUITESPARSE
    solve_serialDirect( m_parameters, mat, sol, rhs, m_result );
#else
    GEOSX_ERROR( "Trilinos direct solver interface: serial direct solver not available (try to compile GEOSX TPLs with SuiteSparse)." );
#endif
  }
}

void TrilinosSolver::solve_krylov( EpetraMatrix & mat,
                                   EpetraVector & sol,
                                   EpetraVector & rhs )
{
  // Time setup and solve
  Stopwatch watch;

  // Create Epetra linear problem.
  Epetra_LinearProblem problem( &mat.unwrapped(),
                                &sol.unwrapped(),
                                &rhs.unwrapped() );

  // Instantiate the AztecOO solver.
  AztecOO solver( problem );

  // Choose the solver type
  CreateTrilinosKrylovSolver( m_parameters, solver );

  // Deal with separate component approximation
  EpetraMatrix separateComponentMatrix;
  if( m_parameters.amg.separateComponents )
  {
    LAIHelperFunctions::SeparateComponentFilter( mat, separateComponentMatrix, m_parameters.dofsPerNode );
  }
  EpetraMatrix & precondMat = m_parameters.amg.separateComponents ? separateComponentMatrix : mat;

  // Create and compute preconditioner
  TrilinosPreconditioner precond( m_parameters );
  precond.compute( precondMat );
  GEOSX_LAI_CHECK_ERROR( solver.SetPrecOperator( &precond.unwrapped() ) );

  m_result.setupTime = watch.elapsedTime();

  // Actually solve
  watch.zero();
  int const result = solver.Iterate( m_parameters.krylov.maxIterations,
                                     m_parameters.krylov.relTolerance );
  m_result.solveTime = watch.elapsedTime();

  m_result.status = result ? LinearSolverResult::Status::NotConverged : LinearSolverResult::Status::Success;
  m_result.numIterations = solver.NumIters();
  m_result.residualReduction = solver.ScaledResidual();
}

} // end geosx namespace
