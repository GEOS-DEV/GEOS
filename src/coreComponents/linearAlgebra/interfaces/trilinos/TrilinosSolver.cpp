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

  if( rhs.norm2() > 0.0 )
  {
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
      solveDirect( mat, sol, rhs );
    }
    else
    {
      solveKrylov( mat, sol, rhs );
    }
  }
  else
  {
    sol.zero();
    m_result.status = LinearSolverResult::Status::Success;
    m_result.setupTime = 0.0;
    m_result.solveTime = 0.0;
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
  GEOSX_LAI_CHECK_ERROR( SLUDData.solve( rhs.extractLocalVector(), sol.extractLocalVector() ) );

  // Save setup and solution times
  result.setupTime = SLUDData.setupTime();
  result.solveTime = SLUDData.solveTime();

  EpetraVector res( rhs );
  mat.gemv( -1.0, sol, 1.0, res );
  result.residualReduction = res.norm2() / rhs.norm2();

  result.status = parameters.direct.checkResidual == 0 ? LinearSolverResult::Status::Success : LinearSolverResult::Status::Breakdown;
  result.numIterations = 1;
  if( parameters.direct.checkResidual )
  {
    if( result.residualReduction < SLUDData.relativeTolerance() )
    {
      result.status = LinearSolverResult::Status::Success;
    }
    else
    {
      real64 const cond = EpetraSuperLU_DistCond( mat, SLUDData );
      if( parameters.logLevel > 0 )
      {
        GEOSX_LOG_RANK_0( "Using a more accurate estimate of the condition number" );
        GEOSX_LOG_RANK_0( "Condition number is " << cond );
      }
      if( result.residualReduction < SLUDData.precisionTolerance() * cond )
      {
        result.status = LinearSolverResult::Status::Success;
      }
    }
  }
}

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

  GEOSX_LAI_CHECK_ERROR( SSData.setup() );
  GEOSX_LAI_CHECK_ERROR( SuiteSparseSolve( SSData, serialMap, importToSerial, rhs, sol ) );

  // Save setup and solution times
  result.setupTime = SSData.setupTime();
  result.solveTime = SSData.solveTime();

  EpetraVector res( rhs );
  mat.gemv( -1.0, sol, 1.0, res );
  result.residualReduction = res.norm2() / rhs.norm2();

  result.status = parameters.direct.checkResidual == 0 ? LinearSolverResult::Status::Success : LinearSolverResult::Status::Breakdown;
  result.numIterations = 1;
  if( parameters.direct.checkResidual )
  {
    if( result.residualReduction < SSData.relativeTolerance() )
    {
      result.status = LinearSolverResult::Status::Success;
    }
    else
    {
      real64 const cond = EpetraSuiteSparseCond( mat, serialMap, importToSerial, SSData );
      if( parameters.logLevel > 0 )
      {
        GEOSX_LOG_RANK_0( "Using a more accurate estimate of the condition number" );
        GEOSX_LOG_RANK_0( "Condition number is " << cond );
      }
      if( result.residualReduction < SSData.precisionTolerance() * cond )
      {
        result.status = LinearSolverResult::Status::Success;
      }
    }
  }

  delete serialMap;
  delete importToSerial;
}

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

void TrilinosSolver::solveDirect( EpetraMatrix & mat,
                                  EpetraVector & sol,
                                  EpetraVector & rhs )
{
  if( m_parameters.direct.parallel )
  {
    solve_parallelDirect( m_parameters, mat, sol, rhs, m_result );
  }
  else
  {
    solve_serialDirect( m_parameters, mat, sol, rhs, m_result );
  }
}

void TrilinosSolver::solveKrylov( EpetraMatrix & mat,
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
