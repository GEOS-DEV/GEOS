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
 * @file TrilinosSolver.cpp
 */

#include "TrilinosSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Epetra_Import.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <ml_MultiLevelPreconditioner.h>

#include <memory>

namespace geosx
{

namespace
{

void createTrilinosKrylovSolver( LinearSolverParameters const & params, AztecOO & solver )
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

TrilinosSolver::TrilinosSolver( LinearSolverParameters parameters )
  : Base( std::move( parameters ) ),
  m_precond( m_params ),
  m_solver( std::make_unique< AztecOO >() )
{
  createTrilinosKrylovSolver( m_params, *m_solver );
}

TrilinosSolver::~TrilinosSolver() = default;

void TrilinosSolver::setup( EpetraMatrix const & mat )
{
  clear();
  Base::setup( mat );
  Stopwatch timer( m_result.setupTime );
  m_precond.setup( mat );

  // HACK: Epetra is not const-correct, so we need the cast. The matrix is not actually modified.
  GEOSX_LAI_CHECK_ERROR( m_solver->SetUserMatrix( &const_cast< Epetra_FECrsMatrix & >( mat.unwrapped() ) ) );
  GEOSX_LAI_CHECK_ERROR( m_solver->SetPrecOperator( &m_precond.unwrapped() ) );
}

int TrilinosSolver::doSolve( EpetraVector const & rhs,
                             EpetraVector & sol ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );

  // HACK: Epetra is not const-correct, so we need the cast. The vector is not actually modified.
  GEOSX_LAI_CHECK_ERROR( m_solver->SetRHS( &const_cast< Epetra_FEVector & >( rhs.unwrapped() ) ) );
  GEOSX_LAI_CHECK_ERROR( m_solver->SetLHS( &sol.unwrapped() ) );
  return m_solver->Iterate( m_params.krylov.maxIterations, m_params.krylov.relTolerance );
}

void TrilinosSolver::apply( EpetraVector const & rhs,
                            EpetraVector & sol ) const
{
  doSolve( rhs, sol );
}

void TrilinosSolver::solve( EpetraVector const & rhs,
                            EpetraVector & sol ) const
{
  if( isZero( rhs.norm2(), 0.0 ) )
  {
    sol.zero();
    m_result.numIterations = 0;
    m_result.residualReduction = 0.0;
    m_result.solveTime = 0.0;
    m_result.status = LinearSolverResult::Status::Success;
    return;
  }

  {
    Stopwatch watch( m_result.solveTime );
    int const result = doSolve( rhs, sol );
    m_result.status = result ? LinearSolverResult::Status::NotConverged : LinearSolverResult::Status::Success;
  }

  m_result.numIterations = m_solver->NumIters();
  m_result.residualReduction = m_solver->ScaledResidual();

  if( m_params.logLevel >= 1 )
  {
    GEOSX_LOG_RANK_0( "\t\tLinear Solver | " << m_result.status <<
                      " | Iterations: " << m_result.numIterations <<
                      " | Final Rel Res: " << m_result.residualReduction <<
                      " | Setup Time: " << m_result.setupTime << " s" <<
                      " | Solve Time: " << m_result.solveTime << " s" );
  }
}

} // end geosx namespace
