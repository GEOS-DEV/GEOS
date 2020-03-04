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
 * @file TrilinosSolver.cpp
 */

// BEGIN_RST_NARRATIVE TrilinosSolver.rst
// ==============================
// Trilinos Solver
// ==============================
// This class implements solvers from the Trilinos library. Iterative solvers come from
// the AztecOO package, and direct solvers from the Amesos package.

// Source inclues
#include "TrilinosSolver.hpp"
#include "linearAlgebra/interfaces/EpetraMatrix.hpp"
#include "linearAlgebra/interfaces/EpetraVector.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

// TPL includes
#include <Epetra_Map.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <ml_MultiLevelPreconditioner.h>

// Put everything under the geosx namespace.
namespace geosx
{

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

TrilinosSolver::TrilinosSolver( LinearSolverParameters const & parameters ) :
  m_parameters( parameters )
{}

TrilinosSolver::~TrilinosSolver()
{
  delete m_solver;
  m_solver = nullptr;
}

// ----------------------------
// Top-Level Solver
// ----------------------------
// We switch between different solverTypes here

void TrilinosSolver::solve( EpetraMatrix & mat,
                            EpetraVector & sol,
                            EpetraVector & rhs )
{
  if( m_parameters.scaling.useRowScaling )
  {
    Epetra_FECrsMatrix * mat_ptr = mat.unwrappedPointer();
    Epetra_MultiVector * rhs_ptr = rhs.unwrappedPointer();

    Epetra_Vector scaling( mat_ptr->RowMap() );
    mat_ptr->InvRowSums( scaling );
    mat_ptr->LeftScale( scaling );

    Epetra_MultiVector tmp( *rhs_ptr );
    rhs_ptr->Multiply( 1.0, scaling, tmp, 0.0 );
  }

  if( m_parameters.solverType == "direct" )
  {
    solve_direct( mat, sol, rhs );
  }
  else
  {
    solve_krylov( mat, sol, rhs );
  }
}

// ----------------------------
// Direct solver
// ----------------------------

void TrilinosSolver::solve_direct( EpetraMatrix & mat,
                                   EpetraVector & sol,
                                   EpetraVector & rhs )
{
  // Create Epetra linear problem and instantiate solver.
  Epetra_LinearProblem problem( mat.unwrappedPointer(),
                                sol.unwrappedPointer(),
                                rhs.unwrappedPointer() );

  if ( m_solver == nullptr )
  {
    // Instantiate the Amesos solver.
    Amesos Factory;

    // Select KLU solver (only one available as of 9/20/2018)
    m_solver = Factory.Create( "Klu", problem );
  }

  // Factorize the matrix
  m_solver->SymbolicFactorization();
  m_solver->NumericFactorization();

  // Solve the system
  m_solver->Solve();

  // Basic output
  if( m_parameters.logLevel > 0 )
  {
    m_solver->PrintStatus();
    m_solver->PrintTiming();
  }
}


// ----------------------------
// Iterative solver
// ----------------------------

void TrilinosSolver::solve_krylov( EpetraMatrix & mat,
                                   EpetraVector & sol,
                                   EpetraVector & rhs )
{
  // Create Epetra linear problem.
  Epetra_LinearProblem problem( mat.unwrappedPointer(),
                                sol.unwrappedPointer(),
                                rhs.unwrappedPointer() );

  // Instantiate the AztecOO solver.
  AztecOO solver( problem );

  // Extra scratch matrix, may or may not be used
  std::unique_ptr<EpetraMatrix> scratch;

  // Choose the solver type
  if( m_parameters.solverType == "gmres" )
  {
    solver.SetAztecOption( AZ_solver, AZ_gmres );
    solver.SetAztecOption( AZ_kspace, m_parameters.krylov.maxRestart );
  }
  else if( m_parameters.solverType == "bicgstab" )
  {
    solver.SetAztecOption( AZ_solver, AZ_bicgstab );
  }
  else if( m_parameters.solverType == "cg" )
  {
    solver.SetAztecOption( AZ_solver, AZ_cg );
  }
  else
    GEOSX_ERROR( "The requested linear solverType doesn't seem to exist" );

  // Create a null pointer to an ML amg preconditioner
  std::unique_ptr<ML_Epetra::MultiLevelPreconditioner> ml_preconditioner;

  // Choose the preconditioner type
  if( m_parameters.preconditionerType == "none" )
  {
    solver.SetAztecOption( AZ_precond, AZ_none );
  }
  else if( m_parameters.preconditionerType == "jacobi" )
  {
    solver.SetAztecOption( AZ_precond, AZ_Jacobi );
  }
  else if( m_parameters.preconditionerType == "ilu" )
  {
    solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
    solver.SetAztecOption( AZ_overlap, m_parameters.dd.overlap );
    solver.SetAztecOption( AZ_subdomain_solve, AZ_ilu );
    solver.SetAztecOption( AZ_graph_fill, m_parameters.ilu.fill );
  }
  else if( m_parameters.preconditionerType == "icc" )
  {
    solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
    solver.SetAztecOption( AZ_overlap, m_parameters.dd.overlap );
    solver.SetAztecOption( AZ_subdomain_solve, AZ_icc );
    solver.SetAztecOption( AZ_graph_fill, m_parameters.ilu.fill );
  }
  else if( m_parameters.preconditionerType == "ilut" )
  {
    solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
    solver.SetAztecOption( AZ_overlap, m_parameters.dd.overlap );
    solver.SetAztecOption( AZ_subdomain_solve, AZ_ilut );
    solver.SetAztecParam( AZ_ilut_fill, (m_parameters.ilu.fill>0 ? real64( m_parameters.ilu.fill ) : 1.0));
  }
  else if( m_parameters.preconditionerType == "amg" )
  {
    Teuchos::ParameterList list;

    if( m_parameters.amg.isSymmetric )
      ML_Epetra::SetDefaults( "SA", list );
    else
      ML_Epetra::SetDefaults( "NSSA", list );

    std::map<string, string> translate; // maps GEOSX to ML syntax

    translate.insert( std::make_pair( "V", "MGV" ));
    translate.insert( std::make_pair( "W", "MGW" ));
    translate.insert( std::make_pair( "direct", "Amesos-KLU" ));
    translate.insert( std::make_pair( "jacobi", "Jacobi" ));
    translate.insert( std::make_pair( "blockJacobi", "block Jacobi" ));
    translate.insert( std::make_pair( "gaussSeidel", "Gauss-Seidel" ));
    translate.insert( std::make_pair( "blockGaussSeidel", "block Gauss-Seidel" ));
    translate.insert( std::make_pair( "chebyshev", "Chebyshev" ));
    translate.insert( std::make_pair( "ilu", "ILU" ));
    translate.insert( std::make_pair( "ilut", "ILUT" ));
    translate.insert( std::make_pair( "icc", "IC" ));

    list.set( "ML output", m_parameters.logLevel );
    list.set( "max levels", m_parameters.amg.maxLevels );
    list.set( "aggregation: type", "Uncoupled" );
    list.set( "PDE equations", m_parameters.dofsPerNode );
    list.set( "smoother: sweeps", m_parameters.amg.numSweeps );
    list.set( "prec type", translate[m_parameters.amg.cycleType] );
    list.set( "smoother: type", translate[m_parameters.amg.smootherType] );
    list.set( "coarse: type", translate[m_parameters.amg.coarseType] );
    //list.set( "aggregation: threshold", 0.0 );
    //list.set( "smoother: pre or post", "post" );

    //TODO: add user-defined null space / rigid body mode support
    //list.set("null space: type","pre-computed");
    //list.set("null space: vectors",&rigid_body_modes[0]);
    //list.set("null space: dimension", n_rbm);

    //TODO: templatization for LAIHelperFunctions needed
    if(m_parameters.amg.separateComponents) // apply separate displacement component filter
    {
      scratch.reset(new EpetraMatrix());
      LAIHelperFunctions::SeparateComponentFilter<TrilinosInterface>(mat,*scratch,m_parameters.dofsPerNode);
      ml_preconditioner.reset( new ML_Epetra::MultiLevelPreconditioner( *scratch->unwrappedPointer(), list ));
    }
    else // just use original matrix to construct amg operator
    {
      ml_preconditioner.reset( new ML_Epetra::MultiLevelPreconditioner( *mat.unwrappedPointer(), list ));
    }

    solver.SetPrecOperator( ml_preconditioner.get() );
  }
  else
    GEOSX_ERROR( "The requested preconditionerType doesn't seem to exist" );

  // Ask for a convergence normalized by the right hand side
  solver.SetAztecOption( AZ_conv, AZ_rhs );

  // Control output
  switch( m_parameters.logLevel )
  {
  case 1:
    solver.SetAztecOption( AZ_output, AZ_summary );
    solver.SetAztecOption( AZ_diagnostics, AZ_all );
    break;
  case 2:
    solver.SetAztecOption( AZ_output, AZ_all );
    solver.SetAztecOption( AZ_diagnostics, AZ_all );
    break;
  default:
    solver.SetAztecOption( AZ_output, AZ_none );
  }

  // Actually solve
  solver.Iterate( m_parameters.krylov.maxIterations,
                  m_parameters.krylov.tolerance );

  //TODO: should we return performance feedback to have GEOSX pretty print details?:
  //      i.e. iterations to convergence, residual reduction, etc.
}

} // end geosx namespace
