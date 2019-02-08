/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * PetscSolver.cpp
 *
 *  Created on: Feb. 8, 2019
 *  Author: Hannah Morgan
 */

// BEGIN_RST_NARRATIVE PetscSolver.rst
// ==============================
// Petsc Solver
// ==============================
// This class implements solvers from the PETSc library. 

// Include the corresponding header file.
#include "PetscSolver.hpp"

// Put everything under the geosx namespace.
namespace geosx
{
/**
 * @brief Empty solver constructor.
 *
 * Create an empty (distributed) vector.
 */

// ----------------------------
// Constructors
// ----------------------------

// """""""""""""""""""""""""""""""""""""""""""""""""""""""""
// Empty constructor
// """""""""""""""""""""""""""""""""""""""""""""""""""""""""

PetscSolver::PetscSolver()
{}

/**
 * @brief Solve system.
 *
 * Solve Ax=b with A an PetscSparseMatrix, x and b PetscVector. Prec is an optional
 * pointer to a preconditioner.
 */

// ----------------------------
// Iterative solver
// ----------------------------

void PetscSolver::solve( MPI_Comm const comm,
                         PetscSparseMatrix &Mat,
                         PetscVector &sol,
                         PetscVector &rhs,
                         integer const max_iter,
                         real64 const newton_tol,
                         PC Prec )
{
  // create linear solver
  KSP ksp;
  KSPCreate(comm, &ksp);

  KSPSetOperators(ksp, M.getMat(), M.getMat());
  KSPSetType(ksp, KSPGMRES);

  if (pc != NULL) 
  {
    // use given preconditioner
    KSPGetPC(ksp, &pc);
  } 

  // solve with default preconditioner: block Jacobi with ILU(0) on each block
  KSPSetTolerances(ksp, newton_tol, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
  KSPSolve(ksp, rhs.getVec(), sol.getVec());
}

/**
 * @brief Solve system using an ML preconditioner. (Testing purposes mostly)
 *
 * Solve Ax=b with A an PetscSparseMatrix, x and b PetscVector.
 */

// ----------------------------
// Iterative solver with ML
// ----------------------------
// Initial test of an ML-based preconditioner.

void PetscSolver::ml_solve( MPI_Comm const comm,
                            PetscSparseMatrix &Mat,
                            PetscVector &sol,
                            PetscVector &rhs,
                            integer const max_iter,
                            real64 const newton_tol,
                            PC MLPrec )
{
  // create linear solver
  KSP ksp;
  KSPCreate(comm, &ksp);

  KSPSetOperators(ksp, M.getMat(), M.getMat());
  KSPSetType(ksp, KSPGMRES);

  if (pc != NULL) 
  {
    // use given preconditioner
    KSPGetPC(ksp, &pc);

  } else {

    // ML preconditioner
    PC prec;
    KSPGetPC(ksp, &prec);
    PCSetType(prec, PCML);

    // or use PETSc geometric algebraic multigrid
    // PCSetType(prec, PCGAMG); 
    // MatSetBlockSize(M.getMat(), M.myRows()); // must be the same on all processes
  }

  // solve system
  KSPSetTolerances(ksp, newton_tol, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
  KSPSolve(ksp, rhs.getVec(), sol.getVec());
}

/**
 * @brief Solve system using a direct solver (sequential!).
 *
 * Solve Ax=b with A an PetscSparseMatrix, x and b PetscVectors.
 */

// ----------------------------
// Direct solver
// ----------------------------

void PetscSolver::dsolve( MPI_Comm const comm,
                          PetscSparseMatrix &Mat,
                          PetscVector &sol,
                          PetscVector &rhs )
{
  // create linear solver
  KSP ksp;
  KSPCreate(comm, &ksp);

  KSPSetOperators(ksp, M.getMat(), M.getMat());
  KSPSetType(ksp, KSPPREONLY);

  // use direct solve preconditioner
  PC prec;
  KSPGetPC(ksp, &prec);
  PCSetType(prec, PCLU);
  PCFactorSetMatSolverType(prec, MATSOLVERMUMPS);
  // more MUMPS parameters can go here

  // solve system
  KSPSolve(ksp, rhs.getVec(), sol.getVec());
}

}
