
#include "PETScSolver.hpp"

/* Empty constructor */
PETScSolver::PETScSolver()
{

}

/* Copy constructor */
// PETScSolver::PETScSolver( const PETScSolver &Solver )
// {

// }

/* Virtual destructor */
// virtual ~PETScSolver() = default;

/* Solve Ax=b with A an PETScSparseMatrix, x and b PETScVector */
void PETScSolver::solve( PETScSparseMatrix &M,
                         PETScVector &rhs,
                         PETScVector &sol,
                         int max_iter,
                         double newton_tol ){

  KSP solver;
  Vec X, B;
  X = sol.getVec();
  B = rhs.getVec();

  VecView(X, PETSC_VIEWER_STDOUT_WORLD);

  KSPCreate(PETSC_COMM_WORLD, &solver);
  KSPSetOperators(solver, M.getMat(), M.getMat());
  KSPSetType(solver, KSPGMRES);
  KSPSetTolerances(solver, newton_tol, newton_tol, NULL, max_iter);

  KSPSetUp(solver);
  KSPSolve(solver, rhs.getVec(), X);

  VecView(X, PETSC_VIEWER_STDOUT_WORLD);

  PETScVector x_(X);
  sol = x_;

}

// /* Solve Ax=b with A an PETScSparseMatrix, x and b PETScVector using ml preconditioner */
// void PETScSolver::ml_solve( PETScSparseMatrix &Mat,
//                PETScVector &rhs,
//                PETScVector &sol,
//                int max_iter,
//                double newton_tol,
//                std::unique_ptr<ML_Epetra::MultiLevelPreconditioner> MLPrec = nullptr )
// {

//   KSP solver;
//   Mat A;
//   A = M.getMat();
//   Vec X, B;
//   X = sol.getVec();
//   B = rhs.getVec();

//   KSPCreate(PETSC_COMM_WORLD, &solver);
//   KSPSetOperators(solver, A, A);

//   // KSPSetType(solver, KSPCG);

//   // set preconditioner

//   KSPSetTolerances(solver, newton_tol, newton_tol, NULL, max_iter);
//   KSPSolve(solver, B, X);

//   // put X in sol

// }

/* Solve Ax=b with A an PETScSparseMatrix, x and b PETScVector, direct solve */
void PETScSolver::dsolve( PETScSparseMatrix &Mat,
                          PETScVector &rhs,
                          PETScVector &sol )
{

  KSP solver;
  PC pc;
  Vec X, B;
  X = sol.getVec();
  B = rhs.getVec();

  VecView(X, PETSC_VIEWER_STDOUT_WORLD);

  KSPCreate(PETSC_COMM_WORLD, &solver);
  KSPSetOperators(solver, Mat.getMat(), Mat.getMat());
  PCSetType(pc,PCLU);
  KSPSetType(solver,KSPPREONLY);

  KSPSetUp(solver);
  
  KSPSolve(solver, rhs.getVec(), X);

  VecView(X, PETSC_VIEWER_STDOUT_WORLD);

  PETScVector x_(X);
  sol = x_;

}


