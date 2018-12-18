
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

  KSP ksp;

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, M.getMat(), M.getMat());
  KSPSetType(ksp, KSPGMRES);
  // KSPSetTolerances(ksp, newton_tol, newton_tol, NULL, max_iter);
  // KSPSetUp(ksp);

  KSPSolve(ksp, rhs.getVec(), sol.getVec());
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
void PETScSolver::dsolve( PETScSparseMatrix &M,
                          PETScVector &rhs,
                          PETScVector &sol )
{
  KSP ksp;
  PC pc;

  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, M.getMat(), M.getMat());
  PCSetType(pc,PCILU);
  KSPSetType(ksp,KSPPREONLY);
  // KSPSetUp(ksp);
  
  KSPSolve(ksp, rhs.getVec(), sol.getVec());
}


