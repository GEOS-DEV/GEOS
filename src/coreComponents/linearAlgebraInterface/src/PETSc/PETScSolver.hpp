
#include "PETScSparseMatrix.hpp"
// #include "PETScVector.hpp"

#include <petscksp.h>

class PETScSolver
{
public:

  /* Empty constructor */
  PETScSolver();

  /* Copy constructor */
  // PETScSolver( const PETScSolver &Solver );

  /* Virtual destructor */
  // virtual ~PETScSolver() = default;

  /* Solve Ax=b with A an PETScSparseMatrix, x and b PETScVector */
  void solve( PETScSparseMatrix &Mat,
              PETScVector &rhs,
              PETScVector &sol,
              int max_iter,
              double newton_tol );
              // preconditioner? Prec = nullptr 

  // /* Solve Ax=b with A an PETScSparseMatrix, x and b PETScVector using ml preconditioner */
  // void ml_solve( PETScSparseMatrix &Mat,
  //                PETScVector &rhs,
  //                PETScVector &sol,
  //                int max_iter,
  //                double newton_tol,
  //                std::unique_ptr<ML_Epetra::MultiLevelPreconditioner> MLPrec = nullptr );

  /* Solve Ax=b with A an PETScSparseMatrix, x and b PETScVector, direct solve */
  void dsolve( PETScSparseMatrix &Mat,
               PETScVector &rhs,
               PETScVector &sol );

protected:

};




