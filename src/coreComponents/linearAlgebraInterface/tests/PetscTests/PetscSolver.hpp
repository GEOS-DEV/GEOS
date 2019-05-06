#include <petscksp.h>
#include <petscsys.h>
#include <map>
#include <string>
#include "PetscSparseMatrix.hpp"
#include "LinearSolverParameters.hpp"

class PetscSolver
{
public:

  // parameter list
  PetscSolver( LinearSolverParameters const & parameters );

  // destructor
  // virtual ~PetscSolver() = default;

  // Solve Ax=b with A an PetscSparseMatrix, x and b PetscVector
  void solve( MPI_Comm const comm,
              PetscSparseMatrix &mat,
              PetscVector &sol,
              PetscVector &rhs );

private:

  LinearSolverParameters const & m_parameters;

  void solve_direct( MPI_Comm const comm,
                     PetscSparseMatrix &mat,
                     PetscVector &sol,
                     PetscVector &rhs );

  void solve_krylov( MPI_Comm const comm,
                     PetscSparseMatrix &mat,
                     PetscVector &sol,
                     PetscVector &rhs );

};


