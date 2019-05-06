#include "PetscSolver.hpp"

// constructor
PetscSolver::PetscSolver( LinearSolverParameters const & parameters )
  :
  m_parameters( parameters )
{}

// top level solver
void PetscSolver::solve( MPI_Comm const comm,
                         PetscSparseMatrix &mat,
                         PetscVector &sol,
                         PetscVector &rhs )
{
  // direct or iterative
  if( m_parameters.solverType == "direct" )
    solve_direct( comm, mat, sol, rhs );
  else
    solve_krylov( comm, mat, sol, rhs );
}

// direct solver
void PetscSolver::solve_direct( MPI_Comm const comm,
                                PetscSparseMatrix &mat,
                                PetscVector &sol,
                                PetscVector &rhs )
{
  // create linear solver
  KSP ksp;
  KSPCreate(comm, &ksp);
  KSPSetOperators(ksp, mat.getMat(), mat.getMat());
  KSPSetType(ksp, KSPPREONLY);

  // use direct solve preconditioner MUMPS
  PC prec;
  KSPGetPC(ksp, &prec);
  PCSetType(prec, PCLU);
  PCFactorSetMatSolverType(prec, MATSOLVERMUMPS);
  // more MUMPS parameters can go here

  // solve system
  KSPSetFromOptions(ksp);
  KSPSolve(ksp, rhs.getVec(), sol.getVec());
}

// iterative solver
void PetscSolver::solve_krylov( MPI_Comm const comm,
                                PetscSparseMatrix &mat,
                                PetscVector &sol,
                                PetscVector &rhs )
{
  // create linear solver
  KSP ksp;
  KSPCreate(comm, &ksp);
  KSPSetOperators(ksp, mat.getMat(), mat.getMat());
  KSPGMRESSetRestart(ksp, m_parameters.krylov.maxRestart);
  KSPSetTolerances(ksp, m_parameters.krylov.tolerance, PETSC_DEFAULT, PETSC_DEFAULT, m_parameters.krylov.maxIterations);

  // pick the solver type
  if( m_parameters.solverType == "gmres" )
  {
    KSPSetType(ksp, KSPGMRES);
  }
  else if( m_parameters.solverType == "bicgstab" )
  {
    KSPSetType(ksp, KSPBCGS);
  }
  else if( m_parameters.solverType == "cg" )
  {
    KSPSetType(ksp, KSPCG);
  }
  else
  {
    printf("The requested linear solverType doesn't seem to exist\n");
    // GEOS_ERROR( "The requested linear solverType doesn't seem to exist" );
  }
  
  // create a preconditioner and pick type
  PC prec;
  KSPGetPC(ksp, &prec);

  if( m_parameters.preconditionerType == "none" )
  {
    PCSetType(prec, PCNONE);
  }
  else if( m_parameters.preconditionerType == "jacobi" )
  {
    PCSetType(prec, PCJACOBI);
  }
  else if( m_parameters.preconditionerType == "ilu" )
  {
    printf("The requested linear preconditionerType isn't availbe in PETSc\n");
    // GEOS_ERROR( "The requested linear preconditionerType isn't available in PETSc" );
  }
  else if( m_parameters.preconditionerType == "icc" )
  {
    printf("The requested linear preconditionerType isn't availbe in PETSc\n");
    // GEOS_ERROR( "The requested linear preconditionerType isn't available in PETSc" );
  }
  else if( m_parameters.preconditionerType == "ilut" )
  {
    PCSetType(prec, PCHYPRE);
    PCHYPRESetType(prec, "pilut");
    // HYPRE Pilut Options
    // -pc_hypre_pilut_maxiter <-2>: Number of iterations (None)
    // -pc_hypre_pilut_tol <-2.>: Drop tolerance (None)
    // -pc_hypre_pilut_factorrowsize <-2>: FactorRowSize (None)
  }
  else if( m_parameters.preconditionerType == "amg" )
  {
    std::map<std::string, std::string> translate; // maps GEOSX to PETSc syntax for Hyper options

    translate.insert( std::make_pair( "jacobi", "Jacobi" ));
    translate.insert( std::make_pair( "sequentialGaussSeidel", "sequential-Gauss-Seidel" ));
    translate.insert( std::make_pair( "seqboundaryGaussSeidel", "seqboundary-Gauss-Seidel" ));
    translate.insert( std::make_pair( "sorJacobi", "SOR/Jacobi" ));
    translate.insert( std::make_pair( "backwardSorJacobi", "backward-SOR/Jacobi" ));
    translate.insert( std::make_pair( "symmetricSorJacobi", "symmetric-SOR/Jacobi" ));
    translate.insert( std::make_pair( "l1scaledSorJacobi", "l1scaled-SOR/Jacobi" ));
    translate.insert( std::make_pair( "direct", "Gaussian-elimination" ));
    translate.insert( std::make_pair( "l1GaussSeidel", "l1-Gauss-Seidel" ));
    translate.insert( std::make_pair( "backwardl1GaussSeidel", "backward-l1-Gauss-Seidel" ));
    translate.insert( std::make_pair( "cg", "CG" ));
    translate.insert( std::make_pair( "chebyshev", "Chebyshev" ));
    translate.insert( std::make_pair( "fcfJacobi", "FCF-Jacobi" ));
    translate.insert( std::make_pair( "l1scaledJacobi", "l1scaled-Jacobi" ));

    PCSetType(prec, PCHYPRE);
    PCHYPRESetType(prec, "boomeramg");
  
    // PETSc needs char[]
    char max_levels[10], cycle_type[10], num_sweeps[10], smoother_type[30], coarse_type[30];
    sprintf(max_levels, "%d", m_parameters.amg.maxLevels);
    sprintf(cycle_type, "%s", m_parameters.amg.cycleType.c_str());
    sprintf(num_sweeps, "%d", m_parameters.amg.maxLevels);
    sprintf(smoother_type, "%s", translate[m_parameters.amg.smootherType].c_str());
    sprintf(coarse_type, "%s", translate[m_parameters.amg.coarseType].c_str());

    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_max_levels", max_levels); 
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_cycle_type", cycle_type); 
    // relaxation method
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_all", smoother_type); // <symmetric-SOR/Jacobi> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination    l1-Gauss-Seidel backward-l1-Gauss-Seidel CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_relax_type_coarse", coarse_type); // <Gaussian-elimination> (choose one of) Jacobi sequential-Gauss-Seidel seqboundary-Gauss-Seidel SOR/Jacobi backward-SOR/Jacobi  symmetric-SOR/Jacobi  l1scaled-SOR/Jacobi Gaussian-elimination    l1-Gauss-Seidel backward-l1-Gauss-Seidel CG Chebyshev FCF-Jacobi l1scaled-Jacobi (None)
    // number of relaxation sweeps
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_grid_sweeps_all", num_sweeps); 
    PetscOptionsSetValue(NULL, "-pc_hypre_boomeramg_grid_sweeps_coarse", num_sweeps); // coarsest grid
  }
  else
  {
    printf("The requested linear preconditionerType isn't availbe in PETSc\n");
    // GEOS_ERROR( "The requested preconditionerType isn't availbe in exist" );
  }

  // HANNAH: in PETSc? Ask for a convergence normalized by the right hand side
  // KSPSetNormType(ksp, KSP_NORM_PRECONDITIONED);

  // display output
  if (m_parameters.verbosity > 0)
  {
    PetscOptionsSetValue(NULL, "-ksp_monitor", NULL); 
    // PetscOptionsSetValue(NULL, "-log_view", "true"); // Hannah: not working?
  }

  // Actually solve
  KSPSetFromOptions(ksp); 
  KSPSolve(ksp, rhs.getVec(), sol.getVec());
}

