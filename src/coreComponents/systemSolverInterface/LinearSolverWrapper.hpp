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
 * LinearSolverWrapper.hpp
 *
 *  Created on: Sep 12, 2017
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_LINEARSOLVERWRAPPER_HPP_
#define SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_LINEARSOLVERWRAPPER_HPP_

#include "LinearSystemRepository.hpp"
#include "systemSolverInterface/SystemSolverParameters.hpp"
#include "linearAlgebraInterface/src/TrilinosInterface.hpp"

namespace geosx
{
namespace systemSolverInterface
{

class LinearSolverWrapper
{
public:
  LinearSolverWrapper();
  virtual ~LinearSolverWrapper();

  template< typename LAI >
  void SolveSingleBlockSystem( LinearSystemRepository * const system,
                               SystemSolverParameters const * const params,
                               string const & blockID );

#ifdef GEOSX_USE_MPI
  Epetra_MpiComm m_epetraComm;
#else
  Epetra_SerialComm m_epetraComm;
#endif

};



template< typename LAI = TrilinosInterface >
void LinearSolverWrapper::SolveSingleBlockSystem( LinearSystemRepository * const blockSystem,
                                                  SystemSolverParameters const * const params,
                                                  string const & blockID)
{

  // initial guess for solver
//  int dummy;
//  double* local_solution = nullptr;


  typename LAI::ParallelMatrix * const matrix = blockSystem->GetMatrix<typename LAI::ParallelMatrix>( blockID,
                                                                                             blockID );
  typename LAI::ParallelVector * const rhs = blockSystem->GetResidualVector<typename LAI::ParallelVector>( blockID );
  typename LAI::ParallelVector * const solution = blockSystem->GetSolutionVector<typename LAI::ParallelVector>( blockID );


//  solution->ExtractView(&local_solution,&dummy);
//  local_solution = solution->getValues();

  if(params->scalingOption())
  {
//    Epetra_Vector scaling(matrix->RowMap());
//    matrix->InvRowSums(scaling);
//    matrix->LeftScale(scaling);
//
//    Epetra_MultiVector tmp (*rhs);
//    rhs->Multiply(1.0,scaling,tmp,0.0);
  }

  Epetra_LinearProblem problem( matrix->getPointer(),
                                solution->getPointer(),
                                rhs->getPointer() );


  if(params->useDirectSolver())
  {
    Amesos_BaseSolver* Solver;
    Amesos Factory;

    std::string SolverType = "Klu";

    if( !(params->solverType().compare("Slu")) )
    {
      SolverType = "Amesos_Superludist";
    }
    Solver = Factory.Create(SolverType, problem);

    int ierr = Solver->SymbolicFactorization();
    if (ierr > 0)
      std::cerr << "ERROR!" << std::endl;

    ierr = Solver->NumericFactorization();
    if (ierr > 0)
      std::cerr << "ERROR!" << std::endl;

    Solver->Solve();

    if( Solver!=nullptr )
      delete Solver;
  }
  else
  {
    AztecOO solver(problem);

//      SetLinearSolverParameters( solver );
    std::unique_ptr<ML_Epetra::MultiLevelPreconditioner> MLPrec;

    if( params->useMLPrecond() )
    {
      Teuchos::ParameterList MLList;
      ML_Epetra::SetDefaults("SA",MLList);
      //MLList.set("aggregation: type", "Uncoupled");
      //MLList.set("aggregation: type", "MIS");
      MLList.set("prec type", "MGW");
      MLList.set("smoother: type","ILU");
      MLList.set("PDE equations",3);
      MLList.set("ML output", params->verbose());
      MLPrec = std::make_unique<ML_Epetra::MultiLevelPreconditioner>( *(matrix->getPointer()), MLList);
      solver.SetPrecOperator(MLPrec.get());




//      Teuchos::ParameterList list;
//                             ML_Epetra::SetDefaults("SA",list);
//                             list.set("ML output",0);
//                             list.set("max levels",parameters.uu.max_levels);
//                             list.set("aggregation: type",parameters.uu.aggregation_type);
//                             list.set("smoother: type",parameters.uu.smoother_type);
//                             list.set("coarse: type",parameters.uu.coarse_type);
//                             list.set("smoother: sweeps",parameters.uu.smoother_sweeps); // Chebyshev polynomial order
//                             list.set("coarse: sweeps",parameters.uu.coarse_sweeps);
//                             list.set("aggregation: threshold",parameters.uu.aggregation_threshold);
//                             list.set("PDE equations",3);
//                             list.set("prec type",parameters.uu.cycle_type);
//                             list.set("coarse: max size",parameters.uu.max_coarse_size);
//                             list.set("smoother: damping factor",parameters.uu.damping);
//                             list.set("coarse: damping factor",parameters.uu.damping);
//
//                             list.set("null space: type","pre-computed");
//                             list.set("null space: vectors",&rigid_body_modes[0]);
//                             list.set("null space: dimension", n_rbm);
    }
    else   // use ILUT preconditioner with domain decomp
    {
      solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
      solver.SetAztecOption(AZ_subdomain_solve,AZ_ilut);
      solver.SetAztecParam(AZ_ilut_fill,params->ilut_fill());
      solver.SetAztecParam(AZ_drop,params->ilut_drop());
    }

    solver.SetAztecOption(AZ_output,params->verbose());
    solver.Iterate(params->numKrylovIter(),
                   params->krylovTol() );

  }
}

}




} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_SYSTEMSOLVERINTERFACE_LINEARSOLVERWRAPPER_HPP_
        */
