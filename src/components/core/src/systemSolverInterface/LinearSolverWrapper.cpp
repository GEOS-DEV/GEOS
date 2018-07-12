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
 * LinearSolverWrapper.cpp
 *
 *  Created on: Sep 12, 2017
 *      Author: settgast
 */

#include "LinearSolverWrapper.hpp"


#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_SolverMap_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "AztecOO.h"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"
#include "ml_RowMatrix.h"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Amesos_BaseSolver.h"
#include "Amesos.h"
#include "Thyra_OperatorVectorClientSupport.hpp"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_MLPreconditionerFactory.hpp"


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace geosx
{
namespace systemSolverInterface
{

LinearSolverWrapper::LinearSolverWrapper():
#if USE_MPI
  m_epetraComm(MPI_COMM_WORLD)
#else
  m_epetraComm()
#endif
{
  // TODO Auto-generated constructor stub

}

LinearSolverWrapper::~LinearSolverWrapper()
{
  // TODO Auto-generated destructor stub
}


void LinearSolverWrapper::SolveSingleBlockSystem( EpetraBlockSystem * const blockSystem,
                                                  SystemSolverParameters const * const params,
                                                  BlockIDs const blockID)
{

  // initial guess for solver
  int dummy;
  double* local_solution = nullptr;

  //  if(m_numerics.m_verbose)
  //  {
  //    EpetraExt::RowMatrixToMatlabFile("umatrix.dat",*m_matrix);
  //    EpetraExt::MultiVectorToMatlabFile("urhs.dat",*m_rhs);
  //  }

  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( blockID,
                                                              blockID );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( blockID );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( blockID );


  solution->ExtractView(&local_solution,&dummy);

//    SetInitialGuess( domain, local_solution  );

  if(params->scalingOption())
  {
    //      printf("Scaling matrix/rhs in place\n");
    Epetra_Vector scaling(matrix->RowMap());
    matrix->InvRowSums(scaling);
    matrix->LeftScale(scaling);

    Epetra_MultiVector tmp (*rhs);
    rhs->Multiply(1.0,scaling,tmp,0.0);
  }
//  matrix->Print(std::cout);
//  rhs->Print(std::cout);

  Epetra_LinearProblem problem( matrix,
                                solution,
                                rhs );


  // @annavarapusr1: Needed to use direct solver without changing it for
  // everyone else
  if(params->useDirectSolver())   // If Chandra's test problems, use direct
                                  // solver
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
      MLList.set("ML output",1);
      MLList.set("PDE equations",3);
      MLList.set("ML output", 0);
      MLPrec = std::make_unique<ML_Epetra::MultiLevelPreconditioner>(*matrix, MLList);
      solver.SetPrecOperator(MLPrec.get());
    }
    else   // use ILUT preconditioner with domain decomp
    {
      solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
      solver.SetAztecOption(AZ_subdomain_solve,AZ_ilut);
      solver.SetAztecParam(AZ_ilut_fill,params->ilut_fill());
      solver.SetAztecParam(AZ_drop,params->ilut_drop());
      solver.SetAztecParam(AZ_output,0);
    }

//    std::cout<<params->numKrylovIter()<<std::endl;
//    std::cout<<params->krylovTol()<<std::endl;

    solver.Iterate(params->numKrylovIter(),
                   params->krylovTol() );

  }

//    // copy vector solution into geos data structures
//
//    realT scalingFactor = CheckSolution( local_solution, domain, 0 );
//    PropagateSolution( local_solution, scalingFactor, domain, 0 );
//
//    // re-sync ghost nodes
//
//    partition.SynchronizeFields(m_syncedFields,
// CommRegistry::lagrangeSolver02);
//
//    // copy vector solution into geos data structures
//
//    PostSyncConsistency( domain, partition );

}



}
} /* namespace geosx */
