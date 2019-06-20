/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
#ifdef GEOSX_USE_MPI
  m_epetraComm(MPI_COMM_GEOSX)
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


  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( blockID,
                                                              blockID );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( blockID );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( blockID );


  solution->ExtractView(&local_solution,&dummy);

  if(params->scalingOption())
  {
    Epetra_Vector scaling(matrix->RowMap());
    matrix->InvRowSums(scaling);
    matrix->LeftScale(scaling);

    Epetra_MultiVector tmp (*rhs);
    rhs->Multiply(1.0,scaling,tmp,0.0);
  }

  Epetra_LinearProblem problem( matrix,
                                solution,
                                rhs );


  if(params->useDirectSolver())
  {
    Amesos_BaseSolver* Solver;
    Amesos Factory;

    std::string SolverType = "Klu";

    if( !(params->solverType().compare("Slu")) )
    {
      SolverType = "Amesos_Superludist";
    }
    else if ( !(params->solverType().compare("Lapack")) )
    {
      SolverType = "Amesos_Lapack";
    }

    //GEOS_LOG( "Solver selected: " << SolverType );
    GEOS_ERROR_IF( !Factory.Query(SolverType), "Requested solver not available: " << SolverType );

    Solver = Factory.Create(SolverType, problem);

    int ierr = Solver->SymbolicFactorization();

    GEOS_WARNING_IF( ierr > 0, "Error during symbolic factorization: " << ierr );

    ierr = Solver->NumericFactorization();

    GEOS_WARNING_IF( ierr > 0, "Error during numeric factorization: " << ierr );

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
      MLPrec = std::make_unique<ML_Epetra::MultiLevelPreconditioner>(*matrix, MLList);
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
