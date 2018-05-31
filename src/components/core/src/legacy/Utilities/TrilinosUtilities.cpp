// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
#include "TrilinosUtilities.h"

#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

/// Convert row column value format to epetra FECrs Matrix
void RCVToEpetraFECrsMatrix(const array<rcv>& K, Epetra_FECrsMatrix& A)
{

  for (size_t i = 0 ; i < K.size() ; ++i)
  {
    int row = K[i].r;
    int col = K[i].c;
    double val = K[i].v;
    A.InsertGlobalValues(1, &row, 1, &col, &val);
    // fixme ideally would be sum - but doesn't work if entry does not exist?
    // should consolidate sparse array first as a result
  }

}

/// Convert array to epetra finite element vector
void ArrayToEpetraFEVector(const array<real64>& a, Epetra_FEVector& v)
{
  if (a.size() >= INT_MAX)
    throw GPException(
            "TrilinosUtilities:ArrayToEpetraFEVector: invalid conversion from unsigned int to int");
  for (unsigned int i = 0 ; i < a.size() ; ++i)
  {
    int ii = i;
    v.ReplaceGlobalValues(1, &ii, &(a[i]));
  }
}

/// Convert epetra finite element vector to real array.
void EpetraFEVectorToArray(Epetra_FEVector& v, array<real64>& a)
{
  a.resize(v.MyLength());
  double ** vp;
  v.ExtractView(&vp);
  for (size_t i = 0 ; i < a.size() ; ++i)
    a[i] = vp[0][i]; // first vector only
}

///
void WriteEpetraFECrsMatrixToMatlabFile(std::string filename, Epetra_FECrsMatrix& A)
{
  EpetraExt::RowMatrixToMatlabFile(filename.c_str(), A);
}

///
void WriteEpetraFEVectorToMatlabFile(const std::string& filename, const Epetra_FEVector& X,
                                     const std::string& descr1, const std::string& descr2)
{
  EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(), X, descr1.c_str(), descr2.c_str(),
                                           false);
}

/// Solve Ax = b with a conjugate gradient solver and jacobi precoditioning
void LinSolve_CG(const array<rcv>& A, array<real64>& x, const array<real64>& b, const Epetra_Comm* commPtr,
                 realT tol, int maxNumIters, bool verboseFlag, bool doDataWrite)
{

  int myPID = commPtr->MyPID();

  const long long int numDofs = (long long int) x.size();
  const long long int tmp = static_cast<long long int>(numDofs);
  int arrayEntriesPerDof = ceil(double(A.size()) / double(numDofs));
  Epetra_Map epMap(tmp, 0, *commPtr);

  int numMyTargetElements = (myPID == 0) ? numDofs : 0;
  Epetra_Map epTargetMap(-1, numMyTargetElements, 0, *commPtr);
  Epetra_Export gatherer(epMap, epTargetMap);

  // Ax=b for solver
  Epetra_FECrsMatrix eA(Copy, epMap, epMap, arrayEntriesPerDof);
  Epetra_FEVector ex(epMap);
  Epetra_FEVector eb(epMap);

  Epetra_FEVector ex0(epTargetMap); // fixme -> ex0 = used to communicate x
                                    // between processor 0 and other pids
                                    // (redundant)

  if (myPID == 0)
  {
    std::cout << "Running matrix solver" << std::endl;
    // copy problem to epetra format
    RCVToEpetraFECrsMatrix(A, eA);
    ArrayToEpetraFEVector(b, eb);
    ArrayToEpetraFEVector(x, ex);
  }

  eA.GlobalAssemble();
  eA.FillComplete(epMap, epMap);

  Epetra_LinearProblem problem(&eA, &ex, &eb);

  // solve problem
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
//  solver.SetAztecOption(AZ_output, AZ_none);
  solver.SetAztecOption(AZ_solver, AZ_cg);

  solver.Iterate(maxNumIters, tol);

  if (verboseFlag)
  {
    std::cout << "Iterations " << solver.NumIters() << std::endl;
    std::cout << "Residual " << solver.TrueResidual() << std::endl;
    std::cout << "Scaled Residual " << solver.ScaledResidual() << std::endl;
    std::cout << "Condition number " << solver.GetAztecOption(AZ_condnum) << std::endl;
  }

  if (doDataWrite)
  {
    WriteEpetraFEVectorToMatlabFile("b.dat", eb, "", "");
    WriteEpetraFEVectorToMatlabFile("x.dat", ex, "", "");
    WriteEpetraFECrsMatrixToMatlabFile("A.dat", eA);
    WriteSparseArrayToFile("K.txt", A);
    exit(0);
  }

  // collect result on PID 0  < ----- fixme should keep on individual processors
  // once parallel version in place
  ex0.Export(ex, gatherer, Add);
  if (myPID == 0)
    EpetraFEVectorToArray(ex0, x);

}

/// Solve Ax = b with a biconjugate gradient stabilized solver and jacobi
// precoditioning
void LinSolve_BICGSTAB(const array<rcv>& A, array<real64>& x, const array<real64>& b,
                       const Epetra_Comm* commPtr, realT tol, int maxNumIters, bool verboseFlag,
                       bool doDataWrite)
{

  int myPID = commPtr->MyPID();

  const long long int numDofs = (long long int) x.size();
  const long long int tmp = static_cast<long long int>(numDofs);
  int arrayEntriesPerDof = ceil(double(A.size()) / double(numDofs));
  Epetra_Map epMap(tmp, 0, *commPtr);

  int numMyTargetElements = (myPID == 0) ? numDofs : 0;
  Epetra_Map epTargetMap(-1, numMyTargetElements, 0, *commPtr);
  Epetra_Export gatherer(epMap, epTargetMap);

  // Ax=b for solver
  Epetra_FECrsMatrix eA(Copy, epMap, epMap, arrayEntriesPerDof);
  Epetra_FEVector ex(epMap);
  Epetra_FEVector eb(epMap);

  Epetra_FEVector ex0(epTargetMap); // fixme -> ex0 = used to communicate x
                                    // between processor 0 and other pids
                                    // (redundant)

  if (myPID == 0)
  {
    std::cout << "Running matrix solver" << std::endl;
    // copy problem to epetra format
    RCVToEpetraFECrsMatrix(A, eA);
    ArrayToEpetraFEVector(b, eb);
    ArrayToEpetraFEVector(x, ex);
  }

  eA.GlobalAssemble();
  eA.FillComplete(epMap, epMap);

  Epetra_LinearProblem problem(&eA, &ex, &eb);

  // solve problem
  AztecOO solver(problem);
  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  // solver.SetAztecOption(AZ_output, AZ_none);
  solver.SetAztecOption(AZ_solver, AZ_bicgstab);

  solver.Iterate(maxNumIters, tol);

  if (verboseFlag)
  {
    std::cout << "Iterations " << solver.NumIters() << std::endl;
    std::cout << "Residual " << solver.TrueResidual() << std::endl;
    std::cout << "Scaled Residual " << solver.ScaledResidual() << std::endl;
    std::cout << "Condition number " << solver.GetAztecOption(AZ_condnum) << std::endl;
  }

  if (doDataWrite)
  {
    WriteEpetraFEVectorToMatlabFile("b.dat", eb, "", "");
    WriteEpetraFEVectorToMatlabFile("x.dat", ex, "", "");
    WriteEpetraFECrsMatrixToMatlabFile("A.dat", eA);
    WriteSparseArrayToFile("K.txt", A);
    exit(0);
  }

  // collect result on PID 0  < ----- fixme should keep on individual processors
  // once parallel version in place
  ex0.Export(ex, gatherer, Add);
  if (myPID == 0)
    EpetraFEVectorToArray(ex0, x);

}

/// FIXME - the following won't work on settgast's machine - some problem with
// Teuchos and compiler
//#ifdef WALSH24
#if 1
int LinSolve_Local(rArray2d& A, array<real64>& x, array<real64>& b)
{

  bool doEquilibrate = false;

  array<real64>::size_type n = x.size();

  Teuchos::SerialDenseMatrix<int, realT> A_Teuch(Teuchos::Copy, A.data(), n, n, n);

  Teuchos::SerialDenseVector<int, realT> x_Teuch(Teuchos::Copy, x.data(), n);

  Teuchos::SerialDenseVector<int, realT> b_Teuch(Teuchos::Copy, b.data(), n);

  Teuchos::SerialDenseSolver<int, realT> Solver;

  int info = 0;
  Solver.setMatrix(Teuchos::rcp(&A_Teuch, false));
  Solver.setVectors(Teuchos::rcp(&x_Teuch, false), Teuchos::rcp(&b_Teuch, false));

  Solver.solveWithTranspose(true);
  if (doEquilibrate)
    Solver.factorWithEquilibration(true); // not working?

  info = Solver.factor();
  if (info != 0)
  {
    std::cout << "Teuchos::SerialDenseSolver::factor() returned : " << info << std::endl;
    return info;
  }
  info = Solver.solve();
  if (info != 0)
    std::cout << "Teuchos::SerialDenseSolver::solve() returned : " << info << std::endl;

  //if(doEquilibrate) Solver.unequilibrateLHS(); // not working? needed?

  // copy back, not required if using view?
  for (rArray2d::size_type i = 0 ; i < n ; ++i)
    x[i] = x_Teuch[i];

  //if ( Solver.shouldEquilibrate()){
  if (false)
  {
    for (rArray2d::size_type i = 0 ; i < A.Dimension(0) ; ++i)
    {
      for (rArray2d::size_type j = 0 ; j < A.Dimension(1) ; ++j)
        std::cout << A(i, j) << " ";
      std::cout << std::endl;
    }
    std::cout << "x" << std::endl;
    for (array<real64>::size_type i = 0 ; i < x.size() ; ++i)
      std::cout << x(i) << " ";
    std::cout << std::endl;
    std::cout << "b" << std::endl;
    for (array<real64>::size_type i = 0 ; i < b.size() ; ++i)
      std::cout << b(i) << " ";
    std::cout << std::endl;
    exit(0);
  }
  return info;
}
#else
int LinSolve_Local(rArray2d& A, array<real64>& x,array<real64>& b)
{
  throw GPException("LinSolve_Local - FIXME - problem with Teuchos compilation.");
  return 1;
}
#endif

void NonlinearSolve(Epetra_Vector& x, NOX::Epetra::Interface::Required& iReq,
                    NOX::Epetra::Interface::Jacobian& iJac, const bool printFlag,
                    const Epetra_MpiComm& comm)
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;

  // Create the top-level parameter list to control NOX.
  //
  // "parameterList" (lowercase initial "p") is a "nonmember
  // constructor" that returns an RCP<ParameterList> with the
  // given name.
  RCP < ParameterList > params = parameterList("NOX");

  // Tell the nonlinear solver to use line search.
  params->set("Nonlinear Solver", "Line Search Based");

  //
  // Set the printing parameters in the "Printing" sublist.
  //
  ParameterList& printParams = params->sublist("Printing");
  printParams.set("MyPID", comm.MyPID());
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);

  // Set verbose=true to see a whole lot of intermediate status
  // output, during both linear and nonlinear iterations.
  const bool verbose = false;
  if (verbose)
  {
    printParams.set(
      "Output Information",
      NOX::Utils::OuterIteration + NOX::Utils::OuterIterationStatusTest + NOX::Utils::InnerIteration + NOX::Utils::Parameters + NOX::Utils::Details +
      NOX::Utils::Warning);
  }
  else
  {
    printParams.set("Output Information", NOX::Utils::Warning);
  }

  //
  // Set the nonlinear solver parameters.
  //

  // Line search parameters.
  ParameterList& searchParams = params->sublist("Line Search");
  searchParams.set("Method", "Full Step");

  // Parameters for picking the search direction.
  ParameterList& dirParams = params->sublist("Direction");
  // Use Newton's method to pick the search direction.
  dirParams.set("Method", "Newton");

  // Parameters for Newton's method.
  ParameterList& newtonParams = dirParams.sublist("Newton");
  newtonParams.set("Forcing Term Method", "Constant");

  //
  // Newton's method invokes a linear solver repeatedly.
  // Set the parameters for the linear solver.
  //
  ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  // Use Aztec's implementation of GMRES, with at most 800
  // iterations, a residual tolerance of 1.0e-4, with output every
  // 50 iterations, and Aztec's native ILU preconditioner.
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 800);
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 50);
  lsParams.set("Aztec Preconditioner", "ilu");

  Epetra_Map Map(x.Map().NumMyElements(), 0, comm);
  //
  // Build the initial Jacobian matrix operator.
  //
  RCP < Epetra_CrsMatrix > J0 = rcp(new Epetra_CrsMatrix(Copy, Map, Map.NumMyElements()));
  iJac.computeJacobian(x, *J0);

  RCP < NOX::Epetra::Interface::Required > rcp_iReq(&iReq);
  RCP < NOX::Epetra::Interface::Jacobian > rcp_iJac(&iJac);

  RCP < NOX::Epetra::LinearSystemAztecOO > linSys = rcp(
    new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, rcp_iReq, rcp_iJac, J0, x));

  // Need a NOX::Epetra::Vector for constructor.
  NOX::Epetra::Vector noxInitGuess(x, NOX::DeepCopy);
  RCP < NOX::Epetra::Group > group = rcp(
    new NOX::Epetra::Group(printParams, rcp_iReq, noxInitGuess, linSys));

  //
  // Set up NOX's iteration stopping criteria ("status tests").
  //

  // ||F(X)||_2 / N < 1.0e-4, where N is the length of F(X).
  //
  // NormF has many options for setting up absolute vs. relative
  // (scaled by the norm of the initial guess) tolerances, scaling
  // or not scaling by the length of F(X), and choosing a
  // different norm (we use the 2-norm here).
  RCP < NOX::StatusTest::NormF > testNormF = rcp(new NOX::StatusTest::NormF(1.0e-4));

  // At most 20 (nonlinear) iterations.
  RCP < NOX::StatusTest::MaxIters > testMaxIters = rcp(new NOX::StatusTest::MaxIters(20));

  // Combine the above two stopping criteria (normwise
  // convergence, and maximum number of nonlinear iterations).
  // The result tells NOX to stop if at least one of them is
  // satisfied.
  RCP < NOX::StatusTest::Combo > combo = rcp(
    new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, testNormF, testMaxIters));

  // Create the NOX nonlinear solver.
  RCP < NOX::Solver::Generic > solver = NOX::Solver::buildSolver(group, combo, params);

  // Solve the nonlinear system.
//  NOX::StatusTest::StatusType status = solver->solve();

  // Get the Epetra_Vector with the final solution from the solver.
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());

  x = dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX()).getEpetraVector();
}


namespace TrilinosUtilities
{

int RowSum( const Epetra_CrsMatrix& A, Epetra_Vector& x)
{
  //
  // Put inverse of the sum of absolute values of the ith row of A in x[i].
  //

  if (!A.Filled())
    EPETRA_CHK_ERR(-1);                // Matrix must be filled.
  int ierr = 0;
  int i, j;
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double * xp = (double*)x.Values();
  if (A.Graph().RangeMap().SameAs(x.Map()) && A.Exporter() != 0)
  {
    Epetra_Vector x_tmp(A.RowMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();
    for (i=0 ; i < A.NumMyRows() ; i++)
    {
      int      NumEntries = A.NumMyEntries(i);
      double * RowValues  = NULL;
      int junk;
      A.ExtractMyRowView(i,junk, RowValues );
      for (j=0 ; j < NumEntries ; j++)
        x_tmp_p[i] += std::abs(RowValues[j]);
    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *A.Exporter(), Add)); //Export partial row
                                                         // sums to x.
  }
  else if (A.Graph().RowMap().SameAs(x.Map()))
  {
    for (i=0 ; i < A.NumMyRows() ; i++)
    {
      int      NumEntries = A.NumMyEntries(i);
//      double * RowValues  = A.Values(i);
      double * RowValues  = NULL;
      int junk;
      A.ExtractMyRowView(i,junk, RowValues );
      double scale = 0.0;
      for (j=0 ; j < NumEntries ; j++)
        scale += std::abs(RowValues[j]);
      xp[i] = scale;
    }
  }
  else   // x.Map different than both Graph().RowMap() and Graph().RangeMap()
  {
    EPETRA_CHK_ERR(-2); // The map of x must be the RowMap or RangeMap of A.
  }
  A.UpdateFlops(A.NumGlobalNonzeros64());
  EPETRA_CHK_ERR(ierr);
  return(0);
}


int RowSum( const Epetra_CrsMatrix& A, const Epetra_CrsMatrix& B, Epetra_Vector& x)
{

  int ierr = 0;
  int i, j;
  x.PutScalar(0.0); // Make sure we sum into a vector of zeros.
  double * xp = (double*)x.Values();

  if (A.Graph().RangeMap().SameAs(x.Map()) && A.Exporter() != 0 &&
      B.Graph().RangeMap().SameAs(x.Map()) && B.Exporter() != 0 )
  {

    Epetra_Vector x_tmp(A.RowMap());
    x_tmp.PutScalar(0.0);
    double * x_tmp_p = (double*)x_tmp.Values();

    for (i=0 ; i < A.NumMyRows() ; i++)
    {
      int      NumEntriesA = A.NumMyEntries(i);
      int      NumEntriesB = B.NumMyEntries(i);
      double * RowValuesA  = NULL;
      double * RowValuesB  = NULL;
      int junk;
      A.ExtractMyRowView(i,junk, RowValuesA );
      B.ExtractMyRowView(i,junk, RowValuesB );
      for (j=0 ; j < NumEntriesA ; j++)
        x_tmp_p[i] += std::abs(RowValuesA[j]);
      for (j=0 ; j < NumEntriesB ; j++)
        x_tmp_p[i] += std::abs(RowValuesB[j]);

    }
    EPETRA_CHK_ERR(x.Export(x_tmp, *A.Exporter(), Add)); //Export partial row
                                                         // sums to x.
  }
  else if ( A.Graph().RowMap().SameAs(x.Map()) &&
            B.Graph().RowMap().SameAs(x.Map()) )
  {
    for (i=0 ; i < A.NumMyRows() ; i++)
    {
      int      NumEntriesA = A.NumMyEntries(i);
      int      NumEntriesB = B.NumMyEntries(i);
      double * RowValuesA  = NULL;
      double * RowValuesB  = NULL;
      int junk;
      A.ExtractMyRowView(i,junk, RowValuesA );
      B.ExtractMyRowView(i,junk, RowValuesB );
      double scale = 0.0;
      for (j=0 ; j < NumEntriesA ; j++)
        scale += std::abs(RowValuesA[j]);
      for (j=0 ; j < NumEntriesB ; j++)
        scale += std::abs(RowValuesB[j]);
      xp[i] = scale;
    }
  }
  else   // x.Map different than both Graph().RowMap() and Graph().RangeMap()
  {
    EPETRA_CHK_ERR(-2); // The map of x must be the RowMap or RangeMap of A.
  }
  A.UpdateFlops(A.NumGlobalNonzeros64());
  B.UpdateFlops(A.NumGlobalNonzeros64());
  EPETRA_CHK_ERR(ierr);
  return(0);
}


}
