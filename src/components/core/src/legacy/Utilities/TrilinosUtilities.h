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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef TRILINOSUTILITIES_H_
#define TRILINOSUTILITIES_H_

// Trilinos libraries
#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_Map.h"

#include "Epetra_Export.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "AztecOO.h"

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"


// GEOS libraries
#include "Common/Common.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "RCVSparse.h"
#include "ArrayT/array.h"
#include "ArrayT/Array2dT.h"

#include <string>

#include "NOX.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_Group.H"

/// Convert row column value format to epetra FECrs Matrix
void RCVToEpetraFECrsMatrix(const array<rcv>& K, Epetra_FECrsMatrix& A);

/// Convert array to epetra finite element vector
void ArrayToEpetraFEVector(const array<real64>& a, Epetra_FEVector& v);

/// Convert epetra finite element vector to real array.
void EpetraFEVectorToArray(Epetra_FEVector& v, array<real64>& a);

///
void WriteEpetraFECrsMatrixToMatlabFile(std::string filename, Epetra_FECrsMatrix& A);
///
void WriteEpetraFEVectorToMatlabFile(const std::string& filename, const Epetra_FEVector& X,const std::string&  descr1="", const std::string& descr2="");

/// Conjugate gradient linear solver
void LinSolve_CG(const array<rcv>& A,
                 array<real64>& x,
                 const array<real64>& b,
                 const Epetra_Comm* commPtr,
                 realT tol=1e-10, int maxNumIters=1000,
                 bool verboseFlag=false, bool doDataWrite=false);

/// Biconjugate gradient linear solver
void LinSolve_BICGSTAB(const array<rcv>& A,
                       array<real64>& x,
                       const array<real64>& b,
                       const Epetra_Comm* commPtr,
                       realT tol=1e-10, int maxNumIters=1000,
                       bool verboseFlag=false, bool doDataWrite=false);

/// Dense local linear solver.
int LinSolve_Local(rArray2d& A, array<real64>& x,array<real64>& b);

void NonlinearSolve( Epetra_Vector& x,
                     NOX::Epetra::Interface::Required& iReq,
                     NOX::Epetra::Interface::Jacobian& iJac,
                     const bool printFlag,
                     const Epetra_MpiComm& comm);

class Epetra_System;
int ConstructSchurEstimate( Epetra_System& epetraBlockSystem );

/// Wrapper for the Trilinos solver
class TrilinosLinearSolver
{
public:

  TrilinosLinearSolver():
    m_tol(1e-10),
    m_maxIters(1000),
    m_solver("bicgstab"),
    m_subdomain_solve("ilu"),
    m_precond("dom_decomp"),
    m_conv("rhs"),
    m_verboseFlag(false),
    m_useMLPrecond(false)
  {
    //empty
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn){
    m_tol = hdn->GetAttributeOrDefault<realT>("tol",m_tol);
    m_maxIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);
    m_verboseFlag = hdn->GetAttributeOrDefault<bool>("verbose",m_verboseFlag);
    m_useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",m_useMLPrecond);

    m_solver = hdn->GetAttributeStringOrDefault("solver",m_solver);
    m_subdomain_solve = hdn->GetAttributeStringOrDefault("subdomain_solve",m_subdomain_solve);
    m_precond = hdn->GetAttributeStringOrDefault("precond",m_precond);
    m_conv = hdn->GetAttributeStringOrDefault("conv",m_conv);

  };

  realT m_tol; // Solver convergence criteria
  int m_maxIters; // Maximum number of solver iterations

  std::string m_solver;
  std::string m_subdomain_solve;
  std::string m_precond;
  std::string m_conv;

  bool m_verboseFlag;
  bool m_useMLPrecond;

  void Solve(Teuchos::RCP<Epetra_FECrsMatrix>& matrix,
             Teuchos::RCP<Epetra_FEVector>&    solution,
             Teuchos::RCP<Epetra_FEVector>& rhs){

    Epetra_LinearProblem problem(&(*matrix),
                                 &(*solution),
                                 &(*rhs));

    // ML preconditioner
    //////////////////////

    ML_Epetra::MultiLevelPreconditioner* MLPrec;
    if(m_useMLPrecond)
    {

      // create a parameter list for ML options
      Teuchos::ParameterList MLList;

      ML_Epetra::SetDefaults("SA",MLList);

      MLPrec = new ML_Epetra::MultiLevelPreconditioner(*matrix, MLList);
    }

    // Aztec solver
    /////////////////
    AztecOO solver(problem);

    // TODO - use parameterlist interface?
    // Multi level preconditioner
    if(m_useMLPrecond) solver.SetPrecOperator(MLPrec);

    // set solver
    if(m_solver == "bicgstab")
    {
      solver.SetAztecOption(AZ_solver,AZ_bicgstab);
    }
    else if(m_solver == "cg")
    {
      solver.SetAztecOption(AZ_solver,AZ_cg);
    }
    else if(m_solver == "gmres")
    {
      solver.SetAztecOption(AZ_solver,AZ_gmres);
    }
    else
    {
      throw GPException("Error TrilinosLinearSolver: Unknown solver " + m_solver );
    }

    // set preconditioner
    if(m_precond == "dom_decomp")
    {
      solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    }

    // set subdomain solver
    if(m_subdomain_solve == "ilu")
    {
      solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
    }
    else if (m_subdomain_solve == "icc")
    {
      solver.SetAztecOption(AZ_subdomain_solve,AZ_icc);
    }
    else if (m_subdomain_solve == "ilut")
    {
      solver.SetAztecOption(AZ_subdomain_solve,AZ_ilut);
    }
    else
    {
      throw GPException("Error TrilinosLinearSolver: Unknown subdomain solver " + m_subdomain_solve );
    }

    // set convergence criteria
    if(m_conv == "rhs")
    {
      solver.SetAztecOption(AZ_conv,AZ_rhs);
    }
    else if(m_conv == "r0")
    {
      solver.SetAztecOption(AZ_conv,AZ_r0);
    }
    else if(m_conv == "Anorm")
    {
      solver.SetAztecOption(AZ_conv,AZ_Anorm);
    }
    else if(m_conv == "noscaled")
    {
      solver.SetAztecOption(AZ_conv,AZ_noscaled);
    }
    else
    {
      throw GPException("Error TrilinosLinearSolver: Unknown convergence criteria " + m_conv );
    }

    // verbosity
    if(!m_verboseFlag) solver.SetAztecOption(AZ_output,AZ_none);

    solver.Iterate(m_maxIters,m_tol);

    // destroy the preconditioner
    if(m_useMLPrecond)
    {
      delete MLPrec;
    }
  }
};


namespace TrilinosUtilities
{
int RowSum( const Epetra_CrsMatrix& A, Epetra_Vector& x);
int RowSum( const Epetra_CrsMatrix& A, const Epetra_CrsMatrix& B, Epetra_Vector& x);
}

#endif
