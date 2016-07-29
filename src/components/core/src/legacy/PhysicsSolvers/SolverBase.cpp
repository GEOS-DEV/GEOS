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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file SolverBase.cpp
 * @author settgast1
 * @date Feb 10, 2011
 */

#include "SolverBase.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "Amesos.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Utilities/TrilinosUtilities.h"

SolverBase::SolverBase( const std::string& name,
                        ProblemManagerT* const pm ):
m_stabledt(),
m_courant(1.0),
m_relaxationCoefficient(0.035),
m_solverErrorFlag(false),
m_syncedFields(),
m_name(name),
m_system(nullptr),
m_trilinosIndexStr(),
m_solvers(pm->m_solvers),
m_commonFields(),
epetra_comm(&(pm->m_epetraComm))
#if USECPP11!=1
,
m_rowMap(nullptr),
m_sparsity(nullptr),
m_matrix(nullptr),
m_solution(nullptr),
m_rhs(nullptr)
#endif

{
  m_stabledt.m_maxdt = std::numeric_limits<double>::max();
}

SolverBase::~SolverBase()
{
#if USECPP11!=1
//  delete m_rowMap;
//  delete m_sparsity;
//  delete m_matrix;
//  delete m_solution;
//  delete m_rhs;
#endif
}


void SolverBase::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{

  m_numerics.krylov_tol        = hdn->GetAttributeOrDefault("tol",1.0e-6);
  m_numerics.newton_tol        = hdn->GetAttributeOrDefault("newtonTol",1.0e-6);
  m_numerics.newton_tol       = hdn->GetAttributeOrDefault<realT>("tolNewton",1e-06);
  m_numerics.m_kspace          = hdn->GetAttributeOrDefault("kspace",300);
  m_numerics.m_ilut_fill       = hdn->GetAttributeOrDefault("ilut_fill",3.0);
  m_numerics.m_ilut_drop       = hdn->GetAttributeOrDefault("ilut_drop",0.0);
  m_numerics.m_maxIters        = hdn->GetAttributeOrDefault<int>("maxSolverIterations",10000);
  m_numerics.m_useMLPrecond    = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);
  m_numerics.m_useInnerSolver  = hdn->GetAttributeOrDefault<bool>("useInnerSolver",false);
  m_numerics.m_scalingOption   = hdn->GetAttributeOrDefault<int>("scalingOption",0);
  m_numerics.m_useBicgstab     = hdn->GetAttributeOrDefault<bool>("useBicgstab",false);
  m_numerics.m_verbose         = hdn->GetAttributeOrDefault<int>("verbose",0);
  m_numerics.m_useDirectSolver = hdn->GetAttributeOrDefault<bool>("useDirectSolver",false);
  m_numerics.m_useNewtonSolve  = hdn->GetAttributeOrDefault<bool>("useNewtonSolve",false);
  m_numerics.m_maxIterNewton   = hdn->GetAttributeOrDefault<int>("maxIterNewton",30);

  m_courant = hdn->GetAttributeOrDefault<realT>("courant",0.5);
  m_relaxationCoefficient = hdn->GetAttributeOrDefault<realT>("relaxationCoefficient", 0.007);
  m_lsp = hdn->GetAttributeStringOrDefault("lsp","Klu");
}

void SolverBase::PostProcess( PhysicalDomainT& domain,
                              SpatialPartition& partition,
                              const sArray1d& namesOfSolverRegions)
{
  //Do nothing.  To be implemented by each solver as needed;
  //For operations that we don't need at each step but need for plotting.
}


void SolverBase::SetMaxStableTimeStep( const realT& time,
                                       PhysicalDomainT& domain,
                                       const sArray1d& namesOfSolverRegions,
                                       SpatialPartition& partition)
{
  TimeStep( time, 0.0, 0, domain, namesOfSolverRegions, partition, nullptr );
}


void SolverBase::WriteSilo( SiloFile& siloFile ) const
{
  siloFile.DBWriteWrapper( "m_stabledt__m_maxdt", m_stabledt.m_maxdt );

  std::map<std::string, sArray1d> syncedFields;
  for( std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d>::const_iterator i=m_syncedFields.begin() ;
      i!= m_syncedFields.end() ; ++i )
  {
    if( !(i->second.empty()) )
    {
      syncedFields[ PhysicalDomainT::GetObjectDataStructureName( i->first ) ] = i->second;
    }
  }

  siloFile.DBWriteWrapper( "syncedFields", syncedFields );

  siloFile.DBWriteWrapper( "m_name", m_name );

  siloFile.DBWriteWrapper( "m_dtInit", m_dtInit );



}

void SolverBase::ReadSilo( const SiloFile& siloFile )
{
  siloFile.DBReadWrapper( "m_stabledt__m_maxdt", m_stabledt.m_maxdt );

  std::map<std::string, sArray1d> syncedFields;
  siloFile.DBReadWrapper( "syncedFields", syncedFields );


  for( std::map<std::string, sArray1d>::const_iterator i=syncedFields.begin() ;
      i!= syncedFields.end() ; ++i )
  {
    m_syncedFields[ PhysicalDomainT::GetObjectDataStructureKey(i->first)] = i->second;
  }


  siloFile.DBReadWrapper( "m_name", m_name );

  siloFile.DBReadWrapper( "m_dtInit", m_dtInit );

}



void SolverBase::Solve ( PhysicalDomainT&  domain,
                         SpatialPartition& partition,
                         const realT time,
                         const realT dt )
{

  // initial guess for solver
  int dummy;
  double* local_solution = nullptr;

//  if(m_numerics.m_verbose)
//  {
//    EpetraExt::RowMatrixToMatlabFile("umatrix.dat",*m_matrix);
//    EpetraExt::MultiVectorToMatlabFile("urhs.dat",*m_rhs);
//  }

  m_solution->ExtractView(&local_solution,&dummy);

  SetInitialGuess( domain, local_solution  );

  if(m_numerics.m_scalingOption)
  {
//      printf("Scaling matrix/rhs in place\n");
      Epetra_Vector scaling(m_matrix->RowMap());
      m_matrix->InvRowSums(scaling); 
      m_matrix->LeftScale(scaling);

      Epetra_MultiVector tmp (*m_rhs);
      m_rhs->Multiply(1.0,scaling,tmp,0.0);
  }

#if USECPP11==1
  Epetra_LinearProblem problem( m_matrix.get(),
                                m_solution.get(),
                                m_rhs.get() );
#else
  Epetra_LinearProblem problem( m_matrix,
                                m_solution,
                                m_rhs );

#endif

  // @annavarapusr1: Needed to use direct solver without changing it for everyone else
  if(m_numerics.m_useDirectSolver) // If Chandra's test problems, use direct solver
  {
    Amesos_BaseSolver* Solver;
    Amesos Factory;

    std::string SolverType = "Klu";

    if( !(m_lsp.compare("Slu")) )
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

    SetLinearSolverParameters( solver );

    ML_Epetra::MultiLevelPreconditioner* MLPrec = nullptr;
    //SetupMLPreconditioner( domain, MLPrec );

    if( m_numerics.m_useMLPrecond )
    {
      Teuchos::ParameterList MLList;
      ML_Epetra::SetDefaults("SA",MLList);
      //MLList.set("aggregation: type", "Uncoupled");
      //MLList.set("aggregation: type", "MIS");
      MLList.set("prec type", "MGW");
      MLList.set("smoother: type","ILU");
      MLList.set("ML output",1);
      MLList.set("PDE equations",3);

#if USECPP11==1
      MLPrec = new ML_Epetra::MultiLevelPreconditioner(*m_matrix.get(), MLList);
#else
      MLPrec = new ML_Epetra::MultiLevelPreconditioner(*m_matrix, MLList);
#endif
      solver.SetPrecOperator(MLPrec);
    }
    else // use ILUT preconditioner with domain decomp
    {
      solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
      solver.SetAztecOption(AZ_subdomain_solve,AZ_ilut);
      solver.SetAztecParam(AZ_ilut_fill,m_numerics.m_ilut_fill);
      solver.SetAztecParam(AZ_drop,m_numerics.m_ilut_drop);
    }

    solver.Iterate(m_numerics.m_maxIters,
                   m_numerics.krylov_tol);

    if(m_numerics.m_useMLPrecond)
    {
      delete MLPrec; MLPrec=nullptr;
    }
  }

  // copy vector solution into geos data structures

  realT scalingFactor = CheckSolution( local_solution, domain, 0 );
  PropagateSolution( local_solution, scalingFactor, domain, 0 );

  // re-sync ghost nodes

  partition.SynchronizeFields(m_syncedFields, CommRegistry::lagrangeSolver02);

  // copy vector solution into geos data structures

  PostSyncConsistency( domain, partition );
}


void SolverBase::SetLinearSolverParameters( AztecOO& solver )
{
    if(m_numerics.m_useBicgstab)
    {
      solver.SetAztecOption(AZ_solver,AZ_bicgstab);
    }
    else
    {
      solver.SetAztecOption(AZ_solver,AZ_gmres);
    }
    solver.SetAztecOption(AZ_conv,AZ_r0);
    solver.SetAztecOption(AZ_kspace,m_numerics.m_kspace);
    if ( m_numerics.m_verbose>=2 )
    {
      solver.SetAztecOption(AZ_output,AZ_all);
    }
    else
    {
      solver.SetAztecOption(AZ_output,AZ_none);
    }
}


		// PRINT NORMS
		// Utility to print matrix scaling

void print_norms(Epetra_System & epetraSystem, std::string nametag)
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   double matnorm[2][2];
   double rhsnorm[2];

   matnorm[0][0] = epetraSystem.m_matrix[0][0]->NormInf();
   matnorm[0][1] = epetraSystem.m_matrix[0][1]->NormInf();
   matnorm[1][0] = epetraSystem.m_matrix[1][0]->NormInf();
   matnorm[1][1] = epetraSystem.m_matrix[1][1]->NormInf();

   epetraSystem.m_rhs[0]->NormInf(&(rhsnorm[0]));
   epetraSystem.m_rhs[1]->NormInf(&(rhsnorm[1]));

   if( rank==0 )
   {
     printf("SolverBase :: Linear system inf-norms (%s)\n",nametag.c_str());
     printf("           ::   | %.1e %.1e | = | %.1e |\n",matnorm[0][0],matnorm[0][1],rhsnorm[0]);
     printf("           ::   | %.1e %.1e |   | %.1e |\n",matnorm[1][0],matnorm[1][1],rhsnorm[1]);
   }
}



// Utility to clear matrix rows
realT SolverBase::ClearRow ( Epetra_FECrsMatrix* const matrix,
                             const unsigned int row,
                             const realT factor )
{
//    Assert (matrix->Filled()==true, ExcMatrixNotCompressed());

                                // Only do this on the rows owned
                                // locally on this processor.
  long long int rowTmp = static_cast<long long int>(row);
  int local_row = matrix->LRID(rowTmp);
  realT LARGE = 0.0;
  
  if (local_row >= 0)
  {
    realT *values = NULL;
    int *col_indices = NULL;
    int num_entries;

    matrix->ExtractMyRowView( local_row, num_entries, values, col_indices);


    if( values!=NULL && col_indices!=NULL && num_entries>0 )
    {
      int* diag_find = std::find(col_indices,col_indices+num_entries-1, local_row);
      int diag_index = (int)(diag_find - col_indices);

      for (int j=0; j<num_entries; ++j)
      {
        if (diag_index != j )
        {
          values[j] = 0.;
        }
      }
      values[diag_index] *= factor;
      LARGE = values[diag_index];
    }
  }
  return (LARGE);
}




		// BLOCK SOLVER
		// Solve a 2x2 linear system using a
		// block-preconditioned Krylov method

void SolverBase::SolveBlock ( PhysicalDomainT&  domain,
                              SpatialPartition& partition,
                              Epetra_System& epetraSystem,
                              const realT time,
                              const realT dt )
{
		// SCHEME CHOICES
		//
		// there are several flags to control solver behavior.
		// these should be compared in a scaling study.
		//
		// 1. whether to use inner solvers or just the 
		//    sub-block preconditioners directly. false
		//    is probably better.
		// 2. whether to use a block diagonal or a full
		//    block triangular preconditioner.  false is
		//    probably better.
		// 3. whether to perform an explicit scaling
		//    of the linear system before solving.  note
		//    that the matrix and rhs are modified in place
		//    by this operation.  true is probably better.
		// 4. whether to use BiCGstab or GMRES for the 
		//    krylov solver.  GMRES is generally more robust,
		//    BiCGstab sometimes shows better parallel performance.
		//    false is probably better.

  const bool use_inner_solver  = m_numerics.m_useInnerSolver;
  const int use_scaling        = m_numerics.m_scalingOption;  // no longer just row
  const bool use_bicgstab      = m_numerics.m_useBicgstab;
  const bool use_diagonal_prec = false;


		// DEBUGGING
		// Write out unscaled linear system to matlab

		// TODO: Josh: I noticed we seem to be storing a lot of
		// zero-valued entries in our sparsity pattern.  We should
		// follow up on this to make sure we are not over-allocating
		// space in our matrices.
  /* 
  {
    EpetraExt::RowMatrixToMatlabFile("umatrix00.dat",*epetraSystem.m_matrix[0][0]);
    EpetraExt::RowMatrixToMatlabFile("umatrix01.dat",*epetraSystem.m_matrix[0][1]);
    EpetraExt::RowMatrixToMatlabFile("umatrix10.dat",*epetraSystem.m_matrix[1][0]);
    EpetraExt::RowMatrixToMatlabFile("umatrix11.dat",*epetraSystem.m_matrix[1][1]);
    EpetraExt::MultiVectorToMatlabFile("urhs0.dat",*epetraSystem.m_rhs[0]);
    EpetraExt::MultiVectorToMatlabFile("urhs1.dat",*epetraSystem.m_rhs[1]);
  }
  */


		// ROW & COLUMN SCALING
		//
		// Scale the linear system with row and column scaling
		// matrices R and C.  The resulting linear system is
		// 	(R.A.C).(Cinv.x) = R.b
		// We use the iterative method of Ruiz (2001) to 
		// repeatedly update R and C until the desired scaling
		// is found. Note also that C must be saved to later 
		// compute the true solution from the temporary solution
		// 	x = C.x' where x' = Cinv.x

		// The diagonal scaling matrices are stored as four
		// vectors, one for each combination of row/column and 
		// block 0/block 1.  We store them in a 2x2 array as
		// [ R0 C0 ; 
		//   R1 C1 ]

		// note that we can extend this methodology to larger
		// block systems by storing a (n_blocks x 2) array:
		// [ R0 C0 ;
		//   R1 C1 ;
		//   .. ..
		//   Rn Cn ]

  const unsigned n_blocks = 2;		       // algorithm *should* work for any block size n 
  enum {ROW,COL};			       // indexing to improve readability (ROW=0,COL=1)

  RCP<Epetra_Vector> scaling   [n_blocks][2];  // complete scaling
  RCP<Epetra_Vector> scaling_k [n_blocks][2];  // scaling at iteration k

  if(use_scaling == 2)
  {
		// first print unscaled norms

    if(m_numerics.m_verbose >= 2)
    {
      print_norms(epetraSystem,"unscaled");
    }

		// allocate storage for our scaling vectors, and initialize
		// them to identity scalings (R=C=I). 

    for(unsigned b=0; b<n_blocks; ++b)
    {
      scaling[b][ROW] = rcp(new Epetra_Vector(epetraSystem.m_matrix[b][b]->RangeMap()));
      scaling[b][COL] = rcp(new Epetra_Vector(epetraSystem.m_matrix[b][b]->DomainMap()));

      scaling[b][ROW]->PutScalar(1.0);
      scaling[b][COL]->PutScalar(1.0);

      scaling_k[b][ROW] = rcp(new Epetra_Vector(epetraSystem.m_matrix[b][b]->RangeMap()));
      scaling_k[b][COL] = rcp(new Epetra_Vector(epetraSystem.m_matrix[b][b]->DomainMap()));
    }

		// begin scaling iterations 

    for(unsigned k=0; k<30; ++k)
    {
		// get row and column max norms for scaling

      for(unsigned a=0; a<n_blocks; ++a)
      {

        scaling_k[a][ROW]->PutScalar(0.0); // clear
        scaling_k[a][COL]->PutScalar(0.0); // clear

        Epetra_Vector tmp_row(epetraSystem.m_matrix[a][a]->RangeMap());
        Epetra_Vector tmp_col(epetraSystem.m_matrix[a][a]->DomainMap());

        for(unsigned b=0; b<n_blocks; ++b)
        {
          epetraSystem.m_matrix[a][b]->InvRowMaxs(tmp_row); // 1/row_norms for block 
          epetraSystem.m_matrix[b][a]->InvColMaxs(tmp_col); // 1/col_norms for block

          tmp_row.Reciprocal(tmp_row); // row_norms for block
          tmp_col.Reciprocal(tmp_col); // col_norms for block

          scaling_k[a][ROW]->Update(1.0,tmp_row,1.0);  // add across blocks (A and B) or (C and D)
          scaling_k[a][COL]->Update(1.0,tmp_col,1.0);  // add across blocks (A and C) or (B and D)
		
		// note this last step defines a weird norm, i.e. the sum inf_norm(A)+inf_norm(B)
		// rather than inf_norm([A B]).  the first is just easier to compute using
		// built in operations.  this should not make much of a difference in terms
		// of actual performance, as we're just trying to get a reasonable scaling.
        }

        for(int i=0; i<scaling_k[a][ROW]->MyLength(); ++i)
           (*scaling_k[a][ROW])[i] = 1./sqrt((*scaling_k[a][ROW])[i]);  // use 1/sqrt(norm) for scaling
        for(int i=0; i<scaling_k[a][COL]->MyLength(); ++i)
           (*scaling_k[a][COL])[i] = 1./sqrt((*scaling_k[a][COL])[i]);  // use 1/sqrt(norm) for scaling

        scaling[a][ROW]->Multiply(1.0,*scaling[a][ROW],*scaling_k[a][ROW],0.0); // save total row scaling over all iterations
        scaling[a][COL]->Multiply(1.0,*scaling[a][COL],*scaling_k[a][COL],0.0); // save total col scaling over all iterations
      }

		// actually scale matrix A(k) = R(k).A(k-1).C(k)
		// also scale rhs b(k) = R(k)*b(k-1)
		// will scale solution x = C*x' after solve

      for(unsigned a=0; a<n_blocks; ++a)
      {
        for(unsigned b=0; b<n_blocks; ++b)
        {
          epetraSystem.m_matrix[a][b]->LeftScale(*scaling_k[a][ROW]);
          epetraSystem.m_matrix[a][b]->RightScale(*scaling_k[b][COL]);
        }
        epetraSystem.m_rhs[a]->Multiply(1.0,*scaling_k[a][ROW],*epetraSystem.m_rhs[a],0.0); 
      }

		// check for convergence in desired row and column norms
		// and print info in verbose mode > 0

      double convergence = 0.0;
      double norm_threshold = 0.2;  

      for(unsigned a=0; a<n_blocks; ++a)
      for(unsigned b=0; b<2; ++b)
      {
         double tmp[1]; 
         scaling_k[a][b]->Reciprocal(*scaling_k[a][b]);
         scaling_k[a][b]->NormInf(&(tmp[0]));
         tmp[0] = abs(1-pow(tmp[0],2));
         convergence = std::max(convergence,tmp[0]);
      }

      if( partition.m_rank == 0 && m_numerics.m_verbose >= 2 )
      {
        if(k==0)
        {
          printf("SolverBase :: Re-scaling matrix \n");
          printf("           ::   %d ... %.1e\n",k,convergence);
        }
        else
          printf("           ::   %d ... %.1e\n",k,convergence);
      }

      if(convergence < norm_threshold && k > 1) break; 
    }

    if(m_numerics.m_verbose >= 2)
    {
      print_norms(epetraSystem,"scaled");
    }
  } // end scaling
  else if( use_scaling==1 )
  {

    // perform an explicit row scaling of the linear system,
    // R*A*x = R*b, where R is a diagonal scaling matrix.
    // we will use inverse row sums for the scaling.

    for(unsigned b=0; b<2; ++b)
    {
      Epetra_Vector scale_one(epetraSystem.m_matrix[b][b]->RowMap());
      Epetra_Vector scale_two(epetraSystem.m_matrix[b][b]->RowMap());

      Epetra_Vector scale_one_inv(epetraSystem.m_matrix[b][b]->RowMap());
      Epetra_Vector scale_two_inv(epetraSystem.m_matrix[b][b]->RowMap());

      epetraSystem.m_matrix[b][0]->InvRowSums(scale_one_inv);
      epetraSystem.m_matrix[b][1]->InvRowSums(scale_two_inv);
      scale_one.Reciprocal(scale_one_inv);
      scale_two.Reciprocal(scale_two_inv);  // not ideal, could choke if 1/0 or 1/NaN appears
      scale_one.Update(1.0,scale_two,1.0);
      scale_one_inv.Reciprocal(scale_one);

      for(unsigned c=0; c<2; ++c)
      {
        epetraSystem.m_matrix[b][c]->LeftScale(scale_one_inv);
      }

      Epetra_MultiVector tmp (*epetraSystem.m_rhs[b]);
      epetraSystem.m_rhs[b]->Multiply(1.0,scale_one_inv,tmp,0.0);
    }
  }


		// set initial guess to zero.  this is not strictly
		// necessary but is good for comparing solver performance.

  epetraSystem.m_solution[0]->PutScalar(0.0);
  epetraSystem.m_solution[1]->PutScalar(0.0);

		// The standard AMG aggregation strategy based on
		// the system matrix A can struggle when using
		// grids with large element aspect ratios.  To fix 
		// this, we can instead build the AMG hierarchy 
		// using an alternative matrix L built using information
		// about nodal positions.  Once the aggregates are 
		// determined, A is then used to construct the actual
		// coarse / fine scale operations. For more details
		// see: ML USER GUIDE V5, sec. 6.4.12, p. 33

		// Here, we simply extract three arrays of nodal
		// positions for the locally owned nodes, for later use.  
		// For vector-valued problems (with multiple dofs per node)
		// ML is going to assume degrees of freedom are ordered as
		// [u_x_0, u_y_0, u_z_0, u_x_1, u_y_1, u_z_1, ... ]
		// where dof components are grouped "node-wise."

#define AGGREGATION 0
#if     AGGREGATION==1

  Array1dT<double> x_coord;
  Array1dT<double> y_coord; // set to to 0 for 1D problems
  Array1dT<double> z_coord; // set to to 0 for 2D problems

  if(m_numerics.m_useMLPrecond)
  {
    const iArray1d & is_ghost = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
    //iArray1d const & trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
    const Array1dT<R1Tensor> & X = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition>();
    x_coord.resize(domain.m_feNodeManager.m_numNodes);
    y_coord.resize(domain.m_feNodeManager.m_numNodes);
    z_coord.resize(domain.m_feNodeManager.m_numNodes);
    localIndex b=0;
    for( auto a=0u ; a<domain.m_feNodeManager.m_numNodes ; ++a )
    {
      if(is_ghost[a] < 0)
      {
        realT const * const X_ref = X[a].Data();

        x_coord[b] = X_ref[0];// + 0.1*((double) rand() / (RAND_MAX));
        y_coord[b] = X_ref[1];// + 0.1*((double) rand() / (RAND_MAX));
        z_coord[b] = X_ref[2];// + 0.1*((double) rand() / (RAND_MAX));
        ++b;
      }
    }
  }

#endif

		// we want to use thyra to wrap epetra operators and vectors
		// for individual blocks.  this is an ugly conversion, but 
		// it is basically just window dressing.
		//
		// note the use of Teuchos::RCP reference counted pointers.
		// The general syntax is usually one of:
		//
		//   RCP<T> Tptr = rcp(new T) 
		//   RCP<T> Tptr = nonMemberConstructor();
		//   RCP<T> Tptr (t_ptr,false)
		//
		// where "false" implies the RCP does not own the object and
		// should not attempt to delete it when finished.


  RCP<const Thyra::LinearOpBase<double> >  matrix_block[2][2];
  RCP<Thyra::MultiVectorBase<double> >     lhs_block[2];
  RCP<Thyra::MultiVectorBase<double> >     rhs_block[2];


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  for(unsigned i=0; i<2; ++i)
  for(unsigned j=0; j<2; ++j)
  {
    RCP<Epetra_Operator> mmm (&*epetraSystem.m_matrix[i][j],false);
    matrix_block[i][j] = Thyra::epetraLinearOp(mmm);
  }
  

  for(unsigned i=0; i<2; ++i)
  { 
    RCP<Epetra_MultiVector> lll (&*epetraSystem.m_solution[i],false);
    RCP<Epetra_MultiVector> rrr (&*epetraSystem.m_rhs[i],false);

    lhs_block[i] = Thyra::create_MultiVector(lll,matrix_block[i][i]->domain());
    rhs_block[i] = Thyra::create_MultiVector(rrr,matrix_block[i][i]->range());
  }

		// now use thyra to create an operator representing
		// the full block 2x2 system

  RCP<const Thyra::LinearOpBase<double> > matrix = Thyra::block2x2(matrix_block[0][0],
                                                                   matrix_block[0][1],
                                                                   matrix_block[1][0],
                                                                   matrix_block[1][1]);

		// creating a representation of the blocked
		// rhs is a little uglier. (todo: check if there is
		// a cleaner way to do this.)

  RCP<Thyra::ProductMultiVectorBase<double> > rhs;
  {
    Teuchos::Array<RCP<Thyra::MultiVectorBase<double> > > mva;
    Teuchos::Array<RCP<const Thyra::VectorSpaceBase<double> > > mvs;

    for(unsigned i=0; i<2; ++i)
    {
      mva.push_back(rhs_block[i]);
      mvs.push_back(rhs_block[i]->range());
    }

    RCP<const Thyra::DefaultProductVectorSpace<double> > vs = Thyra::productVectorSpace<double>(mvs);

    rhs = Thyra::defaultProductMultiVector<double>(vs,mva);
  }

		// do the identical operation for the lhs

  RCP<Thyra::ProductMultiVectorBase<double> > lhs;

  {
    Teuchos::Array<RCP<Thyra::MultiVectorBase<double> > > mva;
    Teuchos::Array<RCP<const Thyra::VectorSpaceBase<double> > > mvs;

    for(unsigned i=0; i<2; ++i)
    {
      mva.push_back(lhs_block[i]);
      mvs.push_back(lhs_block[i]->range());
    }

    RCP<const Thyra::DefaultProductVectorSpace<double> > vs = Thyra::productVectorSpace<double>(mvs);

    lhs = Thyra::defaultProductMultiVector<double>(vs,mva);
  }


		// for the preconditioner, we need two approximate inverses,
		// one for the (0,0) block and one for the approximate
		// schur complement.  for now, we will use the (1,1) block
		// as our schur complement approximation, though we should
		// explore better approaches later. 

		// we store both "sub operators" in a 1x2 array:

  RCP<const Thyra::LinearOpBase<double> > sub_op[2];

		// each implicit "inverse" is based on an inner krylov solver,
		// with their own sub-preconditioners.  this leads to a very
		// accurate approximation of the inverse operator, but can be
		// overly expensive.  the other option is to ditch the inner
		// krylov solver, and just use the sub-preconditioners directly.

		// the implicit inverse for each diagonal block is built in 
		// three steps
		//   1.  define solver parameters
		//   2.  build a solver factory
		//   3.  build the inner solver operator
	

  for(unsigned i=0; i<2; ++i) // loop over diagonal blocks
  {
    RCP<Teuchos::ParameterList> list = rcp(new Teuchos::ParameterList("solver_list"),true);

      list->set("Linear Solver Type","AztecOO");
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Max Iterations",m_numerics.m_maxIters);
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Tolerance",1e-1*m_numerics.krylov_tol);
      if(use_bicgstab)
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","BiCGStab");
      else 
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","GMRES");
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",0);//int(m_numerics.m_verbose));

      if(m_numerics.m_useMLPrecond && i==0 )
      {
        if( m_numerics.m_verbose>=2 && partition.m_rank==0 )
        {
          std::cout<< "SolverBase :: Using ML preconditioner for block " << i << i <<std::endl;
        }

        list->set("Preconditioner Type","ML");
        list->sublist("Preconditioner Types").sublist("ML").set("Base Method Defaults","SA");
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("PDE equations",(i==0?3:1));
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("smoother: type","block Gauss-Seidel");
        //list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("smoother: type","Gauss-Seidel");
        //list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("smoother: type","Chebyshev");
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("ML output", 0);
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("aggregation: type","Uncoupled");
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("smoother: sweeps",3);

#if AGGREGATION==1
          list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("x-coordinates",x_coord.data());
          list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("y-coordinates",y_coord.data());
          list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("z-coordinates",z_coord.data());
          list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("null space: type",(i==0?"elasticity from coordinates":"default vectors"));
#endif

      }
      else
      {
        if( m_numerics.m_verbose>=2 && partition.m_rank==0 )
        {
          std::cout<< "SolverBase :: Using ILU preconditioner for block " << i << i <<std::endl;
        }

        list->set("Preconditioner Type","Ifpack");
        list->sublist("Preconditioner Types").sublist("Ifpack").set("Prec Type","ILU");
      }

    Stratimikos::DefaultLinearSolverBuilder builder;

      builder.setParameterList(list);

    if(use_inner_solver)
    {
      RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > strategy = createLinearSolveStrategy(builder);

      //if(i==0)
        sub_op[i] = Thyra::inverse(*strategy,matrix_block[i][i]);
      //else
      //{
      //  RCP<const Thyra::LinearOpBase<double> > BAinvBt = Thyra::multiply(matrix_block[0][1],sub_op[0],matrix_block[1][0]);
      //  RCP<const Thyra::LinearOpBase<double> > schur = Thyra::add(matrix_block[1][1],Thyra::scale(-1.0,BAinvBt));
      //  sub_op[i] = Thyra::inverse(*strategy,schur);
      //}
    }
    else
    {
      RCP<const Thyra::PreconditionerFactoryBase<double> > strategy = createPreconditioningStrategy(builder);
      RCP<Thyra::PreconditionerBase<double> > tmp;

      //if(i==0) 
        tmp = prec(*strategy,matrix_block[i][i]);
      //else     
      //  tmp = prec(*strategy,SchurEstimate);
 
     sub_op[i] = tmp->getUnspecifiedPrecOp();
    }
  }
 

		// create zero operators for off diagonal blocks

  RCP<const Thyra::LinearOpBase<double> > zero_01
    = rcp(new Thyra::DefaultZeroLinearOp<double>(matrix_block[0][0]->range(),
                                                 matrix_block[1][1]->domain())); 

  RCP<const Thyra::LinearOpBase<double> > zero_10
    = rcp(new Thyra::DefaultZeroLinearOp<double>(matrix_block[1][1]->range(),
                                                 matrix_block[0][0]->domain())); 

		// now build the block preconditioner

  RCP<const Thyra::LinearOpBase<double> > preconditioner;

  if(use_diagonal_prec)
  {
    preconditioner = Thyra::block2x2(sub_op[0],zero_01,zero_10,sub_op[1]);
  }
  else
  {
    RCP<const Thyra::LinearOpBase<double> > eye_00
      = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(matrix_block[0][0]->range())); 
 
    RCP<const Thyra::LinearOpBase<double> > eye_11
      = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(matrix_block[1][1]->range())); 

    RCP<const Thyra::LinearOpBase<double> > mAinvB1, mB2Ainv; 
    
    mAinvB1 = Thyra::scale(-1.0, Thyra::multiply(sub_op[0],matrix_block[0][1]) );
    mB2Ainv = Thyra::scale(-1.0, Thyra::multiply(matrix_block[1][0],sub_op[0]) );

    RCP<const Thyra::LinearOpBase<double> > Linv,Dinv,Uinv,Eye;
    
    //Eye = Thyra::block2x2(eye_00,zero_01,zero_10,eye_11);
    //Linv = Thyra::block2x2(eye_00,zero_01,mB2Ainv,eye_11);
    Dinv = Thyra::block2x2(sub_op[0],zero_01,zero_10,sub_op[1]);
    Uinv = Thyra::block2x2(eye_00,mAinvB1,zero_10,eye_11);
    
    //preconditioner = Eye;
    //preconditioner = Dinv;
    preconditioner = Thyra::multiply(Uinv,Dinv);
    //preconditioner = Thyra::multiply(Dinv,Linv);
    //preconditioner = Thyra::multiply(Uinv,Dinv,Linv);
  }


		// define solver strategy for blocked system. this is
		// similar but slightly different from the sub operator
		// construction, since now we have a user defined preconditioner

  {
    RCP<Teuchos::ParameterList> list = rcp(new Teuchos::ParameterList("list"),true);

      list->set("Linear Solver Type","AztecOO");
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Max Iterations",m_numerics.m_maxIters);
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Tolerance",m_numerics.krylov_tol);
      if(use_bicgstab)
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","BiCGStab");
      else 
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","GMRES");

      if( m_numerics.m_verbose>=3 )
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",1);
      else
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",0);

      list->set("Preconditioner Type","None"); // will use user-defined P

    Stratimikos::DefaultLinearSolverBuilder builder;

      builder.setParameterList(list);

    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > strategy = createLinearSolveStrategy(builder);

    RCP<Thyra::LinearOpWithSolveBase<double> > solver = strategy->createOp();

    Thyra::initializePreconditionedOp<double>(*strategy,
                                               matrix,
                                               Thyra::rightPrec<double>(preconditioner),
                                               solver.ptr());


    		// JAW: check "true" residual before solve.
    		//      should remove after debugging because this is potentially slow
    		//      and should just use iterative residual

    RCP<Thyra::VectorBase<double> > Ax = Thyra::createMember(matrix->range());
    RCP<Thyra::VectorBase<double> > r  = Thyra::createMember(matrix->range());
    {      
      Thyra::apply(*matrix, Thyra::NOTRANS,*lhs,Ax.ptr());
      Thyra::V_VmV<double>(r.ptr(),*rhs,*Ax);
      m_numerics.m_KrylovResidualInit = Thyra::norm(*r);
    }

		// !!!! Actual Solve !!!!

    Thyra::SolveStatus<double> status = solver->solve(Thyra::NOTRANS,*rhs,lhs.ptr()); 
    m_numerics.m_numKrylovIter = status.extraParameters->get<int>("Iteration Count");

    		// JAW: check "true" residual after
    		//      should remove after debugging because this is potentially slow

    {      
      Thyra::apply(*matrix, Thyra::NOTRANS,*lhs,Ax.ptr());
      Thyra::V_VmV<double>(r.ptr(),*rhs,*Ax);
      m_numerics.m_KrylovResidualFinal = Thyra::norm(*r);
    }

		// write a solver profile file

    if( partition.m_rank == 0 && m_numerics.m_verbose>=3 )
    {
      FILE* fp = fopen("solver_profile.txt","a");
      fprintf(fp,"%d %.9e %.9e\n", m_numerics.m_numKrylovIter, m_numerics.m_KrylovResidualInit, m_numerics.m_KrylovResidualFinal);
      fclose(fp);
    }
		
		// apply column scaling C to get true solution x from x' = Cinv*x

    if(use_scaling==2)
    {
      for(unsigned b=0; b<n_blocks; ++b)
        epetraSystem.m_solution[b]->Multiply(1.0,*scaling[b][COL],*epetraSystem.m_solution[b],0.0); 
    }
  }

  		// copy vector solution into geos data structures
		// re-sync ghost nodes

  realT scalingFactor = CheckSolutionBlock( domain );
  PropagateSolutionBlock( scalingFactor, domain );


  partition.SynchronizeFields(m_syncedFields, CommRegistry::lagrangeSolver02);

  PostSyncConsistency( domain, partition );

		// put 00 matrix back to unscaled form

  if(use_scaling==2)
  {
    scaling[0][ROW]->Reciprocal(*scaling[0][ROW]);
    scaling[0][COL]->Reciprocal(*scaling[0][COL]);

    epetraSystem.m_matrix[0][0]->LeftScale(*scaling[0][ROW]);
    epetraSystem.m_matrix[0][0]->RightScale(*scaling[0][COL]);
  }

  //exit(9);
}





/* ................ SCRATCH CODE .....................................

		// make identity and zero operators

  RCP<const Thyra::LinearOpBase<double> > Eye_11
    = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(A_11->range())); 

  RCP<const Thyra::LinearOpBase<double> > Zero_01
    = Teuchos::rcp(new Thyra::DefaultZeroLinearOp<double>(A_00->range(),A_11->domain())); 

  RCP<const Thyra::LinearOpBase<double> > Zero_10
    = Teuchos::rcp(new Thyra::DefaultZeroLinearOp<double>(A_11->range(),A_00->domain())); 

*/

/*
  EpetraExt::RowMatrixToMatlabFile("matrix00.dat",*epetraSystem.m_matrix[0][0]);
  EpetraExt::RowMatrixToMatlabFile("matrix11.dat",*epetraSystem.m_matrix[1][1]);


		// test diagonal blocks

  for(unsigned b=0; b<1; ++b)
  {
    printf("Checking singularity of (%d,%d) block ...\n",b,b);

    Epetra_LinearProblem problem( &*epetraSystem.m_matrix[b][b],
                                  &*epetraSystem.m_solution[b],
                                  &*epetraSystem.m_rhs[b] );


    epetraSystem.m_rhs[b]->PutScalar(1.0);
    epetraSystem.m_solution[b]->PutScalar(0.0);

    AztecOO solver(problem);

      solver.SetAztecOption(AZ_solver,AZ_gmres);
      solver.SetAztecOption(AZ_conv,AZ_r0);
      solver.SetAztecOption(AZ_kspace,m_numerics.m_kspace);
      solver.SetAztecOption(AZ_output,AZ_all);

    SetLinearSolverParameters( solver );
*/
/*
    ML_Epetra::MultiLevelPreconditioner* MLPrec = nullptr;

    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("smoother: type","ILU");
    MLList.set("ML output",0);
    MLList.set("PDE equations",1);//(b==1?1:3));

    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*epetraSystem.m_matrix[b][b], MLList);
    solver.SetPrecOperator(MLPrec);
*/
/*
    solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
    solver.SetAztecParam(AZ_ilut_fill,m_numerics.m_ilut_fill);

    solver.Iterate(100,1e-6);

    std::ofstream sol_output;
    char solNamei[100] = "solution.vec";
    sol_output.open( solNamei );
    epetraSystem.m_solution[0]->Print( sol_output );
    sol_output.close();

//    delete MLPrec;
  }
  return;
*/


/*
		// comms for ML

  ML_Comm * ml_comm;
  ML_Comm_Create(&ml_comm);

		// create ML block operator

  ML_Operator * matrix = ML_Operator_Create(ml_comm);

  ML_Operator_BlkMatInit(matrix, ml_comm, 2, 2, ML_DESTROY_EVERYTHING);

  EpetraExt::CrsMatrix_SolverMap map;

  for(unsigned i=0; i<2; ++i)
  for(unsigned j=0; j<2; ++j)
  {
    Epetra_RowMatrix * epetra_block = ML_Epetra::ModifyEpetraMatrixColMap(*epetraSystem.m_matrix[i][j],map);
    ML_Operator * ml_block = ML_Operator_Create(ml_comm);
    ML_Operator_WrapEpetraMatrix(epetra_block,ml_block);
    ML_Operator_BlkMatInsert(matrix,ml_block,i,j);
  }

  ML_Operator_BlkMatFinalize(matrix);

		// extract diagonal blocks
		// approximate S with A_11 for now

  ML_Operator * matrix_00 = ML_Operator_BlkMatExtract(matrix,0,0);
  ML_Operator * matrix_11 = ML_Operator_BlkMatExtract(matrix,1,1);

  Epetra_RowMatrix* A[2];
  A[0] = dynamic_cast<Epetra_RowMatrix*>((Epetra_CrsMatrix *) matrix_00->data);
  A[1] = dynamic_cast<Epetra_RowMatrix*>((Epetra_CrsMatrix *) matrix_11->data);

		// set up ML parameter lists

  Teuchos::ParameterList TopList;
  Teuchos::ParameterList SubList[2];

  ML_Epetra::SetDefaults("SA",TopList);
  ML_Epetra::SetDefaults("SA",SubList[0]);
  ML_Epetra::SetDefaults("SA",SubList[1]);

  TopList.set("aggregation: type", "Uncoupled");
  TopList.set("smoother: type","ILU");
  
  SubList[0].set("PDE equations",3);
  SubList[0].set("aggregation: type", "Uncoupled");
  SubList[0].set("smoother: type","ILU");

  SubList[1].set("PDE equations",1);
  SubList[1].set("aggregation: type", "Uncoupled");
  SubList[1].set("smoother: type","ILU");

		// create composite AMG preconditioner

  ML_Epetra::MultiLevelPreconditioner * ml_pre = 
      new ML_Epetra::MultiLevelPreconditioner(matrix, TopList);
      //new ML_Epetra::MultiLevelPreconditioner(matrix, TopList, A, SubList, 2);

		// convert ML matrix to Epetra Block Matrix

  ML_Epetra::RowMatrix e_matrix(matrix,epetra_comm);


  RCP<Epetra_Vector> RHS = Teuchos::rcp(new Epetra_Vector(e_matrix.OperatorRangeMap()));
  RHS->PutScalar(7.0);
  RCP<Epetra_Vector> srcRHS = Teuchos::rcp(new Epetra_Vector(*RHS));
  e_matrix.Apply(*srcRHS,*RHS);

  // set initial guess 
  Epetra_Vector LHS(e_matrix.OperatorDomainMap()); 
  LHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(&e_matrix, &LHS, &*RHS);
  AztecOO solver(Problem);
  
  solver.SetPrecOperator(ml_pre);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 1);
  solver.Iterate(100, 1e-8);

  delete ml_pre;
  ML_Operator_Destroy(&matrix);
  ML_Comm_Destroy(&ml_comm);
*/

//const Epetra_RowMatrix & rm = ml_precOp->RowMatrix();


		// ML preconditioners for A00 and A11
/*
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("smoother: type","ILU");
    MLList.set("ML output",0);
    MLList.set("PDE equations",3);

    const RCP<Epetra_Operator> tmp 
      = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*epetraSystem.m_matrix[0][0], MLList));

    const RCP<const Thyra::PreconditionerBase<double> > ml_A00
      = Thyra::epetraLinearOp(tmp);
*/

/*
  const RCP<Thyra::MLPreconditionerFactory> ml_A00_factory
    = Teuchos::rcp(new Thyra::MLPreconditionerFactory());

  const RCP<Thyra::MLPreconditionerFactory> ml_A11_factory
    = Teuchos::rcp(new Thyra::MLPreconditionerFactory());


  const RCP<Teuchos::ParameterList> ml_A00_list 
    = Teuchos::rcp(new Teuchos::ParameterList("ml_A00_list"),true);

  const RCP<Teuchos::ParameterList> ml_A11_list 
    = Teuchos::rcp(new Teuchos::ParameterList("ml_A11_list"),true);


  ML_Epetra::SetDefaults("SA",*ml_A00_list);
    ml_A00_list->set("aggregation: type", "Uncoupled");
    ml_A00_list->set("smoother: type","ILU");
    ml_A00_list->set("ML output",0);
    ml_A00_list->set("PDE equations",3);

  ML_Epetra::SetDefaults("SA",*ml_A11_list);
    ml_A11_list->set("aggregation: type", "Uncoupled");
    ml_A11_list->set("smoother: type","ILU");
    ml_A11_list->set("ML output",0);
    ml_A11_list->set("PDE equations",1);

 
  ml_A00_factory->setParameterList(ml_A00_list);
  ml_A11_factory->setParameterList(ml_A11_list);



  RCP<Thyra::PreconditionerBase<double> > ml_A00 
    = ml_A00_factory->createPrec();

  Thyra::initializePrec(A_00,&*ml_A00,Thyra::ESupportSolveUse()); 
 */

/*
		// old way

  const RCP<Teuchos::ParameterList> inv_A_00_list 
    = Teuchos::rcp(new Teuchos::ParameterList("inv_A_00_list"),true);


    inv_A_00_list->sublist("Forward Solve").set("Max Iterations",100);
    inv_A_00_list->sublist("Forward Solve").set("Tolerance",1e-10);
    inv_A_00_list->sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver", "GMRES");
    inv_A_00_list->sublist("Forward Solve").sublist("AztecOO Settings").set("Size of Krylov Subspace",100);
    inv_A_00_list->sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency", 0);
    inv_A_00_list->sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Preconditioner", "none");

  const RCP<Thyra::LinearOpWithSolveFactoryBase<double> > inv_A_00_factory
    = Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory());

    inv_A_00_factory->setParameterList(inv_A_00_list);

  const RCP<const Thyra::LinearOpBase<double> > inv_A_00
    = Thyra::inverse(*inv_A_00_factory,A_00);

		// do same for A_11 block

  const RCP<const Thyra::LinearOpBase<double> > inv_A_11
    = Thyra::inverse(*inv_A_00_factory,A_11);
*/
/*
#if 0

  Amesos_BaseSolver* Solver;
  Amesos Factory;

  std::string SolverType = "Klu";

//  if( !(m_lsp.compare("Slu")) )
//  {
//    SolverType = "Amesos_Superludist";
//  }

  Solver = Factory.Create(SolverType, problem);
//  if (Solver == 0) {
//      cerr << "Specified solver is not available" << endl;
//  }

//  Amesos_Superludist* Solver = new Amesos_Superludist(problem);

  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  Teuchos::ParameterList List;
  List.set("MaxProcs", Comm.NumProc());
  Solver->SetParameters(List);

  int ierr = Solver->SymbolicFactorization();
  if (ierr > 0)
    cerr << "ERROR!" << endl;

  ierr = Solver->NumericFactorization();
  if (ierr > 0)
    cerr << "ERROR!" << endl;

  Solver->Solve();


  if( Solver!=nullptr )
    delete Solver;
#else
*/

/*
#if 0

  //Amesos_BaseSolver* Solver;
  Amesos Factory;

  std::string SolverType = "Klu";

  if( !(m_lsp.compare("Slu")) )
  {
    SolverType = "Amesos_Superludist";
  }

//  Solver = Factory.Create(SolverType, problem);
//  if (Solver == 0) {
//      cerr << "Specified solver is not available" << endl;
//  }

  Amesos_Superludist* Solver = new Amesos_Superludist(problem);

  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  Teuchos::ParameterList List;
  List.set("MaxProcs", Comm.NumProc());
  Solver->SetParameters(List);

  int ierr = Solver->SymbolicFactorization();
  if (ierr > 0)
    cerr << "ERROR!" << endl;

  ierr = Solver->NumericFactorization();
  if (ierr > 0)
    cerr << "ERROR!" << endl;

  Solver->Solve();


  if( Solver!=nullptr )
    delete Solver;
#else
*/


		// from randy:
		// create diagonal schur complement approximation
		// S = C - B2 diag(A)inv B1

#if 0
  RCP<const Thyra::LinearOpBase<double> > SchurEstimate;
  RCP<const Thyra::LinearOpBase<double> >  invDiag00LinOp;
  Epetra_CrsMatrix SchurEstimate00( Copy, *(epetraSystem.m_rowMap[1]),0 );
  if(1)
  {
    /*
    Epetra_FECrsGraph sparsity_00_diag( Copy, *(epetraSystem.m_rowMap[0]), 1 );

    iArray1d nonEmptyDOF;
    for( int a=0; a<epetraSystem.m_matrix[0][0]->NumMyRows(); ++a )
    {
      nonEmptyDOF.push_back(a);
    }
    sparsity_00_diag.InsertIndices( nonEmptyDOF.size(),
                                    nonEmptyDOF.data(),
                                    nonEmptyDOF.size(),
                                    nonEmptyDOF.data() );
    */

    Epetra_CrsMatrix invDiag00( Copy, *(epetraSystem.m_rowMap[0]), 0 );



    Epetra_Vector diagVec( *(epetraSystem.m_rowMap[0]) );
    Epetra_Vector invDiagVec( *(epetraSystem.m_rowMap[0]) );
    epetraSystem.m_matrix[0][0]->ExtractDiagonalCopy( diagVec );
    invDiagVec.Reciprocal( diagVec );

//    std::cout<<diagVec<<std::endl;
//    std::cout<<invDiagVec<<std::endl;

    invDiag00.ReplaceDiagonalValues( invDiagVec );


    Epetra_CrsMatrix scratch( Copy, *(epetraSystem.m_rowMap[0]),0 );

    EpetraExt::MatrixMatrix::Multiply( invDiag00, false,
                                       *(epetraSystem.m_matrix[1][0]), false,
                                       scratch );

    EpetraExt::MatrixMatrix::Multiply( *(epetraSystem.m_matrix[0][1]), false,
                                       scratch, false,
                                       SchurEstimate00 );

    EpetraExt::MatrixMatrix::Add( *(epetraSystem.m_matrix[1][1]), false, 1.0,
                                  SchurEstimate00, -1.0 );

    SchurEstimate00.FillComplete();
    
    // JAW write for checking
    EpetraExt::RowMatrixToMatlabFile("schur_est.dat",SchurEstimate00);

//    invDiag00LinOp = Thyra::epetraLinearOp( RCP<Epetra_Operator>( &invDiag00, false ) );
    RCP<Epetra_Operator> mmm (&SchurEstimate00,false);
    SchurEstimate = Thyra::epetraLinearOp(mmm);

//    std::cout<<*(epetraSystem.m_matrix[1][1])<<std::endl;
//    std::cout<<SchurEstimate00<<std::endl;
  }
//  SchurEstimate = matrix_block[1][1];

//  SchurEstimate = Thyra::add(matrix_block[1][1],Thyra::scale(-1.0,Thyra::multiply(matrix_block[0][1],invDiag00LinOp,matrix_block[1][0])));
#endif




