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
#include "ImplicitLaplace.h"


/* Constructor */

template <int dim>
ImplicitLaplaceSolver<dim> :: ImplicitLaplaceSolver ( const std::string& name,
                                                      ProblemManagerT* const pm )
  :
  SolverBase(name,pm),
  epetra_comm(pm->m_epetraComm),
  this_mpi_process(epetra_comm.MyPID()),
  n_mpi_processes(epetra_comm.NumProc()),
  verbose(this_mpi_process==0)
{}


/* Read XML */

template <int dim>
void ImplicitLaplaceSolver<dim> :: ReadXML ( TICPP::HierarchicalDataNode* const hdn)
{
  SolverBase::ReadXML(hdn);

  equation_data.diffusion = hdn->GetAttributeOrDefault("diffusion_coeff",1.0);
  equation_data.source    = hdn->GetAttributeOrDefault("source_coeff",1.0);
  numerics.krylov_tol     = hdn->GetAttributeOrDefault("krylov_tol",1.0e-6);
}


/* Destructor */

template <int dim>
ImplicitLaplaceSolver<dim> :: ~ImplicitLaplaceSolver ()
{
  delete row_map;
  delete sparsity;
  delete matrix;
  delete solution;
  delete rhs;
}


/* Solver Name */

template <>
const char* ImplicitLaplaceSolver<2> :: SolverName()
{
  return "ImplicitLaplaceSolver2D";
}

template <>
const char* ImplicitLaplaceSolver<3> :: SolverName()
{
  return "ImplicitLaplaceSolver3D";
}


/* Initialization */

template <int dim>
void ImplicitLaplaceSolver<dim> :: Initialize (PhysicalDomainT * domain, SpatialPartition& partition )
{
  // ... can't do much without partition info ....

  // edit - now that partition has been added could potentially initialize
  // sparsity matrix here (rather than at each timestep)
  // would depend on boundary conditions (and constant DOF interactions).
}


/* Registration */

template <int dim>
void ImplicitLaplaceSolver<dim> :: RegisterFields (PhysicalDomainT * domain)
{
  domain->m_feNodeManager.AddKeyedDataField<FieldInfo::pressure>();
  domain->m_feNodeManager.AddKeylessDataField<int>("TrilinosIndex",true,true);
}


/* Communications */

template <int dim>
void ImplicitLaplaceSolver<dim> :: InitializeCommunications (PartitionBase& partition)
{
  syncedFields.clear();
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::pressure>::Name());
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("TrilinosIndex");
  partition.SetBufferSizes(syncedFields, CommRegistry::implicitLaplaceComm);
}


/* Timestepping */

template <int dim>
double ImplicitLaplaceSolver<dim> :: TimeStep (const realT& time,
                                               const realT& dt,
                                               const int cycleNumber,
                                               PhysicalDomainT * domain,
                                               const array<string>& namesOfSolverRegions,
                                               SpatialPartition& partition,
                                               FractunatorBase* const fractunator )
{
  realT dt_return = dt;

  if(verbose)
    std::cout << std::endl
              << " :: LAPLACE SOLVER | " << dim << "D" << std::endl
              << " :: No. mpi processes ... " << n_mpi_processes << std::endl;

  SetupSystem (domain, partition);
  Assemble    (domain, partition);
  Solve       (domain, partition);

  if(verbose)
    std::cout << std::endl;
  return dt_return;

}


/* Setup System */

template <int dim>
void ImplicitLaplaceSolver<dim> :: SetupSystem (PhysicalDomainT * domain,
                                                SpatialPartition& partition)
{
  // determine the global/local degree of freedom distribution.

  int n_ghost_rows  = domain->m_feNodeManager.GetNumGhosts();
  int n_local_rows  = domain->m_feNodeManager.DataLengths()-n_ghost_rows;
//  int n_hosted_rows = n_local_rows+n_ghost_rows;

  std::vector<int> gather(n_mpi_processes);

  epetra_comm.GatherAll(&n_local_rows,
                        &gather.front(),
                        1);

  int first_local_row = 0;
  int n_global_rows = 0;

  for(int p=0 ; p<n_mpi_processes ; ++p)
  {
    n_global_rows += gather[p];
    if(p<this_mpi_process)
      first_local_row += gather[p];
  }

  // create trilinos dof indexing

  array<integer>& trilinos_index = domain->m_feNodeManager.GetFieldData<int>("TrilinosIndex");
  array<integer>& is_ghost       = domain->m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();

  unsigned local_count = 0;
  for(unsigned r=0 ; r<trilinos_index.size() ; ++r)
  {
    if(is_ghost[r] < 0)
    {
      trilinos_index[r] = first_local_row+local_count;
      local_count++;
    }
    else
    {
      trilinos_index[r] = -INT_MAX;
    }
  }

  assert(static_cast<int>(local_count) == n_local_rows);

  partition.SynchronizeFields(syncedFields, CommRegistry::implicitLaplaceComm);

  // create epetra map

  row_map = new Epetra_Map(n_global_rows,n_local_rows,0,epetra_comm);

  // set up sparsity graph

  sparsity = new Epetra_FECrsGraph(Copy,*row_map,0);

  RegionMap::iterator
    region     = domain->m_feElementManager.m_ElementRegions.begin(),
    end_region = domain->m_feElementManager.m_ElementRegions.end();

  for( ; region != end_region ; ++region)
  {
    const unsigned dofs_per_element = region->second.m_numNodesPerElem;

    std::vector<int> element_index (dofs_per_element);

    for(localIndex element = 0 ; element < region->second.m_numElems ; ++element)
    {
      const localIndex* const local_index = region->second.m_toNodesRelation[element];

      for(unsigned i=0 ; i<dofs_per_element ; ++i)
      {
        const localIndex n = local_index[i];
        element_index[i] = trilinos_index[n];
      }

      sparsity->InsertGlobalIndices(dofs_per_element,
                                    &element_index.front(),
                                    dofs_per_element,
                                    &element_index.front());
    }
  }

  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();
}


/* Assemble */

template <int dim>
void ImplicitLaplaceSolver<dim> :: Assemble (PhysicalDomainT * domain,
                                             SpatialPartition& partition  )

{
  // (re-)init linear system

  matrix   = new Epetra_FECrsMatrix(Copy,*sparsity);
  solution = new Epetra_FEVector(*row_map);
  rhs      = new Epetra_FEVector(*row_map);

  // basic nodal data ( = dof data for our problem)

  array<integer>& trilinos_index = domain->m_feNodeManager.GetFieldData<int>("TrilinosIndex");
//  array<integer>& is_ghost       =
// domain->m_nodeManager.GetFieldData<FieldInfo::ghostRank>();
  const array<integer>& isExternal    = domain->m_feNodeManager.m_isExternal;

  // begin region loop

  RegionMap::iterator
    region     = domain->m_feElementManager.m_ElementRegions.begin(),
    end_region = domain->m_feElementManager.m_ElementRegions.end();

  for( ; region != end_region ; ++region)
  {
    // assume all regions are linear hexes

    LagrangeBasis<dim>    basis(1);
    GaussQuadrature<dim>  quadrature(2);
    FiniteElement<dim>    fe(basis,quadrature);

    const unsigned n_q_points = fe.n_quadrature_points();
    const unsigned dofs_per_element = fe.dofs_per_element();

    // following should hold for scalar node problems

    assert(dofs_per_element == region->second.m_numNodesPerElem);

    // space for element matrix and rhs

    std::vector<R1TensorT<3> > nodal_coord    (dofs_per_element);

    Epetra_IntSerialDenseVector  element_index   (dofs_per_element);
    Epetra_SerialDenseVector     element_rhs     (dofs_per_element);
    Epetra_SerialDenseMatrix     element_matrix  (dofs_per_element,dofs_per_element);

    // determine ghost elements

    array<integer>& elem_is_ghost = region->second.GetFieldData<FieldInfo::ghostRank>();

    // begin element loop, skipping ghost elements

    for(localIndex element = 0 ; element < region->second.m_numElems ; ++element)
      if(elem_is_ghost[element] < 0)
      {
        // get global row indices and element node coordinates

        // TODO: the nodal coordinate handling break dimension
        //       independence because m_refposition assumes 3D

        const localIndex* const local_index = region->second.m_toNodesRelation[element];

        for(unsigned i=0 ; i<dofs_per_element ; ++i)
        {
          const localIndex n = local_index[i];
          element_index[i] = trilinos_index[n];

          nodal_coord[i] = (*domain->m_feNodeManager.m_refposition)[n];
        }

        // now the actual integration loop

        fe.reinit(nodal_coord);

        element_rhs.Scale(0);
        element_matrix.Scale(0);

        for(unsigned q=0 ; q<n_q_points ; ++q)
        {
          for(unsigned i=0 ; i<dofs_per_element ; ++i)
          {
            element_rhs(i) += fe.JxW(q) *
                              equation_data.source *
                              fe.value(i,q);

            for(unsigned j=0 ; j<dofs_per_element ; ++j)
            {
              element_matrix(i,j) += fe.JxW(q) *
                                     equation_data.diffusion *
                                     (fe.gradient(i,q) * fe.gradient(j,q));
            }

          }
        }

        // modify local row if this is a boundary dof

        for(unsigned i=0 ; i<dofs_per_element ; ++i)
          if(isExternal[local_index[i]] > 0) //&& is_ghost[local_index[i]] < 0)
          {
            double row_sum = 0;
            for(unsigned j=0 ; j<dofs_per_element ; ++j)
            {
              row_sum += fabs(element_matrix(i,j));
              element_matrix(i,j) = 0;
            }

            element_matrix(i,i) = row_sum;
            element_rhs(i) = 0;
          }

        // assemble into global system

        matrix->SumIntoGlobalValues(element_index,
                                    element_matrix);

        rhs->SumIntoGlobalValues(element_index,
                                 element_rhs);

      }// end element

  } // end region

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  //rhs->Print(std::cout);
  //std::cout << matrix->NormInf() << std::endl;
  //EpetraExt::RowMatrixToMatlabFile("system-matrix.dat",*matrix);
}


/* Solve */

template <int dim>
void ImplicitLaplaceSolver<dim> :: Solve (PhysicalDomainT * domain,
                                          SpatialPartition& partition)
{
  // krylov solver

  Epetra_LinearProblem problem(&(*matrix),
                               &(*solution),
                               &(*rhs));

  // ml preconditioner

  //Teuchos::ParameterList amg_params;
  //int *options    = new int[AZ_OPTIONS_SIZE];
  //double *params  = new double[AZ_PARAMS_SIZE];
  //ML_Epetra::SetDefaults("SA",amg_params,options,params);

  // note, using asymmetric solver because ILU is asymmetric.
  // system is (mostly) symmetric though.  note that this is
  // a simple, but far from optimal, combination.  should use
  // CG + ML AMG here.

  AztecOO solver(problem);
  solver.SetAztecOption(AZ_solver,AZ_bicgstab);
  solver.SetAztecOption(AZ_precond,AZ_none);
//          solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
//          solver.SetAztecOption(AZ_subdomain_solve,AZ_ilut);
  solver.SetAztecOption(AZ_conv,AZ_rhs);
  //solver.SetAztecOption(AZ_output,AZ_none);

  solver.Iterate(1000,numerics.krylov_tol);
  //solution->Print(std::cout);

  // copy vector solution into geos data structures

  array<integer>& trilinos_index = domain->m_feNodeManager.GetFieldData<int>("TrilinosIndex");
  array<integer>& is_ghost       = domain->m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  array<real64>& geos_pressure  = domain->m_feNodeManager.GetFieldData<FieldInfo::pressure>();

  int dummy;
  double* local_solution = NULL;

  solution->ExtractView(&local_solution,&dummy);

  for(unsigned r=0 ; r<geos_pressure.size() ; ++r)
    if(is_ghost[r] < 0)
    {
      int lid = row_map->LID(trilinos_index[r]);
      geos_pressure[r] = local_solution[lid];
    }

  // re-sync ghost nodes

  partition.SynchronizeFields(syncedFields, CommRegistry::implicitLaplaceComm);
}


/* Explicit Instantiations */

//template class ImplicitLaplaceSolver<2>;
template class ImplicitLaplaceSolver<3>;


/* Register solver in the solver factory */

SolverRegistrator<ImplicitLaplaceSolver<3> >
reg_ImplicitLaplaceSolver3D;

//SolverRegistrator<ImplicitLaplaceSolver<2> >
//  reg_ImplicitLaplaceSolver2D;
