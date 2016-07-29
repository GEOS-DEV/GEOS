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
 * @file ImplicitMechanicsSolver.cpp
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#include "ImplicitMechanics.h"
#include "SolverFactory.h"
#include "PhysicsSolverStrings.h"

#include "Utilities/Utilities.h"
#include "Utilities/Kinematics.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "ElementLibrary/FiniteElement.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"

#include "ObjectManagers/TableManager.h"
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"

template <int dim>
ImplicitMechanicsSolver<dim>::ImplicitMechanicsSolver( const std::string& name,
                                                       ProblemManagerT* const pm ):
SolverBase(name,pm),
epetra_comm(pm->m_epetraComm),
this_mpi_process(epetra_comm.MyPID()),
n_mpi_processes(epetra_comm.NumProc()),
verbose(this_mpi_process==0),
m_useMLPrecond(true),
equation_data(),
numerics(),
row_map(),
sparsity(),
matrix(),
solution(),
rhs(),
m_cfl(1.0),
m_nonContactModulus(0.0),
m_recordIncrementalDisplacement(false),
syncedFields()
{
}

template <int dim>
ImplicitMechanicsSolver<dim>::~ImplicitMechanicsSolver()
{
  // TODO Auto-generated destructor stub
}


template <int dim>
void ImplicitMechanicsSolver<dim>::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML(hdn);

  numerics.krylov_tol     = hdn->GetAttributeOrDefault("tol",1.0e-6);
  
  m_nonContactModulus = hdn->GetAttributeOrDefault("NonContactModulus","2.2 GPa"); // fixme incorrect units
  m_recordIncrementalDisplacement = hdn->GetAttributeOrDefault("RecordIncrementalDisplacement",false); 
  m_doApertureUpdate = hdn->GetAttributeOrDefault("UpdateAperture",false); 
  
  m_useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);

  m_trilinosIndexStr = "IMS_" +  toString<int>(m_instances) + "_GlobalDof";  
  ++m_instances; 
}



template <int dim>
void ImplicitMechanicsSolver<dim>::RegisterFields( PhysicalDomainT& domain )
{

  // register nodal fields
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::displacement>();
  
  if(m_recordIncrementalDisplacement){
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();
  }
  
  domain.m_feNodeManager.AddKeylessDataField<int>("inContact",true,true);

  domain.m_feNodeManager.AddKeylessDataField<int>(m_trilinosIndexStr,true,true);

}


template <int dim>
void ImplicitMechanicsSolver<dim>::Initialize(PhysicalDomainT& domain, SpatialPartition& partition )
{
  using namespace BoundaryConditionFunctions;
  ApplyMultiSetBoundaryCondition<R1Tensor>(this, &ImplicitMechanicsSolver<dim>::NonpenetratingBC_NeighborUpdate, 
                                           domain, domain.m_feNodeManager,
                       NonPenetratingBoundaryCondition::BoundaryConditionName(), 0.0 );
}


template <int dim>
void ImplicitMechanicsSolver<dim>::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_trilinosIndexStr);
  partition.SetBufferSizes(syncedFields, CommRegistry::lagrangeSolver02);
}

template <int dim>
double ImplicitMechanicsSolver<dim>::TimeStep( const realT& time,
                                               const realT& dt ,
                                               const int cycleNumber,
                                               PhysicalDomainT& domain,
                                               const sArray1d& namesOfSolverRegions ,
                                               SpatialPartition& partition,
                                               FractunatorBase* const fractunator )
{

  realT dt_return = dt;

  if(verbose)
  std::cout << std::endl
            << " :: Implicit Mechanics Solver  " << std::endl
            << " :: No. mpi processes ... " << n_mpi_processes << std::endl;

  std::cout << "Setting up system" << std::endl;
  SetupSystem (domain,partition);
  std::cout << "Assembling matrix" << std::endl;
  Assemble    (domain,partition,time);
  std::cout << "Solving system" << std::endl;
  Solve       (domain,partition,time);

  if(verbose)
  std::cout << std::endl;
  return dt_return;


}



template <int dim>
void ImplicitMechanicsSolver<dim> :: SetupSystem (PhysicalDomainT&  domain,
                                                SpatialPartition& partition)
{

  using namespace BoundaryConditionFunctions;

  // determine the global/local degree of freedom distribution.

  int n_ghost_rows  = domain.m_feNodeManager.GetNumGhosts();
  int n_local_rows  = domain.m_feNodeManager.DataLengths()-n_ghost_rows;
//  int n_hosted_rows = n_local_rows+n_ghost_rows;

  std::vector<int> gather(n_mpi_processes);

  epetra_comm.GatherAll(&n_local_rows,
                        &gather.front(),
                        1);

  int first_local_row = 0;
  int n_global_rows = 0;

  for(unsigned int p=0; p<n_mpi_processes; ++p)
  {
    n_global_rows += gather[p];
    if(p<this_mpi_process)
      first_local_row += gather[p];
  }

    // create trilinos dof indexing

  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
  iArray1d& is_ghost       = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();

  unsigned local_count = 0;
  for(unsigned r=0; r<trilinos_index.size(); ++r )
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

  assert(local_count == static_cast<unsigned int>(n_local_rows) );

  partition.SynchronizeFields(syncedFields, CommRegistry::lagrangeSolver02);
 
//  ApplyMultiSetBoundaryCondition<R1Tensor>(this, &ImplicitMechanicsSolver<dim>::NonpenetratingBC_DetectContact,
//                                           domain, domain.m_feNodeManager,
//                                           NonPenetratingBoundaryCondition::BoundaryConditionName(), time );



    // create epetra map

  row_map = Teuchos::rcp(new Epetra_Map(dim*n_global_rows,dim*n_local_rows,0,epetra_comm));

    // set up sparsity graph

  sparsity = Teuchos::rcp(new Epetra_FECrsGraph(Copy,*row_map,0));
  

  dummyDof.resize(0);                                 
//  ApplyMultiSetBoundaryCondition<R1Tensor>(this, &ImplicitMechanicsSolver<dim>::NonpenetratingBC_Sparsity,
//                                           domain, domain.m_feNodeManager,
//                                           NonPenetratingBoundaryCondition::BoundaryConditionName(), time );

  RegionMap::iterator
    region     = domain.m_feElementManager.m_ElementRegions.begin(),
    end_region = domain.m_feElementManager.m_ElementRegions.end();

  for(; region != end_region; ++region)
  {
    const unsigned numNodesPerElement = region->second.m_numNodesPerElem;

    std::vector<int> elementLocalDofIndex (dim*numNodesPerElement);

    for(localIndex element = 0; element < region->second.m_numElems; ++element)
    {
      const localIndex* const localNodeIndices = region->second.m_toNodesRelation[element];

      for(unsigned i=0; i<numNodesPerElement; ++i)
      {
        const localIndex localNodeIndex = localNodeIndices[i];
        for( int d=0 ; d<dim ; ++d )
        {
          elementLocalDofIndex[i*dim+d] = dim*trilinos_index[localNodeIndex]+d;
        }
      }

      sparsity->InsertGlobalIndices(elementLocalDofIndex.size(),
                                    &elementLocalDofIndex.front(),
                                    elementLocalDofIndex.size(),
                                    &elementLocalDofIndex.front());
    }
  }
  
  
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();
}


/* Assemble */

template <int dim>
void ImplicitMechanicsSolver<dim> :: Assemble (PhysicalDomainT&  domain,
                                             SpatialPartition& partition , realT time)

{
  using namespace BoundaryConditionFunctions;
    // (re-)init linear system

  matrix   = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*sparsity));
  solution = Teuchos::rcp(new Epetra_FEVector(*row_map));
  rhs      = Teuchos::rcp(new Epetra_FEVector(*row_map));

    // basic nodal data ( = dof data for our problem)

  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);

    // begin region loop

  RegionMap::iterator
    region     = domain.m_feElementManager.m_ElementRegions.begin(),
    end_region = domain.m_feElementManager.m_ElementRegions.end();

  for(; region != end_region; ++region)
  {
    // assume all regions are linear hexes

    LagrangeBasis<dim>    basis(1);
    GaussQuadrature<dim>  quadrature(2);
    FiniteElement<dim>    fe(basis,quadrature);

    const unsigned n_q_points = fe.n_quadrature_points();
    const unsigned dof_per_fe = fe.dofs_per_element();

    // following should hold for scalar node problems

    assert(dof_per_fe == region->second.m_numNodesPerElem);

    // space for element matrix and rhs

    std::vector<R1TensorT<dim> > nodal_coord    (dof_per_fe);
    std::vector<R1TensorT<dim> > nodal_disp    (dof_per_fe);

    Epetra_IntSerialDenseVector  elementLocalDofIndex   (dim*dof_per_fe);
    Epetra_SerialDenseVector     element_rhs     (dim*dof_per_fe);
    Epetra_SerialDenseMatrix     element_matrix  (dim*dof_per_fe,dim*dof_per_fe);

    Epetra_SerialDenseMatrix B(3*(dim-1),dim*dof_per_fe) ;
    Epetra_SerialDenseMatrix D(3*(dim-1),3*(dim-1)) ;
    Epetra_SerialDenseMatrix D_B_Jxw(3*(dim-1),dim*dof_per_fe) ; // = D_{kl}B_{lj}*Jxw(q)

    

    // determine ghost elements

    iArray1d& elem_is_ghost = region->second.GetFieldData<FieldInfo::ghostRank>();

    // begin element loop, skipping ghost elements
    
    for(localIndex element = 0; element < region->second.m_numElems; ++element){
      if(elem_is_ghost[element] < 0){        
      	
      	// fill D - TODO assumes 1 gauss point per element
      	realT& K = region->second.m_mat->StateData(element,0)->BulkModulus;
        realT& G = region->second.m_mat->StateData(element,0)->ShearModulus;
      	realT lambda = K - 2.0/3.0 * G;
      	
      	D(0,0) = lambda+2*G;
      	D(1,1) = lambda+2*G;
      	D(2,2) = lambda+2*G;
      	D(0,1) = lambda;
      	D(0,2) = lambda;
      	D(1,0) = lambda;
      	D(1,2) = lambda;
      	D(2,0) = lambda;
      	D(2,1) = lambda;
      	D(3,3) = G;
      	D(4,4) = G;
      	D(5,5) = G;
      	
      // get global row indices and element node coordinates

      // TODO: the nodal coordinate handling break dimension
      //       independence because m_refposition assumes 3D

  
        const localIndex* const localNodeIndices = region->second.m_toNodesRelation[element];

        unsigned numNodesPerElement = 8;
        for(unsigned i=0; i<numNodesPerElement; ++i)
        {
          const localIndex localNodeIndex = localNodeIndices[i];
          for( int d=0 ; d<dim ; ++d )
          {
            elementLocalDofIndex[i*dim+d] = dim*trilinos_index[localNodeIndex]+d;
          }
  
          nodal_coord[i] = (*domain.m_feNodeManager.m_refposition)[localNodeIndex];
          nodal_disp[i] = (*domain.m_feNodeManager.m_displacement)[localNodeIndex];

        }

      // now the actual integration loop

        fe.reinit(nodal_coord);
  
        element_rhs.Scale(0);
        element_matrix.Scale(0);
  
        for(unsigned q=0; q<n_q_points; ++q)
        {
          B.Scale(0.0);
          for(unsigned i=0; i<dof_per_fe; ++i)
          {
            int col = i*dim;
            B(0,col+0) = fe.gradient(i,q)[0];
            B(4,col+0) = fe.gradient(i,q)[2];
            B(5,col+0) = fe.gradient(i,q)[1];
            B(1,col+1) = fe.gradient(i,q)[1];
            B(3,col+1) = fe.gradient(i,q)[2];
            B(5,col+1) = fe.gradient(i,q)[0];
            B(2,col+2) = fe.gradient(i,q)[2];
            B(3,col+2) = fe.gradient(i,q)[1];
            B(4,col+2) = fe.gradient(i,q)[0];
          }


          /*
          for(unsigned i=0; i<dim*dof_per_fe; ++i)
          {
            for(unsigned j=0; j<dim*dof_per_fe; ++j)
            {
              for( unsigned k=0 ; k<6 ; ++k  )
                for( unsigned l=0 ; l<6 ; ++l  )
                  element_matrix(i,j) += B(k,i) * D(k,l) * B(l,j) * fe.JxW(q) ;
            }
  
          }
          */
          D_B_Jxw.Scale(0.0);
          realT JxW = fe.JxW(q);
          for(unsigned j=0; j<dim*dof_per_fe; ++j){
            for( unsigned k=0 ; k<6 ; ++k  ){
              for( unsigned l=0 ; l<6 ; ++l  ){
              	D_B_Jxw(k,j) += D(k,l) * B(l,j) ;
              }
              D_B_Jxw(k,j) *= JxW; 
            }
          }
          
          for(unsigned i=0; i<dim*dof_per_fe; ++i)
          {
            for(unsigned j=0; j<dim*dof_per_fe; ++j)
            {
              for( unsigned k=0 ; k<6 ; ++k  )
                  element_matrix(i,j) += B(k,i) * D_B_Jxw(k,j);
            }
          }
          
        }

        // assemble into global system

        matrix->SumIntoGlobalValues(elementLocalDofIndex,
                                    element_matrix);
                                    
        // residual forces
        for(unsigned i=0; i<dof_per_fe; ++i){
          for(unsigned j=0; j<dim; ++j){
            unsigned ii = dim*i+j;
            for(unsigned k=0; k<dim*dof_per_fe; ++k){
              element_rhs(k) -= element_matrix(k,ii)* (nodal_disp[i])[j];
            }
          }
        }

        rhs->SumIntoGlobalValues(elementLocalDofIndex,
                                 element_rhs);
  
      } 
    }  // end element
  } // end region
  

  // Apply boundary conditions
  ApplyBoundaryCondition<R1Tensor>(this, &ImplicitMechanicsSolver<dim>::TractionBC, 
                                   domain, domain.m_feFaceManager, "Traction", time );
  ApplyBoundaryCondition<R1Tensor>(this, &ImplicitMechanicsSolver<dim>::DisplacementBC, 
                                   domain, domain.m_feNodeManager, 
                                   Field<FieldInfo::displacement>::Name(), time );
  ApplyMultiSetBoundaryCondition<R1Tensor>(this, &ImplicitMechanicsSolver<dim>::NonpenetratingBC_Apply, 
                                           domain, domain.m_feNodeManager, 
                                           NonPenetratingBoundaryCondition::BoundaryConditionName(), time );
  ApplyBoundaryCondition<R1Tensor>(this, &ImplicitMechanicsSolver<dim>::PressureBC,
                                   domain, domain.m_feFaceManager, "Pressure", time );

  // Global assemble
  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

 // rhs->Print(std::cout);
 // std::cout << matrix->NormInf() << std::endl;
  EpetraExt::RowMatrixToMatlabFile("system-matrix.dat",*matrix);
 // exit(0);
}


/* Solve */

template <int dim>
void ImplicitMechanicsSolver<dim> :: Solve (PhysicalDomainT&  domain,
                                          SpatialPartition& partition, realT time)
{

  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
  iArray1d& is_ghost       = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  Array1dT<R1Tensor>& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
  Array1dT<R1Tensor>* inc_disp = NULL;
 
  if(m_recordIncrementalDisplacement){ 
    inc_disp = &(domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>() );
  }
  
	
 // initial guess for solver
  int dummy;
  double* local_solution = NULL;

  solution->ExtractView(&local_solution,&dummy);
  
  // copy initial guess from geos data structures
  /*
  if(m_recordIncrementalDisplacement){ 
    for(unsigned r=0; r<disp.size(); ++r)
    if(is_ghost[r] < 0)
    {
      for( int d=0 ; d<dim ; ++d )
      {
        int lid = row_map->LID(dim*trilinos_index[r]+d);
        local_solution[lid] = (*inc_disp)[r][d];
      }
    }
  }
  */
	
  // krylov solver

  Epetra_LinearProblem problem(&(*matrix),
                               &(*solution),
                               &(*rhs));

  // ML preconditioner
  //////////////////////

  // create a parameter list for ML options
  Teuchos::ParameterList MLList;

  ML_Epetra::SetDefaults("SA",MLList);

  // create the preconditioning object.
  ML_Epetra::MultiLevelPreconditioner* MLPrec = NULL;
  if(m_useMLPrecond){
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*matrix, MLList);
  }

  // note, using asymmetric solver because ILU is asymmetric.
  // system is (mostly) symmetric though.  note that this is
  // a simple, but far from optimal, combination.  should use
  // CG + ML AMG here.

  AztecOO solver(problem);
          // tell AztecOO to use the ML preconditioner
          if(m_useMLPrecond) solver.SetPrecOperator(MLPrec);

          solver.SetAztecOption(AZ_solver,AZ_bicgstab);
          solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
          solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
          solver.SetAztecOption(AZ_conv,AZ_rhs);
          //solver.SetAztecOption(AZ_output,AZ_none);

  solver.Iterate(1000,numerics.krylov_tol);
  //solution->Print(std::cout);

  // destroy the preconditioner
  if(m_useMLPrecond){
    delete MLPrec;
  }

  // copy vector solution into geos data structures

  for(unsigned r=0; r<disp.size(); ++r)
  if(is_ghost[r] < 0)
  {
    for( int d=0 ; d<dim ; ++d )
    {
      int lid = row_map->LID(dim*trilinos_index[r]+d);
      disp[r][d] += local_solution[lid];
    }

  }

  if(m_recordIncrementalDisplacement){ 
    for(unsigned r=0; r<disp.size(); ++r)
    if(is_ghost[r] < 0)
    {
      for( int d=0 ; d<dim ; ++d )
      {
        int lid = row_map->LID(dim*trilinos_index[r]+d);
        (*inc_disp)[r][d] = local_solution[lid];
      }
    } 	
  }


  // re-sync ghost nodes

  partition.SynchronizeFields(syncedFields, CommRegistry::lagrangeSolver02);
  
                             
  if(m_doApertureUpdate){
  	using namespace BoundaryConditionFunctions;
    ApplyMultiSetBoundaryCondition<R1Tensor>(this, &ImplicitMechanicsSolver<dim>::NonpenetratingBC_UpdateAperture, 
                                             domain, domain.m_feNodeManager, 
                                             NonPenetratingBoundaryCondition::BoundaryConditionName(), time );
  }
}

// Boundary Conditions
///////////////////////


/**
 * @author walsh24
 * @brief Apply traction boundary condition to a given face set
 * 
 */
template <int dim>
void ImplicitMechanicsSolver<dim>::TractionBC( PhysicalDomainT& domain,
                                               ObjectDataStructureBaseT& object ,
                                               BoundaryConditionBase* bc,
                                               const lSet& set, realT time )
{
   
  TractionBoundaryCondition* trbc 
     = dynamic_cast<TractionBoundaryCondition*> ( bc);
  if( trbc ){
  	iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
    iArray1d& face_is_ghost  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
   
    Epetra_IntSerialDenseVector  node_dof(1);
    Epetra_SerialDenseVector     node_rhs(1);
    
    //const R1Tensor& n = trbc->GetDirection(time); // old way
  
    for( lSet::const_iterator fc=set.begin() ; fc!=set.end() ; ++fc ) {	
      if( face_is_ghost[*fc] < 0 ){
        //std::cout << "time "<< time << " traction " << value << std::endl;
        const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, *fc ); 

        /* old way
        realT value = trbc->GetValue(domain.m_feFaceManager,fc,time);
        R1Tensor traction = n;
        traction *= value * area * 0.25; // TODO hardcoded for 4 node face
        */
        R1Tensor traction = trbc->GetTractionOnFace(domain,fc,time);
        traction *= area*0.25; // / domain.m_feFaceManager.m_toNodesRelation[*fc].size(); // if test fail this may be why (prev hard coded for 4 node face)
        
        for( lArray1d::const_iterator nd=domain.m_feFaceManager.m_toNodesRelation[*fc].begin() ;
                 nd!=domain.m_feFaceManager.m_toNodesRelation[*fc].end() ; ++nd )
        {
      
          for(int ii =0; ii < dim;++ii){
            node_dof(0) = dim*trilinos_index[*nd]+ii; 
            node_rhs(0) = traction(ii); 
            rhs->SumIntoGlobalValues(node_dof, node_rhs);
          }
       }
      } 
    }
  }                     	
}

/**
 * @author walsh24
 * @brief Apply pressure boundary condition to a given face set
 *
 */
template <int dim>
void ImplicitMechanicsSolver<dim>::PressureBC( PhysicalDomainT& domain,
                                               ObjectDataStructureBaseT& ,
                                               BoundaryConditionBase*,
                                               const lSet& set, realT )
{

  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
  iArray1d& face_is_ghost  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  rArray1d& pressures = domain.m_feFaceManager.GetFieldData<realT>(PS_STR::PressureStr);
  rArray1d* apertures;

  Epetra_IntSerialDenseVector  node_dof(1);
  Epetra_SerialDenseVector     node_rhs(1);

  // If aperture is defined and is <= 0, the pressure BC is not enforced.
  bool isApertureControlled = false;
  if( domain.m_feFaceManager.HasField<realT>(PS_STR::ApertureStr) ){
    isApertureControlled = true;
    apertures = &(domain.m_feFaceManager.GetFieldData<realT>(PS_STR::ApertureStr));
  }

  for( lSet::const_iterator fc=set.begin() ; fc!=set.end() ; ++fc ) {
    if( face_is_ghost[*fc] < 0 ){
      realT pressure = pressures[*fc];

      if(isApertureControlled){
        if ( (*apertures)[*fc] <= 0.0) break;
      }

      const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, *fc );

      R1Tensor traction = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, *fc);
      traction *= -pressure * area/ domain.m_feFaceManager.m_toNodesRelation[*fc].size();

      for( lArray1d::const_iterator nd=domain.m_feFaceManager.m_toNodesRelation[*fc].begin() ;
          nd!=domain.m_feFaceManager.m_toNodesRelation[*fc].end() ; ++nd )
      {

        for(int ii =0; ii < dim;++ii){
          node_dof(0) = dim*trilinos_index[*nd]+ii;
          node_rhs(0) = traction(ii);
          rhs->SumIntoGlobalValues(node_dof, node_rhs);
        }
      }
    }
  }


}

/**
 * @author walsh24
 * @brief  Apply displacement boundary condition to a node set.
 * 
 */
template <int dim>
void ImplicitMechanicsSolver<dim> :: 
DisplacementBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
               BoundaryConditionBase* bc, const lSet& set, realT time){

  const realT LARGE = 1e64;

  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);
  iArray1d& node_is_ghost  = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  Array1dT<R1Tensor>& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
  Epetra_IntSerialDenseVector  node_dof(1);
  Epetra_SerialDenseVector     node_rhs(1);
  Epetra_SerialDenseMatrix     node_matrix  (1,1);
  node_matrix(0,0) = LARGE; 
  int component = bc->GetComponent(time);
  if (component == 2 && domain.m_feFaceManager.m_toNodesRelation[0].size() == 2)
    throw GPException("This is a 2D problem but you are applying BC to component 2 (z).");
  for( lSet::const_iterator nd=set.begin() ; nd!=set.end() ; ++nd ) {	 
      if( node_is_ghost[*nd] < 0 ){
   //   int rank = 0;
   //  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   //   std::cout << trilinos_index[*nd] << " " << component << " "<< rank << " " << bc->GetValue(domain.m_faceManager,nd,time) << std::endl;
        
        node_dof(0) = dim*trilinos_index[*nd]+component; 
        node_rhs(0) = LARGE*(bc->GetValue(domain.m_feFaceManager,nd,time) - (disp[*nd])[component] );
        
        rhs->ReplaceGlobalValues(node_dof, node_rhs);     
        matrix->ReplaceGlobalValues(node_dof, node_matrix);
      //  rhs->SumIntoGlobalValues(node_dof, node_rhs);     
      //  matrix->SumIntoGlobalValues(node_dof, node_matrix);
      }
  }
}

/**
 * @author walsh24
 * @brief  Update neighbors for non-penetrating boundary condition. 
 * 
 */
template <int dim>
void ImplicitMechanicsSolver<dim> :: 
NonpenetratingBC_NeighborUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                                BoundaryConditionBase* bc, realT time  ){
               	
  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   	
  if(npbc) npbc->UpdateNearestNeighborMaps(domain);         	
               	
}

/**
 * @author walsh24
 * @brief  Fix displacements of penetrating neighbors - caution - assumes nodes are in close proximity at contact. 
 * 
 */
template <int dim>
void ImplicitMechanicsSolver<dim> :: 
NonpenetratingBC_DetectContact(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
                               BoundaryConditionBase* bc, realT time  ){
  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   
  if(npbc){     
  
  	iArray1d& inContact = domain.m_feNodeManager.GetFieldData<int>("inContact");
  	
    Array1dT<R1Tensor>& disp = object.GetFieldData<FieldInfo::displacement> ();
    const Array1dT<R1Tensor>& X = object.GetFieldData<FieldInfo::referencePosition> ();
    
//    iArray1d& face_is_ghost  = domain.m_faceManager.GetFieldData<FieldInfo::ghostRank>();
    
  	lSet&  contactingNodes = npbc->m_contactingNodes;
  	std::map<localIndex,localIndex>&  nearestNodeNeighborMap = npbc->m_nearestNodeNeighborMap;
  	
  	std::map<localIndex,localIndex>::iterator itr = npbc->m_nearestFaceNeighborMap.begin();
  	std::map<localIndex,localIndex>::iterator iend = npbc->m_nearestFaceNeighborMap.end();
  	for(; itr != iend; ++itr){
  	  localIndex fc  = itr->first;
//  	  localIndex fc_nbr = itr->second;
  	 
  	  // loop over nodes in face
  	  lArray1d& nodeList = domain.m_feFaceManager.m_toNodesRelation[fc];
  	  R1Tensor n = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fc);
  	  for(lArray1d::size_type i =0 ; i< nodeList.size(); ++i){
  	 	
  	    localIndex nd = nodeList[i];
  	    localIndex nbr = nearestNodeNeighborMap[nd];
  	   
  	    R1Tensor pa = disp[nd] + X[nd];
  	    R1Tensor l = disp[nbr] + X[nbr] - pa;  // branch vector
  	       	 
  	    if(Dot(n,l) <=  0.0  && (nearestNodeNeighborMap[nbr] == nd) ){
  	 	  contactingNodes.insert(nd);
  	 	  contactingNodes.insert(nbr);
  	 	    
  	 	  l*=0.5;
  	 	  disp[nd] += l;
  	 	  disp[nbr] -= l;
  	 	    
   	      inContact[nd] = nbr;
   	      inContact[nbr] = nd;
   	    }
  	  }
    }
  }                     	
}

template <int dim>
void ImplicitMechanicsSolver<dim> :: 
NonpenetratingBC_UpdateAperture(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                               BoundaryConditionBase* bc, realT time  ){
  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   
  if(npbc){     
  
  	rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>("Aperture");
  	
  	
  	std::map<localIndex,localIndex>::iterator itr = npbc->m_nearestFaceNeighborMap.begin();
  	std::map<localIndex,localIndex>::iterator iend = npbc->m_nearestFaceNeighborMap.end();
  	for(; itr != iend; ++itr){
  	  localIndex fc  = itr->first;
  	  localIndex fc_nbr = itr->second;
  	  
  	  //if(npbc->m_nearestFaceNeighborMap[fc_nbr] == fc){
  	  	
        R1Tensor norm  = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fc);
        R1Tensor normB = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fc_nbr);
        R1Tensor center; domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager, fc,center);
        R1Tensor centerB; domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager, fc_nbr,centerB);
        norm -= normB;
        norm.Normalize();
        
  	    aperture[fc] = std::max((centerB-center)*norm,0.0);
  	  //}
    }
  }                     	
}

/**
 * @author walsh24
 * @brief  Fix displacements of penetrating neighbors - caution - assumes nodes are in close proximity. 
 * 
 */
template <int dim>
void ImplicitMechanicsSolver<dim> :: 
NonpenetratingBC_Sparsity(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                          BoundaryConditionBase* bc, realT time  ){
  

  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   
  if(npbc){
  	
    iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);	
    iArray1d& node_is_ghost  = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
    
    std::map<localIndex,localIndex>&  nearestNodeNeighborMap = npbc->m_nearestNodeNeighborMap;
    
    std::vector<int>  node_row_dof(1);
    std::vector<int>  node_row_dofB(dim);
    std::vector<int>  node_col_dofB(2*dim);
    
    // give contacting nodes the same trilinos index, 
    // add other index to dummy trilinos indicies.
    lSet&  contactingNodes = npbc->m_contactingNodes;
    lSet::iterator iend = contactingNodes.end();
    for( lSet::iterator itr = contactingNodes.begin(); itr != iend; ++itr ) {	
      localIndex nd =  *itr;
      localIndex nbr = nearestNodeNeighborMap[nd];
      
      if(trilinos_index[nbr] > trilinos_index[nd]){
          
        int triNbr = trilinos_index[nbr];
        trilinos_index[nbr] = trilinos_index[nd];

        for(int component =0; component < dim; ++component){
          node_row_dof[0] = dim*triNbr+component;
          dummyDof.push_back(node_row_dof[0]);
          sparsity->InsertGlobalIndices(node_row_dof.size(),
                                        &node_row_dof.front(),
                                        node_row_dof.size(),
                                        &node_row_dof.front());
      	}
     }
     
    }
    //std::cout << "Number of Contacting Nodes: " <<  contactingNodes.size() <<std::endl;
    
    // add resistance between non contacting nodes.

    std::map<localIndex,localIndex>::iterator iendB = nearestNodeNeighborMap.end();
    for( std::map<localIndex,localIndex>::iterator itr = nearestNodeNeighborMap.begin();
         itr != iendB; ++itr ) {	
      localIndex nd =  itr->first;
      if(node_is_ghost[nd] < 0){
        localIndex nbr = nearestNodeNeighborMap[nd];
        int tri = trilinos_index[nd];
        int triNbr = trilinos_index[nbr];
        if(tri != triNbr){  // ie not in contact
          for(int component =0; component < dim; ++component){
            node_row_dofB[component] = dim*tri+component;
            node_col_dofB[component] = node_row_dofB[component];
            node_col_dofB[dim+component] = dim*triNbr+component;
          }
          sparsity->InsertGlobalIndices(node_row_dofB.size(),
                                        &node_row_dofB.front(),
                                        node_col_dofB.size(),
                                        &node_col_dofB.front());

        }
      }
    }
    	         	
  }             	
}

/**
 * @author walsh24
 * @brief  Fix displacements of penetrating neighbors - caution - assumes neighbors in close proximity. 
 * 
 */
template <int dim>
void ImplicitMechanicsSolver<dim> :: 
NonpenetratingBC_Apply(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
               BoundaryConditionBase* bc, realT time  ){


  NonPenetratingBoundaryCondition* npbc = dynamic_cast<NonPenetratingBoundaryCondition*> (bc);   
  if(npbc){
//    const Array1dT<R1Tensor>& disp = object.GetFieldData<FieldInfo::displacement> ();
  	
    iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);	
//    iArray1d& node_is_ghost  = domain.m_nodeManager.GetFieldData<FieldInfo::ghostRank>();
    iArray1d& face_is_ghost  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  	
    Epetra_IntSerialDenseVector  node_row_dof(1);
    Epetra_IntSerialDenseVector  node_row_dofB(dim);
    Epetra_IntSerialDenseVector  node_col_dofB(2*dim);
    Epetra_SerialDenseVector     node_rhsB(dim);
    Epetra_SerialDenseMatrix     node_matrix  (1,1);
    Epetra_SerialDenseMatrix     node_matrixB  (dim,2*dim);
  	
    // Set diagonal of dummy dofs to 1
    node_matrix(0,0) = 1.0;
    for( iArray1d::size_type i=0; i < dummyDof.size() ; ++i ) {
      node_row_dof(0) = dummyDof[i];	
      matrix->ReplaceGlobalValues(node_row_dof,node_matrix);
    }       	

    // add entries to represent springs between non contacting nodes.
    std::map<localIndex,localIndex>&  nearestFaceNeighborMap = npbc->m_nearestFaceNeighborMap;
    std::map<localIndex,localIndex>&  nearestNodeNeighborMap = npbc->m_nearestNodeNeighborMap;
    std::map<localIndex,localIndex>::iterator iendB = nearestFaceNeighborMap.end();
    for( std::map<localIndex,localIndex>::iterator itr = nearestFaceNeighborMap.begin();
         itr != iendB; ++itr ) {	
      localIndex fc =  itr->first;	
      
      if(face_is_ghost[fc] < 0){
        localIndex fcb =  itr->second;
        R1Tensor norm  = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fc);
        R1Tensor normB = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, fcb);
        R1Tensor center; domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager, fc,center);
        R1Tensor centerB; domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager, fcb,centerB);
        norm -= normB;
        norm.Normalize();

        R1Tensor dl = centerB-center;
 
//        realT l = sqrt(Dot(dl,dl));
        realT a = 0.5*(domain.m_feFaceManager.SurfaceArea(domain.m_feNodeManager, fc)
                       + domain.m_feFaceManager.SurfaceArea(domain.m_feNodeManager, fcb));

        realT k =  0.25*a*m_nonContactModulus; // FIXME /l;    // spring stiffness
   
  	    lArray1d& nodeList = domain.m_feFaceManager.m_toNodesRelation[fc];
  	    for( lArray1d::size_type nn =0 ; nn< nodeList.size(); ++nn){
          localIndex nd = nodeList[nn];
          localIndex nbr = nearestNodeNeighborMap[nd];
          
          int tri = trilinos_index[nd];
          int triNbr = trilinos_index[nbr];
          if(tri != triNbr){  // not in contact

           // R1Tensor du = disp[nd] - disp[nbr];
           // R1Tensor f_rhs = norm; f_rhs *= k*Dot(norm,du);

            for(int i =0; i < dim; ++i){
              node_row_dofB(i) = dim*tri+i;
              node_col_dofB(i) = node_row_dofB(i);
              node_col_dofB(dim+i) = dim*triNbr+i;
              for(int j =0; j < dim; ++j){
                realT val = k*(0.1+0.9*norm[i]*norm[j]);
                node_matrixB(i,j) = val;
                node_matrixB(i,j+dim) = -val;
              }
             // node_rhsB(i) = f_rhs[i];
            }
            matrix->SumIntoGlobalValues(node_row_dofB,node_col_dofB,node_matrixB);
           // rhs->SumIntoGlobalValues(node_row_dofB,node_rhsB);
          }
        }
      }
    }

    if(npbc->m_updatePressure){

      // loop over second set - update pressure based on nearest neighbor from first
      rArray1d& pressures = domain.m_feFaceManager.GetFieldData<realT>(PS_STR::PressureStr);

      const lSet& setB = domain.m_feFaceManager.m_Sets[npbc->m_setNames[1] ];
      for( lSet::const_iterator fcB=setB.begin(); fcB!=setB.end() ; ++fcB ) {
        localIndex nbr = nearestFaceNeighborMap[*fcB];
        pressures[*fcB] = pressures[nbr];
      }
    }

  }             	
}



/* Explicit Instantiations */

//template class ImplicitMechanicsSolver<2>;
template class ImplicitMechanicsSolver<3>;


/* Register solver in the solver factory */

SolverRegistrator<ImplicitMechanicsSolver<3> > reg_ImplicitMechanicsSolver3D;

//SolverRegistrator<ImplicitMechanicsSolver<2> >
//  reg_ImplicitLaplaceSolver2D;



