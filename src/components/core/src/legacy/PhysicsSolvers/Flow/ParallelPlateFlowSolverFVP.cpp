// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/**
 * @file TwoDSteadyStateParallelPlateFlowSolver.cpp
 * @author walsh24
 * @date February 21, 2012
 */

#include "ParallelPlateFlowSolverFVP.h"
#include "PhysicsSolvers/Flow/FractureMatrixFlowSolverFV.h"
#include "../SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"
#include "PhysicsSolvers/Flow/FractureMatrixFlowSolverFV.h"


// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"



#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"


using namespace BoundaryConditionFunctions;
using namespace PS_STR;
using namespace PPFS;


int ParallelPlateFlowSolverFVP::m_instances = 0;



static realT extraMassLoss( const realT extraMass,
                            const realT mass,
                            const realT volume,
                            const realT referenceMass,
                            const realT dt,
                            const realT tref )
{

//  realT const density = mass/volume;
  realT massLoss = 0.0;

//  realT massLoss = referenceMass / 100.0 * dt;


  if( ( mass - referenceMass ) > massLoss )
  {
    massLoss = ( mass - referenceMass );
  }


  if( massLoss > extraMass )
  {
    massLoss = extraMass;
  }

  return massLoss;
}


/**
 * @author settgast
 * @param nodeIndex the local node index relative to the face
 * @param ndofPerNode number of degrees of freedom associated with a node
 *
 * This function gives the dof index for a node-pair associated with an original
 * parent node on a parent face. This function is used for consistency when
 * forming and applying contributions to the system.
 */
static void faceNodePairIndexing( const int nodeIndex,
                                  const int ndofPerNode,
                                  int index[2])
{
  index[0] = ndofPerNode*(nodeIndex*2);
  index[1] = ndofPerNode*(nodeIndex*2+1);
}

ParallelPlateFlowSolverFVP::ParallelPlateFlowSolverFVP(  const std::string& name,
                                                         ProblemManagerT* const pm ):
  ParallelPlateFlowSolverBase(name,pm),
  m_faceSet(),
  m_numFaces(0),
  m_faceDofMap(),
  m_phi(1.0),
  this_mpi_process(pm->m_epetraComm.MyPID()),
  n_mpi_processes(pm->m_epetraComm.NumProc()),
  m_wellboreSolve(NULL)
{
  ++m_instances;
  m_trilinosIndexStr = "ParallelPlateFlowSolverFVP_" +  toString<int>(m_instances) + "_GlobalDof";

  m_dofVariable = dofVariable::pressure;
//  m_dofVariable = dofVariable::mass;
}

ParallelPlateFlowSolverFVP::~ParallelPlateFlowSolverFVP()
{
  // TODO Auto-generated destructor stub
}

void ParallelPlateFlowSolverFVP::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  ParallelPlateFlowSolverBase::ReadXML(hdn);

  // Mixed difference parameter
  m_phi = hdn->GetAttributeOrDefault<realT>("phi",1.0); // backward difference
                                                        // by default.

  // Linear Solver
  m_numerics.krylov_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-10);
  m_numerics.m_maxIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);

  // Flags
  m_doApertureUpdate = hdn->GetAttributeOrDefault<bool>("UpdateAperture",false);

  m_gamma = hdn->GetAttributeOrDefault("gamma","0.0");

  m_negPressAperFix = hdn->GetAttributeOrDefault("negPressAperFix",0);

  m_phiLimiter = hdn->GetAttributeOrDefault("phiLimiter",0);

  m_IncompressibleFlow = hdn->GetAttributeOrDefault("incompressibleFlow",false);

  m_steadyStateFlow = hdn->GetAttributeOrDefault("steadyStateFlow",false);

  m_matrixFlowSolverName = hdn->GetAttributeString("matrixFlowSolverName");

  m_linearizeGoverningEquation = hdn->GetAttributeOrDefault("linearize",0);
  if( m_linearizeGoverningEquation != 0 )
  {
    m_permeabilityUpdateOption = 0;
//    m_numerics.m_maxIterNewton = 1;
  }


  m_leakExtraMass = hdn->GetAttributeOrDefault<int>("leakExtraMass",1);
  m_newVolumeDensityRatio = hdn->GetAttributeOrDefault<realT>("newVolumeDensityRatio",1.0);

  this->m_fluidEOS = new PPFS::BiLinearEOS( this->m_rho_o, this->m_bulk_modulus, this->m_negativeFluidPressureSlopeReduction * this->m_bulk_modulus);

  std::string dofVarString = hdn->GetAttributeStringOrDefault("dofVar","Pressure");
  if( !(dofVarString.compare("Pressure")) )
  {
    m_dofVariable = dofVariable::pressure;
  }
  else if( !(dofVarString.compare("Mass")) )
  {
    m_dofVariable = dofVariable::mass;
  }

  m_subscaleAngle = hdn->GetAttributeOrDefault<int>("subscaleAngle",0);

}


void ParallelPlateFlowSolverFVP::RegisterFields( PhysicalDomainT& domain )
{

  ParallelPlateFlowSolverBase::RegisterFields( domain.m_feFaceManager, domain.m_feEdgeManager );

  const bool plotOldValues = false; // no need to plot old values unless
                                    // debugging - overwritten at end of
                                    // timestep.

  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("lastAper",true,false);

  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FluidVelocityStr,true,true);
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::pressure>();

  domain.m_feFaceManager.AddKeylessDataField<realT>("pressureIncrement",true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("apertureIncrement",true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("massIncrement",true,true);

  domain.m_feFaceManager.AddKeylessDataField<int>("isInitialized",true,true);

  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::density>();
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::mass>();
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::volume>();
//  domain.m_feFaceManager.AddKeylessDataField<realT>("BoundedAperture",true,true);

  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("faceNormal0",true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("extraMass",true,true);

  domain.m_feFaceManager.AddKeylessDataField<int>(m_trilinosIndexStr,true,true);
  domain.m_feEdgeManager.AddKeylessDataField<int>("FlowFaceCount",true,true);// debug
  domain.m_feEdgeManager.AddKeylessDataField<R1Tensor>(EdgeCenterStr,true,plotOldValues);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(EdgeLengthStr,true,true);

  domain.m_feEdgeManager.AddKeylessDataField<realT>("phi",true,true);// debug

  // Wellbore flow solver
  if (!m_wellboreSolverName.empty())
  {
    m_wellboreSolve = dynamic_cast<WellboreFlowSolverImplicit*>(stlMapLookup( m_solvers, m_wellboreSolverName, "m_solvers" ));

    // Override the default pressureCapEOS in the wellbore solver
    m_wellboreSolve->m_fluidEOS = this->m_fluidEOS;
  }
  else
  {
    m_wellboreSolve = NULL;
  }

}

void ParallelPlateFlowSolverFVP::RegisterTemporaryFields( PhysicalDomainT& domain )
{
  domain.m_feFaceManager.AddKeylessDataField<realT>("FluidMass_n",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>("FluidVolume_n",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>("FluidDensity_n",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>("Pressure_n",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>("dPdM",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr+std::string("_n"),false,false);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr+std::string("_n"),false,false);
//  domain.m_feFaceManager.AddKeylessDataField<realT>("BoundedAperture_n",false,false);

  domain.m_feEdgeManager.AddKeylessDataField<R1Tensor>(EdgeCenterStr+std::string("_n"),false,false);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(EdgeLengthStr+std::string("_n"),false,false);


}


void ParallelPlateFlowSolverFVP::DeregisterTemporaryFields( PhysicalDomainT& domain )
{

  domain.m_feFaceManager.RemoveDataField<realT>("FluidMass_n");
  domain.m_feFaceManager.RemoveDataField<realT>("FluidVolume_n");
  domain.m_feFaceManager.RemoveDataField<realT>("FluidDensity_n");
  domain.m_feFaceManager.RemoveDataField<realT>("Pressure_n");
  domain.m_feFaceManager.RemoveDataField<realT>("dPdM");
  domain.m_feFaceManager.RemoveDataField<realT>(ApertureStr+std::string("_n"));
  domain.m_feFaceManager.RemoveDataField<R1Tensor>(FaceCenterStr+std::string("_n"));
//  domain.m_feFaceManager.RemoveDataField<R1Tensor>("BoundedAperture_n");

  domain.m_feEdgeManager.RemoveDataField<R1Tensor>("EdgeCenter_n");
  domain.m_feEdgeManager.RemoveDataField<realT>(EdgeLengthStr+std::string("_n"));
}


void ParallelPlateFlowSolverFVP::FillTemporaryFields( PhysicalDomainT& domain )
{

  const array<real64>& mass_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  array<real64>& mass_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");
  mass_n = mass_np1;

  const array<real64>& volume_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  array<real64>& volume_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidVolume_n");
  volume_n = volume_np1;

  const array<real64>& density_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  array<real64>& density_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");
  density_n = density_np1;

  const array<real64>& pressure_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  array<real64>& pressure_n   = domain.m_feFaceManager.GetFieldData<realT>("Pressure_n");
  pressure_n = pressure_np1;

  const array<real64>& apertures_np1 = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  array<real64>& apertures_n   = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n") );
  apertures_n = apertures_np1;

//  const array<real64>& boundedApertures_np1 =
// domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture" );
//  array<real64>& boundedApertures_n   =
// domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture_n" );
//  boundedApertures_n = boundedApertures_np1;

  const array<R1Tensor>& faceCenters_np1 = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  array<R1Tensor>& faceCenters_n   = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr+std::string("_n") );
  faceCenters_n = faceCenters_np1;

  const array<R1Tensor>& edgeCenters_np1 = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  array<R1Tensor>& edgeCenters_n   = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr+std::string("_n"));
  edgeCenters_n = edgeCenters_np1;

  const array<real64>& edgeLengths_np1 = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  array<real64>& edgeLengths_n   = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr+std::string("_n") );
  edgeLengths_n = edgeLengths_np1;

  array<real64>& phiStep   = domain.m_feEdgeManager.GetFieldData<realT>( "phi" );
  phiStep = m_phi;

  if (m_wellboreSolve)
  {
    m_wellboreSolve->FillTemporaryFields( domain );
  }
}


void ParallelPlateFlowSolverFVP::OverwriteFieldsWithTemporaryFields( PhysicalDomainT& domain )
{
  array<real64>& mass_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const array<real64>& mass_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");
  mass_np1 = mass_n;

  array<real64>& volume_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  const array<real64>& volume_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidVolume_n");
  volume_np1 = volume_n;

  array<real64>& density_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const array<real64>& density_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");
  density_np1 = density_n;

  array<real64>& pressure_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const array<real64>& pressure_n   = domain.m_feFaceManager.GetFieldData<realT>("Pressure_n");
  pressure_np1 = pressure_n;

  array<real64>& apertures_np1 = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const array<real64>& apertures_n   = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n") );
  apertures_np1 = apertures_n;

//  array<real64>& boundedApertures_np1 =
// domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture" );
//  const array<real64>& boundedApertures_n   =
// domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture_n" );
//  boundedApertures_np1 = boundedApertures_n;

  array<R1Tensor>& faceCenters_np1 = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const array<R1Tensor>& faceCenters_n   = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr+std::string("_n") );
  faceCenters_np1 = faceCenters_n;

  array<R1Tensor>& edgeCenters_np1 = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  const array<R1Tensor>& edgeCenters_n   = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr+std::string("_n"));
  edgeCenters_np1 = edgeCenters_n;

  array<real64>& edgeLengths_np1 = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  const array<real64>& edgeLengths_n   = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr+std::string("_n") );
  edgeLengths_np1 = edgeLengths_n;

}


void ParallelPlateFlowSolverFVP::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
  ParallelPlateFlowSolverBase::Initialize( domain, partition );
  DefineFlowSets(domain);

  array<real64>& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");
  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    faceArea[kf] = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, kf, 1 );
  }



  RegisterTemporaryFields( domain );
  GenerateParallelPlateGeometricQuantities( domain,0,0 );
  UpdateEOS(0,0,domain);
  DeregisterTemporaryFields( domain );
//  FillTemporaryFields( domain );

//  InitializeDensity( domain);

  array<real64>& phi   = domain.m_feEdgeManager.GetFieldData<realT>( "phi" );
  phi = m_phi;

  if (m_leakoffModel.GetLeakoffCoefficient() > 0)
  {
    array<real64>& initialSaturatedTime = domain.m_feFaceManager.GetFieldData<realT>("initialSaturatedTime");
    initialSaturatedTime = std::numeric_limits<realT>::max();
  }

}



void ParallelPlateFlowSolverFVP::SetNumRowsAndTrilinosIndices( PhysicalDomainT& domain,
                                                               SpatialPartition& partition,
                                                               int& numLocalRows,
                                                               int& numGlobalRows,
                                                               array<integer>& localIndices,
                                                               int offset )
{

  array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
  array<integer>& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  trilinos_index = -INT_MAX;

  //if(m_doApertureUpdate) UpdateAperture(domain);


  // count local dof
  ///////////////////////////////

  // local rows
  numLocalRows = 0;

  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if( is_ghost[*kf] < 0 )
    {
      ++numLocalRows;
    }
  }


  // determine the global/local degree of freedom distribution.
  ////////////////////////////////////////////////////////////

  std::vector<int> gather(n_mpi_processes);

  epetra_comm->GatherAll(&numLocalRows,
                         &gather.front(),
                         1);

  int first_local_row = 0;
  numGlobalRows = 0;

  for( int p=0 ; p<n_mpi_processes ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<this_mpi_process)
      first_local_row += gather[p];
  }

  // create trilinos dof indexing
  //////////////////////////////////
  unsigned local_count = 0;
  // faces
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      trilinos_index[*kf] = first_local_row+local_count+offset;
      localIndices.push_back(trilinos_index[*kf]);
      local_count++;
    }
    else
    {
      trilinos_index[*kf] = -INT_MAX;
    }
  }



  assert(static_cast<int>(local_count) == numLocalRows);

  partition.SynchronizeFields(m_syncedFields, CommRegistry::parallelPlateFlowSolver);


  if( m_matrixFlowSolver != nullptr )
  {
    int numMatrixLocalRows = 0;
    int numMatrixGlobalRows = 0;
    m_matrixFlowSolver->SetNumRowsAndTrilinosIndices( domain, partition,
                                                      numMatrixLocalRows, numMatrixGlobalRows,
                                                      localIndices,
                                                      numGlobalRows );
    numLocalRows += numMatrixLocalRows;
    numGlobalRows += numMatrixGlobalRows;
  }

  // Record starting trilinos index for wellbore solver
  // Add in dof's for wellbore solver
  if( m_wellboreSolve != nullptr )
  {
    int numWellboreLocalRows = 0;
    int numWellboreGlobalRows = 0;

    m_wellboreSolve->SetNumRowsAndTrilinosIndices( domain, partition,
                                                   numWellboreLocalRows,
                                                   numWellboreGlobalRows,
                                                   localIndices,
                                                   numGlobalRows );
    numLocalRows += numWellboreLocalRows;
    numGlobalRows += numWellboreGlobalRows;

  }
}

void ParallelPlateFlowSolverFVP:: SetupSystem (PhysicalDomainT&  domain,
                                               SpatialPartition& partition )
{
  int numLocalRows;
  int numGlobalRows;
  array<integer> indices;
  SetNumRowsAndTrilinosIndices( domain, partition,
                                numLocalRows,
                                numGlobalRows,
                                indices,
                                0 );

#if 0
  {
    array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
    array<integer>& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

    for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
    {
      if( is_ghost[*kf] < 0 )
      {
        indices.push_back( trilinos_index[*kf] );
      }
    }
    this->m_matrixFlowSolver->GetTrilinosIndices(domain,indices);
  }
#endif

  #if USECPP11==1
  // create epetra map
  ////////////////////

  m_rowMap = std::make_shared<Epetra_Map>(numGlobalRows,numLocalRows,0,*epetra_comm);

  // set up sparsity graph
  ////////////////////////

  m_sparsity = std::make_shared<Epetra_FECrsGraph>(Copy,*m_rowMap,0);
#else
  delete m_rowMap;
  delete m_sparsity;
  m_rowMap = new Epetra_Map(numGlobalRows,numLocalRows,indices.data(),0,*epetra_comm);
  m_sparsity = new Epetra_FECrsGraph(Copy,*m_rowMap,0);

#if 0
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
    const array<integer>& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

    for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
    {
      if(is_ghost[*kf] < 0)
      {
        std::cout<<"rank "<<rank<<", trilinosIndex["<<*kf<<"]= "<<trilinos_index[*kf]<<std::endl;
      }

    }

    for( int row=0 ; row<numGlobalRows ; ++row )
    {
      std::cout<<"rank "<<rank<<", LID("<<row<<")= "<<m_rowMap->LID(row)<<std::endl;
    }
  }
#endif
#endif

  SetSparsityPattern( domain );
  m_sparsity->GlobalAssemble();
  m_sparsity->OptimizeStorage();

#if USECPP11==1
  m_matrix   = std::make_shared<Epetra_FECrsMatrix>(Copy,*m_sparsity);
  m_solution = std::make_shared<Epetra_FEVector>(*m_rowMap);
  m_rhs      = std::make_shared<Epetra_FEVector>(*m_rowMap);
#else
  delete m_matrix;
  delete m_solution;
  delete m_rhs;
  m_matrix   = new Epetra_FECrsMatrix(Copy,*m_sparsity);
  m_solution = new Epetra_FEVector(*m_rowMap);
  m_rhs      = new Epetra_FEVector(*m_rowMap);

#endif
}


void ParallelPlateFlowSolverFVP::SetSparsityPattern( PhysicalDomainT& domain )
{
  array<integer>& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  const array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);

  // loop over edges
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin() ; itr!=itrEnd ; ++itr )
  {
    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0 )
    {
      unsigned int numFaces = itr->second.size();
      if( numFaces > 0)
      {
        std::vector<int> dofIndex (numFaces);
        for(unsigned i=0 ; i<numFaces ; ++i)
        {
          localIndex kf = itr->second[i];
          dofIndex[i] = trilinos_index[kf];
        }

        m_sparsity->InsertGlobalIndices(dofIndex.size(),
                                        &dofIndex.front(),
                                        dofIndex.size(),
                                        &dofIndex.front());
      }
    }
  }



  if( m_matrixFlowSolver != nullptr )
  {
    m_matrixFlowSolver->SetSparsityPtr(m_sparsity);
    m_matrixFlowSolver->SetSparsityPattern(domain);
  }
  // Setup the wellbore sparsity map
  if (m_wellboreSolve)
  {
    m_wellboreSolve->SetupSparsity( domain, m_trilinosIndexStr, m_sparsity);
  }

}



/* Assemble */

realT ParallelPlateFlowSolverFVP::Assemble ( PhysicalDomainT&  domain,
                                             Epetra_System& epetraSystem,
                                             const realT& time,
                                             const realT& dt )
{
  realT maxMassScale = 0.0;
  // (re-)init linear system

  // basic face data ( = dof data for our problem)

  array<integer>& trilinosIndexFace = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);

  array<integer>& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  array<integer>& face_is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

//  const array<R1Tensor>& faceCenters_np1 =
// domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const array<R1Tensor>& faceCenters_n   = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr+std::string("_n") );

//  const array<R1Tensor>& edgeCenters_np1 =
// domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  const array<R1Tensor>& edgeCenters_n   = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr+std::string("_n"));

  const array<real64>& edgeLengths_np1 = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  array<real64>& edgeLengths_n   = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr+std::string("_n") );

  const array<real64>& apertures_np1 = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const array<real64>& apertures_n   = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n") );
//  array<real64> const & apertureIncrement =
// domain.m_feFaceManager.GetFieldData<realT>("apertureIncrement");

//  const array<real64>& mass_np1 =
// domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const array<real64>& mass_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");
  const array<real64>& massIncrement   = domain.m_feFaceManager.GetFieldData<realT>("massIncrement");
  array<real64>& extraMass     = domain.m_feFaceManager.GetFieldData<realT>("extraMass");

  const array<real64>& volume_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  const array<real64>& volume_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidVolume_n");

  const array<real64>& density_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const array<real64>& density_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");

  array<real64>& pressure_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const array<real64>& pressure_n   = domain.m_feFaceManager.GetFieldData<realT>("Pressure_n");
  array<real64> const & pressureIncrement = domain.m_feFaceManager.GetFieldData<realT>("pressureIncrement");

//  const array<real64>& dPdM =
// domain.m_feFaceManager.GetFieldData<realT>("dPdM");

  array<real64>& edgePermeability = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);

  array<real64>& phiStep   = domain.m_feEdgeManager.GetFieldData<realT>( "phi" );

  const array<real64>& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");

  const int dim = 3;

  const array<integer>& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const array<integer>& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");
//  const OrderedVariableOneToManyRelation& childFaceIndex =
// domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );


  const bool displacementCouplingFlag = epetraSystem.HasMatrixBlock( EpetraBlock::fluidBlock, EpetraBlock::solidBlock );

#if USECPP11==1
  std::shared_ptr<Epetra_FECrsMatrix> matrixPD;
#else
  Epetra_FECrsMatrix* matrixPD = nullptr;
#endif

  if (m_leakoffModel.GetLeakoffCoefficient() > 0)
    AssembleLeakoff(domain, time, dt);


  if( displacementCouplingFlag )
  {
    matrixPD = epetraSystem.GetMatrix( EpetraBlock::fluidBlock, EpetraBlock::solidBlock );
  }

  CalculateApertureDerivatives( domain.m_feFaceManager,
                                domain.m_feNodeManager );

  // diagonal terms not involving flow between neighbors
  {
    localIndex numFaces = domain.m_feFaceManager.DataLengths();
    Epetra_IntSerialDenseVector  faceDofIndex (1);
    Epetra_SerialDenseVector     flowRHS     (1);
    Epetra_SerialDenseMatrix     flowMatrix  (1,1);

    Epetra_IntSerialDenseVector  displacementDofIndex;
    Epetra_SerialDenseMatrix     matrix_PD;

    if( !m_steadyStateFlow )
    {
      for( localIndex r=0 ; r<numFaces ; ++r )
      {

        faceDofIndex[0] = trilinosIndexFace[r];
        if( flowFaceType[r] == 1 && face_is_ghost[r] < 0 )
        {

          if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::mass )
          {
            maxMassScale = std::max( maxMassScale, fabs( mass_n[r] + massIncrement[r] ) );
            flowRHS[0] = massIncrement[r];
            flowMatrix(0,0) = 1.0;
          }
          else if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::pressure )
          {
            maxMassScale = std::max( maxMassScale, fabs( volume_np1[r] * density_np1[r] ) );
            const realT dRho_dP = m_IncompressibleFlow ? 0.0 : 1.0 /  m_fluidEOS->dPdRho(density_np1[r]);
            // self contribution to the residual
            flowRHS[0] = (density_np1[r] * volume_np1[r] - density_n[r] * volume_n[r] );

//            std::cout<<density_np1[r]<<" * "<<volume_np1[r]<<" -
// "<<density_n[r]<<" * "<<volume_n[r]<<" = "<<flowRHS[0]<<std::endl;

            /*
                      realT const volumeIncrement = apertureIncrement[r] *
                         faceArea[r];
                      /// TODO this is actually wrong. Using the chain rule to
                         find drho isn't correct unless the EOS is truly linear.
                      realT const densityIncrement = dRho_dP *
                         pressureIncrement[r];

                      realT const Rr = density_n[r] * volumeIncrement
             + volume_n[r] * densityIncrement
             + volumeIncrement * densityIncrement ;
             */
            //flowRHS[0] = -Rr;

            // derivative of the residual wrt pressure
            flowMatrix(0,0) = volume_np1[r] * dRho_dP;

          }

          if( extraMass[r] > 0.0 && m_leakExtraMass > 0 )
          {
            realT massLoss = extraMassLoss( extraMass[r],
                                            mass_n[r],
                                            volume_n[r],
                                            faceArea[r] * ( BoundedAperture(m_zeroApertureOffset) + m_zeroApertureVolume ) * m_rho_o,
                                            dt, 1.0 );
            flowRHS[0] += massLoss;
          }


          m_matrix->SumIntoGlobalValues(faceDofIndex, flowMatrix);
          m_rhs->SumIntoGlobalValues(faceDofIndex,flowRHS);

          if( displacementCouplingFlag && m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::pressure )
          {
            const lArray1d& nodeList= domain.m_feFaceManager.m_toNodesRelation[r];
            const localIndex numNodes = nodeList.size();

            displacementDofIndex.Resize( dim*2*numNodes );
            matrix_PD.Reshape( 1, dim*2*numNodes );
            matrix_PD.Scale(0.0);

            array<real64> dWdU = m_dwdu(r);
            dWdU *= m_dwdw(r);

            for( localIndex a=0 ; a<numNodes ; ++a )
            {

              int rnodeDofIndex[2];
              faceNodePairIndexing( a, dim, rnodeDofIndex);
              for( int i=0 ; i<dim ; ++i )
              {
                for( int side=0 ; side<2 ; ++side )
                {
                  const int ai=rnodeDofIndex[side]+i;
                  displacementDofIndex(ai) = m_dwdu_dof(r)(ai);

                  realT dRrdu = density_np1[r] * faceArea[r] * dWdU(ai);
//                  std::cout<<density_np1[r]<<", "<<faceArea[r]<<",
// "<<dWdU(ai)<<std::endl;
                  matrix_PD(0,ai) += dRrdu;
                }
              }
            }
            matrixPD->SumIntoGlobalValues( faceDofIndex, displacementDofIndex, matrix_PD );
          }
        }
      }
    }
  }


  Epetra_IntSerialDenseVector  faceDofIndex;
  Epetra_SerialDenseVector     flowRHS;
  Epetra_SerialDenseMatrix     flowMatrix;

  realT minKappa = 1.0e99;
  realT maxKappa = 0;

  // loop over edges and process all flow through edges between faces attached
  // to this edge
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin() ; itr!=itrEnd ; ++itr )
  {
    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0  && flowEdgeType[eg] == 1 )
    {

      // number of faces attached to edge.
      const lArray1d& facelist = itr->second;
      const unsigned int numFaces = facelist.size();

      // allocate space for dof map, residual, and dRdM.

      faceDofIndex.Resize(numFaces);
      flowRHS.Resize(numFaces);
      flowMatrix.Reshape(numFaces,numFaces);

      flowRHS.Scale(0.0);
      flowMatrix.Scale(0.0);



      array<real64> dVr_du;
      array<real64> dVs_du;

      array<real64> dKappa_du_rterm;
      array<real64> dKappa_du_sterm;

      // this is the independent face-pair approach. It should be updated to
      // properly handle edges with more than 2 faces.
      for( unsigned int kr=0 ; kr<numFaces ; ++kr )
      {
        const unsigned int r = facelist[kr];
        if( flowFaceType[r] == 1 )
        {

          faceDofIndex[kr] = trilinosIndexFace[r];


          for( unsigned int ks=kr+1 ; ks<numFaces ; ++ks )
          {
            const unsigned int s = facelist[ks];
            if( flowFaceType[s] == 1 )
            {

              localIndex r_eff = r;
              localIndex s_eff = s;

              if( m_negPressAperFix == 1 )
              {
                // if the r face is closed
                if( pressure_n[r] < pressure_n[s] && apertures_n[s] > 3 * apertures_n[r] )
                {
                  r_eff = s;
                }
                // if the s face is closed
                if( pressure_n[s] < pressure_n[r] && apertures_n[r] > 3 * apertures_n[s] )
                {
                  s_eff = r;
                }
              }

              if( edgeLengths_n[eg] == 0.0 )
              {
                edgeLengths_n[eg] = edgeLengths_np1[eg];
              }

              realT dKappa_r_dw = 0;
              realT dKappa_s_dw = 0;

              R1Tensor lengthR_r = edgeCenters_n[eg];
              lengthR_r -= faceCenters_n[r];
              realT length_r = lengthR_r.L2_Norm();

              R1Tensor lengthR_s = edgeCenters_n[eg];
              lengthR_s -= faceCenters_n[s];
              realT length_s = lengthR_s.L2_Norm();

              R1Tensor rs_xdiff;
              rs_xdiff = faceCenters_n[r];
              rs_xdiff -= faceCenters_n[s];
              realT const GravityMag = Dot(rs_xdiff,this->m_gravityVector);


              realT kappa_r = 0.0;
              realT kappa_s = 0.0;

              PPFS::CalculatePermeabilityAndDerivative( length_r,
                                                        (apertures_np1[r_eff]),
                                                        (apertures_n[r_eff]),
                                                        m_permeabilityUpdateOption,
                                                        edgeLengths_np1[eg],
                                                        m_mu,
                                                        kappa_r,
                                                        dKappa_r_dw );

              PPFS::CalculatePermeabilityAndDerivative( length_s,
                                                        (apertures_np1[s_eff]),
                                                        (apertures_n[s_eff]),
                                                        m_permeabilityUpdateOption,
                                                        edgeLengths_np1[eg],
                                                        m_mu,
                                                        kappa_s,
                                                        dKappa_s_dw );

              kappa_r     *= pow( cos( m_subscaleAngle ), 4 );
              dKappa_r_dw *= pow( cos( m_subscaleAngle ), 4 );
              kappa_s     *= pow( cos( m_subscaleAngle ), 4 );
              dKappa_s_dw *= pow( cos( m_subscaleAngle ), 4 );

              realT kappa = kappa_s * kappa_r / ( kappa_s + kappa_r );

              maxKappa = std::max( kappa, maxKappa );
              minKappa = std::min( kappa, minKappa );

              kappa += m_maxKappa * m_gamma;

              realT invDenom = 1 / ( kappa_r + kappa_s );
              realT dKappa_dwr = ( kappa_s - kappa ) * invDenom * dKappa_r_dw;
              realT dKappa_dws = ( kappa_r - kappa ) * invDenom * dKappa_s_dw;

              edgePermeability[eg] = kappa;

              realT rhoBar_n = 0.5 * ( density_n[s] + density_n[r] );
              realT rhoBar_np1 = 0.5 * ( density_np1[s] + density_np1[r] );

              realT const dRho_dPs = m_IncompressibleFlow == true ? 0.0 : 1.0 /  m_fluidEOS->dPdRho(density_np1[s]);
              realT const dRho_dPr = m_IncompressibleFlow == true ? 0.0 : 1.0 /  m_fluidEOS->dPdRho(density_np1[r]);

              realT dRhoBar_dPr = 0.5 * dRho_dPr;
              realT dRhoBar_dPs = 0.5 * dRho_dPs;
              realT const dRhoBar_dMr = 0.5 * ( 1 / volume_np1[r] );
              realT const dRhoBar_dMs = 0.5 * ( 1 / volume_np1[s] );

              if( density_n[s] < 0.0 || density_n[r] < 0.0 )
              {
                rhoBar_n = m_rho_o;
              }
              if( density_np1[s] < 0.0 || density_np1[r] < 0.0 )
              {
                rhoBar_np1 = m_rho_o;
                dRhoBar_dPr = 0;
                dRhoBar_dPs = 0;
              }

              const realT dp_n   = pressure_n[s] - pressure_n[r] + rhoBar_n * GravityMag;
              const realT dp_np1 = pressure_np1[s] - pressure_np1[r] +  rhoBar_np1 * GravityMag;



              if( m_phiLimiter == 1 )
              {
                if( dp_n * dp_np1 < 0.0 )
                {
                  phiStep[eg] = 1.0;
                }
              }
              realT phi = phiStep[eg];



              realT R_r = 0.0;
              realT R_s = 0.0;
              realT dRr_dVARr = 0.0;
              realT dRr_dVARs = 0.0;
              realT dRs_dVARr = 0.0;
              realT dRs_dVARs = 0.0;

              realT const dPr_dRhor = m_fluidEOS->dPdRho(density_np1[r]);
              realT const dPs_dRhos = m_fluidEOS->dPdRho(density_np1[s]);
              realT const dPr_dMr = dPr_dRhor / volume_np1[r];
              realT const dPs_dMs = dPs_dRhos / volume_np1[s];

              realT dPs_dVARs = 0.0;
              realT dPr_dVARr = 0.0;
              realT dRhoBar_dVARs = 0.0;
              realT dRhoBar_dVARr = 0.0;
              if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::mass )
              {
                dPs_dVARs = dPs_dMs;
                dPr_dVARr = dPr_dMr;
                dRhoBar_dVARs = dRhoBar_dMs;
                dRhoBar_dVARr = dRhoBar_dMr;
              }
              else if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::pressure )
              {
                dPs_dVARs = 1.0;
                dPr_dVARr = 1.0;
                dRhoBar_dVARs = dRhoBar_dPs;
                dRhoBar_dVARr = dRhoBar_dPr;
              }

              // while adding contributions, we are doing this for a face pair
              // (s,r) on the residual for r.s
              // add the contributions to the residual

              if( m_linearizeGoverningEquation )
              {
//                R_r = -( ( 1.0 - phi ) * kappa * rhoBar_n * dp_n
//                               + phi   * kappa * rhoBar_n * dp_n ) * dt ;
                // pressure increment is zero for a one iteration linearization
                // approach...so we drop terms in addition to the
                // the second order terms. actual equation looks like
                R_r = -( ( 1.0 - phi ) * kappa * rhoBar_n * dp_n
                         + phi   * kappa * ( rhoBar_n * dp_np1
                                             + 1/2 * ( dRho_dPs * pressureIncrement[s]
                                                       + dRho_dPr * pressureIncrement[r] ) * dp_n ) ) * dt;

                dRr_dVARr = -phi * kappa * ( rhoBar_n * ( -dPr_dVARr + dRhoBar_dVARr * GravityMag ) + dRhoBar_dVARr * dp_n ) * dt;
                dRr_dVARs = -phi * kappa * ( rhoBar_n * ( +dPs_dVARs + dRhoBar_dVARs * GravityMag ) + dRhoBar_dVARs * dp_n ) * dt;

              }
              else
              {
                R_r = -( ( 1.0 - phi ) * kappa * rhoBar_n * dp_n
                         + phi   * kappa * rhoBar_np1 * dp_np1 ) * dt;
                // dp_np1 = pressure_np1[s] - pressure_np1[r] +  rhoBar_np1 *
                // GravityMag;


                dRr_dVARr = -phi * kappa * ( rhoBar_np1 * ( -dPr_dVARr + dRhoBar_dVARr * GravityMag ) + dRhoBar_dVARr * dp_np1 ) * dt;
                dRr_dVARs = -phi * kappa * ( rhoBar_np1 * ( +dPs_dVARs + dRhoBar_dVARs * GravityMag ) + dRhoBar_dVARs * dp_np1 ) * dt;

              }

              R_s = -R_r;
              dRs_dVARr = -dRr_dVARr;
              dRs_dVARs = -dRr_dVARs;


              flowRHS[kr] += R_r;
              flowRHS[ks] += R_s;

              flowMatrix(kr,kr) += dRr_dVARr;
              flowMatrix(kr,ks) += dRr_dVARs;
              flowMatrix(ks,kr) += dRs_dVARr;
              flowMatrix(ks,ks) += dRs_dVARs;


#if 1
              if( displacementCouplingFlag && m_permeabilityUpdateOption!=0 )
              {

                dKappa_du_rterm  = m_dwdu[r];
                dKappa_du_rterm *= dKappa_dwr * m_dwdw(r);

                dKappa_du_sterm  = m_dwdu[s];
                dKappa_du_sterm *= dKappa_dws * m_dwdw(s);


                array<real64> drhoRdu = m_dwdu(r);
                drhoRdu *= -0.5 * density_np1[r]/volume_np1[r]*faceArea[r] * m_dwdw(r);
                array<real64> drhoSdu = m_dwdu(s);
                drhoSdu *= -0.5 * density_np1[s]/volume_np1[s]*faceArea[s] * m_dwdw(s);

                Epetra_IntSerialDenseVector  displacementDofIndex;
                Epetra_SerialDenseMatrix     matrix_10;

                // calculate displacement DOF for these two faces
                const localIndex rnodes = domain.m_feFaceManager.m_toNodesRelation[r].size();
                const localIndex snodes = domain.m_feFaceManager.m_toNodesRelation[s].size();

                Epetra_IntSerialDenseVector  dispSubFaceDofIndex(2);
                dispSubFaceDofIndex(0) = trilinosIndexFace[r];
                dispSubFaceDofIndex(1) = trilinosIndexFace[s];

                displacementDofIndex.Resize( dim*2*(rnodes + snodes ) );
                matrix_10.Reshape( 2, dim*2*(rnodes + snodes ) );
                matrix_10.Scale(0.0);



                // derivatives wrt diplacements on r-face...loop over the nodes
                // on the r-face
                for( localIndex a=0 ; a<rnodes ; ++a )
                {
                  int rnodeDofIndex[2];
                  faceNodePairIndexing( a, dim, rnodeDofIndex);
                  for( int i=0 ; i<dim ; ++i )
                  {
                    for( int side=0 ; side<2 ; ++side )
                    {
                      const int ai=rnodeDofIndex[side]+i;
                      displacementDofIndex(ai) = m_dwdu_dof(r)(ai);

//                        const realT R_r = -( phi * kappa * dp_np1 ) * dt *
// rhoBar;
//                        const realT R_s =  ( phi * kappa * dp_np1 ) * dt *
// rhoBar;

                      // r-face residual derivative
                      if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::mass )
                      {
                        realT dRrdu = -phi * ( rhoBar_np1*dp_np1*dKappa_du_rterm(ai)
                                               + kappa*rhoBar_np1*(-dPr_dRhor*drhoRdu(ai) )
                                               + kappa*drhoRdu(ai) * dp_np1 ) * dt;

                        realT dRsdu = -dRrdu;
                        matrix_10(0,ai) += dRrdu;
                        matrix_10(1,ai) += dRsdu;
                      }
                      else if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::pressure )
                      {
                        realT dRrdu = -phi * dKappa_du_rterm(ai) * dp_np1 * dt  * rhoBar_np1;
                        realT dRsdu = -dRrdu;
                        matrix_10(0,ai) += dRrdu;
                        matrix_10(1,ai) += dRsdu;
                      }
                    }
                  }
                }

                // loop over the nodes on the s-face
                for( localIndex a=0 ; a<snodes ; ++a )
                {
                  int snodeDofIndex[2];
                  faceNodePairIndexing( a, dim, snodeDofIndex);
                  for( int i=0 ; i<dim ; ++i )
                  {
                    for( int side=0 ; side<2 ; ++side )
                    {
                      const int ai=snodeDofIndex[side]+i;
                      const int aiOffset = ai + rnodes*dim*2;
                      displacementDofIndex(aiOffset) = m_dwdu_dof(s)(ai);

                      // r-face residual derivative
                      if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::mass )
                      {
                        realT dRrdu = -phi * ( rhoBar_np1*dp_np1*dKappa_du_sterm(ai)
                                               + kappa*rhoBar_np1*(dPs_dRhos*drhoSdu(ai) )
                                               + kappa*drhoSdu(ai) * dp_np1 ) * dt;
                        realT dRsdu = -dRrdu;
                        matrix_10(0,ai) += dRrdu;
                        matrix_10(1,ai) += dRsdu;

                      }
                      else if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::pressure )
                      {
                        realT dRrdu = -phi * dKappa_du_sterm(ai) * dp_np1 * dt  * rhoBar_np1;
                        realT dRsdu = -dRrdu;

                        matrix_10(0,aiOffset) += dRrdu;
                        matrix_10(1,aiOffset) += dRsdu;
                      }
                    }
                  }
                }
                matrixPD->SumIntoGlobalValues(dispSubFaceDofIndex, displacementDofIndex, matrix_10);
              }
#endif
            }
          }
        }
      }

      m_matrix->SumIntoGlobalValues(faceDofIndex, flowMatrix);
      m_rhs->SumIntoGlobalValues(faceDofIndex,flowRHS);

    } // ghost edge
  } // edge loop

//  m_maxKappa = maxKappa;
  MPI_Allreduce( &minKappa, &m_minKappa, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
  MPI_Allreduce( &maxKappa, &m_maxKappa, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );


//  ZeroNegativePressures( domain );


  if( m_matrixFlowSolver != nullptr )
  {
    m_matrixFlowSolver->SetMatrixPtr( m_matrix );
    m_matrixFlowSolver->SetRhsPtr( m_rhs );
    m_matrixFlowSolver->SetRowMapPtr( m_rowMap );
    m_matrixFlowSolver->UpdateFluidRockProps(domain);
    m_matrixFlowSolver->Assemble( domain, time, dt );
  }

  if (m_wellboreSolve)
  {
    m_wellboreSolve->m_maxKappa = m_maxKappa;
    m_wellboreSolve->m_minKappa = m_minKappa;
    realT maxWellboreMass = m_wellboreSolve->AssembleFVP( domain, time, dt, m_trilinosIndexStr, m_matrix, m_rhs, m_dofVariable );
    maxMassScale = std::max( maxMassScale, maxWellboreMass );
  }

  m_matrix->GlobalAssemble();
  m_rhs->GlobalAssemble();


  // boundary conditions
  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverFVP::PressureBoundaryCondition,
                                domain, domain.m_feFaceManager, std::string(Field<FieldInfo::pressure>::Name()),
                                time + dt );


  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverFVP::VolumeRateBC,
                                domain, domain.m_feFaceManager,
                                "VolumeRate", time + 0.5*dt, dt );

//  ZeroNegativePressures( domain );


  if( displacementCouplingFlag )
  {

#if USECPP11==1
    const std::shared_ptr<Epetra_Map> dispRowMap = epetraSystem.GetRowMap(EpetraBlock::solidBlock);
    const std::shared_ptr<Epetra_Map> flowRowMap = epetraSystem.GetRowMap(EpetraBlock::fluidBlock);
#else
    const Epetra_Map* const dispRowMap = epetraSystem.GetRowMap(EpetraBlock::solidBlock);
    const Epetra_Map* const flowRowMap = epetraSystem.GetRowMap(EpetraBlock::fluidBlock);
#endif

    matrixPD->GlobalAssemble(*dispRowMap,*flowRowMap);
  }

  return maxMassScale;
}


void ParallelPlateFlowSolverFVP::VolumeRateBC( PhysicalDomainT& domain,
                                               ObjectDataStructureBaseT& object,
                                               BoundaryConditionBase* bc,
                                               const lSet& set,
                                               realT time,
                                               realT dt )
{


  array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
  array<integer>& faceGhostRank  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const array<integer>& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const array<real64>& apertures_n   = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n") );
  const array<real64>& edgeLengths_n = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  const array<real64>& pressure_n   = domain.m_feFaceManager.GetFieldData<realT>("Pressure_n");
//  array<real64>& flowRateBC =
// domain.m_feFaceManager.GetFieldData<realT>("flowRateBC");
//  const array<real64>& density_np1 =
// domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const array<real64>& density_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");

  int const option = bc->getOption();

  realT volRate = bc->GetValue(domain.m_feFaceManager,set.end(),time);



  if( option == 0 || option==1 )
  {
    if( option == 1 )
    {
      int numLocalFaces = 0;
      for( lSet::const_iterator faceID=set.begin() ; faceID!=set.end() ; ++faceID )
      {
        if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]== 1 )
        {
          ++numLocalFaces;
        }
      }

      int numFaces = 0;
      MPI_Allreduce( &numLocalFaces, &numFaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      volRate /= numFaces;
//      std::cout<<"volRate = "<<volRate<<std::endl;
    }

    Epetra_IntSerialDenseVector  face_dof(1);
    Epetra_SerialDenseVector     face_rhs(1);

    for( lSet::const_iterator faceID=set.begin() ; faceID!=set.end() ; ++faceID )
    {

      if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]== 1 )
      {
        realT density = density_n[*faceID];
        if( density < m_rho_o )
          density = m_rho_o;
        face_dof(0) = trilinos_index[*faceID];
        face_rhs(0) = -density * volRate * dt;
        m_rhs->SumIntoGlobalValues(face_dof, face_rhs);
      }
    }
  }
  else if ( option==3 )
  {
    int const numLocalFaces = set.size();
    localIndex numFaces = 0;
    realT kP_sum = 0.0;
    array<real64> k_local(numLocalFaces);

    realT k_min = 1.0e99;
    for( lSet::const_iterator faceID=set.begin() ; faceID!=set.end() ; ++faceID )
    {
      if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]== 1 )
      {

        realT length = 0.0;
        for( localIndex ke=0 ; ke<domain.m_feFaceManager.m_toEdgesRelation[*faceID].size() ; ++ke )
        {
          length = std::max(length,edgeLengths_n[domain.m_feFaceManager.m_toEdgesRelation[*faceID][ke]]);
        }

        realT kappa = 0.0;
        realT junk = 0.0;
        PPFS::CalculatePermeabilityAndDerivative( length*0.5,
                                                  apertures_n[*faceID],
                                                  apertures_n[*faceID],
                                                  3,
                                                  length,
                                                  m_mu,
                                                  kappa,
                                                  junk );

        k_local[numFaces] = kappa;
        k_min = std::min( k_min, kappa );
        ++numFaces;
      }
    }

    {
      realT buffer[1] = {  k_min };
      realT rbuffer[1] = { 0.0 };
      MPI_Allreduce(buffer, rbuffer, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      k_min = rbuffer[0];
    }

    int kf=0;
    realT k_sum = 0.0;
    for( lSet::const_iterator faceID=set.begin() ; faceID!=set.end() ; ++faceID )
    {
      if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]== 1 )
      {
        k_local[kf] = k_local[kf] * (0.10*k_min) / ( (0.10*k_min) + k_local[kf] );
        k_sum += k_local[kf];
        kP_sum += k_local[kf]  * pressure_n[*faceID];
        ++kf;
      }
    }

    {
      realT buffer[3] = { static_cast<realT>(numFaces), k_sum, kP_sum };
      realT rbuffer[3] = { 0.0, 0.0, 0.0 };
      MPI_Allreduce(buffer, rbuffer, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      numFaces = static_cast<localIndex>(rbuffer[0]);
      k_sum = rbuffer[1];
      kP_sum = rbuffer[2];
    }

    Epetra_IntSerialDenseVector  face_dof(numFaces);
    Epetra_SerialDenseVector     face_rhs(numFaces);

    realT p_e = ( kP_sum + volRate ) / k_sum;
    kf=0;
    realT q_sum = 0.0;
    for( lSet::const_iterator faceID=set.begin() ; faceID!=set.end() ; ++faceID )
    {
      if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]== 1 )
      {
        face_dof(kf) = trilinos_index[*faceID];

        realT q = k_local[kf] * ( p_e - pressure_n[*faceID] );

        face_rhs(kf) = q * dt;
        q_sum += q;
        ++kf;
      }
    }
    m_rhs->SumIntoGlobalValues(face_dof, face_rhs);
  }
  else
  {
    throw GPException("ParallelPlateFlowSolverFVP::VolumeRateBC. invalid bc option specified.");
  }
}


void ParallelPlateFlowSolverFVP::ZeroNegativePressures( PhysicalDomainT& domain )
{

  realT LARGE;

  array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
  array<integer>& faceGhostRank  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const array<integer>& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const array<real64>& pressure_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  bool displacementCouplingFlag = false;

  if( this->m_system != nullptr )
  {
    displacementCouplingFlag = this->m_system->HasMatrixBlock( EpetraBlock::fluidBlock, EpetraBlock::solidBlock );
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FECrsMatrix> matrixPD;
#else
  Epetra_FECrsMatrix* matrixPD = nullptr;
#endif

  if( displacementCouplingFlag )
  {
    matrixPD = this->m_system->GetMatrix( EpetraBlock::fluidBlock, EpetraBlock::solidBlock );
  }

  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);

  for( lSet::const_iterator faceID=m_negativePressureSet.begin() ; faceID!=m_negativePressureSet.end() ; ++faceID )
  {
    face_dof(0) = trilinos_index[*faceID];

    if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]== 1 )
    {
#if USECPP11==1
      LARGE = ClearRow( m_matrix.get(), face_dof(0), 1.0 );
      if( displacementCouplingFlag )
      {
        ClearRow( matrixPD.get(), face_dof(0),0.0);
      }
#else
      LARGE = ClearRow( m_matrix, face_dof(0), 1.0 );
      if( displacementCouplingFlag )
      {
        ClearRow( matrixPD, face_dof(0),0.0);
      }
#endif


      face_rhs(0) = LARGE*( this->m_cutoffPressure - pressure_np1[*faceID] );
      m_rhs->ReplaceGlobalValues(face_dof, face_rhs);
    }
  }
}


void ParallelPlateFlowSolverFVP::PressureBoundaryCondition( PhysicalDomainT& domain,
                                                            ObjectDataStructureBaseT& object,
                                                            BoundaryConditionBase* bc,
                                                            const lSet& set,
                                                            realT time)
{

  realT LARGE;

  array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
  array<integer>& faceGhostRank  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const array<integer>& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const array<real64>& pressure_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  bool displacementCouplingFlag = false;

  if( this->m_system != nullptr )
  {
    displacementCouplingFlag = this->m_system->HasMatrixBlock( EpetraBlock::fluidBlock, EpetraBlock::solidBlock );
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FECrsMatrix> matrixPD;
#else
  Epetra_FECrsMatrix* matrixPD = nullptr;
#endif

  if( displacementCouplingFlag )
  {
    matrixPD = this->m_system->GetMatrix( EpetraBlock::fluidBlock, EpetraBlock::solidBlock );
  }

  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);

  for( lSet::const_iterator faceID=set.begin() ; faceID!=set.end() ; ++faceID )
  {
    face_dof(0) = trilinos_index[*faceID];

    if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]== 1 )
    {
#if USECPP11==1
      LARGE = ClearRow( m_matrix.get(), face_dof(0), 1.0 );
      if( displacementCouplingFlag )
      {
        ClearRow( matrixPD.get(), face_dof(0),0.0);
      }
#else
      LARGE = ClearRow( m_matrix, face_dof(0), 1.0 );
      if( displacementCouplingFlag )
      {
        ClearRow( matrixPD, face_dof(0),0.0);
      }
#endif


      face_rhs(0) = -LARGE*( bc->GetValue(domain.m_feFaceManager,faceID,time) - pressure_np1[*faceID] );
//      std::cout<<"bc->GetValue(domain.m_feFaceManager,"<<*faceID<<","<<time<<")
// - pressure_np1["<<*faceID<<"] =
// "<<bc->GetValue(domain.m_feFaceManager,faceID,time)<<" -
// "<<pressure_np1[*faceID]<<std::endl;
//      std::cout<<LARGE<<std::endl;
//      std::cout<<face_rhs(0)<<std::endl;

      m_rhs->ReplaceGlobalValues(face_dof, face_rhs);
    }
  }
}


/*
   realT ParallelPlateFlowSolverFVP::TwoFacePermeability(const array<R1Tensor>&
      edgeCenters,
                                                           const array<real64>&
                                                              edgeLengths,
                                                           const
                                                              array<R1Tensor>&
                                                              faceCenters,
                                                           const array<real64>&
                                                              apertures,
                                                           const localIndex eg,
                                                           const localIndex r,
                                                           const localIndex s,
                                                           const
                                                              array<array<real64>>*
                                                              const dwdu,
                                                           array<real64>* const
                                                              dkdu_r,
                                                           array<real64>* const
                                                              dkdu_s)
   {

   R1Tensor edgeCenter = edgeCenters[eg];
   R1Tensor lr, ls;

   lr = edgeCenter;
   lr -= faceCenters[r];

   ls = edgeCenter;
   ls -= faceCenters[s];

   realT l_edge = edgeLengths[eg];

   realT wr = BoundedAperture(apertures[r]);
   realT ws = BoundedAperture(apertures[s]);

   realT kappa = CalculatePermeability(lr.L2_Norm(), ls.L2_Norm(), wr, ws,
      l_edge, m_mu, m_SHP_FCT);


   if( dkdu_r!=NULL && dkdu_s!=NULL && dwdu!=NULL )
   {
    if( dwdu->size() > std::max(r,s) )
    {
      const realT denom = lr.L2_Norm()*ws*ws*ws+ls.L2_Norm()*wr*wr*wr;

 * dkdu_r  = (*dwdu)[r] ;
 * dkdu_r *= lr.L2_Norm()*ws*ws*ws/wr;

 * dkdu_s  = (*dwdu)[s] ;
 * dkdu_s *= ls.L2_Norm()*wr*wr*wr/ws;

 * dkdu_r *= 3 * kappa / denom  * BoundedApertureDerivative(apertures[r]);
 * dkdu_s *= 3 * kappa / denom  * BoundedApertureDerivative(apertures[s]);
    }
   }

   return kappa;
   }

 */


realT ParallelPlateFlowSolverFVP::CheckSolution( const realT* const local_solution,
                                                 const PhysicalDomainT& domain,
                                                 const localIndex dofOffset )
{
  realT rval = 1.0;


  const array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
  const array<integer>& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  array<real64> const & pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      const int intDofOffset = 0; //dofOffset;
      int lid = intDofOffset+m_rowMap->LID(trilinos_index[*kf]);

      realT dP = local_solution[lid];

      if( pressure[*kf] * (pressure[*kf] - dP) < 0.0 )
      {
        realT const e = m_negativeFluidPressureSlopeReduction;
//        dP = ( 1.0 - e )*pressure[*kf] + dP * e;

        realT factor = ( 1.0 - e )*pressure[*kf] / dP + e;
        rval = std::min(rval,factor);
//        std::cout<<"scale = "<<rval<<std::endl;
      }

/*
      realT newRho = this->m_fluidEOS->density( pressure[*kf] - dP );
      if( newRho < 0 )
      {
        realT oldRho = this->m_fluidEOS->density( pressure[*kf] );
        realT minP = this->m_fluidEOS->pressure( oldRho * 0.01 );
        realT factor = ( pressure[*kf] - minP ) / dP;
        rval = std::min(rval,factor);

        std::cout<<"P      = "<<pressure[*kf]<<std::endl;
        std::cout<<"dP     = "<<dP<<std::endl;
        std::cout<<"P-dP   = "<<pressure[*kf] - dP<<std::endl;
        std::cout<<"newRho = "<<newRho<<std::endl;
        std::cout<<"minP   = "<<minP<<std::endl;
      }
 */
    }
  }


#if 0
  // face fields
  const array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
  const array<integer>& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  array<real64> const & pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      const int intDofOffset = dofOffset;
      int lid = m_rowMap->LID(intDofOffset+trilinos_index[*kf]);

      const realT pInc = local_solution[lid];

      if( pressure[*kf] + pInc < m_cutoffPressure * (1.0 + 1.0e-14) )
      {
        double ratio = ( m_cutoffPressure - pressure[*kf] ) / pInc;
        if( ratio < rval )
        {
          rval = ratio;
          this->m_negativePressureSet.insert(*kf);
        }
      }
    }
  }
#endif

  realT localVal = rval;
  MPI_Allreduce( &localVal, &rval,1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

  return rval;
}



void ParallelPlateFlowSolverFVP::PropagateSolution( const realT* const local_solution,
                                                    const realT scalingFactor,
                                                    PhysicalDomainT& domain,
                                                    const localIndex dofOffset  )
{
  // face fields
  const array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
  const array<integer>& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  array<real64>& rho = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();


  if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::mass )
  {
    array<real64>& massIncrement = domain.m_feFaceManager.GetFieldData<realT>("massIncrement");
    array<real64>& mass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
    for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
    {
      if(is_ghost[*kf] < 0)
      {
        const int intDofOffset = dofOffset;
        int lid = m_rowMap->LID(intDofOffset+trilinos_index[*kf]);
        massIncrement[*kf] -= scalingFactor*local_solution[lid];
        mass[*kf] -= scalingFactor*local_solution[lid];
      }
    }
  }
  else if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::pressure )
  {
    array<real64>& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
    array<real64>& pressureIncrement = domain.m_feFaceManager.GetFieldData<realT>("pressureIncrement");

    for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
    {
      if(is_ghost[*kf] < 0)
      {
        const int intDofOffset = 0;//dofOffset;
        int lid = intDofOffset+m_rowMap->LID(trilinos_index[*kf]);

        realT dP = scalingFactor*local_solution[lid];

//      printf("pressure[%4lu] = %e\n",*kf,pressure[*kf]);
//      printf("               = %e\n",pressure[*kf]-dP);
        if( pressure[*kf] * (pressure[*kf] - dP) < 0.0 )
        {
          realT const e = m_negativeFluidPressureSlopeReduction;
          dP = ( 1.0 - e )*pressure[*kf] + dP * e;
        }

        pressure[*kf] -= dP;
        pressureIncrement[*kf]  -= dP;
//      printf("               = %f\n",pressure[*kf]);

//      pressure[*kf] -= scalingFactor*local_solution[lid];
//      pressureIncrement[*kf] -= scalingFactor*local_solution[lid];

        if( this->m_IncompressibleFlow )
        {
          rho[*kf] = m_rho_o;
        }
        else
        {
          rho[*kf] = m_fluidEOS->density( pressure[*kf] );
        }
      }
    }
  }


  if( m_matrixFlowSolver != nullptr )
  {
    m_matrixFlowSolver->PropagateSolution( local_solution, scalingFactor, domain, 0 );
  }
  if (m_wellboreSolve)
  {
    m_wellboreSolve->PropagateSolutionFVP(domain, m_rowMap, local_solution, m_dofVariable, scalingFactor);
  }
}

void ParallelPlateFlowSolverFVP::PostSyncConsistency( PhysicalDomainT& domain,
                                                      SpatialPartition& partition )
{
//  GenerateParallelPlateGeometricQuantities( domain, time,dt );
  UpdateEOS( 0.0, 0.0, domain );
}



void ParallelPlateFlowSolverFVP::InitializeCommunications( PartitionBase& partition )
{
  if( m_solvers.find(m_matrixFlowSolverName) != m_solvers.end() )
  {
    m_matrixFlowSolver = dynamic_cast<FractureMatrixFlowSolverFV*>(stlMapLookup( m_solvers, m_matrixFlowSolverName, "m_solvers" ));
  }

  m_syncedFields.clear();
  m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("massRate");
  m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("Pressure");
//  m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("Density");
  m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(m_trilinosIndexStr);
  m_syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back(VolumetricFluxStr);

  if( m_matrixFlowSolver!=nullptr )
  {
    m_matrixFlowSolver->InitializeCommunications(partition);
    for( auto iter : m_matrixFlowSolver->m_syncedFields )
    {
      PhysicalDomainT::ObjectDataStructureKeys key = iter.first;
      for( auto fieldname : iter.second )
      {
        m_syncedFields[key].push_back(fieldname);
      }
    }
  }

  partition.SetBufferSizes(m_syncedFields, CommRegistry::parallelPlateFlowSolver);
}

void ParallelPlateFlowSolverFVP::TimeStepSetup( const realT& time,
                                                const realT& dt,
                                                PhysicalDomainT& domain,
                                                SpatialPartition& partition,
                                                const bool setupSystem )
{
  DefineFlowSets( domain );


  RegisterTemporaryFields( domain );

  GenerateParallelPlateGeometricQuantities( domain, time,dt );
  UpdateEOS( time, dt, domain );

  FillTemporaryFields( domain );

  if( setupSystem )
  {
    SetupSystem (domain,partition);
  }

  array<real64>& pressureIncrement = domain.m_feFaceManager.GetFieldData<realT>("pressureIncrement");
  array<real64>& massIncrement   = domain.m_feFaceManager.GetFieldData<realT>("massIncrement");
  pressureIncrement = 0.0;


  if( this->m_matrixFlowSolver!=nullptr )
  {
    m_matrixFlowSolver->TimeStepSetup( time, dt, domain, partition, 0 );
  }
  massIncrement = 0.0;
}

double ParallelPlateFlowSolverFVP::TimeStepExecute( const realT& time,
                                                    const realT& dt,
                                                    PhysicalDomainT& domain,
                                                    SpatialPartition& partition )
{
  realT dt_return = dt;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for( int i=0 ; i<m_numerics.m_maxIterNewton ; ++i )
  {
    m_rhs->Scale(0.0);
    m_matrix->Scale(0.0);
    m_solution->Scale(0.0);

    GenerateParallelPlateGeometricQuantities( domain, time, dt );
    UpdateEOS( time, dt, domain );

    Epetra_System junk;
    realT massScale = Assemble(domain, junk, time, dt );

    if( this->m_numerics.m_verbose >=2 )
    {
      std::cout<<"MATRIX"<<std::endl;
      m_matrix->Print(std::cout);
      std::cout<<"RHS"<<std::endl;
      m_rhs->Print(std::cout);
    }

    realT residual = 0;
    m_rhs->Norm2( &residual);

    realT localMaxData = massScale;
    MPI_Allreduce (&localMaxData,&massScale,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);


    if( rank==0 && this->m_numerics.m_verbose >=1 )
    {
      std::cout<<"iter, residual/massScale = "<<i<<", "<<residual<<"/"<<massScale<<" = "<<residual/massScale<<std::endl;
    }
    if( residual/massScale < this->m_numerics.newton_tol )
      break;

    Solve(domain,partition, time, dt );
    UpdateEOS( time, dt, domain );
//    m_solution->Print(std::cout);

    if( this->m_numerics.m_verbose >=2 )
    {
      std::cout<<"SOLN"<<std::endl;
      m_solution->Print(std::cout);
      std::cout<<std::endl<<std::endl<<std::endl;
    }
  }
  if( this->m_numerics.m_verbose >=2 )
  {
    std::cout<<std::endl<<std::endl<<std::endl;
    std::cout<<std::endl<<std::endl<<std::endl;
  }

  UpdateFlux( time, dt, domain,partition );

  return dt_return;
}


void ParallelPlateFlowSolverFVP::TimeStepCleanup( PhysicalDomainT& domain, const realT& dt )
{
  const array<real64>& mass_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");
  array<real64>& extraMass     = domain.m_feFaceManager.GetFieldData<realT>("extraMass");
  const array<real64>& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");
//  const array<real64>& pressure_n   =
// domain.m_feFaceManager.GetFieldData<realT>("Pressure_n");
  array<integer>& face_is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const array<real64>& volume_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidVolume_n");
//  const array<real64>& massIncrement   =
// domain.m_feFaceManager.GetFieldData<realT>("massIncrement");
//  const array<real64>& volume_np1 =
// domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    localIndex const r = *kf;
    if( face_is_ghost[r] < 0 )
    {
      if( extraMass[r] > 0.0  && m_leakExtraMass > 0 )
      {
        realT massLoss = extraMassLoss( extraMass[r],
                                        mass_n[r],
                                        volume_n[r],
                                        faceArea[r] * ( BoundedAperture(m_zeroApertureOffset) + m_zeroApertureVolume ) * m_rho_o,
                                        dt, 1.0 );
        extraMass[r] -= massLoss;
      }
    }
  }
}



double ParallelPlateFlowSolverFVP::TimeStep( const realT& time,
                                             const realT& dt,
                                             const int cycleNumber,
                                             PhysicalDomainT& domain,
                                             const array<string>& namesOfSolverRegions,
                                             SpatialPartition& partition,
                                             FractunatorBase* const fractunator )
{
  realT dt_return = dt;

  RegisterTemporaryFields( domain );

  if( m_matrix==NULL )
  {
    TimeStepSetup( time, dt, domain, partition, true );
  }
  else
  {
    TimeStepSetup( time, dt, domain, partition, false );
  }

  dt_return = TimeStepExecute( time, dt, domain, partition );

  TimeStepCleanup( domain, dt );
  DeregisterTemporaryFields( domain );

  return dt_return;
}

void ParallelPlateFlowSolverFVP::DefineFlowSets( PhysicalDomainT& domain )
{
  FaceManagerT& faceManager = domain.m_feFaceManager;
  EdgeManagerT& edgeManager = domain.m_feEdgeManager;

  array<integer>& flowFaceType = faceManager.GetFieldData<int>("flowFaceType");
  array<integer>& flowEdgeType = edgeManager.GetFieldData<int>("flowEdgeType");

  if( m_flowFaceSetName.empty() )
  {
    const array<lSet>& edgeToFlowFaces = edgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");

    m_faceSet.clear();

    for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
    {
      if( flowFaceType[kf] == 1 )
      {
        m_faceSet.insert( kf );
      }
    }

    m_numFaces = m_faceSet.size();

    m_faceDofMap.clear();
    lSet::const_iterator si=m_faceSet.begin();
    for(localIndex i =0 ; i < m_numFaces ; ++i, ++si)
    {
      localIndex f = *si;
      m_faceDofMap[f] = i;
    }

    m_edgesToFaces.clear();
    array<integer>& ffCount = edgeManager.GetFieldData<int>("FlowFaceCount"); // debug

    for( localIndex ke=0 ; ke < edgeManager.DataLengths() ; ++ke )
    {
      if( flowEdgeType[ke] == 1 && !(edgeToFlowFaces[ke].empty()) )
      {
        m_edgesToFaces[ke].assign( edgeToFlowFaces[ke].begin(), edgeToFlowFaces[ke].end() );
        ffCount[ke] = m_edgesToFaces[ke].size();

      }
    }

  }
  else
  {
    m_faceSet = faceManager.GetSet(m_flowFaceSetName);
    m_numFaces = m_faceSet.size();

    // build face-dof map

    lSet::const_iterator si=m_faceSet.begin();
    for(localIndex i =0 ; i < m_numFaces ; ++i, ++si)
    {
      localIndex f = *si;
      m_faceDofMap[f] = i;
    }

    array<integer>& ffCount = domain.m_feEdgeManager.GetFieldData<int>("FlowFaceCount"); // debug


    for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
    {
      flowFaceType[*kf] = 1;
      localIndex numEdges = faceManager.m_toEdgesRelation[*kf].size();
      for(localIndex a =0 ; a < numEdges ; ++a)
      {
        localIndex eg = faceManager.m_toEdgesRelation[*kf][a];

        flowEdgeType[eg] = 1;
        lSet& edgeFaces = edgeManager.m_toFacesRelation[eg];
        lArray1d edgeList;

        for( lSet::iterator edgeFace=edgeFaces.begin() ; edgeFace!=edgeFaces.end() ; ++edgeFace )
        {
          if(isMember(*edgeFace,m_faceDofMap))
          {
            edgeList.push_back(*edgeFace);
          }
        }
        m_edgesToFaces[eg] = edgeList;
        ffCount[eg] = edgeList.size();
      }
    }
  }
}


void ParallelPlateFlowSolverFVP::GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,realT time,realT dt )
{

  array<R1Tensor>& faceCenter  = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );


  array<real64>& aperture_n   = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n") );
  array<real64>& volume_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidVolume_n");


  array<real64>& aperture = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr );
  array<real64>& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  array<real64>& lastAper = domain.m_feFaceManager.GetFieldData<realT>("lastAper" );

  // array<real64>& mass_n   =
  // domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");
  // array<real64>& mass     =
  // domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  // array<real64>& extraMass     =
  // domain.m_feFaceManager.GetFieldData<realT>("extraMass");

  // array<real64>& density_np1 =
  // domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  // array<real64>& density_n   =
  // domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");

//  array<real64>& boundedAperture = domain.m_feFaceManager.GetFieldData<realT>(
// "BoundedAperture" );

  const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
  array<R1Tensor> const & nodalDisp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
  array<R1Tensor> const & incDisp = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();

  const array<real64>& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");

  array<real64>& apertureIncrement = domain.m_feFaceManager.GetFieldData<realT>("apertureIncrement");
  // array<real64>& pressure =
  // domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  array<real64> const * const contactOffset = domain.m_feFaceManager.GetFieldDataPointer<realT>("contactOffset" );

  m_dwdw.resize( domain.m_feFaceManager.DataLengths() );
  m_dwdw = 1.0;

  // update face quantities
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {

    domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *kf, faceCenter[*kf]);

    if(m_doApertureUpdate)
    {

      R1Tensor Nbar;
      if( childFaceIndex[*kf].size() > 0 )
      {
        const localIndex faceIndex[2] = { *kf, childFaceIndex[*kf][0] };
        const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                                domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

        Nbar = N[0];
        Nbar -= N[1];
        Nbar.Normalize();



        R1Tensor w;
        R1Tensor incw;
        const lArray1d& nodeList= domain.m_feFaceManager.m_toNodesRelation[*kf];
        const localIndex numNodes = nodeList.size();

        for( localIndex a=0 ; a<numNodes ; ++a )
        {
          const localIndex aa = a == 0 ? a : numNodes - a;
          const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[*kf][a];
          const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[childFaceIndex[*kf][0]][aa];

          w += nodalDisp[localNodeIndex2];
          w -= nodalDisp[localNodeIndex1];

          incw += incDisp[localNodeIndex2];
          incw -= incDisp[localNodeIndex1];
        }
        w    /= numNodes;
        incw /= numNodes;

        aperture[*kf] = Dot(w,Nbar) + m_zeroApertureOffset;
        apertureIncrement[*kf] = Dot(incw,Nbar);
      }
      else
      {
        R1Tensor gap;
        Nbar = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *kf );
        gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, *kf );
        aperture[*kf] = Dot(gap,Nbar) + m_zeroApertureOffset;
      }



      if( contactOffset!=nullptr )
      {
        if( (*contactOffset)[*kf] > m_contactOffsetCutoff )
        {
          aperture[*kf] += (*contactOffset)[*kf];
        }
        else
        {
          aperture[*kf] += m_contactOffsetCutoff;
        }
      }

      if( this->m_boundPhysicalAperture )
      {

        m_dwdw(*kf) = BoundedApertureDerivative( aperture[*kf] );
//        std::cout<<"m_dwdw(*kf)   = "<<m_dwdw(*kf)<<std::endl;
//        std::cout<<"aperture[*kf] = "<<aperture[*kf]<<std::endl;
        aperture[*kf] = BoundedAperture(aperture[*kf]);
//        std::cout<<"baperture[*kf] = "<<aperture[*kf]<<std::endl;
      }
      apertureIncrement[*kf] = aperture[*kf] - aperture_n[*kf];

    }


    fluidVolume[*kf] = ( aperture[*kf] + m_zeroApertureVolume) * faceArea[*kf];

    if( aperture_n[*kf] == 0.0 )
    {
      aperture_n[*kf] = aperture[*kf];
    }

    if( volume_n[*kf] == 0.0 )
    {
      volume_n[*kf] = fluidVolume[*kf];
    }


    if( fluidVolume[*kf]<=0.0 )
      throw GPException("you have a negative volume finite volume element!!");

/*
    if( isZero( mass[*kf] ) )
    {
      density_np1[*kf] = 1.0 * m_rho_o;
      density_n[*kf]   = 1.0 * m_rho_o;
      pressure[*kf]    = m_fluidEOS->pressure(density_np1[*kf]);
      mass[*kf]        = fluidVolume[*kf] * density_np1[*kf];
      if( isZero( mass_n[*kf] ) )
      {
        mass_n[*kf] = mass[*kf];
      }
      extraMass[*kf] = mass[*kf];
    }

    if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::mass )
    {
    }
    else if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::pressure
       )
    {
   ///      mass[*kf] = fluidVolume[*kf] * density_np1[*kf];
    }
 */

  }
  lastAper = aperture;
  lastAper += m_zeroApertureVolume;



  // update edge properties
//  array<integer>& edge_is_ghost            =
// domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  array<R1Tensor>& edgeCenter_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  array<realT>& edgeLength_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin() ; itr!=itrEnd ; ++itr )
  {
    localIndex eg = itr->first;
    domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg, edgeCenter_new[eg] );
    edgeLength_new[eg] = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, eg);
  }
}


void ParallelPlateFlowSolverFVP::CalculateMassRate( PhysicalDomainT& domain,
                                                    SpatialPartition& partition,
                                                    realT time, realT dt )
{}

// set initial fluid density based on known pressure;
/*
   void ParallelPlateFlowSolverFVP::InitializeDensity( PhysicalDomainT& domain)
   {
   array<real64>& mass     =
      domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
   array<real64>& density  =
      domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
   const array<real64>& pressure =
      domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
   const array<real64>& fluidVolume  =
      domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

   for( lSet::const_iterator fc=m_faceSet.begin() ; fc!=m_faceSet.end() ; ++fc )
      {
    density[*fc] = rho_EOS(pressure[*fc],m_bulk_modulus,m_rho_o );
    mass[*fc] = density[*fc]*fluidVolume[*fc];
   }
   }*/


void ParallelPlateFlowSolverFVP::UpdateEOS( const realT time, const realT dt, PhysicalDomainT& domain )
{
//  const array<integer>& ghostRank       =
// domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  array<real64>& density = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  array<real64>& density_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");
  array<real64>& mass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  array<real64>& mass_n = domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");

  array<real64>& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
//  array<real64>& pressureIncrement =
// domain.m_feFaceManager.GetFieldData<realT>("pressureIncrement");
  const array<real64>& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  array<real64>& dPdM = domain.m_feFaceManager.GetFieldData<realT>("dPdM");

  array<integer>& isInitialized = domain.m_feFaceManager.GetFieldData<int>("isInitialized");
  array<real64>& extraMass     = domain.m_feFaceManager.GetFieldData<realT>("extraMass");


  if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::pressure )
  {
    for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
    {
      if( isInitialized[*kf] == 0 )
      {
        isInitialized[*kf] = 1;
        density[*kf]     = m_newVolumeDensityRatio * m_rho_o;
        density_n[*kf]   = m_newVolumeDensityRatio * m_rho_o;
        pressure[*kf]    = m_fluidEOS->pressure(density[*kf]);
        mass[*kf]        = fluidVolume[*kf] * density[*kf];
        mass_n[*kf] = mass[*kf];
        extraMass[*kf] = mass[*kf];
      }
      else
      {
        if( this->m_IncompressibleFlow )
        {
          density[*kf] = m_rho_o;
        }
        else
        {
          density[*kf] = m_fluidEOS->density( pressure[*kf] );
        }
        mass[*kf] = density[*kf] * fluidVolume[*kf];
      }
    }
  }
  else if( m_dofVariable == ParallelPlateFlowSolverBase::dofVariable::mass )
  {
    for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
    {
      density[*kf] = mass[*kf] / fluidVolume[*kf];

      pressure[*kf] = m_fluidEOS->pressure(density[*kf]);
      dPdM[*kf] = m_fluidEOS->dPdRho(density[*kf]) / fluidVolume[*kf];


      // propagate pressure to children
      lArray1d& childFaces = domain.m_feFaceManager.m_childIndices[*kf];
      for(unsigned i =0 ; i < childFaces.size() ; ++i)
      {
        pressure[childFaces[i]] =  pressure[*kf];
      }
    }
  }

  if (m_wellboreSolve)
  {
    m_wellboreSolve->UpdateEOSFVP( domain, time, dt, m_dofVariable );
  }

  array<real64>* initialSaturatedTime = domain.m_feFaceManager.GetFieldDataPointer<realT>("initialSaturatedTime");
  if (initialSaturatedTime != NULL)
  {
    MarkFaceSaturationTime(domain, time, dt);
  }
}


void ParallelPlateFlowSolverFVP::UpdateFlux( const realT time,
                                             const realT dt,
                                             PhysicalDomainT& domain,
                                             SpatialPartition& partition )
{


  const array<real64>& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  array<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  volFlux = 0.0;

  const array<real64>& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  array<R1Tensor>& flowVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>( FluidVelocityStr );
  flowVelocity = 0.0;

  const array<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  array<integer>& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  const array<real64>& density_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();

  // loop over edges
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin() ; itr!=itrEnd ; ++itr )
  {
    const int numFaces = itr->second.size();
    if( numFaces > 1)
    {
      localIndex eg = itr->first;
      localIndex kfa = itr->second[0];
      localIndex kfb = itr->second[1];

      realT Pa = pressures[kfa];
      realT Pb = pressures[kfb];

      R1Tensor rs_xdiff;
      rs_xdiff = faceCenters[kfb];
      rs_xdiff -= faceCenters[kfa];
      realT const GravityMag = Dot(rs_xdiff,this->m_gravityVector);


      // will be incorrect at junctions of 3 or more faces
      volFlux[eg] =  edgePermeabilities[eg]*(Pa-Pb + 0.5 * (density_np1[kfa]+density_np1[kfb]) * GravityMag ); // flux
                                                                                                               // A
                                                                                                               // ->
                                                                                                               // B

      R1Tensor vecFlux;

      vecFlux = faceCenters[kfb];
      vecFlux -= faceCenters[kfa];

      vecFlux.Normalize();
      vecFlux *= volFlux[eg];

      if(is_ghost[kfa] < 0)
        flowVelocity[kfa] += vecFlux;
      if(is_ghost[kfb] < 0)
        flowVelocity[kfb] += vecFlux;

    }
  }

//  const array<real64>& faceArea =
// domain.m_feFaceManager.GetFieldData<realT>("faceArea");
//  const array<real64>& apertures_np1 =
// domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  //const array<real64>& apertures_n   =
  // domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n")
  // );
  for( auto kf : this->m_faceSet )
  {
    if( is_ghost[kf] < 0 )
    {
      flowVelocity[kf] *= 0.5;
    }
  }
}



void ParallelPlateFlowSolverFVP::CalculateApertureDerivatives( const FaceManagerT& faceManager,
                                                               const NodeManager& nodeManager )
{
  const array<integer>* const trilinosIndexNode = nodeManager.GetFieldDataPointer<int>("IMS_0_GlobalDof");

  if( trilinosIndexNode!=NULL )
  {
    m_dwdu.resize(faceManager.DataLengths());
    m_dwdu_dof.resize( faceManager.DataLengths() );

    const int dim=3;

    const array<integer>& flowFaceType = faceManager.GetFieldData<int>("flowFaceType");
//    const array<real64>& apertures_np1 = faceManager.GetFieldData<realT>(
// ApertureStr );
    const OrderedVariableOneToManyRelation& childFaceIndex = faceManager.GetVariableOneToManyMap( "childIndices" );

    // set aperture derivatives
    for( localIndex r=0 ; r<faceManager.DataLengths() ; ++r )
    {
      if( /*face_is_ghost[r] < 0 &&*/ flowFaceType[r] == 1 )
      {
        m_dwdu(r).resize( faceManager.m_toNodesRelation[r].size() * dim * 2 );
        m_dwdu_dof(r).resize( faceManager.m_toNodesRelation[r].size() * dim * 2 );

        const localIndex numNodes = faceManager.m_toNodesRelation[r].size();
        for( localIndex a=0 ; a<numNodes ; ++a )
        {
          const localIndex faceIndex[2] = { r, childFaceIndex[r][0] };

          const R1Tensor N[2] = { faceManager.FaceNormal( nodeManager, faceIndex[0] ),
                                  faceManager.FaceNormal( nodeManager, faceIndex[1] )};

          R1Tensor Nbar = N[0];
          Nbar -= N[1];
          Nbar.Normalize();

          const localIndex aa = a == 0 ? a : numNodes - a;
          const localIndex node0 = faceManager.m_toNodesRelation[faceIndex[0]][a];
          const localIndex node1 = faceManager.m_toNodesRelation[faceIndex[1]][aa];


          int nodeDofIndex[2];
          faceNodePairIndexing( a, dim, nodeDofIndex);

          for( int i=0 ; i<dim ; ++i )
          {
            m_dwdu(r)(nodeDofIndex[0]+i) = -Nbar[i]/faceManager.m_toNodesRelation[r].size();
            m_dwdu(r)(nodeDofIndex[1]+i) =  Nbar[i]/faceManager.m_toNodesRelation[r].size();

            m_dwdu_dof(r)(nodeDofIndex[0]+i) = dim*(*trilinosIndexNode)[node0]+i;
            m_dwdu_dof(r)(nodeDofIndex[1]+i) = dim*(*trilinosIndexNode)[node1]+i;
          }
        }
      }
    }
  }
}


void ParallelPlateFlowSolverFVP::AssembleLeakoff    (PhysicalDomainT& domain, const realT& time, const realT& dt)
{
  // leakoff
  array<real64>& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  array<integer>& face_is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  array<real64>& faceToMatrixLeakOffRate = domain.m_feFaceManager.GetFieldData<realT>("faceToMatrixLeakOffRate");
  array<real64> &totalLeakedThickness = domain.m_feFaceManager.GetFieldData<realT>("totalLeakedThickness");
  const array<integer>& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const array<real64>&  faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea" );
  array<real64>& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);

  if (m_leakoffModel.GetLeakoffCoefficient() > 0)
  {
    array<real64>& initialSaturatedTime = domain.m_feFaceManager.GetFieldData<realT>("initialSaturatedTime");
    for( localIndex k=0 ; k<domain.m_feFaceManager.DataLengths() ; ++k )
    {
      if (flowFaceType[k] == 1)
      {
        // calculate leakoff
        realT P = faceFluidPressure[k];
        if (P>0.1)
        {
          faceToMatrixLeakOffRate[k] = m_leakoffModel.CalculateFlux(P,time+dt/2,initialSaturatedTime[k],this->m_mu); // flux
                                                                                                                     // out.
                                                                                                                     // Don't
                                                                                                                     // leak
                                                                                                                     // if
                                                                                                                     // pressure
                                                                                                                     // is
                                                                                                                     // too
                                                                                                                     // low.
          faceToMatrixLeakOffRate[domain.m_feFaceManager.m_childIndices[k][0]] = faceToMatrixLeakOffRate[k];
          totalLeakedThickness[k] += faceToMatrixLeakOffRate[k] * dt;
          totalLeakedThickness[domain.m_feFaceManager.m_childIndices[k][0]] += faceToMatrixLeakOffRate[domain.m_feFaceManager.m_childIndices[k][0]] * dt;
        }
      }
    }
  }

  Epetra_IntSerialDenseVector  faceDofIndex (1);
  Epetra_SerialDenseVector     flowRHS     (1);

  for( localIndex k=0 ; k<domain.m_feFaceManager.DataLengths() ; ++k )
  {
    if (flowFaceType[k] == 1 && face_is_ghost[k] < 0)
    {
      faceDofIndex[0] = trilinos_index[k];

      realT A = faceArea[k];
      flowRHS(0) =std::min (2*A*m_rho_o*faceToMatrixLeakOffRate[k]*dt, faceFluidMass[k] * m_maxLeakOffRatio);

      m_rhs->SumIntoGlobalValues(faceDofIndex,flowRHS);
    }
  }
}



void ParallelPlateFlowSolverFVP::SetInitialGuess( PhysicalDomainT& domain )
{
#if 0
  array<real64>& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const array<real64>& pressureIncrement = domain.m_feFaceManager.GetFieldData<realT>("pressureIncrement");

  pressure += pressureIncrement;
#endif
}

void ParallelPlateFlowSolverFVP::SetInitialGuess( const PhysicalDomainT& domain,
                                                  realT* const local_solution )
{
#if 0
  const array<integer>& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_trilinosIndexStr);
  const array<integer>& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const array<real64>& pressureIncrement = domain.m_feFaceManager.GetFieldData<realT>("pressureIncrement");


  std::cout<<"ParallelPlateFlowSolverFVP::SetInitialGuess"<<std::endl;

  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
//      const int intDofOffset = dofOffset;
      int lid = m_rowMap->LID( /*intDofOffset+*/ trilinos_index[*kf]);
      local_solution[lid] = pressureIncrement[*kf];
//      std::cout<<local_solution[lid]<<std::endl;
    }
  }
#endif
}

/// Register solver in the solver factory
REGISTER_SOLVER( ParallelPlateFlowSolverFVP )
