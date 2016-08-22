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

#include "LagrangeSolverBase.h"
#include "LagrangeHelperFunctions.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "PhysicsSolvers/PhysicsSolverStrings.h"

#include "Utilities/Utilities.h"
#include "Utilities/Kinematics.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "ElementLibrary/FiniteElement.h"
#include "ElementLibrary/FiniteElementUtilities.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"

#include "ObjectManagers/TableManager.h"
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"

#include "Constitutive/CohesiveZone/CohesiveZoneFactory.h"

#include "Constitutive/Material/MaterialFactory.h"
#include "MeshUtilities/MesquiteRelaxer.h"
/*
#include "SurfaceGeneration/Fractunator2.h"
#include "SurfaceGeneration/Fractunator3.h"
#include "SurfaceGeneration/Fractunator2D.h"
 */

#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
//#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_LinearProblem.h"

//#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "EpetraExt_MatrixMatrix.h"

#include <iomanip>


int LagrangeSolverBase::m_instances = 0;

static void ML_Coord2RBM( int Nnodes,
                          const Array1dT<R1Tensor>& x,
                          rArray1d& rbm,
                          const int Ndof );


LagrangeSolverBase::LagrangeSolverBase( const std::string& name,
                                        ProblemManagerT* const pm ):
    SolverBase(name,pm),
    m_timeIntegrationOption(ExplicitDynamic),
    this_mpi_process(epetra_comm->MyPID()),
    n_mpi_processes(epetra_comm->NumProc()),
    verbose(0),
    m_cfl(1.0),
    m_2dOption(PlaneStrain),
    m_bulkQLinear(0.0),
    m_bulkQQuadratic(0.0),
    m_gravityVector(0.0)
{
  for( int a=FirstTimeIntegrationEnum ; a<numTimeIntegrationEnums ; ++a )
  {
    m_enableTimeIntegrationOption[a] = 0;
  }
}


LagrangeSolverBase::~LagrangeSolverBase()
{}



void LagrangeSolverBase::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML(hdn);
  m_numerics.krylov_tol     = hdn->GetAttributeOrDefault("tol",1.0e-6);

  m_numerics.m_useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);

  m_trilinosIndexStr = "IMS_" +  toString<int>(m_instances) + "_GlobalDof";  
  ++m_instances;


  m_courant = hdn->GetAttributeOrDefault<realT>("courant",0.5);
  m_dampingM = hdn->GetAttributeOrDefault<realT>("dampingM",0.0);
  m_dampingK = hdn->GetAttributeOrDefault<realT>("dampingK",0.0);

  m_bulkQLinear = hdn->GetAttributeOrDefault<realT>("bulkQ1",0.0);
  m_bulkQQuadratic = hdn->GetAttributeOrDefault<realT>("bulkQ2",0.0);
  m_newmarkGamma = hdn->GetAttributeOrDefault<realT>("newmarkGamma",0.5);
  m_newmarkBeta  = hdn->GetAttributeOrDefault<realT>("newmarkBeta", pow(m_newmarkGamma+0.5,2.0)/4.0 );
  m_useVelocityQS  = hdn->GetAttributeOrDefault<int>("useVelocityQS", 0 );
  m_displacementPenalty = hdn->GetAttributeOrDefault<realT>("displacementPenalty",0.0);


  m_tiedNodesFlag = hdn->GetAttributeOrDefault<int> ("tiedNodesFlag", 0);
  m_tiedNodeNormalRuptureStress = hdn->GetAttributeOrDefault<realT> ("tiedNodeNormalRuptureStress", 1.0e99);
  m_tiedNodeShearRuptureStress = hdn->GetAttributeOrDefault<realT> ("tiedNodeShearRuptureStress", 1.0e99);
  m_tiedNodeTolerance = hdn->GetAttributeOrDefault<realT> ("tiedNodeTolerance", 1.0e-8);

  m_enableTimeIntegrationOption[QuasiStatic]= hdn->GetAttributeOrDefault<int> ("QS", 0);
  m_enableTimeIntegrationOption[ImplicitDynamic]= hdn->GetAttributeOrDefault<int> ("ID", 0);
  m_enableTimeIntegrationOption[ExplicitDynamic]= hdn->GetAttributeOrDefault<int> ("ED", 0);

  {
    int temp = hdn->GetAttributeOrDefault<int> ("timeIntegrationOption", 0);
    m_timeIntegrationOption = IntToTimeIntegration( temp );
    m_enableTimeIntegrationOption[m_timeIntegrationOption] = 1;
  }

  {
    int temp = hdn->GetAttributeOrDefault<int> ("twoD_option", 0);
    m_2dOption = IntToTwoDOptions( temp );
  }

  m_timeToSnapshotDisp = hdn->GetAttributeOrDefault<realT> ("timeToSnapshotDisp",  std::numeric_limits<realT>::max());

  m_refTemperature = hdn->GetAttributeOrDefault<realT> ("referenceTemperature", 0.0);
  m_defaultCTE = hdn->GetAttributeOrDefault<realT> ("defaultCTE", 0.8e-5);
  m_useNodalTemperature = hdn->GetAttributeOrDefault<int> ("useNodalTemperature", 0);
  m_thermalRegionNames = hdn->GetStringVector("thermalRegionNames");

  R1Tensor zeroVector;
  zeroVector *= 0.0;
  m_gravityVector = hdn->GetAttributeOrDefault<R1Tensor>("gravityVector", zeroVector);

  m_relax = hdn->GetAttributeOrDefault<int> ("relax", 0);
  m_writeNodalStress =  hdn->GetAttributeOrDefault<bool> ("writeNodalStress", false);


}




void LagrangeSolverBase::RegisterFields( PhysicalDomainT& domain )
{
  //FIXME: Just noticed that the small strain solvers does not call this function.  Need to fix it.

  // register nodal fields
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::displacement>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();
  domain.m_feNodeManager.AddKeyedDataField<FieldInfo::hgforce>();


  if( m_displacementPenalty > 0.0 )
  {
    domain.m_feNodeManager.AddKeylessDataField<realT>("displacementPenaltyForce", true, true);
  }
  //  for( std::map< std::string, ElementRegionT >::iterator i=domain.m_feElementManager.m_ElementRegions.begin() ;
  //      i != domain.m_feElementManager.m_ElementRegions.end() ; ++i )
  //  {
  //    i->second.AddKeylessDataField<int>("partitionColor", false, true) ;
  //  }

  for(localIndex i = 0; i < m_thermalRegionNames.size(); ++i)
  {
    std::map<std::string, ElementRegionT>::iterator it = domain.m_feElementManager.m_ElementRegions.find(m_thermalRegionNames[i]);
    if (it == domain.m_feElementManager.m_ElementRegions.end() )
      throw GPException("Cannot find the thermal region specified!");

    ElementRegionT& elemRegion = it -> second;

    elemRegion.AddKeylessDataField<realT>("temperature", true, true);
    elemRegion.AddKeylessDataField<realT>("linearCTE", true, false);
    // This is the linear coefficient of thermal expansion, not the volumetric one.
    elemRegion.AddKeylessDataField<realT>("antiThermalStress", false, true);
    // This is the hydrostatic stress that would cause the same amount of volumetric strain as the thermal expantion/contraction would

    rArray1d& temperature = elemRegion.GetFieldData<realT>("temperature");
    rArray1d& CTE = elemRegion.GetFieldData<realT>("linearCTE");
    temperature = m_refTemperature;
    CTE = m_defaultCTE;

    //We don't register a field for reference temperature.  If the user needs varying ref temp in the domain, he/she will initialize the field in initial conditions, which will automatically create the field.

  }

  if( m_enableTimeIntegrationOption[QuasiStatic] )
  {
    domain.m_feNodeManager.AddKeylessDataField<int>(m_trilinosIndexStr,true,true);
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::volume>();
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::contactForce>();
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::velocity>();
  }

  if( m_enableTimeIntegrationOption[ImplicitDynamic] || m_enableTimeIntegrationOption[ExplicitDynamic] )
  {
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::acceleration>();
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::mass>();
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::velocity>();

    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::force>();
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::contactForce>();

    domain.m_feNodeManager.AddKeylessDataField<realT>("work", true, false);
    domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("dampingForce", false, true);

  }

  if( m_enableTimeIntegrationOption[ImplicitDynamic] )
  {
    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::volume>();
    //    domain.m_feNodeManager.AddKeylessDataField<int>(m_trilinosIndexStr,true,true);
    domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("vel_tilde",false,false);
    domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("inc_disp_tilde",false,false);
  }

  if( m_enableTimeIntegrationOption[ExplicitDynamic] )
  {
    // register discrete element fields - all of the nodal fields + currentPosition + rotation
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::displacement>();
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::velocity>();
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::acceleration>();
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::force>();
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::mass>();
    domain.m_discreteElementManager.AddKeylessDataField<realT>("work", true, false);
    //rotation states
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::rotationalAxisIncrement>();
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::rotationalMagnitudeIncrement>();
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::rotationalVelocity>();
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::rotationalAcceleration>();
    domain.m_discreteElementManager.AddKeyedDataField<FieldInfo::moment>();


    // trying out ellipsoidal discrete elements
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::referencePosition>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::displacement>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::velocity>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::acceleration>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::force>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::mass>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::currentPosition>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeylessDataField<realT>("work", true, false);

    //rotation states
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationAxis>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationMagnitude>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalAxisIncrement>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalMagnitudeIncrement>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalVelocity>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalAcceleration>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::moment>();
    domain.m_ellipsoidalDiscreteElementManager.AddKeyedDataField<FieldInfo::rotationalInertia>();

    if  (m_timeToSnapshotDisp < 0.99*std::numeric_limits<realT>::max())
    {
      domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("refDisplacement", true, true);
    }
  }


  if (m_writeNodalStress)
  {
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_x", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_y", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_z", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_xy", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_yz", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_zx", false, true);
  }

  // Define common visualization fields
  m_commonFields.push_back("sigma_x");
  m_commonFields.push_back("sigma_y");
  m_commonFields.push_back("sigma_z");
  m_commonFields.push_back("sigma_xy");
  m_commonFields.push_back("sigma_xz");
  m_commonFields.push_back("sigma_yz");
  m_commonFields.push_back("Displacement");
  m_commonFields.push_back("Velocity");

}




void LagrangeSolverBase::RegisterTemporaryFields( PhysicalDomainT& domain )
{
  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("displacement_n",false,false);
}



void LagrangeSolverBase::DeregisterTemporaryFields( PhysicalDomainT& domain )
{
  domain.m_feNodeManager.RemoveDataField<R1Tensor>("displacement_n");
}



void LagrangeSolverBase::FillTemporaryFields( PhysicalDomainT& domain  )
{
  const Array1dT<R1Tensor>& disp_np1 = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
  Array1dT<R1Tensor>& disp_n   = domain.m_feNodeManager.GetFieldData<R1Tensor>("displacement_n");
  disp_n = disp_np1;
}



void LagrangeSolverBase::OverwriteFieldsWithTemporaryFields( PhysicalDomainT& domain )
{
  Array1dT<R1Tensor>& disp_np1 = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
  const Array1dT<R1Tensor>& disp_n   = domain.m_feNodeManager.GetFieldData<R1Tensor>("displacement_n");
  disp_np1 = disp_n;
}



void LagrangeSolverBase::Initialize(PhysicalDomainT& domain, SpatialPartition& partition )
{
  CalculateAllNodalMassAndVolume( domain, partition );

  if( m_tiedNodesFlag )
  {
    domain.m_feFaceManager.AddKeylessDataField<int>("numKCBCnodesInFace",true,true);
    BoundaryConditionFunctions::BuildKinematicConstraintBoundaryCondition( domain.m_feNodeManager,
                                                                           domain.m_feFaceManager,
                                                                           m_KinematicConstraintNodes,
                                                                           m_tiedNodeTolerance );
  }
  dim = domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;

}



void LagrangeSolverBase::InitializeCommunications( PartitionBase& partition )
{
  m_syncedFields.clear();

  if( m_enableTimeIntegrationOption[QuasiStatic] )
  {
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::incrementalDisplacement>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_trilinosIndexStr);

  }

  if( m_enableTimeIntegrationOption[ImplicitDynamic] || m_enableTimeIntegrationOption[ExplicitDynamic] )
  {
    //    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_trilinosIndexStr);
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::incrementalDisplacement>::Name());


    m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("cohesiveTraction");
    m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("ruptureState");

  }

  if( m_enableTimeIntegrationOption[ExplicitDynamic] )
  {
    m_syncedFields[PhysicalDomainT::DiscreteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    m_syncedFields[PhysicalDomainT::DiscreteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    m_syncedFields[PhysicalDomainT::DiscreteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());

  }

  partition.SetBufferSizes(m_syncedFields, CommRegistry::lagrangeSolver01);
  partition.SetBufferSizes(m_syncedFields, CommRegistry::lagrangeSolver02);

}


double LagrangeSolverBase::TimeStep( const realT& time,
                                     const realT& dt,
                                     const int cycleNumber,
                                     PhysicalDomainT& domain,
                                     const sArray1d& namesOfSolverRegions,
                                     SpatialPartition& partition,
                                     FractunatorBase* const fractunator)
{
  realT dt_return = dt;
  dim = domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;

#ifdef SRC_INTERNAL2
  if(domain.m_xfemManager != NULL)
  {
    TimeStepXFEM(time, dt, cycleNumber, domain, namesOfSolverRegions, partition);
  }
  else
#endif
  {
    if( m_timeIntegrationOption==ExplicitDynamic)
    {
      if( fractunator!=NULL )
      {
        if( cycleNumber%fractunator->m_checkInterval==0 )
        {
          fractunator->SeparationDriver( domain.m_feNodeManager,
                                         domain.m_feEdgeManager,
                                         domain.m_feFaceManager,
                                         domain.m_externalFaces,
                                         domain.m_feElementManager,
                                         partition, false , time);
        }
      }

      TimeStepExplicitDynamic( time, dt , domain, namesOfSolverRegions, partition );
    }
    else if( m_timeIntegrationOption==ImplicitDynamic ||
        m_timeIntegrationOption==QuasiStatic )
    {
      dt_return = TimeStepImplicit( time, dt , domain, namesOfSolverRegions, partition, fractunator );
    }
  }
  return dt_return;
}

//double LagrangeSolverBase::TimeStep( const realT& time,
//                                        const realT& dt,
//                                        const int cycleNumber,
//                                        PhysicalDomainT& domain,
//                                        const sArray1d& namesOfSolverRegions,
//                                        SpatialPartition& partition,
//                                        FractunatorBase* const fractunator,
//                                        XfemManager* elementSplitter)
//{
//
//  realT dt_return = dt;
//
//  if( m_timeIntegrationOption==QuasiStatic)
//  {
//
//    Array1dT<R1Tensor>& incdisp  = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
//
//    int resolve = 1; int count = 0;
//    while( resolve != 0 )
//    {
//
//      incdisp = 0; domain.m_splitElemThisIter = false;
//      if(verbose)
//        std::cout << std::endl
//        << " :: Implicit Mechanics Solver  " << std::endl
//        << " :: No. mpi processes ... " << n_mpi_processes << std::endl;
//
//      std::cout << "Setting up system" << std::endl;
//      SetupSystem (domain,partition,time+dt, dt);
//      std::cout << "Assembling matrix" << std::endl;
//
//
//#if USECPP11==1
//      m_matrix   = std::make_shared<Epetra_FECrsMatrix>(Copy,*m_sparsity);
//      m_solution = std::make_shared<Epetra_FEVector>(*m_rowMap);
//      m_rhs      = std::make_shared<Epetra_FEVector>(*m_rowMap);
//#else
//      if( m_matrix!=NULL )
//        delete m_matrix;
//      m_matrix   = new Epetra_FECrsMatrix(Copy,*m_sparsity);
//
//      if( m_solution!=NULL )
//        delete m_solution;
//      m_solution = new Epetra_FEVector(*m_rowMap);
//
//      if( m_rhs!=NULL )
//        delete m_rhs;
//      m_rhs      = new Epetra_FEVector(*m_rowMap);
//#endif
//
//      m_matrix->Scale(0.0);
//      m_rhs->Scale(0.0);
//      m_solution->Scale(0.0);
//
//      Assemble    (domain,partition,time+dt);
//
////      std::string sMatrix;
////      sMatrix = "systemMatrix.dat";
////      EpetraExt::RowMatrixToMatlabFile(sMatrix.c_str(),*m_matrix);
////
////      std::string sRhs;
////      sRhs = "system-RHS-branch.dat";
////      EpetraExt::MultiVectorToMatlabFile(sRhs.c_str(),*m_rhs);
////
//      std::cout << "Solving system" << std::endl;
//      Solve       (domain,partition,time+dt, dt);
//
//      auto& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
//      auto& kinematicSibling = domain.m_feNodeManager.GetOneToOneMap("NodeToKinematicSibling");
//
//      for(localIndex iNod=0; iNod<domain.m_feNodeManager.m_numNodes; iNod++)
//      {
//        if(kinematicSibling(iNod)!=-1)
//        {
//          disp[iNod] = disp[kinematicSibling(iNod)];
//        }
//      }
//      // Evaluate displacement values at vertex nodes
//      // This is a post-processing step
//      //FIXME: Does not work for quadrilaterals
//      std::cout<<"Calculating displacement values at vertices"<<std::endl;
//      EvaluateDispAtVertex(domain,partition,time+dt);
//
//      std::cout<<"Calculating nodal forces"<<std::endl;
//      GetNodalForces(domain, partition, dt);
//
//      if(verbose)
//        std::cout << std::endl;
//
//      // Uncomment this to split elements that are not prefractured
//      elementSplitter->SeparationDriver(false,time,0,domain.m_feNodeManager, domain.m_feEdgeManager, domain.m_feFaceManager, domain.m_externalFaces, domain.m_feElementManager, partition,
//                                        domain.m_crackManager, domain.m_crackSurfaceVertex, domain.m_splitElemThisIter);
//
//      if(domain.m_splitElemThisIter == false)
//        resolve = 0;
//      else
//      {
//        if(elementSplitter->m_fixed_crack_increment)
//        {
//          count = count+1;
//          if(count==10)
//            resolve = 0;
//        }
//      }
//
//      int localResolve = resolve;
//      MPI_Allreduce (&localResolve,&resolve,1,MPI_INT,MPI_MAX ,MPI_COMM_WORLD);
//
//      if(domain.m_WriteXFEM==false)
//      {
//        if(domain.m_splitElemThisIter == true)
//          domain.m_WriteXFEM = true;
//      }
//    }
//    //this->m_stabledt.m_maxdt = dt;
//
//    this->m_stabledt.m_maxdt = std::numeric_limits<double>::max()* 0.95;
//
////    WriteSilo( false );
//  }
//  else
//  {
//    throw GPException("Chosen time integration option not implemented yet!");
//  }
//
//  return dt_return;
//}
#ifdef SRC_INTERNAL2
void LagrangeSolverBase::TimeStepXFEM( const realT& time,
                                       const realT& dt,
                                       const int cycleNumber,
                                       PhysicalDomainT& domain,
                                       const sArray1d& namesOfSolverRegions,
                                       SpatialPartition& partition)
{
  if( m_timeIntegrationOption==QuasiStatic)
  {

    Array1dT<R1Tensor>& incdisp  = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();

    int resolve = 1; int count = 0;
    while( resolve != 0 )
    {

      incdisp = 0; domain.m_splitElemThisIter = false;
      if(verbose)
        std::cout << std::endl
        << " :: Implicit Mechanics Solver  " << std::endl
        << " :: No. mpi processes ... " << n_mpi_processes << std::endl;

      std::cout << "Setting up system" << std::endl;
      SetupSystem (domain,partition);
      std::cout << "Assembling matrix" << std::endl;


#if USECPP11==1
      m_matrix   = std::make_shared<Epetra_FECrsMatrix>(Copy,*m_sparsity);
      m_solution = std::make_shared<Epetra_FEVector>(*m_rowMap);
      m_rhs      = std::make_shared<Epetra_FEVector>(*m_rowMap);
#else
      if( m_matrix!=NULL )
        delete m_matrix;
      m_matrix   = new Epetra_FECrsMatrix(Copy,*m_sparsity);

      if( m_solution!=NULL )
        delete m_solution;
      m_solution = new Epetra_FEVector(*m_rowMap);

      if( m_rhs!=NULL )
        delete m_rhs;
      m_rhs      = new Epetra_FEVector(*m_rowMap);
#endif

      m_matrix->Scale(0.0);
      m_rhs->Scale(0.0);
      m_solution->Scale(0.0);

      Assemble    (domain,partition,time+dt, dt);

      //      std::string sMatrix;
      //      sMatrix = "systemMatrix.dat";
      //      EpetraExt::RowMatrixToMatlabFile(sMatrix.c_str(),*m_matrix);
      //
      //      std::string sRhs;
      //      sRhs = "system-RHS-branch.dat";
      //      EpetraExt::MultiVectorToMatlabFile(sRhs.c_str(),*m_rhs);
      //
      std::cout << "Solving system" << std::endl;
      Solve       (domain,partition,time+dt, dt);

      auto& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
      auto& kinematicSibling = domain.m_feNodeManager.GetOneToOneMap("NodeToKinematicSibling");

      for(localIndex iNod=0; iNod<domain.m_feNodeManager.m_numNodes; iNod++)
      {
        if( kinematicSibling(iNod)!=static_cast<localIndex>(-1) )
        {
          disp[iNod] = disp[kinematicSibling(iNod)];
        }
      }
      // Evaluate displacement values at vertex nodes
      // This is a post-processing step
      //FIXME: Does not work for quadrilaterals
      std::cout<<"Calculating displacement values at vertices"<<std::endl;
      EvaluateDispAtVertex(domain,partition,time+dt);

      std::cout<<"Calculating nodal forces"<<std::endl;
      GetNodalForces(domain, partition, dt);

      if(verbose)
        std::cout << std::endl;

      // Uncomment this to split elements that are not prefractured
      domain.m_xfemManager->SeparationDriver(false,time,0,domain.m_feNodeManager, domain.m_feEdgeManager, domain.m_feFaceManager, domain.m_externalFaces, domain.m_feElementManager, partition,
                                             *(domain.m_crackManager), *(domain.m_crackSurfaceVertex), domain.m_splitElemThisIter);

      if(domain.m_splitElemThisIter == false)
        resolve = 0;
      else
      {
        if(domain.m_xfemManager->m_fixed_crack_increment)
        {
          count = count+1;
          if(count==10)
            resolve = 0;
        }
      }

      int localResolve = resolve;
      MPI_Allreduce (&localResolve,&resolve,1,MPI_INT,MPI_MAX ,MPI_COMM_WORLD);

      if(domain.m_WriteXFEM==false)
      {
        if(domain.m_splitElemThisIter == true)
          domain.m_WriteXFEM = true;
      }
    }
    //this->m_stabledt.m_maxdt = dt;

    this->m_stabledt.m_maxdt = std::numeric_limits<double>::max()* 0.95;

    //    WriteSilo( false );
  }
  else
  {
    throw GPException("Chosen time integration option not implemented yet!");
  }
}
#endif




void LagrangeSolverBase::TimeStepImplicitSetup( const realT& time,
                                                const realT& dt ,
                                                PhysicalDomainT& domain )
{

  if( this->m_timeIntegrationOption == ImplicitDynamic )
  {
    Array1dT<R1Tensor> const & v_n = domain.m_feNodeManager.GetFieldData<FieldInfo::velocity>();
    Array1dT<R1Tensor> const & a_n = domain.m_feNodeManager.GetFieldData<FieldInfo::acceleration>();
    Array1dT<R1Tensor> & vtilde   = domain.m_feNodeManager.GetFieldData<R1Tensor>("vel_tilde");
    Array1dT<R1Tensor> & uhatTilde   = domain.m_feNodeManager.GetFieldData<R1Tensor>("inc_disp_tilde");

    Array1dT<R1Tensor>& uhat  = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
    Array1dT<R1Tensor>& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();

    for( auto a = 0u ; a < domain.m_feNodeManager.m_numNodes ; ++a )
    {
      for( int i=0 ; i<dim ; ++i )
      {
        vtilde[a][i] = v_n[a][i] + (1.0-m_newmarkGamma) * a_n[a][i] * dt;
        uhatTilde[a][i] = ( v_n[a][i] + 0.5 * ( 1.0 - 2.0*m_newmarkBeta ) * a_n[a][i] * dt ) *dt;
      }
    }
    uhat = uhatTilde;
    disp += uhatTilde;

  }
  else if( this->m_timeIntegrationOption == QuasiStatic  )
  {

    Array1dT<R1Tensor>& uhat  = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
    if( m_useVelocityQS==1 )
    {
      Array1dT<R1Tensor> const & v_n = domain.m_feNodeManager.GetFieldData<FieldInfo::velocity>();
      Array1dT<R1Tensor>& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();

      for( auto a = 0u ; a < domain.m_feNodeManager.m_numNodes ; ++a )
      {
        for( int i=0 ; i<dim ; ++i )
        {
          uhat[a][i] = v_n[a][i] * dt;
          disp[a][i] += uhat[a][i];
        }
      }
    }
    else
    {
      uhat = 0.0;
    }
  }
}

void LagrangeSolverBase::TimeStepImplicitComplete( const realT& time,
                                                   const realT& dt ,
                                                   PhysicalDomainT& domain )
{
  if( this->m_timeIntegrationOption == ImplicitDynamic )
  {
    Array1dT<R1Tensor>& v_n = domain.m_feNodeManager.GetFieldData<FieldInfo::velocity>();
    Array1dT<R1Tensor>& uhat  = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
    Array1dT<R1Tensor>& a_n = domain.m_feNodeManager.GetFieldData<FieldInfo::acceleration>();

    Array1dT<R1Tensor> & vtilde   = domain.m_feNodeManager.GetFieldData<R1Tensor>("vel_tilde");
    Array1dT<R1Tensor> & uhatTilde   = domain.m_feNodeManager.GetFieldData<R1Tensor>("inc_disp_tilde");

    for( auto a = 0u ; a < domain.m_feNodeManager.m_numNodes ; ++a )
    {
      for( int i=0 ; i<dim ; ++i )
      {
        //        realT a_np1 = 4.0 / (dt*dt) * ( uhat[a][i] - uhatTilde[a][i] );
        a_n[a][i] = 1.0 / ( m_newmarkBeta * dt*dt) * ( uhat[a][i] - uhatTilde[a][i] );
        v_n[a][i] = vtilde[a][i] + m_newmarkGamma * a_n[a][i] * dt;
      }
    }
  }
  else if( this->m_timeIntegrationOption == QuasiStatic && dt > 0.0)
  {
    Array1dT<R1Tensor>& v_n = domain.m_feNodeManager.GetFieldData<FieldInfo::velocity>();
    Array1dT<R1Tensor>& uhat  = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
    for( auto a = 0u ; a < domain.m_feNodeManager.m_numNodes ; ++a )
    {
      for( int i=0 ; i<dim ; ++i )
      {
        v_n[a][i] = uhat[a][i] / dt;
      }
    }
  }
}


realT LagrangeSolverBase::TimeStepImplicit( const realT& time,
                                            const realT& dt ,
                                            PhysicalDomainT& domain,
                                            const sArray1d& namesOfSolverRegions,
                                            SpatialPartition& partition,
                                            FractunatorBase* const fractunator )
{

//  CalculateAllNodalMassAndVolume( domain, partition );

  realT dt_return = dt;
//  Array1dT<R1Tensor>& uhat  = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();


  TimeStepImplicitSetup( time, dt, domain );


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(m_numerics.m_useNewtonSolve == false)
  {

    int resolve = 1;
    while( resolve != 0 )
    {
      if(verbose)
      {
        std::cout << std::endl
            << " :: Implicit Mechanics Solver  " << std::endl
            << " :: No. mpi processes ... " << n_mpi_processes << std::endl;

        std::cout << "Setting up system" << std::endl;
      }
      if(domain.m_externalFaces.m_contactActive && domain.m_contactManager.m_use_contact_search )
      {
        UpdateContactDataStructures(domain, true);
      }

      SetupSystem (domain,partition);

      if(verbose)
        std::cout << "Assembling matrix" << std::endl;


#if USECPP11==1
      m_matrix   = std::make_shared<Epetra_FECrsMatrix>(Copy,*m_sparsity);
      m_solution = std::make_shared<Epetra_FEVector>(*m_rowMap);
      m_rhs      = std::make_shared<Epetra_FEVector>(*m_rowMap);
#else
      if( m_matrix!=NULL )
        delete m_matrix;
      m_matrix   = new Epetra_FECrsMatrix(Copy,*m_sparsity);

      if( m_solution!=NULL )
        delete m_solution;
      m_solution = new Epetra_FEVector(*m_rowMap);

      if( m_rhs!=NULL )
        delete m_rhs;
      m_rhs      = new Epetra_FEVector(*m_rowMap);
#endif

      m_matrix->Scale(0.0);
      m_rhs->Scale(0.0);
      m_solution->Scale(0.0);

      Assemble    (domain,partition,time+dt, dt);
      if(verbose)
        std::cout << "Solving system" << std::endl;

      //      std::cout<<*m_matrix<<std::endl;
      Solve(domain,partition,time+dt, dt);

      // Calling this function just to visualize contact tractions
      if(domain.m_contactManager.m_implicitContactActive)
        PostProcessFieldsForVisualizationAndConsistency(domain);

      if(verbose)
        std::cout << std::endl;

      if( fractunator!=NULL)
      {
        int localResolve = 0;
        localResolve=fractunator->SeparationDriver( domain.m_feNodeManager,
                                                    domain.m_feEdgeManager,
                                                    domain.m_feFaceManager,
                                                    domain.m_externalFaces,
                                                    domain.m_feElementManager,
                                                    partition, false, time );

        MPI_Allreduce (&localResolve,&resolve,1,MPI_INT,MPI_MAX ,MPI_COMM_WORLD);

        if( m_relax == 1 )
        {
          MesquiteRelaxer dummyObject;
          dummyObject.RelaxMesh( domain );
        }
      }
      else
      {
        resolve = 0;
      }
    }
    //this->m_stabledt.m_maxdt = dt;

    this->m_stabledt.m_maxdt = std::numeric_limits<double>::max()* 0.95;
    // Fu: If we don't give a hard dt in the solver application (i.e. letting other solvers dictate), maxdt=dt will make the dt be stuck at a small value.
    // Since this is quasi-static and we usually give a hard dt, this m_maxdt thing does not do anything.
    // So setting it to a large number will enable other solvers in the same solver application to dictate the time stepping.
  }
  else
  {
    int resolve = 1;
    while( resolve != 0 )
    {

      if(verbose)
      {
        std::cout << std::endl
            << " :: Implicit Mechanics Solver  " << std::endl
            << " :: No. mpi processes ... " << n_mpi_processes << std::endl;

        std::cout << "Setting up system" << std::endl;
      }

      if(domain.m_externalFaces.m_contactActive && domain.m_contactManager.m_use_contact_search )
      {
        UpdateContactDataStructures(domain, time==0 ? true : false );
      }
      SetupSystem( domain, partition );
      RegisterTemporaryFields( domain );

      //        int rank, size;
      //        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      //        MPI_Comm_size(MPI_COMM_WORLD, &size);

      realT stepDt = dt;
      int converged = 0;
      realT thisNorm;
      realT thisNorm0 = 0;
      realT energyNorm = 0;

      FillTemporaryFields( domain );

      if( rank== 0 )
      {

        std::cout<<std::endl;
        std::cout<<std::endl;
        std::cout<<"*************** Begin Newton loop ******************************"<<std::endl;
      }

      while( converged == 0)
      {
        for( int iter=0 ; iter<m_numerics.m_maxIterNewton ; ++iter )
        {
          //          uhat = 0;
          thisNorm = 0;
          partition.SynchronizeFields(m_syncedFields, CommRegistry::lagrangeSolver02);


#if USECPP11==1
          m_matrix   = std::make_shared<Epetra_FECrsMatrix>(Copy,*m_sparsity);
          m_solution = std::make_shared<Epetra_FEVector>(*m_rowMap);
          m_rhs      = std::make_shared<Epetra_FEVector>(*m_rowMap);
#else
          if( m_matrix!=NULL )
            delete m_matrix;
          m_matrix   = new Epetra_FECrsMatrix(Copy,*m_sparsity);

          if( m_solution!=NULL )
            delete m_solution;
          m_solution = new Epetra_FEVector(*m_rowMap);

          if( m_rhs!=NULL )
            delete m_rhs;
          m_rhs      = new Epetra_FEVector(*m_rowMap);
#endif

          m_matrix->Scale(0.0);
          m_rhs->Scale(0.0);
          m_solution->Scale(0.0);

          Assemble(domain,partition,time+stepDt, stepDt);
          Solve(domain,partition,time+stepDt, stepDt);

          int size_res;
          double* global_residual = NULL;
          m_rhs->ExtractView(&global_residual,&size_res);

          int size_deltaU;
          double* deltaU = NULL;
          m_solution->ExtractView(&deltaU,&size_deltaU);

          for (int i = 0; i<size_deltaU; ++i)
          {
            thisNorm += fabs(deltaU[i]*global_residual[i]);
          }

          realT globalThisNorm = 0.0;
          MPI_Allreduce (&thisNorm,&globalThisNorm,1,MPI_DOUBLE,MPI_SUM ,MPI_COMM_WORLD);

          thisNorm = globalThisNorm;

          if(iter==0)
          {
            thisNorm0 = thisNorm;
            if( thisNorm < m_numerics.newton_tol * 1.0e-4 )
            {
              converged = 1;
              energyNorm = thisNorm;
              break;
            }
          }

          energyNorm = thisNorm/thisNorm0;

          //            incdisp; disp;

          if( rank== 0 )
          {
            std::cout << "  Iteration " <<iter<<":    "<<"Energy norm: "<<std::setprecision(3)<<energyNorm<< std::endl;
          }

          if( energyNorm < m_numerics.newton_tol)
          {
            converged = 1;
            break;
          }
        }
        if( converged==0 )
        {
          stepDt *= 0.5;
          this->m_stabledt.m_maxdt = stepDt * 1.1;
          //Fu: Why do this?  It has not effect since this is quasi-static and we have a hard dt.  If we remove this line we will have to rebaseline a lot of tests.
          OverwriteFieldsWithTemporaryFields(domain);
          std::cout<<"No newton convergence... substepping at dt="<<stepDt<<" for this step "<<std::endl;
          dt_return = stepDt;

        }
        else
        {
          if(domain.m_contactManager.m_sliding_law != 0)
          {
            if(rank==0 && verbose)
            {
              std::cout<<"  Updating plastic slip variables " <<std::endl;
            }
          }
          PostProcessFieldsForVisualizationAndConsistency(domain);

          if(domain.m_contactManager.m_sliding_law != 0)
          {
            /* Scratch code
             * uJumpPl(n) = uJumpPl(n+1);
             * uJumpPl(n+1) = 0;*/
            StoreHistoryVariablesForCurrentLoadStepAndResetTheField(domain);
          }

          if(rank==0)
          {
            std::cout<<"  Converged:      Energy norm at convergence = "<<energyNorm<<std::endl;

            std::cout<<"*************** End Newton loop ******************************"<<std::endl;
          }

          this->m_stabledt.m_maxdt = stepDt * 1.1;  //Fu: Why do this?  It has not effect since this is quasi-static and we have a hard dt.
        }
      }
      DeregisterTemporaryFields( domain );

      if(verbose)
        std::cout << std::endl;

      if( fractunator!=NULL )
      {
        int localResolve = 0;
        localResolve=fractunator->SeparationDriver( domain.m_feNodeManager,
                                                    domain.m_feEdgeManager,
                                                    domain.m_feFaceManager,
                                                    domain.m_externalFaces,
                                                    domain.m_feElementManager,
                                                    partition, false, time );

        MPI_Allreduce (&localResolve,&resolve,1,MPI_INT,MPI_MAX ,MPI_COMM_WORLD);
      }
      else
      {
        resolve = 0;
      }

    }

  }

  TimeStepImplicitComplete( time, dt ,  domain );

  //    partition.SynchronizeFields( m_syncedFields, CommRegistry::lagrangeSolver01 );
  return dt_return;
}

void LagrangeSolverBase::TimeStepExplicitDynamic( const realT& time,
                                                  const realT& dt ,
                                                  PhysicalDomainT& domain,
                                                  const sArray1d& namesOfSolverRegions,
                                                  SpatialPartition& partition )
{
  m_stabledt.m_maxdt = std::numeric_limits<double>::max();
  Array1dT<R1Tensor>& hgforce = domain.m_feNodeManager.GetFieldData<FieldInfo::hgforce> ();
  hgforce = 0.0;

  // update nodes
  LagrangeHelperFunctions::LinearPointUpdatePart1( domain.m_feNodeManager, time, dt);

  LagrangeHelperFunctions::LinearPointUpdatePart1( domain.m_discreteElementManager, time, dt, false );

  //We only need to zero the velocity of detached node here to prevent a huge damping force.
  //We will update their location and velocity in the postprocessing part.
  domain.m_feNodeManager.ZeroDetachedNodeVelocity();

  LagrangeHelperFunctions::RotationalPointUpdatePart1( domain.m_discreteElementManager, time, dt );
  LagrangeHelperFunctions::RotationalPointUpdatePart1b( domain.m_discreteElementManager);

  LagrangeHelperFunctions::RotationalPointUpdatePart1( domain.m_ellipsoidalDiscreteElementManager, time, dt );
  LagrangeHelperFunctions::LinearPointUpdatePart1( domain.m_ellipsoidalDiscreteElementManager, time, dt, false);

  ApplyTemperatureBC(domain, time);

  ProcessElementRegions( domain.m_feNodeManager,
                         domain.m_feElementManager,
                         namesOfSolverRegions, dt ) ;


  BoundaryConditionFunctions::ApplyTractionBoundaryCondition( domain,time );

  // -- CONTACT --
  if (domain.m_externalFaces.m_contactActive)
  {
    ApplyForcesFromContact(domain, this->m_stabledt, dt);
  }

  domain.m_ellipsoidalDiscreteElementManager.UpdateCylindricalBoundary(m_stabledt, dt, time);

  if( m_tiedNodesFlag )
  {
    BoundaryConditionFunctions::ApplyKinematicConstraintBoundaryCondition( domain.m_feFaceManager,
                                                                           domain.m_feNodeManager,
                                                                           this->m_KinematicConstraintNodes,
                                                                           m_tiedNodeNormalRuptureStress,
                                                                           m_tiedNodeShearRuptureStress );
  }


  ProcessCohesiveZones( domain.m_feNodeManager,
                        domain.m_feFaceManager, dt );


#ifdef SRC_INTERNAL
  //FIXME: replace the following call with domain.m_externalFaces.GeodynCouplingParallel( domain.m_feNodeManager );
  domain.m_externalFaces.GeodynCoupling( domain.m_feNodeManager );
#endif
  R1Tensor zero(0);
  LagrangeHelperFunctions::LinearPointUpdatePart2( domain.m_feNodeManager, time, dt, zero, m_dampingM, m_dampingK );
  LagrangeHelperFunctions::LinearPointUpdatePart2( domain.m_discreteElementManager, time, dt, m_gravityVector );
  LagrangeHelperFunctions::RotationalPointUpdatePart2( domain.m_discreteElementManager, time, dt );

  LagrangeHelperFunctions::LinearPointUpdatePart2( domain.m_ellipsoidalDiscreteElementManager, time, dt, m_gravityVector );
  LagrangeHelperFunctions::RotationalPointUpdatePart2( domain.m_ellipsoidalDiscreteElementManager, time, dt );


  partition.SynchronizeFields( m_syncedFields, CommRegistry::lagrangeSolver01 );

  m_stabledt.m_maxdt *= this->m_courant;

  if (time > m_timeToSnapshotDisp)
  {
    SnapshotNodalDisplacement(domain.m_feNodeManager);
    m_timeToSnapshotDisp = std::numeric_limits<realT>::max();
  }


}

void LagrangeSolverBase::ProcessElementRegions( NodeManagerT& nodeManager,
                                                ElementManagerT& elementManager,
                                                const sArray1d& namesOfSolverRegions,
                                                const realT dt )
{

  std::set<std::string> usedRegionNames;
  for( sArray1d::const_iterator regionName = namesOfSolverRegions.begin() ;
      regionName != namesOfSolverRegions.end() ; ++regionName )
  {
    //this conditional supports DE, since an element region in DE may not exist in FE; this case should not
    //throw an exception, which occurs in the absence of the following block
    std::map<std::string, ElementRegionT>::iterator iter = elementManager.m_ElementRegions.find(*regionName);
    if( iter != elementManager.m_ElementRegions.end() && usedRegionNames.count(*regionName)==0 )
    {
      ElementRegionT& elementRegion = iter->second; // stlMapLookup( domain.m_elementManager.m_ElementRegions, *regionName );
      // We cannot apply lagrange solvers to flow face regions, but sometime we have to apply the solver to them (because we have applied the hydrofrac solver to all regions
      if (!elementRegion.m_elementGeometryID.compare(0,4,"STRI") ||  !elementRegion.m_elementGeometryID.compare(0,4,"CPE4") ||  !elementRegion.m_elementGeometryID.compare(0,2,"C3"))
      {
        ProcessElementRegion( nodeManager, elementRegion, dt );
        R1Tensor zeroVector;
        zeroVector *= 0.0;
        ApplyGravity( nodeManager, elementRegion, dt );
        usedRegionNames.insert( *regionName );
      }
    }



  }

  for(localIndex i = 0; i < m_thermalRegionNames.size(); ++i)
  {
    Epetra_SerialDenseVector * nullpointer = nullptr;
    std::map<std::string, ElementRegionT>::iterator it = elementManager.m_ElementRegions.find(m_thermalRegionNames[i]);

    ElementRegionT& elemRegion = it -> second;

    for(localIndex elemID = 0; elemID < elemRegion.m_numElems; ++elemID)
    {
      ApplyThermalStress( elemRegion, nodeManager, elemID, nullpointer);
    }
  }
}





void LagrangeSolverBase::ProcessCohesiveZones( NodeManagerT& nodeManager,
                                               FaceManagerT& faceManager,
                                               const realT dt )
{
  if( nodeManager.HasField<R1Tensor>("cohesiveForce") &&
      faceManager.HasField<int>("ruptureState") &&
      faceManager.m_cohesiveZone != NULL )
  {
    // apply cohesive forces
    const OrderedVariableOneToManyRelation& childFaceIndex = faceManager.GetVariableOneToManyMap( "childIndices" );
    Array1dT<R1Tensor>& nodalForce = nodeManager.GetFieldData<FieldInfo::force>();
    Array1dT<R1Tensor>& cohesiveForce = nodeManager.GetFieldData<R1Tensor>("cohesiveForce");

    iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");

    cohesiveForce = 0.0;

    for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
    {
      if( ruptureState[kf] == 2 && !(childFaceIndex[kf].empty()) )
      {


        const localIndex faceIndex[2] = { kf, childFaceIndex[kf][0] };


        const R1Tensor N[2] = { faceManager.FaceNormal( nodeManager, faceIndex[0] ),
                                faceManager.FaceNormal( nodeManager, faceIndex[1] )};

        R1Tensor Nbar = N[0];
        Nbar -= N[1];
        Nbar.Normalize();

        R1Tensor gap = faceManager.CalculateGapVector( nodeManager, kf );
        R2Tensor cStiffness;
        R1Tensor cohesiveTraction;


        ruptureState[kf] = faceManager.m_cohesiveZone->UpdateCohesiveZone( kf, gap, Nbar,
                                                                           faceManager.m_toElementsRelation[faceIndex[0]][0],
                                                                           faceManager.m_toElementsRelation[faceIndex[1]][0],
                                                                           cohesiveTraction, cStiffness );


        for( int side=0 ; side<2 ; ++side )
        {
          const localIndex faceID = faceIndex[side];
          const realT area = faceManager.SurfaceArea( nodeManager, faceID );

          R1Tensor cForce;
          int direction = -1;
          if(side==0)
            direction = 1;

          cForce = cohesiveTraction;
          cForce *= direction * area / faceManager.m_toNodesRelation[faceID].size();


          for( lArray1d::const_iterator nodeID=faceManager.m_toNodesRelation[faceID].begin() ;
              nodeID!=faceManager.m_toNodesRelation[faceID].end() ; ++nodeID )
          {
            nodalForce[*nodeID] += cForce;
            cohesiveForce[*nodeID] += cForce;
          }

        }
      }
    }
  }
}


void LagrangeSolverBase::SetNumRowsAndTrilinosIndices( PhysicalDomainT& domain,
                                                       SpatialPartition& partition,
                                                       int& numLocalRows,
                                                       int& numGlobalRows,
                                                       iArray1d& localIndices,
                                                       int offset )
{

  int n_ghost_rows  = domain.m_feNodeManager.GetNumGhosts();
  numLocalRows  = domain.m_feNodeManager.DataLengths()-n_ghost_rows;

  dim = domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;

  std::vector<int> gather(n_mpi_processes);

  epetra_comm->GatherAll(&numLocalRows,
                         &gather.front(),
                         1);

  int first_local_row = 0;
  numGlobalRows = 0;

  for(unsigned int p=0; p<n_mpi_processes; ++p)
  {
    numGlobalRows += gather[p];
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
      trilinos_index[r] = first_local_row+local_count+offset;
      localIndices.push_back(trilinos_index[r]);
      local_count++;
    }
    else
    {
      trilinos_index[r] = -INT_MAX;
    }
  }

  assert(local_count == static_cast<unsigned int>(numLocalRows) );

  partition.SynchronizeFields(m_syncedFields, CommRegistry::lagrangeSolver02);

  /*
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout<<"Trilinos Indices for rank "<<rank<<std::endl;

  for( localIndex a=0 ; a<trilinos_index.size() ; ++a )
  {
    if( is_ghost[a] < 0 )
    {
      std::cout<<rank<<", "<<a<<", "<<trilinos_index[a]<<std::endl;
    }
  }
   */
}


void LagrangeSolverBase :: SetupSystem ( PhysicalDomainT&  domain,
                                         SpatialPartition& partition )
{

  using namespace BoundaryConditionFunctions;

  // determine the global/local degree of freedom distribution.

  //  const auto& kinematicSibling = domain.m_feNodeManager.GetOneToOneMap("NodeToKinematicSibling");
  //  const auto& isNodeDead = domain.m_feNodeManager.GetFieldData<int>("isDead");
  //  iArray1d kinematicSibling(domain.m_feNodeManager.m_numNodes);

  dim = domain.m_feElementManager.m_ElementRegions.begin()->second.m_ElementDimension;
  int n_ghost_rows  = domain.m_feNodeManager.GetNumGhosts();
  int n_local_rows  = domain.m_feNodeManager.DataLengths()-n_ghost_rows;
  //  int n_hosted_rows = n_local_rows+n_ghost_rows;
  int n_global_rows = 0;

  iArray1d displacementIndices;
  SetNumRowsAndTrilinosIndices( domain,
                                partition,
                                n_local_rows,
                                n_global_rows,
                                displacementIndices,
                                0 );
  // create epetra map





#if USECPP11==1
  m_rowMap = std::make_shared<Epetra_Map>(dim*n_global_rows,dim*n_local_rows,0,*epetra_comm);
  m_sparsity = std::make_shared<Epetra_FECrsGraph>(Copy,*m_rowMap,0);
#else
  if( m_rowMap!=NULL )
  {
    delete m_rowMap;
  }
  m_rowMap = new Epetra_Map(dim*n_global_rows,dim*n_local_rows,0,*epetra_comm);

  if( m_sparsity!=NULL )
  {
    delete m_sparsity;
  }
  m_sparsity = new Epetra_FECrsGraph(Copy,*m_rowMap,100);

#endif



  dummyDof.resize(0);                                 

  RegionMap::iterator
  region     = domain.m_feElementManager.m_ElementRegions.begin(),
  end_region = domain.m_feElementManager.m_ElementRegions.end();

//  const Array1dT<int>* isDetachedFromSolidMesh = domain.m_feNodeManager.GetFieldDataPointer<int> ("isDetachedFromSolidMesh");

//  if(domain.m_externalFaces.m_contactActive && domain.m_contactManager.m_use_contact_search)
//  {
//    UpdateContactDataStructures(domain, false);
//  }

  if(domain.m_externalFaces.m_contactActive)
  {
    InsertGlobalIndices( domain);
    const bool planeStress =  LagrangeSolverBase::m_2dOption==LagrangeSolverBase::PlaneStress ;
    if(domain.m_contactManager.m_nitsche_active)
      domain.m_externalFaces.GetProjectionTensorAndWeightingAndStabilizationParameters(dim, planeStress, domain);
  }

#ifdef SRC_INTERNAL2
  if(domain.m_xfemManager != NULL)
  {
    SetSparsityPatternXFEM(domain);
  }
  else
#endif
  {
    SetSparsityPattern( domain );
  }

  m_sparsity->GlobalAssemble();
  m_sparsity->OptimizeStorage();

//  std::cout<<"m_sparsity->NumGlobalRows()     = "<<m_sparsity->NumGlobalRows()<<std::endl;
//  std::cout<<"m_sparsity->NumGlobalNonzeros() = "<<m_sparsity->NumGlobalNonzeros()<<std::endl;
}

void LagrangeSolverBase::SetSparsityPattern( PhysicalDomainT& domain )
{
  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);

  RegionMap::iterator
  region     = domain.m_feElementManager.m_ElementRegions.begin(),
  end_region = domain.m_feElementManager.m_ElementRegions.end();

  const Array1dT<int>* isDetachedFromSolidMesh = domain.m_feNodeManager.GetFieldDataPointer<int> ("isDetachedFromSolidMesh");

  for(; region != end_region; ++region)
  {
    ElementRegionT& elemRegion = region->second;
    const unsigned numNodesPerElement = elemRegion.m_numNodesPerElem;

    iArray1d elementLocalDofIndex (dim*numNodesPerElement);

    for(localIndex element = 0; element < elemRegion.m_numElems; ++element)
    {
      const localIndex* const localNodeIndices = elemRegion.m_toNodesRelation[element];

      for(unsigned i=0; i<numNodesPerElement; ++i)
      {
        const localIndex localNodeIndex = localNodeIndices[i];
        for( int d=0 ; d<dim ; ++d )
        {
          elementLocalDofIndex[i*dim+d] = dim*trilinos_index[localNodeIndex]+d;
        }
      }

      m_sparsity->InsertGlobalIndices(elementLocalDofIndex.size(),
                                    &elementLocalDofIndex.front(),
                                    elementLocalDofIndex.size(),
                                    &elementLocalDofIndex.front());
    }

    //HACK: This is to fix the all-zero rows associated with dead nodes.
    if (isDetachedFromSolidMesh != NULL)
    {
      for (localIndex iNd = 0; iNd != domain.m_feNodeManager.DataLengths(); ++iNd)
      {
        if ( (*isDetachedFromSolidMesh)[iNd] == 1 )
        {
          for(unsigned i=0; i<numNodesPerElement; ++i)
          {
            for( int d=0 ; d<dim ; ++d )
            {
              elementLocalDofIndex[i*dim+d] = dim*trilinos_index[iNd]+d;
            }
          }
          m_sparsity->InsertGlobalIndices(elementLocalDofIndex.size(),
                                        &elementLocalDofIndex.front(),
                                        elementLocalDofIndex.size(),
                                        &elementLocalDofIndex.front());

        }
      }
    }
  }



//  m_sparsity->GlobalAssemble();
//  m_sparsity->OptimizeStorage();

}


void LagrangeSolverBase::SetSparsityPatternXFEM ( PhysicalDomainT& domain)
{
  const auto& kinematicSibling = domain.m_feNodeManager.GetOneToOneMap("NodeToKinematicSibling");
  const auto& isNodeDead = domain.m_feNodeManager.GetFieldData<int>("isDead");
  const Array1dT<int>* isDetachedFromSolidMesh = domain.m_feNodeManager.GetFieldDataPointer<int> ("isDetachedFromSolidMesh");
  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);

  RegionMap::iterator
  region     = domain.m_feElementManager.m_ElementRegions.begin(),
  end_region = domain.m_feElementManager.m_ElementRegions.end();

  for(; region != end_region; ++region)
  {
    ElementRegionT& elemRegion = region->second;
    const unsigned numNodesPerElement = elemRegion.m_numNodesPerElem;

    iArray1d elementLocalDofIndex (dim*numNodesPerElement);
    iArray1d diagIndices(dim);

    for(localIndex element = 0; element < elemRegion.m_numElems; ++element)
    {
      const auto elem_is_physical = elemRegion.GetFieldDataPointer<int>("isPhysical");

      if((*elem_is_physical)[element])
      {
        const localIndex* const localNodeIndices = elemRegion.m_toNodesRelation[element];

        for(unsigned i=0; i<numNodesPerElement; ++i)
        {
          localIndex localNodeIndex = localNodeIndices[i];

          if( kinematicSibling(localNodeIndex)!=static_cast<localIndex>(-1) )
          {
            for (int d = 0; d<dim; ++d)
            {
              diagIndices[d] = dim*trilinos_index[localNodeIndex]+d;

              m_sparsity->InsertGlobalIndices(1,
                                              &diagIndices[d],
                                              1,
                                              &diagIndices[d]);
            }

            localNodeIndex = kinematicSibling(localNodeIndex);
          }
          for( int d=0 ; d<dim ; ++d )
          {
            elementLocalDofIndex[i*dim+d] = dim*trilinos_index[localNodeIndex]+d;
          }
        }

        m_sparsity->InsertGlobalIndices(elementLocalDofIndex.size(),
                                        &elementLocalDofIndex.front(),
                                        elementLocalDofIndex.size(),
                                        &elementLocalDofIndex.front());
      }
    }

    //HACK: This is to fix the all-zero rows associated with dead nodes.
    if (isDetachedFromSolidMesh != NULL)
    {
      for (localIndex iNd = 0; iNd != domain.m_feNodeManager.DataLengths(); ++iNd)
      {
        if ( (*isDetachedFromSolidMesh)[iNd] == 1 )
        {
          for(unsigned i=0; i<numNodesPerElement; ++i)
          {
            for( int d=0 ; d<dim ; ++d )
            {
              elementLocalDofIndex[i*dim+d] = dim*trilinos_index[iNd]+d;
            }
          }
          m_sparsity->InsertGlobalIndices(elementLocalDofIndex.size(),
                                          &elementLocalDofIndex.front(),
                                          elementLocalDofIndex.size(),
                                          &elementLocalDofIndex.front());
        }
      }
    }

    for (localIndex iNd = 0; iNd != domain.m_feNodeManager.DataLengths(); ++iNd)
    {
      if(isNodeDead[iNd] == 1)
      {
        for (int d = 0; d<dim; ++d)
        {
          diagIndices[d] = dim*trilinos_index[iNd]+d;

          m_sparsity->InsertGlobalIndices(1,
                                          &diagIndices[d],
                                          1,
                                          &diagIndices[d]);
        }
      }
    }
  }
}


realT LagrangeSolverBase :: Assemble ( PhysicalDomainT&  domain,
                                       const SpatialPartition& partition,
                                       const realT time,
                                       realT const dt )
{
  realT maxForce = 0.0;
  const Array1dT<int>* isDetachedFromSolidMesh = domain.m_feNodeManager.GetFieldDataPointer<int> ("isDetachedFromSolidMesh");

  using namespace BoundaryConditionFunctions;
  // (re-)init linear system

  const Array1dT<R1Tensor>& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
  const Array1dT<R1Tensor>& uhat = domain.m_feNodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
  //  const Array1dT<R1Tensor>& vel = domain.m_feNodeManager.GetFieldData<FieldInfo::velocity>();

  Array1dT<R1Tensor> const * uhattilde = nullptr;
  Array1dT<R1Tensor> const * vtilde = nullptr;

  if( this->m_timeIntegrationOption == ImplicitDynamic )
  {
    uhattilde = domain.m_feNodeManager.GetFieldDataPointer<R1Tensor>("inc_disp_tilde");
    vtilde = domain.m_feNodeManager.GetFieldDataPointer<R1Tensor>("vel_tilde");
  }

  // basic nodal data ( = dof data for our problem)

  const iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);

  if(domain.m_contactManager.m_contact && domain.m_contactManager.m_implicitContactActive)
  {
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // @author: Chandra Annavarapu
    // @brief: Stiffness contributions to tie fractured faces together
    // Penalty/Nitsche's method to enforce contact - choice through a boolean flag in the xml file
    // For penalty method, a user-defined penalty parameter needs to be specified in the xml file
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    GetContactStiffnessContribution(domain);
  }




  static Array1dT< R1Tensor > u_local(8);
  static Array1dT< R1Tensor > uhat_local(8);
  static Array1dT< R1Tensor > vtilde_local(8);
  static Array1dT< R1Tensor > uhattilde_local(8);



  // begin region loop

  RegionMap::const_iterator
  region     = domain.m_feElementManager.m_ElementRegions.begin(),
  end_region = domain.m_feElementManager.m_ElementRegions.end();


  for(; region != end_region; ++region)
  {
    const ElementRegionT& elemRegion = region->second;

    Array1dT<R2SymTensor> const * const refStress = elemRegion.GetFieldDataPointer<R2SymTensor>("referenceStress");

    if (elemRegion.m_elementType.compare("flow_only"))
    {

      const unsigned int dof_per_fe = elemRegion.m_finiteElement->dofs_per_element();

      u_local.resize(dof_per_fe);
      uhat_local.resize(dof_per_fe);
      vtilde_local.resize(dof_per_fe);
      uhattilde_local.resize(dof_per_fe);

      // following should hold for scalar node problems

      assert(dof_per_fe == elemRegion.m_numNodesPerElem);

      // space for element matrix and rhs

      Epetra_IntSerialDenseVector  elementLocalDofIndex   (dim*dof_per_fe);
      Epetra_SerialDenseVector     element_rhs     (dim*dof_per_fe);
      Epetra_SerialDenseMatrix     element_matrix  (dim*dof_per_fe,dim*dof_per_fe);

      Epetra_SerialDenseVector     element_dof_np1  (dim*dof_per_fe);

      realT reasonableDiagValue = 0.0 ;
      // determine ghost elements

      const iArray1d& elem_is_ghost = elemRegion.GetFieldData<FieldInfo::ghostRank>();
      const auto elem_is_physical = elemRegion.GetFieldDataPointer<int>("isPhysical");
      //      const iArray1d& elem_is_physical = elemRegion.GetFieldData<int>("isPhysical");

      // begin element loop, skipping ghost elements

      for(localIndex element = 0; element < elemRegion.m_numElems; ++element)
      {
        bool isElemPhysical;
        if(elem_is_physical !=NULL)
          isElemPhysical = (*elem_is_physical)[element];
        else
          isElemPhysical = 1;

        if(elem_is_ghost[element] < 0 && isElemPhysical)
        {

          const localIndex* const localNodeIndices = elemRegion.m_toNodesRelation[element];

          for(unsigned i=0; i<elemRegion.m_numNodesPerElem; ++i)
          {

            localIndex localNodeIndex = localNodeIndices[i];

#ifdef SRC_INTERNAL2
            if(domain.m_xfemManager != NULL)
            {
              const auto& kinematicSibling = domain.m_feNodeManager.GetOneToOneMap("NodeToKinematicSibling");

              if( kinematicSibling(localNodeIndex)!=static_cast<localIndex>(-1) )
              {
                localNodeIndex = kinematicSibling(localNodeIndex);
              }
            }
#endif

            for( int d=0 ; d<dim ; ++d )
            {
              elementLocalDofIndex[i*dim+d] = dim*trilinos_index[localNodeIndex]+d;

              // TODO must add last solution estimate for this to be valid
              element_dof_np1(i*dim+d) = disp[localNodeIndex][d];
            }
          }

          if( this->m_timeIntegrationOption == ImplicitDynamic )
          {
            CopyGlobalToLocal( localNodeIndices,
                               disp, uhat, *vtilde, *uhattilde,
                               u_local, uhat_local, vtilde_local, uhattilde_local );
          }
          else
          {
            CopyGlobalToLocal( localNodeIndices,
                               disp, uhat,
                               u_local, uhat_local );
          }

          // assemble into global system
          const localIndex paramIndex = elemRegion.m_mat->NumParameterIndex0() > 1 ? element : 0 ;

          R2SymTensor const * const referenceStress = refStress==nullptr ? nullptr : &((*refStress)[element]);
          realT maxElemForce = CalculateElementResidualAndDerivative( *(elemRegion.m_mat->ParameterData(paramIndex)),
                                                                      *(elemRegion.m_finiteElement),
                                                                      elemRegion.m_dNdX[element],
                                                                      elemRegion.m_detJ[element],
                                                                      referenceStress,
                                                                      u_local,
                                                                      uhat_local,
                                                                      uhattilde_local,
                                                                      vtilde_local,
                                                                      dt,
                                                                      element_matrix,
                                                                      element_rhs);


          if( maxElemForce > maxForce )
            maxForce = maxElemForce;

          //          element_rhs.Scale( -1 );
          reasonableDiagValue += fabs(element_matrix(0,0));

          m_matrix->SumIntoGlobalValues(elementLocalDofIndex,
                                        element_matrix);


          m_rhs->SumIntoGlobalValues(elementLocalDofIndex,
                                     element_rhs);

          //          std::string sMatrix;
          //          std::stringstream out;
          //          out << element;
          //
          //          sMatrix = "system-matrix.dat" + out.str();
          //          EpetraExt::RowMatrixToMatlabFile(sMatrix.c_str(),*m_matrix);
          //
          //          std::string sRhs;
          //          sRhs = "system-RHS.dat" + out.str();
          //          EpetraExt::MultiVectorToMatlabFile(sRhs.c_str(),*m_rhs);

          //          std::string sMatrix;
          //          sMatrix = "systemMatrix.dat";
          //          EpetraExt::RowMatrixToMatlabFile(sMatrix.c_str(),*m_matrix);

          //          std::string sRhs;
          //          sRhs = "system-RHS-branch.dat";
          //          EpetraExt::MultiVectorToMatlabFile(sRhs.c_str(),*m_rhs);
        }
      }  // end element

      //      std::string sMatBefDN = "system-matrixBeforeBC.dat";
      //      EpetraExt::RowMatrixToMatlabFile(sMatBefDN.c_str(),*m_matrix);
      //HACK: This is to fix the all-zero rows associated with dead nodes.
      if (region->second.m_numElems > 0)  reasonableDiagValue /= region->second.m_numElems;

      if (isDetachedFromSolidMesh != NULL)
      {
        for (localIndex iNd = 0; iNd != domain.m_feNodeManager.DataLengths(); ++iNd)
        {
          if ( (*isDetachedFromSolidMesh)[iNd] == 1 )
          {
            element_matrix.Scale(0.0);
            element_matrix(0,0) = reasonableDiagValue;
            element_matrix(1,1) = reasonableDiagValue;

            for(unsigned i=0; i<elemRegion.m_numNodesPerElem; ++i)
            {
              for( int d=0 ; d<dim ; ++d )
              {
                elementLocalDofIndex[i*dim+d] = dim*trilinos_index[iNd]+d;
              }
            }
            m_matrix->SumIntoGlobalValues(elementLocalDofIndex,
                                          element_matrix);


          }
        }
      }

#ifdef SRC_INTERNAL2
      if(domain.m_xfemManager != NULL)
      {
        const auto& isNodeDead = domain.m_feNodeManager.GetFieldData<int>("isDead");
        const auto& kinematicSibling = domain.m_feNodeManager.GetOneToOneMap("NodeToKinematicSibling");

        for (localIndex iNd = 0; iNd != domain.m_feNodeManager.DataLengths(); ++iNd)
        {
          if ( kinematicSibling(iNd)!=static_cast<localIndex>(-1) || isNodeDead(iNd) == 1)
          {
            element_matrix.Scale(0.0);
            element_matrix(0,0) = reasonableDiagValue;
            element_matrix(1,1) = reasonableDiagValue;

            for(unsigned i=0; i<elemRegion.m_numNodesPerElem; ++i)
            {
              for( int d=0 ; d<dim ; ++d )
              {
                elementLocalDofIndex[i*dim+d] = dim*trilinos_index[iNd]+d;
              }
            }
            m_matrix->SumIntoGlobalValues(elementLocalDofIndex,
                                          element_matrix);


          }
        }
      }
#endif

      //      std::string sMat = "system-matrixBeforeBC.dat";
      //      EpetraExt::RowMatrixToMatlabFile(sMat.c_str(),*m_matrix);
      //
      //      std::cout<<"reasonableDiagValue = "<<reasonableDiagValue<<std::endl;
    }
  } // end region


  if( m_gravityVector.L2_Norm() > 0.0 )
  {// Now we apply the gravity on nodes regardless whether the elements making the contribution are alive or dead.
    // Need to rewrite this into an element loop when we need to deal with inactive or dead elements.
    Epetra_IntSerialDenseVector  localDofIndex   (dim*domain.m_feNodeManager.DataLengths());
    Epetra_SerialDenseVector     localRHS     (dim*domain.m_feNodeManager.DataLengths());

    const realT* const gravityVector = m_gravityVector.Data();
    const iArray1d& nodeIsGhost = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
    const Array1dT<realT>& mass = domain.m_feNodeManager.GetFieldData<FieldInfo::mass>();


    for(unsigned i=0; i<domain.m_feNodeManager.DataLengths(); ++i)
    {
      if( nodeIsGhost[i] < 0 )
      {
        for( int d=0 ; d<dim ; ++d )
        {

          localDofIndex[i*dim+d] = dim*trilinos_index[i]+d;
          localRHS[i*dim+d] = gravityVector[d] * mass[i];
        }
      }
    }
    m_rhs->SumIntoGlobalValues(localDofIndex,
                               localRHS);

  }



  // We follow the procedure outlined by section 2.10 in
  // Cook, R. D., Malkus, D.S., Plesha, M.E., Witt, R.J., 2001. Concepts and Applications of Finite Element Analysis. John Wiley & Sons, In
  // First we apply the nodal forces corresponding to that develops in a zero-deformation body subjected to the prescribed temperature field as boundary condition
  // Then we subtract this thermal stress in the calculation of stresses.
  // Adding nodal forces caused by thermal stress

  LagrangeSolverBase::ApplyTemperatureBC(domain, time);
  for(localIndex i = 0; i < m_thermalRegionNames.size(); ++i)
  {
    std::map<std::string, ElementRegionT>::iterator it = domain.m_feElementManager.m_ElementRegions.find(m_thermalRegionNames[i]);

    ElementRegionT& elemRegion = it -> second;

    const unsigned int dof_per_fe = elemRegion.m_finiteElement->dofs_per_element();
    assert(dof_per_fe == elemRegion.m_numNodesPerElem);

    const iArray1d& elem_is_ghost = elemRegion.GetFieldData<FieldInfo::ghostRank>();

    Epetra_SerialDenseVector     ele_rhs(dim*dof_per_fe);

    //EpetraExt::MultiVectorToMatlabFile("RHSPreThermal",*m_rhs);

    for(localIndex elemID = 0; elemID < elemRegion.m_numElems; ++elemID)
    {
      if(elem_is_ghost[elemID] < 0)
      {
        const localIndex* const localNodeIndices = elemRegion.m_toNodesRelation[elemID];

        ApplyThermalStress( elemRegion, domain.m_feNodeManager, elemID, &ele_rhs);

        Epetra_IntSerialDenseVector  elementLocalDofIndex(dim*dof_per_fe);

        for(localIndex j=0; j<elemRegion.m_numNodesPerElem; ++j)
        {
          const localIndex localNodeIndex = localNodeIndices[j];
          for( int d=0 ; d<dim ; ++d )
            elementLocalDofIndex[j*dim+d] = dim*trilinos_index[localNodeIndex]+d;
        }
        m_rhs->SumIntoGlobalValues(elementLocalDofIndex, ele_rhs);
      }
      else
      {
        //We need this is update the element temperature from nodal temperature (if applicable), because we need the temperature on ghosts for the final correction for thermal stress.
        ApplyThermalStress( elemRegion, domain.m_feNodeManager, elemID, &ele_rhs);
      }
    }
  }




  //  const rArray1d* const ppFlowPressure = domain.m_feFaceManager.GetFieldDataPointer<FieldInfo::pressure>();
  //  const rArray1d& dPdM = domain.m_feFaceManager.GetFieldData<realT>("dPdM");
  //  const rArray1d& density  = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  //  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  if(false) // Randy turned this off in his branch.  I am turned it on temporarily to pass the tests.
    if( domain.m_feFaceManager.HasField<int>("ruptureState") &&
        domain.m_feFaceManager.m_cohesiveZone )
    {



      Epetra_IntSerialDenseVector  faceLocalDofIndex;
      Epetra_SerialDenseVector     face_rhs;
      Epetra_SerialDenseMatrix     face_matrix;

      // determine ghost elements

      const iArray1d& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();


      // apply cohesive forces
      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
      Array1dT<R1Tensor>& cohesiveForce = domain.m_feNodeManager.GetFieldData<R1Tensor>("cohesiveForce");

      R1Tensor cohesiveTraction;
      iArray1d& ruptureState = domain.m_feFaceManager.GetFieldData<int>("ruptureState");

      cohesiveForce = 0.0;



      // begin element loop, skipping ghost elements

      for(localIndex kf = 0 ; kf < domain.m_feFaceManager.m_numFaces ; ++kf)
      {
        if(isGhost[kf] < 0)
          if( ruptureState[kf] >= 2 && !(childFaceIndex[kf].empty()) )
          {


            const localIndex faceIndex[2] = { kf, childFaceIndex[kf][0] };

            const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                                    domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

            R1Tensor Nbar = N[0];
            Nbar -= N[1];
            Nbar.Normalize();

            const unsigned int numNodes = domain.m_feFaceManager.m_toNodesRelation[kf].size();
            {
              faceLocalDofIndex.Resize(dim*domain.m_feFaceManager.m_toNodesRelation[kf].size()*2);
              face_rhs.Resize(dim*domain.m_feFaceManager.m_toNodesRelation[kf].size()*2);
              face_matrix.Reshape(dim*domain.m_feFaceManager.m_toNodesRelation[kf].size()*2,dim*domain.m_feFaceManager.m_toNodesRelation[kf].size()*2);

              for( localIndex a=0 ; a<numNodes ; ++a )
              {
                const localIndex aa = a == 0 ? a : numNodes - a;
                const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[0]][a];
                const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[1]][aa];
                for( int d=0 ; d<dim ; ++d )
                {
                  faceLocalDofIndex[a*dim*2+d]       = dim*trilinos_index[localNodeIndex1]+d;
                  faceLocalDofIndex[a*dim*2+dim+d]   = dim*trilinos_index[localNodeIndex2]+d;
                }
              }


              R1Tensor gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, kf );
              R2Tensor cStiffness;


              //        ruptureState[kf] = domain.m_feFaceManager.m_cohesiveZone->UpdateCohesiveZone( kf, gap, Nbar, cohesiveTraction, cStiffness );
              ruptureState[kf] = domain.m_feFaceManager.m_cohesiveZone->UpdateCohesiveZone( kf, gap, Nbar,
                                                                                            domain.m_feFaceManager.m_toElementsRelation[faceIndex[0]][0],
                                                                                            domain.m_feFaceManager.m_toElementsRelation[faceIndex[1]][0],
                                                                                            cohesiveTraction, cStiffness );


              const localIndex faceID = faceIndex[0];
              const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, faceID );

              //        cohesiveTraction = 0;
              //        cStiffness = 0;

              IntegrateCohesiveZoneContributions( domain.m_feFaceManager.m_toNodesRelation[kf].size(),
                                                  area,
                                                  cohesiveTraction,
                                                  cStiffness,
                                                  face_matrix,
                                                  face_rhs );

              // assemble into global system


              m_matrix->SumIntoGlobalValues(faceLocalDofIndex,
                                            face_matrix);

              m_rhs->SumIntoGlobalValues(faceLocalDofIndex,
                                         face_rhs);

              for( localIndex a=0 ; a<numNodes ; ++a )
              {
                const localIndex aa = a == 0 ? a : numNodes - a;
                const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[0]][a];
                const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[1]][aa];

                for( int i=0 ; i<dim ; ++i )
                {
                  cohesiveForce[localNodeIndex1][i] -= face_rhs[a*dim*2+i];
                  cohesiveForce[localNodeIndex2][i] -= face_rhs[a*dim*2+dim+i];
                }
              }
            }
          }
      }
    }  // end element



  ApplyDisplacementPenalty( domain.m_feNodeManager );




  // Global assemble
  m_matrix->GlobalAssemble(true);
  m_rhs->GlobalAssemble();



  /*
  Array1dT<R1Tensor>& force = domain.m_feNodeManager.GetFieldData<FieldInfo::force>();
  const iArray1d& is_ghost       = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();

  int dummy;
  double* rhs = NULL;

  m_rhs->ExtractView(&rhs,&dummy);

  for(unsigned r=0; r<force.size(); ++r)
  {
    if(is_ghost[r] < 0)
    {
      for( int d=0 ; d<dim ; ++d )
      {
        int intDofOffset = 0;
        int lid = m_rowMap->LID(intDofOffset+dim*trilinos_index[r]+d);
        force[r][d] = rhs[lid];
      }
    }
  }
   */

  // Apply kinematic constraints
  //  ApplyBoundaryCondition<R1Tensor>(this, &LagrangeSolverBase::kinematicBC,
  //                                   domain, domain.m_feNodeManager,
  //                                   Field<FieldInfo::displacement>::Name(), time );

  // Apply boundary conditions
  ApplyBoundaryCondition<R1Tensor>(this, &LagrangeSolverBase::TractionBC,
                                   domain, domain.m_feFaceManager, "Traction", time );
  //  ApplyBoundaryCondition<R1Tensor>(this, &LagrangeSolverBase::PressureBC,
  //                                   domain, domain.m_feFaceManager, "Pressure", time );
  ApplyBoundaryCondition<R1Tensor>(this, &LagrangeSolverBase::DisplacementBC,
                                   domain, domain.m_feNodeManager,
                                   Field<FieldInfo::displacement>::Name(), time );

  //m_rhs->Print(std::cout);
  //std::cout << matrix->NormInf() << std::endl;
  //EpetraExt::RowMatrixToMatlabFile("system-matrix.dat",*matrix);

  // exit(0);

  return maxForce;
}









/*

realT LagrangeSolverBase::AssembleFluidPressureContributions( PhysicalDomainT& domain,
                                                                  const iArray1d& deformationTrilinosIndex,
                                                                  const iArray1d& flowTrilinosIndex,
                                                                  const Array1dT< rArray1d >& dwdu,
                                                                  const rArray1d& dwdw,
                                                                  const int flowDofOffset  )
{

  realT maxForce = 0.0 ;

  const rArray1d* const ppFlowPressure = domain.m_feFaceManager.GetFieldDataPointer<FieldInfo::pressure>();

  if( ppFlowPressure!=NULL )
  {
    Array1dT<R1Tensor>& hydroForce    = domain.m_feNodeManager.GetFieldData<R1Tensor>("hydroForce");
    Array1dT<R1Tensor>& hgForce    = domain.m_feNodeManager.GetFieldData<FieldInfo::hgforce>();
    hydroForce = 0;
    hgForce = 0;
    const rArray1d& density_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();

    Epetra_IntSerialDenseVector  rowDofIndex;
    Epetra_IntSerialDenseVector  colDofIndex;


    Epetra_SerialDenseVector     face_rhs;
    Epetra_SerialDenseMatrix     face_matrix;

    // determine ghost elements

    const iArray1d& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();


    // apply cohesive forces
    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );


    const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");


    const rArray1d& dPdM = domain.m_feFaceManager.GetFieldData<realT>("dPdM");

//    const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();


//    const Array1dT<R1Tensor>& effectiveNodalNormal = domain.m_feNodeManager.GetFieldData<R1Tensor>( "effectiveNormal");


    for(localIndex kf = 0 ; kf < domain.m_feFaceManager.m_numFaces ; ++kf)
    {
      if( isGhost[kf] < 0 && flowFaceType[kf]==0 )
      {



        const localIndex faceIndex[2] = { kf, childFaceIndex[kf][0] };

        const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                                domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

        R1Tensor Nbar = N[0];
        Nbar -= N[1];
        Nbar.Normalize();


        const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf].size();

        {
          const int nrows = 2*dim*numNodes;
          const int ncols = 2*dim*numNodes + 1;
          rowDofIndex.Resize(nrows);
          colDofIndex.Resize(ncols);
          face_rhs.Resize(nrows);
          face_matrix.Reshape(nrows,ncols);
          face_rhs.Scale(0.0);
          face_matrix.Scale(0.0);

          for( localIndex a=0 ; a<numNodes ; ++a )
          {
            const localIndex aa = a == 0 ? a : numNodes - a;
            const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[0]][a];
            const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[1]][aa];

            const int aDof = 2*a*dim;
            for( int i=0 ; i<dim ; ++i )
            {

              rowDofIndex[aDof+i]       = dim*deformationTrilinosIndex[localNodeIndex1]+i;
              rowDofIndex[aDof+dim+i]   = dim*deformationTrilinosIndex[localNodeIndex2]+i;

              colDofIndex[aDof+i]       = dim*deformationTrilinosIndex[localNodeIndex1]+i;
              colDofIndex[aDof+dim+i]   = dim*deformationTrilinosIndex[localNodeIndex2]+i;

            }
          }
          colDofIndex[ncols-1] = flowDofOffset + flowTrilinosIndex[kf];

          const realT dPdV = -dPdM[kf]*density_np1[kf];

          // todo use correct area
          const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, kf );

          const realT hydroForceMag = (*ppFlowPressure)[kf] * area;
          const realT hydroStiffnessMag = dPdM[kf] * area;


          for( localIndex a=0 ; a<numNodes ; ++a )
          {
            const localIndex aa = a == 0 ? a : numNodes - a;
            const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[0]][a];
            const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[1]][aa];




            // rhs and mass derivative

            const int aDof = 2*a*dim;
            for( int i=0 ; i<dim ; ++i )
            {
              face_rhs[aDof+i]     +=  hydroForceMag * Nbar[i] / numNodes ;
              face_rhs[aDof+dim+i] += -hydroForceMag * Nbar[i] / numNodes ;

              maxForce = max( maxForce, fabs(hydroForceMag * Nbar[i] / numNodes) );

              hydroForce[localNodeIndex1][i] += -hydroForceMag * Nbar[i] / numNodes;
              hydroForce[localNodeIndex2][i] +=  hydroForceMag * Nbar[i] / numNodes;

              face_matrix(aDof+i,ncols-1)     += -hydroStiffnessMag * Nbar[i] / numNodes;
              face_matrix(aDof+dim+i,ncols-1) +=  hydroStiffnessMag * Nbar[i] / numNodes;
            }



            // displacement derivative
            for( localIndex b=0 ; b<numNodes ; ++b)
            {
              const int bDof = 2*b*dim;
              for( int i=0 ; i<dim ; ++i )
              {
                for( int j=0 ; j<dim ; ++j )
                {
                  face_matrix(aDof+i,bDof+j)          -= dPdV*area*dwdw(kf)*dwdu(kf)(bDof+j)*Nbar[i]     * ( area / numNodes );
                  face_matrix(aDof+i,bDof+dim+j)      -= dPdV*area*dwdw(kf)*dwdu(kf)(bDof+dim+j)*Nbar[i] * ( area / numNodes );
                  face_matrix(aDof+dim+i,bDof+j)      += dPdV*area*dwdw(kf)*dwdu(kf)(bDof+j)*Nbar[i]     * ( area / numNodes );
                  face_matrix(aDof+dim+i,bDof+dim+j)  += dPdV*area*dwdw(kf)*dwdu(kf)(bDof+dim+j)*Nbar[i] * ( area / numNodes );

                }
              }
            }


          }

          m_matrix->SumIntoGlobalValues(rowDofIndex,
                                        colDofIndex,
                                        face_matrix);

          m_rhs->SumIntoGlobalValues(rowDofIndex,
                                     face_rhs);

        }
      }
    }
  }  // end element

  m_matrix->GlobalAssemble(true);
  m_rhs->GlobalAssemble();

  return maxForce;
}
 */





void LagrangeSolverBase::IntegrateCohesiveZoneContributions( const int numNodesInFace,
                                                             const realT area,
                                                             const R1Tensor& traction,
                                                             const R2Tensor& stiffness,
                                                             Epetra_SerialDenseMatrix& face_matrix,
                                                             Epetra_SerialDenseVector& face_rhs )
{

  R1Tensor cohesiveNodalForce;
  cohesiveNodalForce.cA( area/numNodesInFace, traction );

  R2Tensor cohesiveNodalStiffness;
  cohesiveNodalStiffness.cA( area/numNodesInFace, stiffness );

  face_matrix.Scale(0.0);
  face_rhs.Scale(0.0);

  for( int a=0 ; a<numNodesInFace ; ++a )
  {
    for( int i=0 ; i<3 ; ++i )
    {
      for( int j=0 ; j<3 ; ++j )
      {
        face_matrix( a*dim*2+i     , a*dim*2+j     )  = -cohesiveNodalStiffness(i,j);
        face_matrix( a*dim*2+i     , a*dim*2+dim+j )  = cohesiveNodalStiffness(i,j);
        face_matrix( a*dim*2+dim+i , a*dim*2+j     )  = cohesiveNodalStiffness(i,j);
        face_matrix( a*dim*2+dim+i , a*dim*2+dim+j )  = -cohesiveNodalStiffness(i,j);
      }

      face_rhs( a*dim*2+i     ) = -cohesiveNodalForce[i];
      face_rhs( a*dim*2+dim+i ) = cohesiveNodalForce[i];
    }
  }
}


/* Solve */
void LagrangeSolverBase::SetupMLPreconditioner( const PhysicalDomainT& domain, ML_Epetra::MultiLevelPreconditioner* MLPrec )
{

  Teuchos::ParameterList MLList;

  ML_Epetra::SetDefaults("SA",MLList);
  //  MLList.set("aggregation: type", "Uncoupled");
  //  MLList.set("smoother: type","Chebyshev");
  //  MLList.set("smoother: sweeps",2);
  //  MLList.set("smoother: ifpack overlap",0);
  MLList.set("prec type","MGW");
  //  MLList.set("cycle applications",1);
  MLList.set("ML output", 5);
  MLList.set("PDE equations",3);

  // create the preconditioning object.
  rArray1d rigid;


  if(m_numerics.m_useMLPrecond)
  {
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*m_matrix, MLList);

    MLPrec->PrintUnused();

    const Array1dT<R1Tensor>& x = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition>();



    int num_PDE_eqns = 3;
    int Nrigid = 6;
    int N_update = domain.m_feNodeManager.DataLengths();
    rigid.resize( Nrigid * N_update * num_PDE_eqns );
    ML_Coord2RBM( domain.m_feNodeManager.DataLengths(),
                  x, rigid, num_PDE_eqns );

    MLList.set("null space: type","pre-computed");
    MLList.set("null space: vectors", rigid.data() );
    MLList.set("null space: dimension", Nrigid );

    //    ML_Aggregate_Set_NullSpace(ag, num_PDE_eqns, Nrigid, rigid.data(), N_update);


  }
}


static void ML_Coord2RBM( int Nnodes,
                          const Array1dT<R1Tensor>& x,
                          rArray1d& rbm,
                          const int Ndof )
{

  int vec_leng, ii, jj, offset, node, dof;

  vec_leng = Nnodes*Ndof;

  for( node = 0 ; node < Nnodes; node++ )
  {
    dof = node*Ndof;
    for(ii=0;ii<3;ii++)
    { /* upper left = [ I ] */
      for(jj=0;jj<3;jj++)
      {
        offset = dof+ii+jj*vec_leng;
        rbm[offset] = (ii==jj) ? 1.0 : 0.0;
      }
    }
    for(ii=0;ii<3;ii++)
    { /* upper right = [ Q ] */
      for(jj=3;jj<6;jj++)
      {
        offset = dof+ii+jj*vec_leng;
        if( ii == jj-3 )
        {
          rbm[offset] = 0.0;
        }
        else
        {
          if (ii+jj == 4) rbm[offset] = x[node][2];
          else if ( ii+jj == 5 ) rbm[offset] = x[node][1];
          else if ( ii+jj == 6 ) rbm[offset] = x[node][0];
          else rbm[offset] = 0.0;
        }
      }
    }
    ii = 0; jj = 5;
    offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
    ii = 1; jj = 3;
    offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;
    ii = 2; jj = 4;
    offset = dof+ii+jj*vec_leng; rbm[offset] *= -1.0;

  } /*for( node = 0 ; node < Nnodes; node++ )*/

}



void LagrangeSolverBase::ApplyForcesFromContact( PhysicalDomainT& domain,
                                                 StableTimeStep& timeStep,
                                                 const realT dt )
{
  //-----------------------------
  //external faces
  //-----------------------------
  {
    Array1dT<R1Tensor>& contactForce = domain.m_feNodeManager.GetFieldData<FieldInfo::contactForce>();
    contactForce = 0.0;
    Array1dT<R1Tensor>& decontactForce = domain.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::contactForce>();
    decontactForce = 0.0;

    //update nodal positions, velocities, and accelerations before updating face geometry
    //also, reset rotational and translational accelerations as well as forces and moments
    domain.m_discreteElementManager.UpdateNodalStatesZeroForcesAndAccelerations();

    //update face geometry and sort faces if necessary
    const bool resort = domain.m_externalFaces.RecalculateNeighborList(domain.m_feNodeManager,
                                                                       domain.m_discreteElementSurfaceNodes,
                                                                       domain.m_discreteElementManager,
                                                                       dt);

    //if a resort has been triggered, then you also need to update the contact manager
    if (resort)
      domain.m_contactManager.Update(domain.m_externalFaces.m_neighborList);

    {
      Array1dT<Array1dT<R1Tensor> > xs;
      xs.resize(domain.m_externalFaces.DataLengths());
      domain.m_externalFaces.UpdateGeometricContactProperties(dt, domain, xs);
#ifdef STATES_ON_CONTACTS
      domain.m_externalFaces.UpdateAndApplyContactForcesFromStatesOnCommonPlanes(dt, domain);
#else
      domain.m_externalFaces.UpdateAndApplyContactStresses(timeStep, dt, domain, xs);
#endif
    }

    //for parallel: do NOT need to synchronize nodal force fields across processes for DE nodes
    //before moving to centroid, since each process will do that calculation (redundantly) itself
    // ... there's therefore no need for explicit synchrony
    //as long as nodal states remain synchronous, everything will be fine!
#ifndef DEEFC
    domain.m_discreteElementManager.ApplyDiscreteElementContactForces(45, 0.2);
#else
    domain.m_discreteElementManager.ApplyNodalForces();
#endif
  }


  //-----------------------------
  //ellipsoidal discrete elements
  //-----------------------------
  {
    //update nodal positions, velocities, and accelerations before updating face geometry
    //also, reset rotational and translational accelerations as well as forces and moments
    domain.m_ellipsoidalDiscreteElementManager.UpdateNodalStatesZeroForcesAndAccelerations();

    //update ellipsoidal discrete elements
    bool resort = domain.m_ellipsoidalDiscreteElementManager.RecalculateNeighborList(dt);

    //if a resort has been triggered, then you also need to update the contact manager
    if(resort)
      domain.m_ellipsoidalContactManager.Update(domain.m_ellipsoidalDiscreteElementManager.m_neighborList);

    domain.m_ellipsoidalDiscreteElementManager.UpdateAndApplyContactStresses(m_stabledt, dt, domain.m_ellipsoidalContactManager);
  }

}


void LagrangeSolverBase::ApplyTemperatureBC( PhysicalDomainT& domain, realT time )
{
  for(localIndex i = 0; i < m_thermalRegionNames.size(); ++i)
  {
    std::map<std::string, ElementRegionT>::iterator it = domain.m_feElementManager.m_ElementRegions.find(m_thermalRegionNames[i]);
    ElementRegionT& elemRegion = it -> second;
    BoundaryConditionFunctions::ApplyDirichletBoundaryCondition<realT>(elemRegion,
                                                                       "temperature",
                                                                       time);
  }
}


void LagrangeSolverBase::PostProcess (PhysicalDomainT& domain,
                                      SpatialPartition& partition,
                                      const sArray1d& namesOfSolverRegions)
{

  if (domain.m_feNodeManager.HasField<int>("isDetachedFromSolidMesh"))
  {
    domain.m_feNodeManager.UpdateDetachedNodeLocationAndVelocity();
  }

  if (m_writeNodalStress)
  {
    rArray1d& sigma_x = domain.m_feNodeManager.GetFieldData<realT>("nsigma_x");
    rArray1d& sigma_y = domain.m_feNodeManager.GetFieldData<realT>("nsigma_y");
    rArray1d& sigma_z = domain.m_feNodeManager.GetFieldData<realT>("nsigma_z");
    rArray1d& sigma_xy = domain.m_feNodeManager.GetFieldData<realT>("nsigma_xy");
    rArray1d& sigma_yz = domain.m_feNodeManager.GetFieldData<realT>("nsigma_yz");
    rArray1d& sigma_zx = domain.m_feNodeManager.GetFieldData<realT>("nsigma_zx");

    sigma_x = 0.0;
    sigma_y = 0.0;
    sigma_z = 0.0;
    sigma_xy = 0.0;
    sigma_yz = 0.0;
    sigma_zx = 0.0;

    for( sArray1d::const_iterator regionName = namesOfSolverRegions.begin() ;
        regionName != namesOfSolverRegions.end() ; ++regionName )
    {
      std::map<std::string, ElementRegionT>::iterator iter = domain.m_feElementManager.m_ElementRegions.find(*regionName);
      if(iter == domain.m_feElementManager.m_ElementRegions.end())
        continue;

      ElementRegionT& elementRegion = iter->second;


      for( localIndex k=0 ; k<elementRegion.m_numElems ; ++k )
      {
        R2SymTensor s;
        realT pressure = 0.0;
        realT ex, ey, ez, exy, eyz, ezx;
        for( localIndex a=0 ; a<elementRegion.m_numIntegrationPointsPerElem ; ++a )
        {
          const MaterialBaseStateData& state = *(elementRegion.m_mat->StateData(k,a));
          s += state.devStress;
          pressure += state.pressure;
        }
        s /= elementRegion.m_numIntegrationPointsPerElem;
        pressure /= elementRegion.m_numIntegrationPointsPerElem;
        ex =  s(0,0) + pressure;
        ey =  s(1,1) + pressure;
        ez =  s(2,2) + pressure;
        exy =  s(0,1);
        eyz =  s(1,2);
        ezx =  s(0,2);

        for (localIndex j = 0; j < elementRegion.m_toNodesRelation.Dimension(1); ++j)
        {
          localIndex ind = elementRegion.m_toNodesRelation[k][j];
          sigma_x[ind] +=ex;
          sigma_y[ind] +=ey;
          sigma_z[ind] +=ez;
          sigma_xy[ind] +=exy;
          sigma_yz[ind] +=eyz;
          sigma_zx[ind] +=ezx;
        }
      }

      for (localIndex i=0; i<domain.m_feNodeManager.DataLengths(); ++i)
      {
        if (domain.m_feNodeManager.m_toElementsRelation[i].size() > 0)
        {
          sigma_x[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
          sigma_y[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
          sigma_z[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
          sigma_xy[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
          sigma_yz[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
          sigma_zx[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
        }
      }

      {
        std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;

        if (m_writeNodalStress)
        {
          syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_x");
          syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_y");
          syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_z");
          syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_xy");
          syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_yz");
          syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_zx");
        }

        partition.SynchronizeFields( syncedFields, CommRegistry::lagrangeSolver01);
      }
    }
  }
}

void LagrangeSolverBase::SnapshotNodalDisplacement( NodeManagerT& nodeManager)
{
  Array1dT<R1Tensor>& refDisplacement = nodeManager.GetFieldData<R1Tensor>("refDisplacement");
  Array1dT<R1Tensor>& u = nodeManager.GetFieldData<FieldInfo::displacement>();
  for (localIndex i = 0; i < nodeManager.DataLengths(); ++i)
  {
    refDisplacement[i] = u[i];
  }
}

// Make it more general later. For now works only for triangles
void LagrangeSolverBase::EvaluateDispAtVertex(PhysicalDomainT&  domain,
                                              const SpatialPartition& partition,
                                              const realT time )
{
  RegionMap::const_iterator
  region     = domain.m_feElementManager.m_ElementRegions.begin(),
  end_region = domain.m_feElementManager.m_ElementRegions.end();

  auto& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();
  const auto& refPos = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition>();
  const auto& isNodeVertex = domain.m_feNodeManager.GetFieldData<int>("isVertex");
  iArray1d vertexAccounted(domain.m_feNodeManager.m_numNodes);

  for(; region != end_region; ++region)
  {
    const ElementRegionT& elemRegion = region->second;
    //    auto& elemToCracks = elemRegion.GetOneToOneMap("ElementToCracks");
    auto& elemToVertexNodes = elemRegion.GetVariableOneToManyToManyMap("ElementToCrackToVertexNodes");
    auto& elementToPhysicalNodes = elemRegion.GetVariableOneToManyMap("ElementToPhysicalNodes");

    if (!elemRegion.m_elementGeometryID.compare(0, 4, "STRI"))
    {
      if (elemRegion.m_elementType.compare("flow_only"))
      {
        const iArray1d& elem_is_ghost = elemRegion.GetFieldData<FieldInfo::ghostRank>();
        const iArray1d& elem_is_physical = elemRegion.GetFieldData<int>("isPhysical");

        // begin element loop, skipping ghost elements

        for(localIndex element = 0; element < elemRegion.m_numElems; ++element)
        {
          if(elem_is_ghost[element] < 0 && elem_is_physical[element]==1 && !(elemToVertexNodes[element].empty()))
          {
            rArray2d Amat(3, 3); Amat(0,0) = 1; Amat(0,1) = 1; Amat(0,2) = 1;
            rArray2d invAmat(3,3);

            const localIndex* const localNodeIndices = elemRegion.m_toNodesRelation[element];

            for(unsigned i=0; i<elemRegion.m_numNodesPerElem; ++i)
            {
              const localIndex localNodeIndex = localNodeIndices[i];
              for( int d=0 ; d<dim ; ++d )
              {
                Amat(d+1,i) = refPos[localNodeIndex][d];
              }
            }

            GetInverseOfA3By3Mat(Amat,invAmat);

            for (localIndex iPNod = 0; iPNod < elementToPhysicalNodes[element].size(); iPNod++)
            {
              localIndex currentNod = elementToPhysicalNodes[element][iPNod];
              if(isNodeVertex[currentNod] && vertexAccounted[currentNod]==0)
              {
                rArray1d vertPos(3), shapeFuncPos(3); vertPos(0) = 1;
                for( int d=0 ; d<dim ; ++d )
                {
                  vertPos(d+1) = refPos[currentNod][d];
                }
                MultiplyArray(invAmat, vertPos, shapeFuncPos);

                disp[currentNod] = 0;
                for (localIndex iSF = 0; iSF<shapeFuncPos.size(); iSF++)
                {
                  const localIndex localNodeIndex = localNodeIndices[iSF];
                  for (int d = 0; d< dim; d++)
                  {
                    disp[currentNod][d] = disp[currentNod][d] + shapeFuncPos[iSF]*disp[localNodeIndex][d];
                  }
                }
                vertexAccounted[currentNod] = 1;
              }
            }
          }
        }  // end element
      }
    }
  } // end region
}

void LagrangeSolverBase::GetInverseOfA3By3Mat(const rArray2d& mat,
                                              rArray2d& invMat)
{
  realT o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11;

  o1 = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1);
  o2 = mat(0,2)*mat(2,1) - mat(0,1)*mat(2,2);
  o3 = mat(0,1)*mat(1,2) - mat(0,2)*mat(1,1);
  o4 = mat(1,2)*mat(2,0) - mat(1,0)*mat(2,2);
  o5 = mat(0,0)*mat(2,2) - mat(0,2)*mat(2,0);
  o6 = mat(0,2)*mat(1,0) - mat(0,0)*mat(1,2);
  o7 = mat(1,0)*mat(2,1) - mat(1,1)*mat(2,0);
  o8 = mat(0,1)*mat(2,0) - mat(0,0)*mat(2,1);
  o9 = mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0);

  // Get determinant
  o10 = mat(0,0)*o1 + mat(1,0)*o2 + mat(2,0)*o3;

  realT MaxVal = *std::max_element(mat.begin(),mat.end());

  const realT tol = 1.0e-14 * MaxVal;

  if( o10<=tol && o10>=-tol )
  {
    throw std::exception();
  }
  else
  {
    o11 = 1.0/o10;
  }

  // Get inverse
  invMat(0,0) = o1*o11;
  invMat(0,1) = o2*o11;
  invMat(0,2) = o3*o11;
  invMat(1,0) = o4*o11;
  invMat(1,1) = o5*o11;
  invMat(1,2) = o6*o11;
  invMat(2,0) = o7*o11;
  invMat(2,1) = o8*o11;
  invMat(2,2) = o9*o11;
}

void LagrangeSolverBase::MultiplyArray(const rArray2d& A,
                                       const rArray1d& B,
                                       rArray1d& C)
{
  const int m = A.Dimension(0), n = A.Dimension(1);
  C.resize(m);

  realT Temp=0;
  for( int iRow=0 ; iRow<m; iRow++)
  {
    Temp = 0.;
    for( int s=0; s< n; s++)
    {
      Temp += A(iRow,s)*B(s);
    }
    C(iRow) = Temp;
  }
}

void LagrangeSolverBase::GetNodalForces( PhysicalDomainT&  domain,
                                         const SpatialPartition& partition,
                                         const realT dt )
{
  RegionMap::iterator
  region     = domain.m_feElementManager.m_ElementRegions.begin(),
  end_region = domain.m_feElementManager.m_ElementRegions.end();

  //  auto& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();

  for(; region != end_region; ++region)
  {
    ElementRegionT& elemRegion = region->second;
    //    auto& elemToCracks = elemRegion.GetOneToOneMap("ElementToCracks");

    if (elemRegion.m_elementType.compare("flow_only"))
    {
      const iArray1d& elem_is_ghost = elemRegion.GetFieldData<FieldInfo::ghostRank>();
      const iArray1d& elem_is_physical = elemRegion.GetFieldData<int>("isPhysical");

      // begin element loop, skipping ghost elements

      for(localIndex element = 0; element < elemRegion.m_numElems; ++element)
      {
        if(elem_is_ghost[element] < 0 && elem_is_physical[element]==1)
        {
          R2SymTensor stress; realT pressureElm;
          elemRegion.m_mat->MeanPressureDevStress(element,pressureElm, stress);
          stress.PlusIdentity(pressureElm);

          CalculateNodalForces(domain.m_feNodeManager, elemRegion, element, stress);
        }
      }
    }
  }
}

void LagrangeSolverBase::CalculateNodalForces ( NodeManagerT& nodeManager,
                                                ElementRegionT& elemRegion,
                                                const localIndex elementID,
                                                R2SymTensor& stress)
{
  R2Tensor F;
  R2Tensor Finv;
  R2Tensor dUhatdX;


  static Array1dT< R1Tensor > u_local;
  static Array1dT< R1Tensor > uhat_local;
  static Array1dT<R1Tensor> f_local;

  rArray2d& detJ = elemRegion.m_detJ;

  Array2dT< R2Tensor >&  dUdX = elemRegion.m_dUdX;

  const unsigned int numNodesPerElem = elemRegion.m_numNodesPerElem;

  if( u_local.size() != numNodesPerElem )
  {
    u_local.resize(numNodesPerElem);
    uhat_local.resize(numNodesPerElem);
    f_local.resize(numNodesPerElem);
  }


  const Array1dT<R1Tensor>& incrementalDisplacement = nodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
  const Array1dT<R1Tensor>& totalDisplacement = nodeManager.GetFieldData<FieldInfo::displacement>();


  Array1dT<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force>();

  const localIndex* const elemToNodeMap = elemRegion.m_toNodesRelation[elementID];


  CopyGlobalToLocal( elemToNodeMap,
                     incrementalDisplacement, totalDisplacement,
                     uhat_local, u_local );

  f_local = 0.0;

  for( unsigned int a=0 ; a<elemRegion.m_numIntegrationPointsPerElem ; ++a )
  {

    R1Tensor* const dNdX = elemRegion.m_dNdX(elementID)[a];

    // Velocity Gradient

    // calculate dUhat/dX at beginning of step
    CalculateGradient( dUhatdX ,uhat_local, dNdX );


    // calculate gradient (end of step)
    dUdX(elementID,a) += dUhatdX;
    F = dUdX(elementID,a);
    F.PlusIdentity(1.0);
    realT detF = F.Det();


    Finv.Inverse(F);
    FiniteElementUtilities::Integrate( stress,
                                       elemRegion.m_dNdX(elementID)[a],
                                       detJ(elementID,a),
                                       detF,
                                       Finv,
                                       f_local.size(),
                                       f_local.data() );
  }


  AddLocalToGlobal( elemToNodeMap, f_local, force);
}

// Rui Wang, apply gravity on elements
void LagrangeSolverBase::ApplyGravity( NodeManagerT& nodeManager,
                                       ElementRegionT& elemRegion,
                                       const realT dt )
{
  Array1dT<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force>();
  Array1dT<realT>& mass = elemRegion.GetFieldData<FieldInfo::mass> ();
  Array1dT<R1Tensor> gravityForce;
  gravityForce.resize(elemRegion.m_numNodesPerElem);
  gravityForce *= 0.0;
  R1Tensor gravityForceElem;

  for( localIndex k=0 ; k < elemRegion.m_numElems ; ++k )
  {
    const localIndex* const elemToNodeMap = elemRegion.m_toNodesRelation[k];
    gravityForceElem = m_gravityVector * mass[k];
    for (localIndex i=0 ; i < elemRegion.m_numNodesPerElem ; ++i)
    {
      gravityForce[i] = gravityForceElem / elemRegion.m_numNodesPerElem;
    }
    AddLocalToGlobal( elemToNodeMap, gravityForce, force);
  }
}

void LagrangeSolverBase::ApplyDisplacementPenalty( NodeManagerT& nodeManager )
{
  if( m_displacementPenalty > 0.0 )
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Epetra_FECrsMatrix junk(*m_matrix);
    junk.PutScalar(0.0);

    Array1dT<R1Tensor> const & disp = nodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
//    Array1dT<R1Tensor> const & disp = nodeManager.GetFieldData<FieldInfo::displacement>();
    rArray1d const & volume = nodeManager.GetFieldData<FieldInfo::volume>();
    iArray1d const & trilinos_index = nodeManager.GetFieldData<int>(m_trilinosIndexStr);
    iArray1d const & ghostRank       = nodeManager.GetFieldData<FieldInfo::ghostRank>();

    Array1dT<realT>& penaltyForce = nodeManager.GetFieldData<realT>("displacementPenaltyForce");
    penaltyForce = 0.0;
    localIndex const numNodes = nodeManager.m_numNodes;
    Epetra_IntSerialDenseVector  elementLocalDofIndex   (dim);
    Epetra_SerialDenseVector     element_rhs     (dim);
    Epetra_SerialDenseMatrix     element_matrix  (dim,dim);

    element_matrix.Scale(0.0);

    realT maxForce = 0;
    for( localIndex a=0 ; a<numNodes ; ++a )
    {
      if( ghostRank[a] < 0 )
      {
        realT const stiffness = - m_displacementPenalty * volume[a];
        for( int d=0 ; d<dim ; ++d )
        {
          elementLocalDofIndex[d] = dim*trilinos_index[a]+d;
          element_matrix(d,d) = stiffness;
          element_rhs(d) = disp[a][d] * stiffness;

          maxForce = std::max( fabs(element_rhs(d)), maxForce );
          penaltyForce[a] = std::max( fabs(element_rhs(d)), penaltyForce[a] );
        }
        m_matrix->SumIntoGlobalValues(elementLocalDofIndex,
                                      element_matrix);

        junk.SumIntoGlobalValues(elementLocalDofIndex,
                                 element_matrix);

        m_rhs->SumIntoGlobalValues(elementLocalDofIndex,
                                   element_rhs);
      }
    }

    if( m_numerics.m_verbose )
    {
      realT normF = junk.NormInf()/m_matrix->NormInf();

      if( rank==0 )
      {
        std::cout<<"Disp Penalty force, stiffness ratio = "<<maxForce<<", "<<normF<<std::endl;
      }
    }
  }

}

void LagrangeSolverBase::CalculateAllNodalMassAndVolume( PhysicalDomainT& domain, SpatialPartition& partition )
{
  rArray1d& mass = domain.m_feNodeManager.GetFieldData<FieldInfo::mass>();
  rArray1d * const volume = domain.m_feNodeManager.GetFieldDataPointer<FieldInfo::volume>();
  iArray1d * isDetachedFromSolidMesh = domain.m_feNodeManager.GetFieldDataPointer<int> ("isDetachedFromSolidMesh");

  // Zero out values
  mass = 0.0;
  if( volume!=nullptr )
  {
    (*volume) = 0.0;
  }

  for( std::map< std::string, ElementRegionT >::iterator i=domain.m_feElementManager.m_ElementRegions.begin() ;
      i != domain.m_feElementManager.m_ElementRegions.end() ; ++i )
  {
    i->second.CalculateNodalMasses( domain.m_feNodeManager ) ;

//    iArray1d& partitionColor = i->second.GetFieldData<int>("partitionColor");
//    partitionColor = partition.Color();
  }

  // Set detached nodal mass to a small value
  if( isDetachedFromSolidMesh!=nullptr )
  {
    for( localIndex ii=0 ; ii<mass.size() ; ++ii )
    {
      if ((*isDetachedFromSolidMesh)[ii] == 1)
      {
        mass[ii] = 1e-99;
      }
    }
  }


  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;

  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::mass>::Name());
  if( volume!=nullptr )
  {
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::volume>::Name());
  }

  partition.SynchronizeFields( syncedFields, CommRegistry::lagrangeSolver01 );
}
