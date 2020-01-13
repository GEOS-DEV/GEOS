/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HydrofractureSolver.cpp
 *
 */


#include "HydrofractureSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

HydrofractureSolver::HydrofractureSolver( const std::string& name,
                                      Group * const parent ):
  SolverBase(name,parent),
  m_solidSolverName(),
  m_flowSolverName(),
  m_couplingTypeOptionString("FixedStress"),
  m_couplingTypeOption(),
  m_solidSolver(nullptr),
  m_flowSolver(nullptr),
  m_maxNumResolves(10)
{
  registerWrapper(viewKeyStruct::solidSolverNameString, &m_solidSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the solid mechanics solver to use in the poroelastic solver");

  registerWrapper(viewKeyStruct::fluidSolverNameString, &m_flowSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the fluid mechanics solver to use in the poroelastic solver");

  registerWrapper(viewKeyStruct::couplingTypeOptionStringString, &m_couplingTypeOptionString, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Coupling option: (FixedStress, TightlyCoupled)");

  registerWrapper(viewKeyStruct::contactRelationNameString, &m_contactRelationName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of contact relation to enforce constraints on fracture boundary.");

  registerWrapper(viewKeyStruct::maxNumResolvesString, &m_maxNumResolves, 0)->
    setApplyDefaultValue(10)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value to indicate how many resolves may be executed to perform surface generation after the execution of flow and mechanics solver. ");

  m_numResolves[0] = 0;
}

void HydrofractureSolver::RegisterDataOnMesh( dataRepository::Group * const GEOSX_UNUSED_ARG( MeshBodies ) )
{

}

void HydrofractureSolver::ImplicitStepSetup( real64 const & time_n,
                                             real64 const & dt,
                                             DomainPartition * const domain,
                                             DofManager & GEOSX_UNUSED_ARG( dofManager ),
                                             ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                             ParallelVector & GEOSX_UNUSED_ARG( rhs ),
                                             ParallelVector & GEOSX_UNUSED_ARG( solution ) )
{
  m_solidSolver = this->getParent()->GetGroup<SolidMechanicsLagrangianFEM>(m_solidSolverName);
  m_flowSolver = this->getParent()->GetGroup<FlowSolverBase>(m_flowSolverName);

  m_solidSolver->ImplicitStepSetup( time_n, dt, domain,
                                    m_solidSolver->getDofManager(),
                                    m_solidSolver->getSystemMatrix(),
                                    m_solidSolver->getSystemRhs(),
                                    m_solidSolver->getSystemSolution() );

  m_flowSolver->ImplicitStepSetup( time_n, dt, domain,
                                   m_flowSolver->getDofManager(),
                                   m_flowSolver->getSystemMatrix(),
                                   m_flowSolver->getSystemRhs(),
                                   m_flowSolver->getSystemSolution() );
}

void HydrofractureSolver::ImplicitStepComplete( real64 const& time_n,
                                                real64 const& dt,
                                                DomainPartition * const domain)
{
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
}

void HydrofractureSolver::PostProcessInput()
{
  string ctOption = this->getReference<string>(viewKeyStruct::couplingTypeOptionStringString);

  if( ctOption == "FixedStress" )
  {
    this->m_couplingTypeOption = couplingTypeOption::FixedStress;
  }
  else if( ctOption == "TightlyCoupled" )
  {
    this->m_couplingTypeOption = couplingTypeOption::TightlyCoupled;
  }
  else
  {
    GEOSX_ERROR("invalid coupling type option");
  }
}

void HydrofractureSolver::InitializePostInitialConditions_PreSubGroups(Group * const GEOSX_UNUSED_ARG( problemManager ) )
{

}

HydrofractureSolver::~HydrofractureSolver()
{
  // TODO Auto-generated destructor stub
#ifdef GEOSX_USE_HYPRE_MGR
  if (IJ_matrix != nullptr)
  {
    HYPRE_IJMatrixDestroy(IJ_matrix);
    IJ_matrix = nullptr;
  }
  if (IJ_matrix_uu != nullptr)
  {
    HYPRE_IJMatrixDestroy(IJ_matrix_uu);
    IJ_matrix_uu = nullptr;
  }
  if (IJ_rhs != nullptr)
  {
    HYPRE_IJVectorDestroy(IJ_rhs);
    IJ_rhs = nullptr;
  }
  if (IJ_lhs != nullptr)
  {
    HYPRE_IJVectorDestroy(IJ_lhs);
    IJ_lhs = nullptr;
  }
#endif
}

void HydrofractureSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  m_flowSolver->ResetStateToBeginningOfStep(domain);
  m_solidSolver->ResetStateToBeginningOfStep(domain);
}

real64 HydrofractureSolver::SolverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition * const domain )
{
  real64 dtReturn = dt;

  SolverBase * const surfaceGenerator =  this->getParent()->GetGroup<SolverBase>("SurfaceGen");

  if( m_couplingTypeOption == couplingTypeOption::FixedStress )
  {
    dtReturn = SplitOperatorStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>() );
  }
  else if( m_couplingTypeOption == couplingTypeOption::TightlyCoupled )
  {

    ImplicitStepSetup( time_n,
                       dt,
                       domain,
                       m_dofManager,
                       m_matrix,
                       m_rhs,
                       m_solution );

    int const maxIter = m_maxNumResolves + 1;
    m_numResolves[1] = m_numResolves[0];
    int solveIter;
    for( solveIter=0 ; solveIter<maxIter ; ++solveIter )
    {
      int locallyFractured = 0;
      int globallyFractured = 0;

      SetupSystem( domain,
                   m_dofManager,
                   m_matrix,
                   m_rhs,
                   m_solution  );

      if( solveIter>0 )
      {
        m_solidSolver->ResetStressToBeginningOfStep( domain );
      }

      // currently the only method is implicit time integration
      dtReturn = this->NonlinearImplicitStep( time_n,
                                              dt,
                                              cycleNumber,
                                              domain,
                                              m_dofManager,
                                              m_matrix,
                                              m_rhs,
                                              m_solution );

      m_solidSolver->updateStress( domain );

      if( surfaceGenerator!=nullptr )
      {
        if( surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain ) > 0 )
        {
          locallyFractured = 1;
        }
        MpiWrapper::allReduce( &locallyFractured,
                               &globallyFractured,
                               1,
                               MPI_MAX,
                               MPI_COMM_GEOSX );
      }
      if( globallyFractured == 0 )
      {
        // surfaceGenerator->group_cast<SurfaceGenerator*>()->getSurfaceElementsRupturedThisSolve().clear();
        // set < localIndex > surfaceElemsRupturedThisSolve = surfaceGenerator->group_cast<surfaceGenerator*>()->getSurfaceElementsRupturedThisSolve();
        // surfaceElemsRupturedThisSolve.clear();
        break;
      }
      else
      {
        std::map<string, string_array > fieldNames;
        fieldNames["node"].push_back( keys::IncrementalDisplacement );
        fieldNames["node"].push_back( keys::TotalDisplacement );
        fieldNames["elems"].push_back( FlowSolverBase::viewKeyStruct::pressureString );
        fieldNames["elems"].push_back( "elementAperture" );

        CommunicationTools::SynchronizeFields( fieldNames,
                                               domain->getMeshBody(0)->getMeshLevel(0),
                                               domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );


        if( getLogLevel() >= 1 )
        {
          GEOSX_LOG_RANK_0("++ Fracture propagation. Re-entering Newton Solve.");
        }
        m_flowSolver->ResetViews(domain);
        this->UpdateDeformationForCoupling(domain);
      }
    }

    // final step for completion of timestep. typically secondary variable updates and cleanup.
    ImplicitStepComplete( time_n, dtReturn, domain );
    n_cycles++;
    m_numResolves[1] = solveIter;
  }

  return dtReturn;
}

void HydrofractureSolver::UpdateDeformationForCoupling( DomainPartition * const domain )
{
  MeshLevel * const meshLevel = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();

  arrayView1d<R1Tensor> const & u = nodeManager->getReference< array1d<R1Tensor> >( keys::TotalDisplacement );
  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();
  // arrayView1d<real64 const> const & faceArea = faceManager->faceArea();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  ConstitutiveManager const * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  ContactRelationBase const * const
  contactRelation = constitutiveManager->GetGroup<ContactRelationBase>(m_contactRelationName);

  elemManager->forElementRegions<FaceElementRegion>([&]( FaceElementRegion * const faceElemRegion )
  {
    faceElemRegion->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion * const subRegion )
    {
      arrayView1d<real64> const & aperture = subRegion->getElementAperture();
      arrayView1d<real64> const & volume = subRegion->getElementVolume();
      arrayView1d<real64> const & deltaVolume = subRegion->getReference<array1d<real64> >(FlowSolverBase::viewKeyStruct::deltaVolumeString);
      arrayView1d<real64 const> const & area = subRegion->getElementArea();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();

      for( localIndex kfe=0 ; kfe<subRegion->size() ; ++kfe )
      {
        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const kf1 = elemsToFaces[kfe][1];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray(kf0);
        R1Tensor temp;
        for( localIndex a=0 ; a<numNodesPerFace ; ++a )
        {
          temp += u[faceToNodeMap(kf0, a)];
          temp -= u[faceToNodeMap(kf1, a)];
        }

        // TODO this needs a proper contact based strategy for aperture
        aperture[kfe] = -Dot(temp,faceNormal[kf0]) / numNodesPerFace;
        aperture[kfe] = contactRelation->effectiveAperture( aperture[kfe] );

        deltaVolume[kfe] = aperture[kfe] * area[kfe] - volume[kfe];
      }

    });
  });
}

real64 HydrofractureSolver::SplitOperatorStep( real64 const & GEOSX_UNUSED_ARG( time_n ),
                                               real64 const & dt,
                                               integer const GEOSX_UNUSED_ARG( cycleNumber ),
                                               DomainPartition * const GEOSX_UNUSED_ARG( domain ) )
{
  real64 dtReturn = dt;
//  real64 dtReturnTemporary = dtReturn;
//
//  m_flowSolver->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
//  m_solidSolver->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
//  this->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
//
//
//
//  fluidSolver.ImplicitStepSetup( time_n, dt, domain,
//                                 fluidSolver.getDofManager(),
//                                 fluidSolver.getSystemMatrix(),
//                                 fluidSolver.getSystemRhs(),
//                                 fluidSolver.getSystemSolution() );
//
//  solidSolver.ImplicitStepSetup( time_n, dt, domain,
//                                 solidSolver.getDofManager(),
//                                 solidSolver.getSystemMatrix(),
//                                 solidSolver.getSystemRhs(),
//                                 solidSolver.getSystemSolution() );
//
//  this->UpdateDeformationForCoupling(domain);
//
//  int iter = 0;
//  while (iter < solverParams->maxIterNewton() )
//  {
//    if (iter == 0)
//    {
//      // reset the states of all slave solvers if any of them has been reset
//      m_flowSolver->ResetStateToBeginningOfStep( domain );
//      m_solidSolver->ResetStateToBeginningOfStep( domain );
//      ResetStateToBeginningOfStep( domain );
//    }
//    LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );
//
//    // call assemble to fill the matrix and the rhs
//    m_flowSolver->AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );
//
//    // apply boundary conditions to system
//    m_flowSolver->ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
//
//    // call the default linear solver on the system
//    m_flowSolver->SolveSystem( getLinearSystemRepository(),
//                 getSystemSolverParameters() );
//
//    // apply the system solution to the fields/variables
//    m_flowSolver->ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );
//
//    if (dtReturnTemporary < dtReturn)
//    {
//      iter = 0;
//      dtReturn = dtReturnTemporary;
//      continue;
//    }
//
////    if (m_fluidSolver->getSystemSolverParameters()->numNewtonIterations() == 0 && iter > 0 && getLogLevel() >= 1)
////    {
////      GEOSX_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
////      break;
////    }
//
//    if (getLogLevel() >= 1)
//    {
//      GEOSX_LOG_RANK_0( "\tIteration: " << iter+1  << ", MechanicsSolver: " );
//    }
//
//    // call assemble to fill the matrix and the rhs
//    m_solidSolver->AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );
//
//
//    ApplyFractureFluidCoupling( domain, *getLinearSystemRepository() );
//
//    // apply boundary conditions to system
//    m_solidSolver->ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
//
//    // call the default linear solver on the system
//    m_solidSolver->SolveSystem( getLinearSystemRepository(),
//                 getSystemSolverParameters() );
//
//    // apply the system solution to the fields/variables
//    m_solidSolver->ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );
//
//    if( m_flowSolver->CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() &&
//        m_solidSolver->CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() )
//    {
//      GEOSX_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
//      break;
//    }
//
//    if (dtReturnTemporary < dtReturn)
//    {
//      iter = 0;
//      dtReturn = dtReturnTemporary;
//      continue;
//    }
////    if (m_solidSolver->getSystemSolverParameters()->numNewtonIterations() > 0)
//    {
//      this->UpdateDeformationForCoupling(domain);
////      m_fluidSolver->UpdateState(domain);
//    }
//    ++iter;
//  }
//
//  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

real64 HydrofractureSolver::ExplicitStep( real64 const& time_n,
                                          real64 const& dt,
                                          const int cycleNumber,
                                          DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ExplicitStep( time_n, dt, cycleNumber, domain );
  m_flowSolver->SolverStep( time_n, dt, cycleNumber, domain );

  return dt;
}


void HydrofractureSolver::SetupDofs( DomainPartition const * const domain,
                                     DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );
  m_flowSolver->SetupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString,
                          DofManager::Connectivity::Elem );
}

void HydrofractureSolver::SetupSystem( DomainPartition * const domain,
                                       DofManager & GEOSX_UNUSED_ARG( dofManager ),
                                       ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                       ParallelVector & GEOSX_UNUSED_ARG( rhs ),
                                       ParallelVector & GEOSX_UNUSED_ARG( solution ) )
{
  GEOSX_MARK_FUNCTION;
  m_flowSolver->ResetViews( domain );

  m_solidSolver->SetupSystem( domain,
                           m_solidSolver->getDofManager(),
                           m_solidSolver->getSystemMatrix(),
                           m_solidSolver->getSystemRhs(),
                           m_solidSolver->getSystemSolution() );

  m_flowSolver->SetupSystem( domain,
                           m_flowSolver->getDofManager(),
                           m_flowSolver->getSystemMatrix(),
                           m_flowSolver->getSystemRhs(),
                           m_flowSolver->getSystemSolution() );

  // TODO: once we move to a monolithic matrix, we can just use SolverBase implementation

//  dofManager.setSparsityPattern( m_matrix01,
//                                 keys::TotalDisplacement,
//                                 FlowSolverBase::viewKeyStruct::pressureString );
//
//  dofManager.setSparsityPattern( m_matrix10,
//                                 FlowSolverBase::viewKeyStruct::pressureString,
//                                 keys::TotalDisplacement );

  m_matrix01.createWithLocalSize( m_solidSolver->getSystemMatrix().localRows(),
                                  m_flowSolver->getSystemMatrix().localCols(),
                                  9,
                                  MPI_COMM_GEOSX);
  m_matrix10.createWithLocalSize( m_flowSolver->getSystemMatrix().localCols(),
                                  m_solidSolver->getSystemMatrix().localRows(),
                                  24,
                                  MPI_COMM_GEOSX);

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  std::unique_ptr<CRSMatrix<real64,localIndex,localIndex> > &
  derivativeFluxResidual_dAperture = m_flowSolver->getRefDerivativeFluxResidual_dAperture();
  {

    localIndex numRows = 0;
    localIndex numCols = 0;
    string_array const & flowRegions = m_flowSolver->getTargetRegions();
    elemManager->forElementSubRegions( flowRegions, [&]( ElementSubRegionBase const * const elementSubRegion )
    {
      numRows += elementSubRegion->size();
      numCols += elementSubRegion->size();
    });

    derivativeFluxResidual_dAperture = std::make_unique<CRSMatrix<real64,localIndex,localIndex>>( numRows, numCols );

    derivativeFluxResidual_dAperture->reserveNonZeros( m_flowSolver->getSystemMatrix().localNonzeros() );
    localIndex maxRowSize = -1;
    for( localIndex row=0 ; row<m_flowSolver->getSystemMatrix().localRows() ; ++row )
    {
      localIndex const rowSize = m_flowSolver->getSystemMatrix().getLocalRowGlobalLength( row );
      maxRowSize = maxRowSize > rowSize ? maxRowSize : rowSize;

      derivativeFluxResidual_dAperture->reserveNonZeros( row,
                                                         rowSize );
    }
    for( localIndex row=m_flowSolver->getSystemMatrix().localRows() ; row<numRows ; ++row )
    {
      derivativeFluxResidual_dAperture->reserveNonZeros( row,
                                                         maxRowSize );
    }


  }
//  CRSMatrixView<real64,localIndex,localIndex const> const &
//  derivativeFluxResidual_dAperture = m_flowSolver->getDerivativeFluxResidual_dAperture();

  string const presDofKey = m_flowSolver->getDofManager().getKey( FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  arrayView1d<globalIndex> const &
  dispDofNumber =  nodeManager->getReference<globalIndex_array>( dispDofKey );

  elemManager->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion const * const elementSubRegion )
  {
    localIndex const numElems = elementSubRegion->size();
    array1d<array1d<localIndex > > const & elemsToNodes = elementSubRegion->nodeList();
    arrayView1d<globalIndex> const &
    faceElementDofNumber = elementSubRegion->getReference< array1d<globalIndex> >( presDofKey );

    for( localIndex k=0 ; k<numElems ; ++k )
    {
      globalIndex const activeFlowDOF = faceElementDofNumber[k];
      localIndex const numNodesPerElement = elemsToNodes[k].size();
      array1d<globalIndex> activeDisplacementDOF(3 * numNodesPerElement);
      array1d<real64> values( 3*numNodesPerElement );
      values = 1;

      for( localIndex a=0 ; a<numNodesPerElement ; ++a )
      {
        for( int d=0 ; d<3 ; ++d )
        {
          activeDisplacementDOF[a * 3 + d] = dispDofNumber[elemsToNodes[k][a]] + d;
        }
      }

      m_matrix01.insert( activeDisplacementDOF.data(),
                         &activeFlowDOF,
                         values.data(),
                         activeDisplacementDOF.size(),
                         1 );

      m_matrix10.insert( &activeFlowDOF,
                         activeDisplacementDOF.data(),
                         values.data(),
                         1,
                         activeDisplacementDOF.size() );
    }
  });

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_flowSolver->getDiscretization() );


  fluxApprox->forStencils<FaceElementStencil>( [&]( FaceElementStencil const & stencil )
  {
//    forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
    for( localIndex iconn=0 ; iconn<stencil.size() ; ++iconn)
    {
      localIndex const numFluxElems = stencil.stencilSize(iconn);
      typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      FaceElementSubRegion const * const
      elementSubRegion = elemManager->GetRegion(seri[iconn][0])->GetSubRegion<FaceElementSubRegion>(sesri[iconn][0]);

//      GEOS_LOG_RANK("connector, numLocal, numGhost: "<<iconn<<", "<<elementSubRegion->size()-elementSubRegion->GetNumberOfGhosts()<<", "<<elementSubRegion->GetNumberOfGhosts());
//      GEOS_LOG_RANK("connector, numRows, numCols: "<<iconn<<", "<<derivativeFluxResidual_dAperture->numRows()<<", "<<derivativeFluxResidual_dAperture->numColumns());

      array1d<array1d<localIndex > > const & elemsToNodes = elementSubRegion->nodeList();

//      arrayView1d<integer const> const & ghostRank = elementSubRegion->GhostRank();

      arrayView1d<globalIndex> const &
      faceElementDofNumber = elementSubRegion->getReference< array1d<globalIndex> >( presDofKey );
      for( localIndex k0=0 ; k0<numFluxElems ; ++k0 )
      {
        globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];

        for( localIndex k1=0 ; k1<numFluxElems ; ++k1 )
        {
//          GEOS_LOG_RANK("ei0, ei1, nonZeroCapacitys: "<<sei[iconn][k0]<<", "<<sei[iconn][k1]<<", "<<derivativeFluxResidual_dAperture->nonZeroCapacity(sei[iconn][k0]));
          derivativeFluxResidual_dAperture->insertNonZero( sei[iconn][k0],sei[iconn][k1], 0.0 );

          localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();
          array1d<globalIndex> activeDisplacementDOF(3 * numNodesPerElement);
          array1d<real64> values( 3*numNodesPerElement );
          values = 1;

          for( localIndex a=0 ; a<numNodesPerElement ; ++a )
          {
            for( int d=0 ; d<3 ; ++d )
            {
              activeDisplacementDOF[a * 3 + d] = dispDofNumber[elemsToNodes[sei[iconn][k1]][a]] + d;
            }
          }

          m_matrix10.insert( &activeFlowDOF,
                             activeDisplacementDOF.data(),
                             values.data(),
                             1,
                             activeDisplacementDOF.size() );
        }
      }
    }//);
  });

  m_matrix01.close();
  m_matrix10.close();
}

void HydrofractureSolver::AssembleSystem( real64 const time,
                                          real64 const dt,
                                          DomainPartition * const domain,
                                          DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                                          ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                          ParallelVector & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->AssembleSystem( time,
                                 dt,
                                 domain,
                                 m_solidSolver->getDofManager(),
                                 m_solidSolver->getSystemMatrix(),
                                 m_solidSolver->getSystemRhs() );

  m_flowSolver->AssembleSystem( time,
                                dt,
                                domain,
                                m_flowSolver->getDofManager(),
                                m_flowSolver->getSystemMatrix(),
                                m_flowSolver->getSystemRhs() );



  AssembleForceResidualDerivativeWrtPressure( domain, &m_matrix01, &(m_solidSolver->getSystemRhs()) );
  AssembleFluidMassResidualDerivativeWrtDisplacement( domain, &m_matrix10, &(m_flowSolver->getSystemRhs()) );
}

void HydrofractureSolver::ApplyBoundaryConditions( real64 const time,
                                                   real64 const dt,
                                                   DomainPartition * const domain,
                                                   DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                                                   ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                                   ParallelVector & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          m_solidSolver->getDofManager(),
                                          m_solidSolver->getSystemMatrix(),
                                          m_solidSolver->getSystemRhs() );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );
  NodeManager const * const nodeManager = mesh->getNodeManager();
  arrayView1d<globalIndex const> const & dispDofNumber = nodeManager->getReference<globalIndex_array>( dispDofKey );
  arrayView1d<integer const> const & nodeGhostRank = nodeManager->GhostRank();

  fsManager.Apply( time + dt,
                   domain,
                   "nodeManager",
                   keys::TotalDisplacement,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        set<localIndex> const & targetSet,
                        Group * const ,
                        string const )
  {
    set<localIndex> localSet;
    for( auto const & a : targetSet )
    {
      if( nodeGhostRank[a]<0 )
      {
        localSet.insert(a);
      }
    }
    bc->ZeroSystemRowsForBoundaryCondition<LAInterface>( localSet,
                                                         dispDofNumber,
                                                         m_matrix01 );
  } );


  m_flowSolver->ApplyBoundaryConditions( time,
                                         dt,
                                         domain,
                                         m_flowSolver->getDofManager(),
                                         m_flowSolver->getSystemMatrix(),
                                         m_flowSolver->getSystemRhs() );

  string const presDofKey = m_flowSolver->getDofManager().getKey( FlowSolverBase::viewKeyStruct::pressureString );

  fsManager.Apply( time + dt,
                    domain,
                    "ElementRegions",
                    FlowSolverBase::viewKeyStruct::pressureString,
                    [&]( FieldSpecificationBase const * const fs,
                         string const &,
                         set<localIndex> const & lset,
                         Group * subRegion,
                         string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( presDofKey );
    arrayView1d<integer const> const & ghostRank = subRegion->group_cast<ObjectManagerBase*>()->GhostRank();

    set<localIndex> localSet;
    for( auto const & a : lset )
    {
      if( ghostRank[a]<0 )
      {
        localSet.insert(a);
      }
    }

    fs->ZeroSystemRowsForBoundaryCondition<LAInterface>( localSet,
                                                         dofNumber,
                                                         m_matrix10 );
  });

  // debugging info.  can be trimmed once everything is working.
  if( getLogLevel()>=10 )
  {
    // Before outputting anything generate permuation matrix and permute.
    ElementRegionManager * const elemManager = mesh->getElemManager();

    LAIHelperFunctions::CreatePermutationMatrix(nodeManager,
                                                m_solidSolver->getSystemMatrix().globalRows(),
                                                m_solidSolver->getSystemMatrix().globalCols(),
                                                3,
                                                m_solidSolver->getDofManager().getKey( keys::TotalDisplacement ),
                                                m_permutationMatrix0);

    LAIHelperFunctions::CreatePermutationMatrix(elemManager,
                                                m_flowSolver->getSystemMatrix().globalRows(),
                                                m_flowSolver->getSystemMatrix().globalCols(),
                                                1,
                                                m_flowSolver->getDofManager().getKey( FlowSolverBase::viewKeyStruct::pressureString ),
                                                m_permutationMatrix1);

    GEOSX_LOG_RANK_0("***********************************************************");
    GEOSX_LOG_RANK_0("matrix00");
    GEOSX_LOG_RANK_0("***********************************************************");
//    LAIHelperFunctions::PrintPermutedMatrix(m_solidSolver->getSystemMatrix(), m_permutationMatrix0, std::cout);
    m_solidSolver->getSystemMatrix().print(std::cout);
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0("***********************************************************");
    GEOSX_LOG_RANK_0("matrix01");
    GEOSX_LOG_RANK_0("***********************************************************");
//    LAIHelperFunctions::PrintPermutedMatrix(m_matrix01, m_permutationMatrix0, m_permutationMatrix1, std::cout);
    m_matrix01.print(std::cout);
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0("***********************************************************");
    GEOSX_LOG_RANK_0("matrix10");
    GEOSX_LOG_RANK_0("***********************************************************");
//    LAIHelperFunctions::PrintPermutedMatrix(m_matrix10, m_permutationMatrix1, m_permutationMatrix0, std::cout);
    m_matrix10.print(std::cout);
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0("***********************************************************");
    GEOSX_LOG_RANK_0("matrix11");
    GEOSX_LOG_RANK_0("***********************************************************");
//    LAIHelperFunctions::PrintPermutedMatrix(m_flowSolver->getSystemMatrix(), m_permutationMatrix1, std::cout);
    m_flowSolver->getSystemMatrix().print(std::cout);
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0("***********************************************************");
    GEOSX_LOG_RANK_0("residual0");
    GEOSX_LOG_RANK_0("***********************************************************");
//    LAIHelperFunctions::PrintPermutedVector(m_solidSolver->getSystemRhs(), m_permutationMatrix0, std::cout);
    m_solidSolver->getSystemRhs().print(std::cout);
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0("***********************************************************");
    GEOSX_LOG_RANK_0("residual1");
    GEOSX_LOG_RANK_0("***********************************************************");
//    LAIHelperFunctions::PrintPermutedVector(m_flowSolver->getSystemRhs(), m_permutationMatrix1, std::cout);
    m_flowSolver->getSystemRhs().print(std::cout);
    MpiWrapper::Barrier();
  }

  if( getLogLevel() >= 10 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    {
      string filename = "matrix00_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemMatrix().write( filename, true );
      GEOSX_LOG_RANK_0( "matrix00: written to " << filename );
    }
    {
      string filename = "matrix01_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_matrix01.write( filename, true );
      GEOSX_LOG_RANK_0( "matrix01: written to " << filename );
    }
    {
      string filename = "matrix10_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_matrix10.write( filename, true );
      GEOSX_LOG_RANK_0( "matrix10: written to " << filename );
    }
    {
      string filename = "matrix11_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_flowSolver->getSystemMatrix().write( filename, true );
      GEOSX_LOG_RANK_0( "matrix11: written to " << filename );
    }
    {
      string filename = "residual0_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemRhs().write( filename, true );
      GEOSX_LOG_RANK_0( "residual0: written to " << filename );
    }
    {
      string filename = "residual1_" + std::to_string( time ) + "_" + std::to_string( newtonIter ) + ".mtx";
      m_flowSolver->getSystemRhs().write( filename, true );
      GEOSX_LOG_RANK_0( "residual1: written to " << filename );
    }
  }
}

real64
HydrofractureSolver::
CalculateResidualNorm( DomainPartition const * const domain,
                       DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                       ParallelVector const & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  /*
  real64 const fluidResidual = m_flowSolver->getSystemRhs().norm2();
  real64 const solidResidual = m_solidSolver->getSystemRhs().norm2();
  */

  real64 const fluidResidual = m_flowSolver->CalculateResidualNorm( domain,
                                                                    m_flowSolver->getDofManager(),
                                                                    m_flowSolver->getSystemRhs() );

  real64 const solidResidual = m_solidSolver->CalculateResidualNorm( domain,
                                                                     m_solidSolver->getDofManager(),
                                                                     m_solidSolver->getSystemRhs() );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output,
             "( Rfluid, Rsolid ) = (%4.2e, %4.2e) ; ",
             fluidResidual,
             solidResidual);
    std::cout<<output;
  }

  return fluidResidual + solidResidual;
}



void
HydrofractureSolver::
AssembleForceResidualDerivativeWrtPressure( DomainPartition * const domain,
                                            ParallelMatrix * const matrix01,
                                            ParallelVector * const rhs0 )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

  arrayView1d<R1Tensor> const &
  fext = nodeManager->getReference< array1d<R1Tensor> >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternal );
  fext = {0,0,0};

  string const presDofKey = m_flowSolver->getDofManager().getKey( FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  arrayView1d<globalIndex> const &
  dispDofNumber =  nodeManager->getReference<globalIndex_array>( dispDofKey );


  matrix01->open();
  matrix01->zero();
  rhs0->open();

  elemManager->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion * const subRegion )->void
  {

    arrayView1d<globalIndex> const &
    faceElementDofNumber = subRegion->getReference< array1d<globalIndex> >( presDofKey );

    if( subRegion->hasWrapper( "pressure" ) )
    {
      arrayView1d<real64 const> const & fluidPressure = subRegion->getReference<array1d<real64> >("pressure");
      arrayView1d<real64 const> const & deltaFluidPressure = subRegion->getReference<array1d<real64> >("deltaPressure");
      arrayView1d<integer const> const & ghostRank = subRegion->GhostRank();
      arrayView1d<real64> const & area = subRegion->getElementArea();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();

      forall_in_range<serialPolicy>( 0,
                                   subRegion->size(),
                                   GEOSX_LAMBDA ( localIndex const kfe )
      {
        R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
        Nbar -= faceNormal[elemsToFaces[kfe][1]];
        Nbar.Normalize();

        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray(kf0);

        globalIndex rowDOF[24];
        real64 nodeRHS[24];
        stackArray2d<real64, 12*12> dRdP(numNodesPerFace*3, 1);
        globalIndex colDOF = faceElementDofNumber[kfe];

        real64 const Ja = area[kfe] / numNodesPerFace;

        //          std::cout<<"fluidPressure["<<kfe<<"] = "<<fluidPressure[kfe]+deltaFluidPressure[kfe]<<std::endl;
        real64 nodalForceMag = ( fluidPressure[kfe]+deltaFluidPressure[kfe] ) * Ja;
        R1Tensor nodalForce(Nbar);
        nodalForce *= nodalForceMag;

        //          std::cout << "    rank " << MpiWrapper::Comm_rank(MPI_COMM_GEOSX) << ", faceElement " << kfe << std::endl;
        //          std::cout << "    fluid pressure " << fluidPressure[kfe]+deltaFluidPressure[kfe] << std::endl;
        //          std::cout << "    nodalForce " << nodalForce << std::endl;
        for( localIndex kf=0 ; kf<2 ; ++kf )
        {
          localIndex const faceIndex = elemsToFaces[kfe][kf];


          for( localIndex a=0 ; a<numNodesPerFace ; ++a )
          {

            for( int i=0 ; i<3 ; ++i )
            {
              rowDOF[3*a+i] = dispDofNumber[faceToNodeMap(faceIndex, a)] + i;
              nodeRHS[3*a+i] = - nodalForce[i] * pow(-1,kf);
              fext[faceToNodeMap(faceIndex, a)][i] += - nodalForce[i] * pow(-1,kf);

              dRdP(3*a+i,0) = - Ja * Nbar[i] * pow(-1,kf);
              // this is for debugging
              //                if (dispDofNumber[faceToNodeMap(faceIndex, a)] == 0 || dispDofNumber[faceToNodeMap(faceIndex, a)] == 6 || dispDofNumber[faceToNodeMap(faceIndex, a)] == 12 || dispDofNumber[faceToNodeMap(faceIndex, a)] == 18)
              //                  std::cout << "rank " << MpiWrapper::Comm_rank(MPI_COMM_GEOSX) << "DOF index " << dispDofNumber[faceToNodeMap(faceIndex, a)] + i << " contribution " << nodeRHS[3*a+i] << std::endl;

            }
          }
          if( ghostRank[kfe] < 0 )
          {

            rhs0->add( rowDOF,
                       nodeRHS,
                       numNodesPerFace*3 );


            matrix01->add( rowDOF,
                           &colDOF,
                           dRdP.data(),
                           numNodesPerFace * 3,
                           1 );
          }
        }
      });
    }
  });

  rhs0->close();
  matrix01->close();
  rhs0->close();
}

void
HydrofractureSolver::
AssembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const * const domain,
                                                    ParallelMatrix * const matrix10,
                                                    ParallelVector * const GEOSX_UNUSED_ARG( rhs0 ) )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = mesh->getElemManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager const * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager const * const constitutiveManager = domain->getConstitutiveManager();

  string const constitutiveName = constitutiveManager->GetGroup(m_flowSolver->fluidIndex())->getName();
  string const presDofKey = m_flowSolver->getDofManager().getKey( FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey = m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  CRSMatrixView<real64 const,localIndex const,localIndex const> const &
  dFluxResidual_dAperture = m_flowSolver->getDerivativeFluxResidual_dAperture();

  ContactRelationBase const * const
  contactRelation = constitutiveManager->GetGroup<ContactRelationBase>( m_contactRelationName );

  matrix10->open();
  matrix10->zero();

  elemManager->forElementSubRegionsComplete<FaceElementSubRegion>( this->m_targetRegions,
                                                                   [&] ( localIndex GEOSX_UNUSED_ARG( er ),
                                                                         localIndex GEOSX_UNUSED_ARG( esr ),
                                                                         ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                                                         FaceElementSubRegion const * const subRegion )
  {


    dataRepository::Group const * const constitutiveGroup = subRegion->GetConstitutiveModels();
    dataRepository::Group const * const constitutiveRelation = constitutiveGroup->GetGroup(constitutiveName);

    arrayView1d<integer const>     const & elemGhostRank = subRegion->GhostRank();
    arrayView1d<globalIndex const> const & presDofNumber = subRegion->getReference<array1d<globalIndex>>( presDofKey );
    arrayView1d<globalIndex const> const & dispDofNumber = nodeManager->getReference<array1d<globalIndex>>( dispDofKey );

    arrayView2d<real64 const> const &
    dens = constitutiveRelation->getReference<array2d<real64>>(SingleFluidBase::viewKeyStruct::densityString);

    arrayView1d<real64 const> const & aperture  = subRegion->getElementAperture();
    arrayView1d<real64 const> const & area      = subRegion->getElementArea();

    arrayView2d<localIndex const> const & elemsToFaces = subRegion->faceList();
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager->nodeList();

    arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();


    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      //if (elemGhostRank[ei] < 0)
      {
        globalIndex const elemDOF = presDofNumber[ei];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray(elemsToFaces[ei][0]);
        real64 const dAccumulationResidualdAperture = dens[ei][0] * area[ei];


        globalIndex nodeDOF[8*3];

        R1Tensor Nbar = faceNormal[elemsToFaces[ei][0]];
        Nbar -= faceNormal[elemsToFaces[ei][1]];
        Nbar.Normalize();

        stackArray1d<real64, 24> dRdU(2*numNodesPerFace*3);

        // Accumulation derivative
        if (elemGhostRank[ei] < 0)
        {
          //GEOS_LOG_RANK( "dAccumulationResidualdAperture("<<ei<<") = "<<dAccumulationResidualdAperture );
          for( localIndex kf=0 ; kf<2 ; ++kf )
          {
            for( localIndex a=0 ; a<numNodesPerFace ; ++a )
            {
              for( int i=0 ; i<3 ; ++i )
              {
                nodeDOF[ kf*3*numNodesPerFace + 3*a+i] = dispDofNumber[faceToNodeMap(elemsToFaces[ei][kf],a)] +i;
                real64 const dGap_dU = - pow(-1,kf) * Nbar[i] / numNodesPerFace;
                real64 const dAper_dU = contactRelation->dEffectiveAperture_dAperture( aperture[ei] ) * dGap_dU;
                dRdU(kf*3*numNodesPerFace + 3*a+i) = dAccumulationResidualdAperture * dAper_dU;
              }
            }
          }
          matrix10->add( elemDOF,
                         nodeDOF,
                         dRdU.data(),
                         2*numNodesPerFace*3 );
        }

        // flux derivative
        localIndex const numColumns = dFluxResidual_dAperture.numNonZeros(ei);
        arraySlice1d<localIndex const> const & columns = dFluxResidual_dAperture.getColumns( ei );
        arraySlice1d<real64 const> const & values = dFluxResidual_dAperture.getEntries( ei );

        for( localIndex kfe2=0 ; kfe2<numColumns ; ++kfe2 )
        {
          real64 dRdAper = values[kfe2];
          localIndex const ei2 = columns[kfe2];
//          GEOS_LOG_RANK( "dRdAper("<<ei<<", "<<ei2<<") = "<<dRdAper );

          for( localIndex kf=0 ; kf<2 ; ++kf )
          {
            for( localIndex a=0 ; a<numNodesPerFace ; ++a )
            {
              for( int i=0 ; i<3 ; ++i )
              {
                nodeDOF[ kf*3*numNodesPerFace + 3*a+i] = dispDofNumber[faceToNodeMap(elemsToFaces[ei2][kf],a)] +i;
                real64 const dGap_dU = - pow(-1,kf) * Nbar[i] / numNodesPerFace;
                real64 const dAper_dU = contactRelation->dEffectiveAperture_dAperture( aperture[ei2] ) * dGap_dU;
                dRdU(kf*3*numNodesPerFace + 3*a+i) = dRdAper * dAper_dU;
              }
            }
          }
          matrix10->add( elemDOF,
                         nodeDOF,
                         dRdU.data(),
                         2*numNodesPerFace*3 );

        }
      }
    });
  });

  matrix10->close();
}

void
HydrofractureSolver::
ApplySystemSolution( DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                     ParallelVector const & GEOSX_UNUSED_ARG( solution ),
                     real64 const scalingFactor,
                     DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplySystemSolution( m_solidSolver->getDofManager(),
                                      m_solidSolver->getSystemSolution(),
                                      scalingFactor,
                                      domain );
  m_flowSolver->ApplySystemSolution( m_flowSolver->getDofManager(),
                                     m_flowSolver->getSystemSolution(),
                                     -scalingFactor,
                                     domain );

  this->UpdateDeformationForCoupling(domain);
}

}

#include "EpetraExt_MatrixMatrix.h"
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
#include "Teuchos_Time.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

namespace geosx
{
  
void HydrofractureSolver::SolveSystem( DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                                       ParallelMatrix & ,
                                       ParallelVector & ,
                                       ParallelVector &  )
{
  GEOSX_MARK_FUNCTION;

  /*
  globalIndex numU = m_solidSolver->getSystemRhs().globalSize();
  globalIndex numP = m_flowSolver->getSystemRhs().globalSize();
  GEOSX_LOG_RANK_0("size = " << numU << " + " << numP);
  */

  SystemSolverParameters * const params = &m_systemSolverParameters;
  double setupTime, solveTime, auxTime;
  setupTime = 0.0;
  solveTime = 0.0;
  auxTime = 0.0;

  integer const newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

  using namespace Teuchos;
  using namespace Thyra;

  Teuchos::Time clock("solveClock");  

  GEOSX_MARK_BEGIN(Setup);
  Epetra_FECrsMatrix * p_matrix[2][2];
  Epetra_FEVector * p_rhs[2];
  Epetra_FEVector * p_solution[2];

  p_rhs[0] = m_solidSolver->getSystemRhs().unwrappedPointer();
  p_rhs[1] = m_flowSolver->getSystemRhs().unwrappedPointer();

  p_solution[0] = m_solidSolver->getSystemSolution().unwrappedPointer();
  p_solution[1] = m_flowSolver->getSystemSolution().unwrappedPointer();

  // set initial guess to zero
  p_solution[0]->PutScalar(0.0);
  p_solution[1]->PutScalar(0.0);

  p_matrix[0][0] = m_solidSolver->getSystemMatrix().unwrappedPointer();
  p_matrix[0][1] = m_matrix01.unwrappedPointer();
  p_matrix[1][0] = m_matrix10.unwrappedPointer();
  p_matrix[1][1] = m_flowSolver->getSystemMatrix().unwrappedPointer();

  // scale and symmetrize

  m_densityScaling = 1e-3;
  m_pressureScaling = 1e9;

  p_matrix[0][1]->Scale(m_pressureScaling);
  p_matrix[1][0]->Scale(m_pressureScaling*m_densityScaling);
  p_matrix[1][1]->Scale(m_pressureScaling*m_pressureScaling*m_densityScaling);
  p_rhs[1]->Scale(m_pressureScaling*m_densityScaling);

#ifdef GEOSX_USE_HYPRE_MGR
  // ordering = 0: reduce U first; 1: reduce P first.
  integer ordering = 0;
  timeval tim;
  gettimeofday(&tim, nullptr);
  const real64 t_start = tim.tv_sec + (tim.tv_usec / 1000000.0);

  int my_id = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
  globalIndex nrows_u, nrows_p;
  nrows_u = p_matrix[0][0]->NumGlobalRows64();
  nrows_p = p_matrix[1][1]->NumGlobalRows64();

  unsigned int num_processes = MpiWrapper::Comm_size(MPI_COMM_GEOSX);
  const std::vector<int> n_local_rows = {p_matrix[0][0]->NumMyRows(),
                                         p_matrix[1][0]->NumMyRows()};
  const int n_local_rows_all = n_local_rows[0] + n_local_rows[1];

  if (my_id == 0)
  {
    printf("Global matrix size solid = %lld, fluid = %lld\n", nrows_u, nrows_p);
  }
  //printf("my_id = %d, local matrix size solid = %d, fluid = %d\n", my_id, n_local_rows[0], n_local_rows[1]);

  int nnz_max = p_matrix[0][0]->MaxNumEntries();
  std::vector<double> col_values(nnz_max);
  std::vector<globalIndex>    col_indices(nnz_max);

  //############################################# BEGIN COPYING HYPRE MATRIX ###############################################
  clock.start(true);
  std::vector<globalIndex> offset(2*num_processes, 0);
  std::vector<globalIndex> iDOF_offset(num_processes, 0);

  if (newtonIter == 0)
  {
    GEOSX_MARK_BEGIN(COPY_HYPRE_MATRIX);
    if (IJ_matrix != nullptr)
    {
      HYPRE_IJMatrixDestroy(IJ_matrix);
      IJ_matrix = nullptr;
    }
    if (IJ_rhs != nullptr)
    {
      HYPRE_IJVectorDestroy(IJ_rhs);
      IJ_rhs = nullptr;
    }
    if (IJ_lhs != nullptr)
    {
      HYPRE_IJVectorDestroy(IJ_lhs);
      IJ_lhs = nullptr;
    }

    for(int iDOF=0; iDOF<2; ++iDOF)
    {
      MpiWrapper::Allgather<int,globalIndex>(&n_local_rows[iDOF], 1, iDOF_offset.data(), 1, MPI_COMM_GEOSX);

      for (unsigned i=0; i<num_processes; ++i)
      {
        offset[i*2+iDOF] = iDOF_offset[i];
        // if (my_id == 0) printf("Proc %d, iDOF %d, offset %lld,\n", i, iDOF, iDOF_offset[i]);
      }
    }

    int this_offset = 0;
    int next_offset = 0;
    for (unsigned iproc = 0; iproc < num_processes; ++iproc)
    {
      for (int iDOF=0; iDOF<2; ++iDOF)
      {
        next_offset += offset[iproc*2 + iDOF];
        offset[iproc*2+iDOF] = this_offset;

        //if (my_id == 0) printf("Proc %d, iDOF %d, this offset %d, next offset %d\n", iproc, iDOF, this_offset, next_offset);
        this_offset = next_offset;
      }
    }

    std::map<globalIndex, globalIndex> GID_trilinos_to_hypre_U;
    std::map<globalIndex, globalIndex> GID_trilinos_to_hypre_P;
    std::vector<int>   rank_IDs(nnz_max);
    std::vector<int>   local_IDs(nnz_max);
    int n_entries;

    // .... .... displacement mapping
    for (int iDOF=0; iDOF < 2; iDOF++)
    {
      for(int row=0; row<n_local_rows[iDOF]; ++row)
      {
        const globalIndex global_row = p_matrix[iDOF][0]->GRID64(row);
        p_matrix[iDOF][0]->ExtractGlobalRowCopy(global_row, nnz_max, n_entries, col_values.data(), col_indices.data());
        std::sort(col_indices.begin(), col_indices.begin()+n_entries);
        p_matrix[iDOF][0]->DomainMap().RemoteIDList(n_entries, col_indices.data(), rank_IDs.data(), local_IDs.data());

        for (int i = 0; i < n_entries; ++i)
        {
          GID_trilinos_to_hypre_U[col_indices[i]] = local_IDs[i] + offset[rank_IDs[i]*2];
          //if (my_id == 0 && (row < 10 || row > (n_local_rows[0] - 10)))
          //{
            //printf("row = %lld, column index = %lld, local_id = %d, offset = %lld, new GID_U = %lld\n", global_row, col_indices[i], local_IDs[i], offset[rank_IDs[i]*2], local_IDs[i] + offset[rank_IDs[i]*2]);
          //}
        }
      }
    }
    GID_trilinos_to_hypre[0] = std::move(GID_trilinos_to_hypre_U);

    // .... .... pressure mapping
    for (int iDOF = 0; iDOF < 2; iDOF++)
    {
      for(int row=0; row<n_local_rows[iDOF]; ++row)
      {
        const globalIndex global_row = p_matrix[iDOF][1]->GRID64(row);
        p_matrix[iDOF][1]->ExtractGlobalRowCopy(global_row, nnz_max, n_entries, col_values.data(), col_indices.data());
        std::sort(col_indices.begin(), col_indices.begin()+n_entries);
        p_matrix[iDOF][1]->DomainMap().RemoteIDList(n_entries, col_indices.data(), rank_IDs.data(), local_IDs.data());

        for (int i = 0; i < n_entries; ++i)
        {
          GID_trilinos_to_hypre_P[col_indices[i]] = local_IDs[i] + offset[rank_IDs[i]*2 + 1];
        }
      }
    }
    GID_trilinos_to_hypre[1] = std::move(GID_trilinos_to_hypre_P);

    // .... Create the matrix
    HYPRE_Int ilower = GID_trilinos_to_hypre[0].at(p_matrix[0][0]->GRID64(0));
    HYPRE_Int iupper = ilower + n_local_rows[0] + n_local_rows[1] - 1;

    std::vector<HYPRE_Int> nnz_local_rows(n_local_rows_all);
    int shift = n_local_rows[0];
    for (int iDOF=0; iDOF<2; ++iDOF)
    {
      for (int row=0; row<n_local_rows[iDOF]; ++row)
      {
        for (int jDOF=0; jDOF<2; ++jDOF)
        {
          int row_idx = row + (iDOF == 0 ? 0 : 1) * shift;
          nnz_local_rows[row_idx] += p_matrix[iDOF][jDOF]->NumMyEntries(row);
        }
      }
    }

    HYPRE_IJMatrixCreate(MPI_COMM_GEOSX, ilower, iupper, ilower, iupper, &IJ_matrix);
    HYPRE_IJMatrixSetObjectType(IJ_matrix, HYPRE_PARCSR);
    HYPRE_IJMatrixSetRowSizes(IJ_matrix, nnz_local_rows.data());
    //printf("Done creating matrix\n");

    // .... Create rhs and lhs
    HYPRE_IJVectorCreate(MPI_COMM_GEOSX, ilower, iupper,&IJ_rhs);
    HYPRE_IJVectorSetObjectType(IJ_rhs, HYPRE_PARCSR);

    HYPRE_IJVectorCreate(MPI_COMM_GEOSX, ilower, iupper,&IJ_lhs);
    HYPRE_IJVectorSetObjectType(IJ_lhs, HYPRE_PARCSR);
    //printf("Done creating lhs and rhs\n");
  }

  // copy entries
  // .... matrix and vector
  HYPRE_IJMatrixInitialize(IJ_matrix);
  HYPRE_IJVectorInitialize(IJ_rhs);
  HYPRE_IJVectorInitialize(IJ_lhs);

  std::vector<HYPRE_BigInt>    row_hypre_GIDs(n_local_rows_all);
  std::vector<double> rhs_values(n_local_rows_all);

  for(int iDOF=0; iDOF<2; ++iDOF)
  {
    for(int row=0; row<n_local_rows[iDOF]; ++row)
    {
      const globalIndex global_row = p_matrix[iDOF][iDOF]->GRID64(row);
      const HYPRE_Int hypre_row = GID_trilinos_to_hypre[iDOF].at(global_row);
      int n_entries;

      // fill hypre matrix row-by-row
      for (int jDOF=0; jDOF<2; ++jDOF)
      {
        p_matrix[iDOF][jDOF]->ExtractGlobalRowCopy(global_row, nnz_max, n_entries, col_values.data(), col_indices.data());
        for (int col = 0; col < n_entries; ++col)
        {
          col_indices[col] = GID_trilinos_to_hypre[jDOF].at(col_indices[col]);
        }
        HYPRE_Int num_entries = n_entries;
        HYPRE_IJMatrixSetValues(IJ_matrix, 1, &num_entries, &hypre_row, col_indices.data(), col_values.data());
      }
    }

    // fill the rhs row-by-row
    int shift = n_local_rows[0];
    for (int row=0; row<n_local_rows[iDOF]; ++row)
    {
      const int global_row = p_matrix[iDOF][iDOF]->GRID64(row);
      const globalIndex hypre_row = GID_trilinos_to_hypre[iDOF].at(global_row);
      int row_idx = row + (iDOF == 0 ? 0 : 1) * shift;
      row_hypre_GIDs[row_idx] = hypre_row;

      *(rhs_values.begin() + row_idx) = p_rhs[iDOF]->Values()[row];
    }
  }

  HYPRE_IJVectorSetValues(IJ_rhs, n_local_rows_all, row_hypre_GIDs.data(), rhs_values.data());
  hypre_IJVectorZeroValues(IJ_lhs);
  //printf("Done copy matrix and rhs values\n");

  // finalize matrix and make it ready to use
  HYPRE_IJMatrixAssemble(IJ_matrix);
  HYPRE_IJMatrixGetObject(IJ_matrix, (void**) &parcsr_matrix); // Get the parcsr matrix object to use
  HYPRE_IJVectorGetObject(IJ_rhs, (void **) &par_rhs);
  HYPRE_IJVectorGetObject(IJ_lhs, (void **) &par_lhs);

  /*
  if (print_matrix < 1)
  {
    print_matrix++;
    //hypre_ParCSRMatrixPrintIJ(parcsr_matrix,0,0,"full_mat");
    //hypre_ParVectorPrintIJ(par_rhs,0,"full_rhs");
    for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
    {
      char fname[256];
      sprintf(fname,"%s%s_block.%05d",i==0?"U":"P",j==0?"U":"P",my_id);
      std::ofstream myfile(fname, std::ios::out);
      p_matrix[i][j]->Print(myfile);
      myfile.close();
      //EpetraExt::RowMatrixToMatrixMarketFile(fname,*p_matrix[i][j]);
    }
  }
  */
  //######################################### END COPYING HYPRE MATRIX #########################################

  //########################################## BEGIN COPYING UU MATRIX #########################################
  if (ordering == 0)
  {
    if(newtonIter==0)
    {
      if (IJ_matrix_uu != nullptr)
      {
        HYPRE_IJMatrixDestroy(IJ_matrix_uu);
        IJ_matrix_uu = nullptr;
      }
      m_blockDiagUU.reset(new ParallelMatrix());
      LAIHelperFunctions::SeparateComponentFilter(m_solidSolver->getSystemMatrix(),*m_blockDiagUU,3);
      const Epetra_FECrsMatrix *blockDiagUU = m_blockDiagUU->unwrappedPointer();
      
      HYPRE_Int ilower = blockDiagUU->RowMatrixRowMap().MinMyGID64();
      HYPRE_Int iupper = blockDiagUU->RowMatrixRowMap().MaxMyGID64();
      HYPRE_IJMatrixCreate(MPI_COMM_GEOSX, ilower, iupper, ilower, iupper, &IJ_matrix_uu);
      HYPRE_IJMatrixSetObjectType(IJ_matrix_uu, HYPRE_PARCSR);
      HYPRE_IJMatrixInitialize(IJ_matrix_uu);
      for(int i = 0; i < blockDiagUU->NumMyRows(); i++)
      {
        int numElements;
        blockDiagUU->NumMyRowEntries(i,numElements);
        globalIndex global_row = blockDiagUU->GRID64(i);
        std::vector<HYPRE_Int> indices; indices.resize(numElements);
        std::vector<double> values; values.resize(numElements);
        int numEntries;
        blockDiagUU->ExtractGlobalRowCopy(global_row, numElements, numEntries, values.data(), indices.data());
        HYPRE_Int n_entries = numEntries;
        HYPRE_IJMatrixSetValues(IJ_matrix_uu, 1, &n_entries, &global_row, indices.data(), values.data());
      }
      HYPRE_IJMatrixAssemble(IJ_matrix_uu);
      HYPRE_IJMatrixGetObject(IJ_matrix_uu, (void**) &parcsr_uu);
    } 
  }
  //########################################## END COPYING UU MATRIX #########################################
  auxTime = clock.stop();
  GEOSX_MARK_END(COPY_HYPRE_MATRIX);

  // ############# MGR OPTIONS ******************
  /* mgr options */
  HYPRE_Int mgr_bsize = 2;
  HYPRE_Int mgr_nlevels = 1;
  HYPRE_Int mgr_non_c_to_f = 1;
  HYPRE_Int *mgr_idx_array = NULL;
  HYPRE_Int *mgr_num_cindexes = NULL; 
  HYPRE_Int **mgr_cindexes = NULL;
  HYPRE_Int mgr_relax_type = 0;
  HYPRE_Int mgr_num_relax_sweeps = 1;
  HYPRE_Int mgr_num_interp_sweeps = 0;
  HYPRE_Int mgr_gsmooth_type = 16;
  HYPRE_Int mgr_num_gsmooth_sweeps = 0;
  HYPRE_Int mgr_num_restrict_sweeps = 0;
  HYPRE_Int *mgr_level_interp_type = hypre_CTAlloc(HYPRE_Int, mgr_nlevels, HYPRE_MEMORY_HOST);
  mgr_level_interp_type[0] = 2;
  HYPRE_Int *mgr_level_restrict_type = hypre_CTAlloc(HYPRE_Int, mgr_nlevels, HYPRE_MEMORY_HOST);
  mgr_level_restrict_type[0] = 0;
  HYPRE_Int *mgr_coarse_grid_method = hypre_CTAlloc(HYPRE_Int, mgr_nlevels, HYPRE_MEMORY_HOST);
  mgr_coarse_grid_method[0] = 0;

  HYPRE_Int *lv1 = hypre_CTAlloc(HYPRE_Int, mgr_bsize, HYPRE_MEMORY_HOST);
  lv1[0] = ordering == 0 ? 1 : 0;
  mgr_cindexes = hypre_CTAlloc(HYPRE_Int*, mgr_nlevels, HYPRE_MEMORY_HOST);
  mgr_cindexes[0] = lv1;
  mgr_num_cindexes = hypre_CTAlloc(HYPRE_Int, mgr_nlevels, HYPRE_MEMORY_HOST);
  mgr_num_cindexes[0] = 1;

  mgr_idx_array = hypre_CTAlloc(HYPRE_Int, mgr_bsize, HYPRE_MEMORY_HOST);
  HYPRE_Int ilower = GID_trilinos_to_hypre[0].at(p_matrix[0][0]->GRID64(0));
  mgr_idx_array[0] = ilower;
  mgr_idx_array[1] = n_local_rows[0] + ilower;

  HYPRE_Int *mgr_level_frelax_method = hypre_CTAlloc(HYPRE_Int, mgr_nlevels, HYPRE_MEMORY_HOST);
  mgr_level_frelax_method[0] = 99;

  HYPRE_ParCSRGMRESCreate(MPI_COMM_GEOSX, &pgmres_solver);
  HYPRE_GMRESSetKDim(pgmres_solver, params->m_maxIters);
  HYPRE_GMRESSetMaxIter(pgmres_solver, params->m_maxIters);
  HYPRE_GMRESSetTol(pgmres_solver, params->m_krylovTol);
  HYPRE_GMRESSetLogging(pgmres_solver, 1);
  HYPRE_GMRESSetPrintLevel(pgmres_solver, 0);

  HYPRE_MGRCreate(&mgr_precond);
  /* set MGR data by block */
  HYPRE_MGRSetCpointsByContiguousBlock( mgr_precond, mgr_bsize, mgr_nlevels, mgr_idx_array, mgr_num_cindexes, mgr_cindexes);
  /* set intermediate coarse grid strategy */
  HYPRE_MGRSetNonCpointsToFpoints(mgr_precond, mgr_non_c_to_f);
  /* set F relaxation strategy */
  HYPRE_MGRSetLevelFRelaxMethod(mgr_precond, mgr_level_frelax_method);
  /* set relax type for single level F-relaxation and post-relaxation */
  HYPRE_MGRSetRelaxType(mgr_precond, mgr_relax_type);
  HYPRE_MGRSetNumRelaxSweeps(mgr_precond, mgr_num_relax_sweeps);
  /* set restrict type */
  HYPRE_MGRSetLevelRestrictType(mgr_precond, mgr_level_restrict_type);
  HYPRE_MGRSetNumRestrictSweeps(mgr_precond, mgr_num_restrict_sweeps);
  /* set interpolation type */
  HYPRE_MGRSetLevelInterpType(mgr_precond, mgr_level_interp_type);
  HYPRE_MGRSetNumInterpSweeps(mgr_precond, mgr_num_interp_sweeps);
  /* set P_max_elmts for coarse grid */
  //HYPRE_MGRSetPMaxElmts(mgr_precond, P_max_elmts);
  /* set print level */
  HYPRE_MGRSetPrintLevel(mgr_precond, 0);
  /* set max iterations */
  HYPRE_MGRSetMaxIter(mgr_precond, 1);
  HYPRE_MGRSetTol(mgr_precond, 0.0);
  //HYPRE_MGRSetCoarseGridMethod(mgr_precond, mgr_coarse_grid_method);

  HYPRE_MGRSetGlobalsmoothType(mgr_precond, mgr_gsmooth_type);
  HYPRE_MGRSetMaxGlobalsmoothIters( mgr_precond, mgr_num_gsmooth_sweeps );   

  /* create AMG coarse grid solver */
  HYPRE_BoomerAMGCreate(&cg_amg_solver); 
  if (ordering == 0)
  {
    HYPRE_BoomerAMGSetPrintLevel(cg_amg_solver, 0);
    HYPRE_BoomerAMGSetRelaxOrder(cg_amg_solver, 1);
    HYPRE_BoomerAMGSetMaxIter(cg_amg_solver, 1);
    HYPRE_BoomerAMGSetNumFunctions(cg_amg_solver, 1);
    HYPRE_BoomerAMGSetNumSweeps(cg_amg_solver, 3);
  }
  else
  {
    HYPRE_BoomerAMGSetPrintLevel(cg_amg_solver, 0);
    HYPRE_BoomerAMGSetMaxIter(cg_amg_solver, 1);
    HYPRE_BoomerAMGSetRelaxOrder(cg_amg_solver, 1);
    HYPRE_BoomerAMGSetAggNumLevels(cg_amg_solver, 1);
    HYPRE_BoomerAMGSetNumFunctions(cg_amg_solver, 3);
    HYPRE_BoomerAMGSetNumSweeps(cg_amg_solver, 3);
    //HYPRE_BoomerAMGSetCoarsenType(cg_amg_solver, 6);
    //HYPRE_BoomerAMGSetRelaxType(cg_amg_solver, 3);
    //HYPRE_BoomerAMGSetInterpType(cg_amg_solver, 0);
    //HYPRE_BoomerAMGSetPMaxElmts(cg_amg_solver, 0);
    /*
    HYPRE_BoomerAMGSetSmoothType(cg_amg_solver, 6);
    //HYPRE_BoomerAMGSetEuLevel(cg_amg_solver, 5);
    HYPRE_BoomerAMGSetSmoothNumLevels(cg_amg_solver, 5);
    HYPRE_BoomerAMGSetSmoothNumSweeps(cg_amg_solver, 1);
    HYPRE_BoomerAMGSetSchwarzUseNonSymm(cg_amg_solver, 1);
    */
  }
  /* set the MGR coarse solver. Comment out to use default CG solver in MGR */
  HYPRE_MGRSetCoarseSolver(mgr_precond, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, cg_amg_solver);

  // set fine grid solver
  if (ordering == 0)
  {
    clock.start(true);
    /* create AMG solver for F-relaxation */
    HYPRE_BoomerAMGCreate(&uu_amg_solver); 
    HYPRE_BoomerAMGSetPrintLevel(uu_amg_solver, 0);
    HYPRE_BoomerAMGSetRelaxOrder(uu_amg_solver, 1);
    HYPRE_BoomerAMGSetMaxIter(uu_amg_solver, 1);
    HYPRE_BoomerAMGSetNumFunctions(uu_amg_solver, 3);
    //HYPRE_BoomerAMGSetStrongThreshold(uu_amg_solver, 0.25);
    HYPRE_BoomerAMGSetAggNumLevels(uu_amg_solver, 1);
    HYPRE_BoomerAMGSetNumSweeps(uu_amg_solver, 1);
    //HYPRE_BoomerAMGSetRelaxType(uu_amg_solver, 3);

    HYPRE_BoomerAMGSetup
      (uu_amg_solver, parcsr_uu, par_rhs_uu, par_lhs_uu);
    setupTime = clock.stop();

    HYPRE_MGRSetFSolver(mgr_precond, HYPRE_BoomerAMGSolve, HYPRE_BoomerAMGSetup, uu_amg_solver);
  }

  /* setup MGR-PCG solver */
  HYPRE_GMRESSetPrecond(pgmres_solver,
                        (HYPRE_PtrToSolverFcn) HYPRE_MGRSolve,
                        (HYPRE_PtrToSolverFcn) HYPRE_MGRSetup,
                        mgr_precond);

  GEOSX_MARK_BEGIN(MGR_SETUP);
  clock.start(true);
  HYPRE_GMRESSetup(pgmres_solver, (HYPRE_Matrix)parcsr_matrix, (HYPRE_Vector)par_rhs, (HYPRE_Vector)par_lhs);
  setupTime += clock.stop();
  GEOSX_MARK_END(MGR_SETUP);

  GEOSX_MARK_BEGIN(MGR_SOLVE);
  clock.start(true);
  HYPRE_GMRESSolve
    (pgmres_solver, (HYPRE_Matrix)parcsr_matrix, (HYPRE_Vector)par_rhs, (HYPRE_Vector)par_lhs);
  solveTime = clock.stop();
  GEOSX_MARK_END(MGR_SOLVE);

  HYPRE_Int num_iterations;
  HYPRE_Real final_res_norm;
  HYPRE_GMRESGetNumIterations(pgmres_solver, &num_iterations);
  HYPRE_GMRESGetFinalRelativeResidualNorm(pgmres_solver, &final_res_norm);

  if (my_id == 0)
    printf("Using hypreMGR, cycle = %d, iters = %lld, Final Residual = %1.5e\n", n_cycles, num_iterations, final_res_norm);

  params->m_numKrylovIter = num_iterations;

  HYPRE_MGRDestroy(mgr_precond);
  HYPRE_BoomerAMGDestroy(cg_amg_solver);
  HYPRE_ParCSRGMRESDestroy(pgmres_solver);
  hypre_TFree(lv1, HYPRE_MEMORY_HOST);
  hypre_TFree(mgr_cindexes, HYPRE_MEMORY_HOST);
  hypre_TFree(mgr_num_cindexes, HYPRE_MEMORY_HOST);
  hypre_TFree(mgr_idx_array, HYPRE_MEMORY_HOST);

  //############################################ BEGIN COPYING SOLUTION #####################################
  // copy the hypre solution mapping back to the original Trilinos ordering
  HYPRE_IJVectorGetValues(IJ_lhs, n_local_rows_all, row_hypre_GIDs.data(), rhs_values.data());

  for(int iDOF=0; iDOF<2; ++iDOF)
  {
    // fill the rhs row-by-row
    for (int row=0; row<n_local_rows[iDOF]; ++row)
    {
      int row_idx = (iDOF == 0) ? row : row + n_local_rows[0];
      p_solution[iDOF]->Values()[row] = *(rhs_values.begin() + row_idx);
    }
  }

  gettimeofday(&tim, nullptr);
  const real64 t_end = tim.tv_sec + (tim.tv_usec / 1000000.0);
  if( getLogLevel()>=2 )
  {
    GEOSX_LOG_RANK_0("\t\tLinear Solver | Iter = " << params->m_numKrylovIter <<
                    " | TargetReduction " << params->m_krylovTol <<
                    " | AuxTime " << auxTime <<
                    " | SetupTime " << setupTime <<
                    " | SolveTime " << solveTime <<
                    " | TotalTime " << t_end - t_start);
  }

#else

    // SCHEME CHOICES
    //
    // there are several flags to control solver behavior.
    // these should be compared in a scaling study.
    //
    // -- whether to use a block diagonal or a 
    //    block triangular preconditioner.
    // -- whether to use BiCGstab or GMRES for the
    //    krylov solver.  GMRES is generally more robust,
    //    BiCGstab sometimes shows better parallel performance.
    //    false is probably better.

  timeval tim;
  gettimeofday(&tim, nullptr);
  const real64 t_start = tim.tv_sec + (tim.tv_usec / 1000000.0);

  const bool use_diagonal_prec = true;
  const bool use_bicgstab      = params->m_useBicgstab;

  // create separate displacement component matrix

  clock.start(true);
  if(newtonIter==0)
  {
    m_blockDiagUU.reset(new ParallelMatrix());
    LAIHelperFunctions::SeparateComponentFilter(m_solidSolver->getSystemMatrix(),*m_blockDiagUU,3);
  }

    // create schur complement approximation matrix
  Epetra_CrsMatrix* schurApproxPP = NULL; // confirm we delete this at end of function!
  {
    Epetra_Vector diag(p_matrix[0][0]->RowMap());
    Epetra_Vector diagInv(p_matrix[0][0]->RowMap());
 
    p_matrix[0][0]->ExtractDiagonalCopy(diag); 
    diagInv.Reciprocal(diag);
 
    Epetra_FECrsMatrix DB(*p_matrix[0][1]);
    DB.LeftScale(diagInv);
    DB.FillComplete();

    Epetra_FECrsMatrix BtDB(Epetra_DataAccess::Copy,p_matrix[1][1]->RowMap(),1); 
    EpetraExt::MatrixMatrix::Multiply(*p_matrix[1][0],false,DB,false,BtDB);
    EpetraExt::MatrixMatrix::Add(BtDB,false,-1.0,*p_matrix[1][1],false,1.0,schurApproxPP);

    schurApproxPP->FillComplete();
  }
  auxTime = clock.stop();
  GEOSX_MARK_END(Setup);

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

  GEOSX_MARK_BEGIN(THYRA_SETUP);

  RCP<const Thyra::LinearOpBase<double> >  matrix_block[2][2];
  RCP<Thyra::MultiVectorBase<double> >     lhs_block[2];
  RCP<Thyra::MultiVectorBase<double> >     rhs_block[2];

  for(unsigned i=0; i<2; ++i)
  for(unsigned j=0; j<2; ++j)
  {
    RCP<Epetra_Operator> mmm (&*p_matrix[i][j],false);
    matrix_block[i][j] = Thyra::epetraLinearOp(mmm);
  }

  RCP<Epetra_Operator> bbb(m_blockDiagUU->unwrappedPointer(),false);
  RCP<Epetra_Operator> ppp(schurApproxPP,false);

  RCP<const Thyra::LinearOpBase<double> >  blockDiagOp = Thyra::epetraLinearOp(bbb);
  RCP<const Thyra::LinearOpBase<double> >  schurOp = Thyra::epetraLinearOp(ppp);

  for(unsigned i=0; i<2; ++i)
  {
    RCP<Epetra_MultiVector> lll (&*p_solution[i],false);
    RCP<Epetra_MultiVector> rrr (&*p_rhs[i],false);

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
    // rhs and lhs is a little uglier. 

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

  GEOSX_MARK_END(THYRA_SETUP);

    // for the preconditioner, we need two approximate inverses,
    // we store both "sub operators" in a 1x2 array:

  RCP<const Thyra::LinearOpBase<double> > sub_op[2];

  clock.start(true);
  GEOSX_MARK_BEGIN(PRECONDITIONER);

  for(unsigned i=0; i<2; ++i) // loop over diagonal blocks
  {
    RCP<Teuchos::ParameterList> list = rcp(new Teuchos::ParameterList("precond_list"),true);

    if(params->m_useMLPrecond)
    {
      list->set("Preconditioner Type","ML");
      list->sublist("Preconditioner Types").sublist("ML").set("Base Method Defaults","SA");
      list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("PDE equations",(i==0?3:1));
      list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("ML output", 0);
      list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("aggregation: type","Uncoupled");
      list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("aggregation: threshold",1e-3);

      if(i==0) // smoother for mechanics block
      {
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("smoother: type","Chebyshev");
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("smoother: sweeps",3);
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("coarse: type","Chebyshev");
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("coarse: sweeps",3);
      }
      else // smoother for flow block
      {
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("smoother: type","Chebyshev");
        list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("smoother: sweeps",3);
      }

    }
    else // use ILU for both blocks
    {
      list->set("Preconditioner Type","Ifpack");
      list->sublist("Preconditioner Types").sublist("Ifpack").set("Prec Type","ILU");
    }

    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList(list);

    RCP<const Thyra::PreconditionerFactoryBase<double> > strategy = createPreconditioningStrategy(builder);
    RCP<Thyra::PreconditionerBase<double> > tmp;

    if(i==0)
      tmp = prec(*strategy,blockDiagOp);
    else
      tmp = prec(*strategy,schurOp);
      //tmp = prec(*strategy,matrix_block[i][i]);

    sub_op[i] = tmp->getUnspecifiedPrecOp();
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

    Linv = Thyra::block2x2(eye_00,zero_01,mB2Ainv,eye_11);
    Dinv = Thyra::block2x2(sub_op[0],zero_01,zero_10,sub_op[1]);
    Uinv = Thyra::block2x2(eye_00,mAinvB1,zero_10,eye_11);

    //preconditioner = Thyra::multiply(Uinv,Dinv);
    //preconditioner = Thyra::multiply(Dinv,Linv);
    preconditioner = Thyra::multiply(Uinv,Dinv,Linv);
  }

  GEOSX_MARK_END(PRECONDITIONER);
  setupTime = clock.stop();

    // define solver strategy for blocked system. this is
    // similar but slightly different from the sub operator
    // construction, since now we have a user defined preconditioner

  {
    RCP<Teuchos::ParameterList> list = rcp(new Teuchos::ParameterList("list"));
    
      list->set("Linear Solver Type","AztecOO");
      list->set("Preconditioner Type","None"); // will use user-defined P
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Max Iterations",params->m_maxIters);
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Tolerance",params->m_krylovTol);

      if(use_bicgstab)
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","BiCGStab");
      else
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","GMRES");

      if( params->getLogLevel()>=2 )
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",1);

    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList(list);

    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > strategy = createLinearSolveStrategy(builder);
    RCP<Thyra::LinearOpWithSolveBase<double> > solver = strategy->createOp();

    Thyra::initializePreconditionedOp<double>(*strategy,
                                               matrix,
                                               Thyra::rightPrec<double>(preconditioner),
                                               solver.ptr());

    clock.start(true);
    GEOSX_MARK_BEGIN(SOLVER);

      // !!!! Actual Solve !!!!
      Thyra::SolveStatus<double> status = solver->solve(Thyra::NOTRANS,*rhs,lhs.ptr());

    GEOSX_MARK_END(SOLVER);
    solveTime = clock.stop();
    params->m_numKrylovIter = status.extraParameters->get<int>("Iteration Count");

    gettimeofday(&tim, nullptr);
    const real64 t_end = tim.tv_sec + (tim.tv_usec / 1000000.0);

    if( getLogLevel()>=2 )
    {
      GEOSX_LOG_RANK_0("\t\tLinear Solver | Iter = " << params->m_numKrylovIter <<
                      " | TargetReduction " << params->m_krylovTol <<
                      " | AuxTime " << auxTime <<
                      " | SetupTime " << setupTime <<
                      " | SolveTime " << solveTime <<
                      " | TotalTime" << t_end - t_start);
    }
  }
  delete schurApproxPP;
#endif

  p_solution[1]->Scale(m_pressureScaling);
  p_rhs[1]->Scale(1/(m_pressureScaling*m_densityScaling));

  /*
  if ((n_cycles % 10) == 0 && newtonIter == 0)
  {
    char fname[256];
    sprintf(fname, "solution_u_%03d.txt", n_cycles);
    EpetraExt::MultiVectorToMatrixMarketFile(fname, *p_solution[0]);
    sprintf(fname, "solution_p_%03d.txt", n_cycles);
    EpetraExt::MultiVectorToMatrixMarketFile(fname, *p_solution[1]);
  }
  */

  //TODO: remove all this once everything is working
  if( getLogLevel() == 2 )
  {
    /*
    ParallelVector permutedSol;
    ParallelVector const & solution = m_solidSolver->getSystemSolution();
    permutedSol.createWithLocalSize(m_solidSolver->getSystemMatrix().localRows(), MPI_COMM_GEOSX);
    m_permutationMatrix0.multiply(solution, permutedSol);
    permutedSol.close();
    */

    /*
    GEOSX_LOG_RANK_0("***********************************************************");
    GEOSX_LOG_RANK_0("solution0");
    GEOSX_LOG_RANK_0("***********************************************************");
    solution.print(std::cout);
    std::cout<<std::endl;
    MPI_Barrier(MPI_COMM_GEOSX);

    GEOSX_LOG_RANK_0("***********************************************************");
    GEOSX_LOG_RANK_0("solution0");
    GEOSX_LOG_RANK_0("***********************************************************");
    permutedSol.print(std::cout);

    GEOSX_LOG_RANK_0("***********************************************************");
    GEOSX_LOG_RANK_0("solution1");
    GEOSX_LOG_RANK_0("***********************************************************");
    p_solution[1]->Print(std::cout);
    */
  }
}

real64
HydrofractureSolver::ScalingForSystemSolution( DomainPartition const * const domain,
                                                 DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                                                 ParallelVector const & GEOSX_UNUSED_ARG( solution ) )
{
  return m_solidSolver->ScalingForSystemSolution( domain,
                                                  m_solidSolver->getDofManager(),
                                                  m_solidSolver->getSystemSolution() );
}

void HydrofractureSolver::SetNextDt( real64 const & currentDt ,
                                     real64 & nextDt )
{
  SolverBase * const surfaceGenerator =  this->getParent()->GetGroup<SolverBase>("SurfaceGen");

  if (m_numResolves[0] == 0 & m_numResolves[1] == 0)
  {
    this->SetNextDtBasedOnNewtonIter(currentDt, nextDt);
  } else
  {
    nextDt = surfaceGenerator->GetTimestepRequest() < 1e99 ? surfaceGenerator->GetTimestepRequest() : currentDt;
  }
  GEOSX_LOG_RANK_0(this->getName() << ": nextDt request is "  << nextDt);
}

void HydrofractureSolver::initializeNewFaceElements( DomainPartition const &  )
{
//  m_flowSolver->
}

REGISTER_CATALOG_ENTRY( SolverBase, HydrofractureSolver, std::string const &, Group * const )
} /* namespace geosx */
