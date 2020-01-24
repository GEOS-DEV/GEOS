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

/*
 * SolidMechanicsEmbeddedFractures.cpp
 */

#include "SolidMechanicsEmbeddedFractures.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SolidMechanicsEmbeddedFractures::SolidMechanicsEmbeddedFractures( const std::string& name,
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

  registerWrapper(viewKeyStruct::contactRelationNameString, &m_contactRelationName, 0)->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Name of contact relation to enforce constraints on fracture boundary.");

  registerWrapper(viewKeyStruct::maxNumResolvesString, &m_maxNumResolves, 0)->
      setApplyDefaultValue(10)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Value to indicate how many resolves may be executed to perform surface generation after the execution of flow and mechanics solver. ");

}

SolidMechanicsEmbeddedFractures::~SolidMechanicsEmbeddedFractures()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsEmbeddedFractures::RegisterDataOnMesh( dataRepository::Group * const GEOSX_UNUSED_ARG( MeshBodies ) )
{

}

void SolidMechanicsEmbeddedFractures::ImplicitStepSetup( real64 const & time_n,
                                                         real64 const & dt,
                                                         DomainPartition * const domain,
                                                         DofManager & GEOSX_UNUSED_ARG( dofManager ),
                                                         ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                                         ParallelVector & GEOSX_UNUSED_ARG( rhs ),
                                                         ParallelVector & GEOSX_UNUSED_ARG( solution ) )
{
  m_solidSolver = this->getParent()->GetGroup<SolidMechanicsLagrangianFEM>(m_solidSolverName);

  m_solidSolver->ImplicitStepSetup( time_n, dt, domain,
                                    m_solidSolver->getDofManager(),
                                    m_solidSolver->getSystemMatrix(),
                                    m_solidSolver->getSystemRhs(),
                                    m_solidSolver->getSystemSolution() );
}

void SolidMechanicsEmbeddedFractures::ImplicitStepComplete( real64 const& time_n,
                                                            real64 const& dt,
                                                            DomainPartition * const domain)
{
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
}

void SolidMechanicsEmbeddedFractures::PostProcessInput()
{

}

void SolidMechanicsEmbeddedFractures::InitializePostInitialConditions_PreSubGroups(Group * const GEOSX_UNUSED_ARG( problemManager ) )
{

}

void SolidMechanicsEmbeddedFractures::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  m_solidSolver->ResetStateToBeginningOfStep(domain);
}

real64 SolidMechanicsEmbeddedFractures::SolverStep( real64 const & time_n,
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

    int const maxNumResolves = m_maxNumResolves;
    for( int solveIter=0 ; solveIter<maxNumResolves ; ++solveIter )
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
        break;
      }
    }

    // final step for completion of timestep. typically secondary variable updates and cleanup.
    ImplicitStepComplete( time_n, dtReturn, domain );
  }
  return dtReturn;
}

void SolidMechanicsEmbeddedFractures::SetupDofs( DomainPartition const * const domain,
                                                 DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString,
                          DofManager::Connectivity::Elem );
}

void SolidMechanicsEmbeddedFractures::SetupSystem( DomainPartition * const domain,
                                                   DofManager & GEOSX_UNUSED_ARG( dofManager ),
                                                   ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                                   ParallelVector & GEOSX_UNUSED_ARG( rhs ),
                                                   ParallelVector & GEOSX_UNUSED_ARG( solution ) )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->SetupSystem( domain,
                              m_solidSolver->getDofManager(),
                              m_solidSolver->getSystemMatrix(),
                              m_solidSolver->getSystemRhs(),
                              m_solidSolver->getSystemSolution() );



  // TODO: once we move to a monolithic matrix, we can just use SolverBase implementation

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

  m_matrix01.close();
  m_matrix10.close();

}

void SolidMechanicsEmbeddedFractures::AssembleSystem( real64 const time,
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


  AssembleForceResidualDerivativeWrtPressure( domain, &m_matrix01, &(m_solidSolver->getSystemRhs()) );

  AssembleFluidMassResidualDerivativeWrtDisplacement( domain, &m_matrix10, &(m_flowSolver->getSystemRhs()) );

}

void SolidMechanicsEmbeddedFractures::ApplyBoundaryConditions( real64 const time,
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
}

real64
SolidMechanicsEmbeddedFractures::
CalculateResidualNorm( DomainPartition const * const domain,
                       DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                       ParallelVector const & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  real64 const solidResidual = m_solidSolver->CalculateResidualNorm( domain,
                                                                     m_solidSolver->getDofManager(),
                                                                     m_solidSolver->getSystemRhs() );

  GEOSX_LOG_RANK_0("residuals for fluid, solid: "<<fluidResidual<<", "<<solidResidual);

  return fluidResidual + solidResidual;
}



void
SolidMechanicsEmbeddedFractures::
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

}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsEmbeddedFractures, std::string const &, Group * const )
} /* namespace geosx */



