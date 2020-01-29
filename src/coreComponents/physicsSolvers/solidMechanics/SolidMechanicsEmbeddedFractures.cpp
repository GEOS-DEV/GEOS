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
      m_solidSolver(nullptr)
{
  registerWrapper(viewKeyStruct::solidSolverNameString, &m_solidSolverName, 0)->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Name of the solid mechanics solver in the rock matrix");

  registerWrapper(viewKeyStruct::contactRelationNameString, &m_contactRelationName, 0)->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Name of contact relation to enforce constraints on fracture boundary.");

}

SolidMechanicsEmbeddedFractures::~SolidMechanicsEmbeddedFractures()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsEmbeddedFractures::RegisterDataOnMesh( dataRepository::Group * const  MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {

    elemManager->forElementRegions<EmbeddedSurfaceRegion>( [&] ( EmbeddedSurfaceRegion * const region )
       {
         region->forElementSubRegions<EmbeddedSurfaceSubRegion>( [&]( EmbeddedSurfaceSubRegion * const subRegion )
         {
           subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::dispJumpString )->setPlotLevel(PlotLevel::LEVEL_0)->
               setDefaultValue( region->getDefaultAperture() );
           subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::deltaDispJumpString );
         });
       });
  }
}

void SolidMechanicsEmbeddedFractures::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  m_solidSolver->ResetStateToBeginningOfStep(domain);
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

real64 SolidMechanicsEmbeddedFractures::SolverStep( real64 const & time_n,
                                                    real64 const & dt,
                                                    int const cycleNumber,
                                                    DomainPartition * const domain )
{
  real64 dtReturn = dt;

  ImplicitStepSetup( time_n,
                     dt,
                     domain,
                     m_dofManager,
                     m_matrix,
                     m_rhs,
                     m_solution );

  SetupSystem( domain,
               m_dofManager,
               m_matrix,
               m_rhs,
               m_solution  );

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

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dtReturn, domain );

  return dtReturn;
}

void SolidMechanicsEmbeddedFractures::SetupDofs( DomainPartition const * const domain,
                                                 DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );

  MeshLevel const * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  array1d<string> regions;
  elemManager->forElementRegions<EmbeddeSurfaceRegion>( [&]( EmbeddeSurfaceRegion const * const region ) {
    regions.push_back( region->getName() );
  } );

  dofManager.addField( keys::DispJump,
                       DofManager::Location::Elem,
                       3,
                       regions );

  dofManager.addCoupling( viewKeyStruct::dispJumpString,
                          viewKeyStruct::dispJumpString,
                          DofManager::Connectivity::Elem,
                          regions );
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


  // setup coupled DofManager
  m_dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );

  // By not calling dofManager.reorderByRank(), we keep separate dof numbering for each field,
  // which allows constructing separate sparsity patterns for off-diagonal blocks of the matrix.
  // Once the solver moves to monolithic matrix, we can remove this method and just use SolverBase::SetupSystem.

  m_matrix11.createWithLocalSize( m_dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  m_dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  1,
                                  MPI_COMM_GEOSX);

  m_matrix01.createWithLocalSize( m_solidSolver->getSystemMatrix().localRows(),
                                  m_dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  9,
                                  MPI_COMM_GEOSX);

  m_matrix10.createWithLocalSize( m_dofManager.numLocalDofs(viewKeyStruct::dispJumpString),
                                  m_solidSolver->getSystemMatrix().localRows(),
                                  24,
                                  MPI_COMM_GEOSX);

  //dofManager.setSparsityPattern( m_matrix01, keys::TotalDisplacement, keys::DispJump ); I am guessing that this won't work coz coupling has not been created.
  //dofManager.setSparsityPattern( m_matrix10, keys::DispJump, keys::TotalDisplacement );
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


  if( getLogLevel() == 2 )
    {
      // Before outputting anything generate permuation matrix and permute.
      MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
      NodeManager * const nodeManager = mesh->getNodeManager();

      LAIHelperFunctions::CreatePermutationMatrix(nodeManager,
                                                  m_solidSolver->getSystemMatrix().localRows(),
                                                  m_solidSolver->getSystemMatrix().localCols(),
                                                  3,
                                                  m_solidSolver->getDofManager().getKey( keys::TotalDisplacement ),
                                                  m_permutationMatrix0);

      GEOSX_LOG_RANK_0("***********************************************************");
      GEOSX_LOG_RANK_0("matrix00");
      GEOSX_LOG_RANK_0("***********************************************************");
      LAIHelperFunctions::PrintPermutedMatrix(m_solidSolver->getSystemMatrix(), m_permutationMatrix0, std::cout);
      MpiWrapper::Barrier();

      GEOSX_LOG_RANK_0("***********************************************************");
      GEOSX_LOG_RANK_0("residual0");
      GEOSX_LOG_RANK_0("***********************************************************");
      LAIHelperFunctions::PrintPermutedVector(m_solidSolver->getSystemRhs(), m_permutationMatrix0, std::cout);
      MpiWrapper::Barrier();
    }
}

real64
SolidMechanicsEmbeddedFractures::
CalculateResidualNorm( DomainPartition const * const domain,
                       DofManager const & GEOSX_UNUSED_ARG( dofManager ),
                       ParallelVector const & GEOSX_UNUSED_ARG( rhs ) )
{
  GEOSX_MARK_FUNCTION;

  real64 const solidResidualNorm = m_solidSolver->CalculateResidualNorm( domain,
                                                                     m_solidSolver->getDofManager(),
                                                                     m_solidSolver->getSystemRhs() );

  GEOSX_LOG_RANK_0("residual = "<< solidResidualNorm);

  return solidResidualNorm;
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



