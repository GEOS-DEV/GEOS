/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PhaseFieldFractureSolver.cpp
 *
 */

#include "PhaseFieldFractureSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/simplePDE/PhaseFieldDamageFEM.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

PhaseFieldFractureSolver::PhaseFieldFractureSolver( const std::string & name,
                                                    Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_damageSolverName(),
  m_couplingTypeOption( CouplingTypeOption::FixedStress )

{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription(
    "Name of the solid mechanics solver to use in the PhaseFieldFracture solver" );

  registerWrapper( viewKeyStruct::damageSolverNameString, &m_damageSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription(
    "Name of the damage mechanics solver to use in the PhaseFieldFracture solver" );

  registerWrapper( viewKeyStruct::couplingTypeOptionString, &m_couplingTypeOption )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Coupling option. Valid options:\n* " + EnumStrings< CouplingTypeOption >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::subcyclingOptionString, &m_subcyclingOption )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "turn on subcycling on each load step" );

}

void PhaseFieldFractureSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    ElementRegionManager * const elemManager = mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getElemManager();

    elemManager->forElementSubRegions< CellElementSubRegion,
                                       FaceElementSubRegion >( [ &]( auto & elementSubRegion ) -> void
    {
      elementSubRegion.template registerWrapper< array1d< real64 > >( viewKeyStruct::totalMeanStressString )->
        setDescription( "Total Mean Stress" );
      elementSubRegion.template registerWrapper< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString )->
        setDescription( "Total Mean Stress" );
    } );
  }
}

void PhaseFieldFractureSolver::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                  DomainPartition & domain )
{
  MeshLevel * const mesh = domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const totalMeanStress =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::totalMeanStressString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > oldTotalMeanStress =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::oldTotalMeanStressString );

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [ &]( localIndex const er,
                                 localIndex const esr,
                                 localIndex const k ) -> void
  {
    oldTotalMeanStress[er][esr][k] = totalMeanStress[er][esr][k];
  } );
}

void PhaseFieldFractureSolver::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                     real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                     DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void PhaseFieldFractureSolver::PostProcessInput()
{
  if( m_couplingTypeOption == CouplingTypeOption::FixedStress )
  {
    // For this coupled solver the minimum number of Newton Iter should be 0 for both flow and solid solver otherwise it
    // will never converge.
    SolidMechanicsLagrangianFEM &
    solidSolver = *( this->getParent()->GetGroup( m_solidSolverName )->group_cast< SolidMechanicsLagrangianFEM * >() );
    integer & minNewtonIterSolid = solidSolver.getNonlinearSolverParameters().m_minIterNewton;

    PhaseFieldDamageFEM &
    damageSolver = *( this->getParent()->GetGroup( m_damageSolverName )->group_cast< PhaseFieldDamageFEM * >() );
    integer & minNewtonIterFluid = damageSolver.getNonlinearSolverParameters().m_minIterNewton;

    minNewtonIterSolid = 0;
    minNewtonIterFluid = 0;
  }
}

void PhaseFieldFractureSolver::InitializePostInitialConditions_PreSubGroups( Group * const )
{}

PhaseFieldFractureSolver::~PhaseFieldFractureSolver()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldFractureSolver::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel * const mesh = domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const totalMeanStress =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::totalMeanStressString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > oldTotalMeanStress =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::oldTotalMeanStressString );

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [ &]( localIndex const er,
                                 localIndex const esr,
                                 localIndex const k ) -> void
  {
    totalMeanStress[er][esr][k] = oldTotalMeanStress[er][esr][k];
  } );
}

real64 PhaseFieldFractureSolver::SolverStep( real64 const & time_n,
                                             real64 const & dt,
                                             int const cycleNumber,
                                             DomainPartition & domain )
{
  real64 dtReturn = dt;
  if( m_couplingTypeOption == CouplingTypeOption::FixedStress )
  {
    dtReturn = SplitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::TightlyCoupled )
  {
    GEOSX_ERROR( "CouplingTypeOption::FullyImplicit not yet implemented" );
  }
  return dtReturn;
}

real64 PhaseFieldFractureSolver::SplitOperatorStep( real64 const & time_n,
                                                    real64 const & dt,
                                                    integer const cycleNumber,
                                                    DomainPartition & domain )
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsLagrangianFEM &
  solidSolver = *( this->getParent()->GetGroup( m_solidSolverName )->group_cast< SolidMechanicsLagrangianFEM * >() );

  PhaseFieldDamageFEM &
  damageSolver = *( this->getParent()->GetGroup( m_damageSolverName )->group_cast< PhaseFieldDamageFEM * >() );

  damageSolver.SetupSystem( domain,
                            damageSolver.getDofManager(),
                            damageSolver.getLocalMatrix(),
                            damageSolver.getLocalRhs(),
                            damageSolver.getLocalSolution(),
                            true );

  solidSolver.SetupSystem( domain,
                           solidSolver.getDofManager(),
                           solidSolver.getLocalMatrix(),
                           solidSolver.getLocalRhs(),
                           solidSolver.getLocalSolution() );

  damageSolver.ImplicitStepSetup( time_n, dt, domain );

  solidSolver.ImplicitStepSetup( time_n, dt, domain );

  this->ImplicitStepSetup( time_n, dt, domain );

  NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
  integer & iter = solverParams.m_numNewtonIterations;
  iter = 0;
  bool isConverged = false;
  while( iter < solverParams.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      damageSolver.ResetStateToBeginningOfStep( domain );
      solidSolver.ResetStateToBeginningOfStep( domain );
      ResetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", MechanicsSolver: " );

    solidSolver.ResetStressToBeginningOfStep( domain );
    dtReturnTemporary = solidSolver.NonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

    std::cout << dtReturnTemporary << std::endl;

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( solidSolver.getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter << " iterations! *****\n" );
      isConverged = true;
      break;
    }
    else if( m_subcyclingOption == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
      isConverged = true;
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", DamageSolver: " );

    dtReturnTemporary = damageSolver.NonlinearImplicitStep( time_n,
                                                            dtReturn,
                                                            cycleNumber,
                                                            domain );

    mapDamageToQuadrature( domain );

    std::cout << dtReturnTemporary << std::endl;

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }
    ++iter;
  }

  GEOSX_ERROR_IF( !isConverged, "PhaseFieldFractureSolver::SplitOperatorStep() did not converge" );

  damageSolver.ImplicitStepComplete( time_n, dt, domain );
  solidSolver.ImplicitStepComplete( time_n, dt, domain );
  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

void PhaseFieldFractureSolver::mapDamageToQuadrature( DomainPartition & domain )
{

  MeshLevel * const mesh = domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = mesh->getNodeManager();

  SolidMechanicsLagrangianFEM &
  solidSolver = *( this->getParent()->GetGroup( m_solidSolverName )->group_cast< SolidMechanicsLagrangianFEM * >() );

  PhaseFieldDamageFEM const &
  damageSolver = *( this->getParent()->GetGroup( m_damageSolverName )->group_cast< PhaseFieldDamageFEM * >() );

  string const & damageFieldName = damageSolver.getFieldName();

  //should get reference to damage field here.
  arrayView1d< real64 const > const nodalDamage = nodeManager->getReference< array1d< real64 > >( damageFieldName );

  ElementRegionManager * const elemManager = mesh->getElemManager();

  ConstitutiveManager * const constitutiveManager = domain.GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );

  ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase >
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor< ConstitutiveBase >( constitutiveManager );


  // begin region loop
  forTargetSubRegionsComplete< CellElementSubRegion >( *mesh, [this, &solidSolver, nodalDamage]
                                                         ( localIndex const targetIndex, localIndex, localIndex, ElementRegionBase &, CellElementSubRegion & elementSubRegion )
  {
    constitutive::ConstitutiveBase * const
    solidModel = elementSubRegion.getConstitutiveModel< constitutive::ConstitutiveBase >( solidSolver.solidMaterialNames()[targetIndex] );

    ConstitutivePassThru< DamageBase >::Execute( solidModel, [this, &elementSubRegion, nodalDamage]( auto * const damageModel )
    {
      using CONSTITUTIVE_TYPE = TYPEOFPTR( damageModel );
      typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel->createKernelUpdates();

      arrayView2d< real64 > const damageFieldOnMaterial = constitutiveUpdate.m_damage;
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemNodes = elementSubRegion.nodeList();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( m_discretizationName );

      finiteElement::dispatch3D( fe, [nodalDamage, &elementSubRegion, damageFieldOnMaterial, elemNodes]( auto & finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );
        constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
        constexpr localIndex n_q_points = FE_TYPE::numQuadraturePoints;

        forAll< serialPolicy >( elementSubRegion.size(), [nodalDamage, damageFieldOnMaterial, elemNodes] ( localIndex const k )
        {
          for( localIndex q = 0; q < n_q_points; ++q )
          {
            real64 N[ numNodesPerElement ];
            FE_TYPE::calcN( q, N );

            damageFieldOnMaterial( k, q ) = 0;
            for( localIndex a = 0; a < numNodesPerElement; ++a )
            {
              damageFieldOnMaterial( k, q ) += N[a] * nodalDamage[elemNodes( k, a )];
              //solution is probably not going to work because the solution of the coupled solver
              //has both damage and displacements. Using the damageResult field from the Damage solver
              //is probably better
              //            std::cout<<"q, N, Dnode = "<<q<<", "<<feDiscretization->m_finiteElement->value(a, q)<<",
              // "<<nodalDamage[elemNodes(k, a)]<<std::endl;
            }
            //          std::cout<<"damage("<<k<<","<<q<<") = "<<damageFieldOnMaterial(k,q)<<std::endl;
          }
        } );
      } );
    } );
  } );


}

REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldFractureSolver, std::string const &, Group * const )

} /* namespace geosx */
