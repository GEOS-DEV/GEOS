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

PhaseFieldFractureSolver::PhaseFieldFractureSolver( const std::string& name,
                                      Group * const parent ):
  SolverBase(name,parent),
  m_solidSolverName(),
  m_damageSolverName(),
  m_couplingTypeOptionString("FixedStress"),
  m_couplingTypeOption()

{
  registerWrapper(viewKeyStruct::solidSolverNameString, &m_solidSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the solid mechanics solver to use in the PhaseFieldFracture solver");

  registerWrapper(viewKeyStruct::damageSolverNameString, &m_damageSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the damage mechanics solver to use in the PhaseFieldFracture solver");

  registerWrapper(viewKeyStruct::couplingTypeOptionStringString, &m_couplingTypeOptionString, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Coupling option: (FixedStress, TightlyCoupled)");

}

void PhaseFieldFractureSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    ElementRegionManager * const elemManager = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0)->getElemManager();


    elemManager->forElementSubRegions<CellElementSubRegion,
                                      FaceElementSubRegion>( [&]( auto * const elementSubRegion ) -> void
      {
        elementSubRegion->template registerWrapper< array1d<real64> >( viewKeyStruct::totalMeanStressString )->
          setDescription("Total Mean Stress");
        elementSubRegion->template registerWrapper< array1d<real64> >( viewKeyStruct::oldTotalMeanStressString )->
          setDescription("Total Mean Stress");
      });
  }
}

void PhaseFieldFractureSolver::ImplicitStepSetup( real64 const & GEOSX_UNUSED_ARG( time_n ),
                                           real64 const & GEOSX_UNUSED_ARG( dt ),
                                           DomainPartition * const domain,
                                           DofManager & GEOSX_UNUSED_ARG( dofManager ),
                                           ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                           ParallelVector & GEOSX_UNUSED_ARG( rhs ),
                                           ParallelVector & GEOSX_UNUSED_ARG( solution ) )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const totalMeanStress =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(viewKeyStruct::totalMeanStressString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> oldTotalMeanStress =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(viewKeyStruct::oldTotalMeanStressString);

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    oldTotalMeanStress[er][esr][k] = totalMeanStress[er][esr][k];
  });
}

void PhaseFieldFractureSolver::ImplicitStepComplete( real64 const& GEOSX_UNUSED_ARG( time_n ),
                                              real64 const& GEOSX_UNUSED_ARG( dt ),
                                              DomainPartition * const GEOSX_UNUSED_ARG( domain ) )
{
}

void PhaseFieldFractureSolver::PostProcessInput()
{
  string ctOption = this->getReference<string>(viewKeyStruct::couplingTypeOptionStringString);

  if( ctOption == "FixedStress" )
  {
    this->m_couplingTypeOption = couplingTypeOption::FixedStress;

    // For this coupled solver the minimum number of Newton Iter should be 0 for both flow and solid solver otherwise it will never converge.
    SolidMechanicsLagrangianFEM &
      solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolidMechanicsLagrangianFEM*>());
      integer & minNewtonIterSolid = solidSolver.getSystemSolverParameters()->minIterNewton();

    PhaseFieldDamageFEM &
      damageSolver = *(this->getParent()->GetGroup(m_damageSolverName)->group_cast<PhaseFieldDamageFEM*>());
    integer & minNewtonIterFluid = damageSolver.getSystemSolverParameters()->minIterNewton();

    minNewtonIterSolid = 0;
    minNewtonIterFluid = 0;
  }
  else if( ctOption == "TightlyCoupled" )
  {
    this->m_couplingTypeOption = couplingTypeOption::TightlyCoupled;
  }
  else
  {
    GEOS_ERROR("invalid coupling type option");
  }

}

void PhaseFieldFractureSolver::InitializePostInitialConditions_PreSubGroups(Group * const  )
{

}

PhaseFieldFractureSolver::~PhaseFieldFractureSolver()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldFractureSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const totalMeanStress =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(viewKeyStruct::totalMeanStressString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> oldTotalMeanStress =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(viewKeyStruct::oldTotalMeanStressString);

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    totalMeanStress[er][esr][k] = oldTotalMeanStress[er][esr][k];
  });
}

real64 PhaseFieldFractureSolver::SolverStep( real64 const & time_n,
                                      real64 const & dt,
                                      int const cycleNumber,
                                      DomainPartition * const domain )
{
  real64 dtReturn = dt;
  if( m_couplingTypeOption == couplingTypeOption::FixedStress )
  {
    dtReturn = SplitOperatorStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>() );
  }
  else if( m_couplingTypeOption == couplingTypeOption::TightlyCoupled )
  {
    GEOS_ERROR( "couplingTypeOption::FullyImplicit not yet implemented");
  }
  return dtReturn;
}



real64 PhaseFieldFractureSolver::SplitOperatorStep( real64 const& time_n,
                                             real64 const& dt,
                                             integer const cycleNumber,
                                             DomainPartition * const domain)
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsLagrangianFEM &
  solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolidMechanicsLagrangianFEM*>());

  PhaseFieldDamageFEM &
  damageSolver = *(this->getParent()->GetGroup(m_damageSolverName)->group_cast<PhaseFieldDamageFEM*>());

  damageSolver.SetupSystem( domain,
                           damageSolver.getDofManager(),
                           damageSolver.getSystemMatrix(),
                           damageSolver.getSystemRhs(),
                           damageSolver.getSystemSolution() );

  solidSolver.SetupSystem( domain,
                           solidSolver.getDofManager(),
                           solidSolver.getSystemMatrix(),
                           solidSolver.getSystemRhs(),
                           solidSolver.getSystemSolution() );

  damageSolver.ImplicitStepSetup( time_n, dt, domain,
                                 damageSolver.getDofManager(),
                                 damageSolver.getSystemMatrix(),
                                 damageSolver.getSystemRhs(),
                                 damageSolver.getSystemSolution() );

  solidSolver.ImplicitStepSetup( time_n, dt, domain,
                                 solidSolver.getDofManager(),
                                 solidSolver.getSystemMatrix(),
                                 solidSolver.getSystemRhs(),
                                 solidSolver.getSystemSolution() );

  this->ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  int iter = 0;
  while (iter < (*(this->getSystemSolverParameters())).maxIterNewton() )
  {
    if (iter == 0)
    {
      // reset the states of all slave solvers if any of them has been reset
      damageSolver.ResetStateToBeginningOfStep( domain );
      solidSolver.ResetStateToBeginningOfStep( domain );
      ResetStateToBeginningOfStep( domain );
    }


    GEOS_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", MechanicsSolver: " );

    solidSolver.ResetStressToBeginningOfStep(domain);
    dtReturnTemporary = solidSolver.NonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain,
                                                           solidSolver.getDofManager(),
                                                           solidSolver.getSystemMatrix(),
                                                           solidSolver.getSystemRhs(),
                                                           solidSolver.getSystemSolution() );

    solidSolver.updateStress( domain );

    if (dtReturnTemporary < dtReturn)
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if (solidSolver.getSystemSolverParameters()->numNewtonIterations() == 0 && iter > 0)
    {
      GEOS_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
      break;
    }


    GEOS_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", DamageSolver: " );

    dtReturnTemporary = damageSolver.LinearImplicitStep( time_n,
                                                         dtReturn,
                                                         cycleNumber,
                                                         domain,
                                                         damageSolver.getDofManager(),
                                                         damageSolver.getSystemMatrix(),
                                                         damageSolver.getSystemRhs(),
                                                         damageSolver.getSystemSolution() );




    if (dtReturnTemporary < dtReturn)
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }
    ++iter;
  }

  damageSolver.ImplicitStepComplete( time_n, dt, domain );
  solidSolver.ImplicitStepComplete( time_n, dt, domain );
  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

void PhaseFieldFractureSolver::ApplySystemSolution(DofManager const & GEOSX_UNUSED_ARG(dofManager),
                                               ParallelVector const & GEOSX_UNUSED_ARG(solution),
                                               real64 const GEOSX_UNUSED_ARG(scalingFactor),
                                               DomainPartition *const domain)
{

  MeshLevel *const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  NodeManager *const nodeManager = mesh->getNodeManager();

  PhaseFieldDamageFEM const &
  damageSolver = *(this->getParent()->GetGroup(m_damageSolverName)->group_cast<PhaseFieldDamageFEM*>());

  string const & damageFieldName = damageSolver.getFieldName();

  //should get reference to damage field here.
  arrayView1d<real64> const & nodalDamage = nodeManager->getReference<array1d<real64>>(damageFieldName);

  ElementRegionManager *const elemManager = mesh->getElemManager();
  ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  NumericalMethodsManager const * const
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const * const
  feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  // begin region loop
  for (localIndex er = 0; er < elemManager->numRegions(); ++er)
  {
    ElementRegionBase * const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const * const
    feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>( [&]( localIndex const esr,
                                                                         CellElementSubRegion const *const elementSubRegion)
    {
      localIndex m_solidMaterialFullIndex = 0;
      SolidBase * solidModel =  constitutiveRelations[er][esr][m_solidMaterialFullIndex]->group_cast<SolidBase*>();

      arrayView2d<real64> const & damageFieldOnMaterial = solidModel->getDamage();

      localIndex const numNodesPerElement = elementSubRegion->numNodesPerElement();

      arrayView2d<localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM> const &
      elemNodes = elementSubRegion->nodeList();

      globalIndex_array elemDofIndex(numNodesPerElement);

      localIndex const n_q_points = feDiscretization->m_finiteElement->n_quadrature_points();

      for (localIndex k = 0; k < elementSubRegion->size(); ++k)
      {
        for (localIndex q = 0; q < n_q_points; ++q)
        {
          damageFieldOnMaterial(k,q) = 0;
          for (localIndex a = 0; a < numNodesPerElement; ++a)
          {
            damageFieldOnMaterial(k,q) +=
                //solution is probably not going to work because the solution of the coupled solver
                //has both damage and displacements. Using the damageResult field from the Damage solver
                //is probably better
                feDiscretization->m_finiteElement->value(a, q) * nodalDamage[elemNodes(k, a)];
          }
        }
      }
    });
  }

}



REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldFractureSolver, std::string const &, Group * const )

} /* namespace geosx */
