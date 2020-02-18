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
 * @file PoroelasticSolver.cpp
 *
 */


#include "PoroelasticSolver.hpp"

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

PoroelasticSolver::PoroelasticSolver( const std::string& name,
                                      Group * const parent ):
  SolverBase(name,parent),
  m_solidSolverName(),
  m_flowSolverName(),
  m_couplingTypeOptionString("FixedStress"),
  m_couplingTypeOption()

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

}

void PoroelasticSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
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

void PoroelasticSolver::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const & GEOSX_UNUSED_PARAM( dt ),
                                           DomainPartition * const domain,
                                           DofManager & GEOSX_UNUSED_PARAM( dofManager ),
                                           ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                           ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                                           ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
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

void PoroelasticSolver::ImplicitStepComplete( real64 const& GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const& GEOSX_UNUSED_PARAM( dt ),
                                              DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{
}

void PoroelasticSolver::PostProcessInput()
{
  string ctOption = this->getReference<string>(viewKeyStruct::couplingTypeOptionStringString);

  if( ctOption == "FixedStress" )
  {
    this->m_couplingTypeOption = couplingTypeOption::FixedStress;

    // For this coupled solver the minimum number of Newton Iter should be 0 for both flow and solid solver otherwise it will never converge.
    SolidMechanicsLagrangianFEM &
      solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolidMechanicsLagrangianFEM*>());
      integer & minNewtonIterSolid = solidSolver.getNonlinearSolverParameters().m_minIterNewton;


    SinglePhaseBase &
      fluidSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseBase*>());
    integer & minNewtonIterFluid = fluidSolver.getNonlinearSolverParameters().m_minIterNewton;

    minNewtonIterSolid = 0;
    minNewtonIterFluid = 0;
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

void PoroelasticSolver::InitializePostInitialConditions_PreSubGroups(Group * const problemManager)
{
  this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseBase*>()->setPoroElasticCoupling();
  // Calculate initial total mean stress
  this->UpdateDeformationForCoupling(problemManager->GetGroup<DomainPartition>(keys::domain));
}

PoroelasticSolver::~PoroelasticSolver()
{
  // TODO Auto-generated destructor stub
}

void PoroelasticSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
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

real64 PoroelasticSolver::SolverStep( real64 const & time_n,
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
    GEOSX_ERROR( "couplingTypeOption::FullyImplicit not yet implemented");
  }
  return dtReturn;
}

void PoroelasticSolver::UpdateDeformationForCoupling( DomainPartition * const domain )
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  
  ElementRegionManager * const elemManager = mesh->getElemManager();
  
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteElementDiscretizationManager const * const feDiscretizationManager =
    numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & X = nodeManager->referencePosition();
  arrayView2d<real64 const, nodes::TOTAL_DISPLACEMENT_USD> const & u = nodeManager->totalDisplacement();


  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegionBase const * const elemRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellElementSubRegion const * const cellElementSubRegion = elemRegion->GetSubRegion<CellElementSubRegion>(esr);

      arrayView2d<localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = cellElementSubRegion->nodeList();

      arrayView1d<real64> const &
      totalMeanStress = cellElementSubRegion->getReference< array1d<real64> >( viewKeyStruct::totalMeanStressString );

      arrayView1d<real64> const &
      oldTotalMeanStress = cellElementSubRegion->getReference< array1d<real64> >( viewKeyStruct::oldTotalMeanStressString );

      arrayView1d<real64 const> const &
      pres = cellElementSubRegion->getReference< array1d<real64> >( FlowSolverBase::viewKeyStruct::pressureString );

      arrayView1d<real64 const > const &
      dPres = cellElementSubRegion->getReference< array1d<real64> >( FlowSolverBase::viewKeyStruct::deltaPressureString );

      arrayView1d<real64 > const &
      poro = cellElementSubRegion->getReference< array1d<real64> >( SinglePhaseBase::viewKeyStruct::porosityString );

      arrayView1d<real64 const > const &
      poroOld = cellElementSubRegion->getReference< array1d<real64> >(SinglePhaseBase::viewKeyStruct::porosityOldString);

      arrayView1d<real64 const > const &
      volume = cellElementSubRegion->getReference< array1d<real64> >(CellBlock::viewKeyStruct::elementVolumeString);

      arrayView1d<real64 > const &
      dVol = cellElementSubRegion->getReference< array1d<real64> >(SinglePhaseBase::viewKeyStruct::deltaVolumeString);

      string const solidModelName = this->getSolidSolver()->getSolidMaterialName();

      arrayView1d<real64 const > const &
      bulkModulus = cellElementSubRegion->GetConstitutiveModels()->GetGroup(solidModelName)->getReference< array1d<real64> >( "BulkModulus" );

      real64 const
      biotCoefficient = cellElementSubRegion->GetConstitutiveModels()->GetGroup(solidModelName)->getReference<real64>( "BiotCoefficient");

      arrayView3d<real64 const, solid::STRESS_USD > const &
      stress = cellElementSubRegion->GetConstitutiveModels()->GetGroup(solidModelName)->
               getReference< array3d<real64,solid::STRESS_PERMUTATION> >( SolidBase::viewKeyStruct::stressString);


      localIndex const numNodesPerElement = elemsToNodes.size(1);
      localIndex const numQuadraturePoints = feDiscretization->m_finiteElement->n_quadrature_points() ;


      forall_in_range< parallelHostPolicy >( 0, cellElementSubRegion->size(),
                                             GEOSX_HOST_DEVICE_LAMBDA ( localIndex const ei )
      {

        R1Tensor u_local[10];

        for ( localIndex i = 0; i < numNodesPerElement; ++i )
        {
          localIndex const nodeIndex = elemsToNodes( ei, i );
          u_local[ i ] = u[ nodeIndex ];
        }

        real64 effectiveMeanStress = 0.0;
        for( localIndex q=0 ; q<numQuadraturePoints; ++q )
        {
          effectiveMeanStress += ( stress(ei,q,0) + stress(ei,q,2) + stress(ei,q,5) );
        }
        effectiveMeanStress /= ( 3 * numQuadraturePoints );

        totalMeanStress[ei] = effectiveMeanStress - biotCoefficient * (pres[ei] + dPres[ei]);

        poro[ei] = poroOld[ei] + (biotCoefficient - poroOld[ei]) / bulkModulus[ei]
                                                 * (totalMeanStress[ei] - oldTotalMeanStress[ei] + dPres[ei]);

        // update element volume
        R1Tensor Xlocal[ElementRegionManager::maxNumNodesPerElem];
        for (localIndex a = 0; a < elemsToNodes.size(1); ++a)
        {
          Xlocal[a] = X[elemsToNodes[ei][a]];
          Xlocal[a] += u[elemsToNodes[ei][a]] ;
        }

        dVol[ei] = computationalGeometry::HexVolume(Xlocal) - volume[ei];
      } );
    }
  }

}

real64 PoroelasticSolver::SplitOperatorStep( real64 const& time_n,
                                             real64 const& dt,
                                             integer const cycleNumber,
                                             DomainPartition * const domain)
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsLagrangianFEM &
  solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolidMechanicsLagrangianFEM*>());

  SinglePhaseBase &
  fluidSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseBase*>());

  fluidSolver.SetupSystem( domain,
                           fluidSolver.getDofManager(),
                           fluidSolver.getSystemMatrix(),
                           fluidSolver.getSystemRhs(),
                           fluidSolver.getSystemSolution() );

  solidSolver.SetupSystem( domain,
                           solidSolver.getDofManager(),
                           solidSolver.getSystemMatrix(),
                           solidSolver.getSystemRhs(),
                           solidSolver.getSystemSolution() );

  fluidSolver.ImplicitStepSetup( time_n, dt, domain,
                                 fluidSolver.getDofManager(),
                                 fluidSolver.getSystemMatrix(),
                                 fluidSolver.getSystemRhs(),
                                 fluidSolver.getSystemSolution() );

  solidSolver.ImplicitStepSetup( time_n, dt, domain,
                                 solidSolver.getDofManager(),
                                 solidSolver.getSystemMatrix(),
                                 solidSolver.getSystemRhs(),
                                 solidSolver.getSystemSolution() );

  this->ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  int iter = 0;
  while (iter < m_nonlinearSolverParameters.m_maxIterNewton )
  {
    if (iter == 0)
    {
      // reset the states of all slave solvers if any of them has been reset
      fluidSolver.ResetStateToBeginningOfStep( domain );
      solidSolver.ResetStateToBeginningOfStep( domain );
      ResetStateToBeginningOfStep( domain );
    }
    
    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );
    
    dtReturnTemporary = fluidSolver.NonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain,
                                                           fluidSolver.getDofManager(),
                                                           fluidSolver.getSystemMatrix(),
                                                           fluidSolver.getSystemRhs(),
                                                           fluidSolver.getSystemSolution() );

    if (dtReturnTemporary < dtReturn)
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if ( fluidSolver.getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0)
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", MechanicsSolver: " );

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
    if (solidSolver.getNonlinearSolverParameters().m_numNewtonIterations > 0)
    {
      this->UpdateDeformationForCoupling(domain);
    }
    ++iter;
  }

  fluidSolver.ImplicitStepComplete( time_n, dt, domain );
  solidSolver.ImplicitStepComplete( time_n, dt, domain );
  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}


REGISTER_CATALOG_ENTRY( SolverBase, PoroelasticSolver, std::string const &, Group * const )

} /* namespace geosx */
