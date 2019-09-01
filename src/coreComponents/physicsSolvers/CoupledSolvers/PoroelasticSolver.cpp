/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file PoroelasticSolver.cpp
 *
 */


#include "PoroelasticSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "../FiniteVolume/SinglePhaseFlow.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"


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

void PoroelasticSolver::ImplicitStepSetup( real64 const & time_n,
                                           real64 const & dt,
                                           DomainPartition * const domain,
                                           DofManager & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution )
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

void PoroelasticSolver::ImplicitStepComplete( real64 const& time_n,
                                              real64 const& dt,
                                              DomainPartition * const domain)
{
}

void PoroelasticSolver::PostProcessInput()
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
    GEOS_ERROR("invalid coupling type option");
  }

}

void PoroelasticSolver::InitializePostInitialConditions_PreSubGroups(Group * const problemManager)
{
  this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseFlow*>()->setPoroElasticCoupling();
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
    GEOS_ERROR( "couplingTypeOption::FullyImplicit not yet implemented");
  }
  return dtReturn;
}

void PoroelasticSolver::UpdateDeformationForCoupling( DomainPartition * const domain )
{
//  SolverBase & solidSolver =
//    *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolverBase*>());

  SinglePhaseFlow & fluidSolver = 
    *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseFlow*>());

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  
  ElementRegionManager * const elemManager = mesh->getElemManager();
  
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteElementDiscretizationManager const * const feDiscretizationManager =
    numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  ConstitutiveManager * const constitutiveManager =
    domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  arrayView1d<R1Tensor> const & X = nodeManager->getReference<r1_array>(nodeManager->viewKeys.referencePosition);
  arrayView1d<R1Tensor> const & u = nodeManager->getReference<r1_array>(keys::TotalDisplacement);
  arrayView1d<R1Tensor> const & uhat = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);

  ElementRegionManager::ElementViewAccessor<arrayView2d<localIndex>> const elemsToNodes = 
    elemManager->ConstructViewAccessor<FixedOneToManyRelation, arrayView2d<localIndex>>( CellElementSubRegion::viewKeyStruct::nodeListString );

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> totalMeanStress =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(viewKeyStruct::totalMeanStressString);
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const oldTotalMeanStress =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(viewKeyStruct::oldTotalMeanStressString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const pres =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(FlowSolverBase::viewKeyStruct::pressureString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const dPres =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(FlowSolverBase::viewKeyStruct::deltaPressureString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> poro =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(SinglePhaseFlow::viewKeyStruct::porosityString);
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const poroOld =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(SinglePhaseFlow::viewKeyStruct::porosityOldString);
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const volume =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(CellBlock::viewKeyStruct::elementVolumeString);
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> dVol =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>(SinglePhaseFlow::viewKeyStruct::deltaVolumeString);

  ElementRegionManager::MaterialViewAccessor< arrayView1d<real64> > const bulkModulus =
    elemManager->ConstructFullMaterialViewAccessor< array1d<real64>, arrayView1d<real64> >( "BulkModulus", constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const pvmult =
    elemManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                        constitutiveManager );
  ElementRegionManager::MaterialViewAccessor<real64> const biotCoefficient =
    elemManager->ConstructFullMaterialViewAccessor<real64>( "BiotCoefficient", constitutiveManager);

  localIndex const solidIndex = domain->getConstitutiveManager()->GetConstitutiveRelation( fluidSolver.solidIndex() )->getIndexInParent();

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegionBase const * const elemRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellElementSubRegion const * const cellElementSubRegion = elemRegion->GetSubRegion<CellElementSubRegion>(esr);

      arrayView3d<R1Tensor> const & dNdX = cellElementSubRegion->getReference< array3d<R1Tensor> >(keys::dNdX);

      localIndex const numNodesPerElement = elemsToNodes[er][esr].size(1);
      r1_array u_local( numNodesPerElement );
      r1_array uhat_local( numNodesPerElement );

      for( localIndex ei=0 ; ei<cellElementSubRegion->size() ; ++ei )
      {
        CopyGlobalToLocal<R1Tensor>( elemsToNodes[er][esr][ei], u, uhat, u_local, uhat_local, numNodesPerElement );

        real64 volumetricStrain = 0.0;
        localIndex const numQuadraturePoints = feDiscretization->m_finiteElement->n_quadrature_points() ;
        for( localIndex q=0 ; q<numQuadraturePoints; ++q )
        {
          R2Tensor dUdX;
          CalculateGradient( dUdX, u_local, dNdX[ei][q], numNodesPerElement );
          volumetricStrain += dUdX.Trace();
        }
        volumetricStrain /= numQuadraturePoints;
        totalMeanStress[er][esr][ei] = volumetricStrain * bulkModulus[er][esr][solidIndex][ei] - biotCoefficient[er][esr][solidIndex] * (pres[er][esr][ei] + dPres[er][esr][ei]);

        poro[er][esr][ei] = poroOld[er][esr][ei] + (biotCoefficient[er][esr][solidIndex] - poroOld[er][esr][ei]) / bulkModulus[er][esr][solidIndex][ei]
                                                 * (totalMeanStress[er][esr][ei] - oldTotalMeanStress[er][esr][ei] + dPres[er][esr][ei]);

        // update element volume
        R1Tensor Xlocal[ElementRegionManager::maxNumNodesPerElem];
        for (localIndex a = 0; a < elemsToNodes[er][esr].size(1); ++a)
        {
          Xlocal[a] = X[elemsToNodes[er][esr][ei][a]];
          Xlocal[a] += u[elemsToNodes[er][esr][ei][a]] ;
        }

        dVol[er][esr][ei] = computationalGeometry::HexVolume(Xlocal) - volume[er][esr][ei];
      }
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

  SolverBase &
  solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolverBase*>());

  SinglePhaseFlow &
  fluidSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseFlow*>());

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
  while (iter < (*(this->getSystemSolverParameters())).maxIterNewton() )
  {
    if (iter == 0)
    {
      // reset the states of all slave solvers if any of them has been reset
      fluidSolver.ResetStateToBeginningOfStep( domain );
      solidSolver.ResetStateToBeginningOfStep( domain );
      ResetStateToBeginningOfStep( domain );
    }
    if (this->verboseLevel() >= 1)
    {
      GEOS_LOG_RANK_0( "\tIteration: " << iter+1  << ", FlowSolver: " );
    }
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

    if (fluidSolver.getSystemSolverParameters()->numNewtonIterations() == 0 && iter > 0 && this->verboseLevel() >= 1)
    {
      GEOS_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
      break;
    }

    if (this->verboseLevel() >= 1)
    {
      GEOS_LOG_RANK_0( "\tIteration: " << iter+1  << ", MechanicsSolver: " );
    }
    dtReturnTemporary = solidSolver.NonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain,
                                                           solidSolver.getDofManager(),
                                                           solidSolver.getSystemMatrix(),
                                                           solidSolver.getSystemRhs(),
                                                           solidSolver.getSystemSolution() );
    if (dtReturnTemporary < dtReturn)
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }
    if (solidSolver.getSystemSolverParameters()->numNewtonIterations() > 0)
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
