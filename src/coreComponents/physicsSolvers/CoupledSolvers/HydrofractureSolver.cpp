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
 * @file HydrofractureSolver.cpp
 *
 */


#include "HydrofractureSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "../FiniteVolume/SinglePhaseFlow.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

HydrofractureSolver::HydrofractureSolver( const std::string& name,
                                      ManagedGroup * const parent ):
  SolverBase(name,parent),
  m_solidSolverName(),
  m_flowSolverName(),
  m_couplingTypeOptionString("FixedStress"),
  m_couplingTypeOption()

{
  RegisterViewWrapper(viewKeyStruct::solidSolverNameString, &m_solidSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the solid mechanics solver to use in the poroelastic solver");

  RegisterViewWrapper(viewKeyStruct::fluidSolverNameString, &m_flowSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the fluid mechanics solver to use in the poroelastic solver");

  RegisterViewWrapper(viewKeyStruct::couplingTypeOptionStringString, &m_couplingTypeOptionString, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Coupling option: (FixedStress, TightlyCoupled)");

}

void HydrofractureSolver::RegisterDataOnMesh( dataRepository::ManagedGroup * const MeshBodies )
{

}

void HydrofractureSolver::ImplicitStepSetup( real64 const& time_n,
                                             real64 const& dt,
                                             DomainPartition * const domain,
                                             systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  SolverBase & solidSolver =
    *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolverBase*>());

  SinglePhaseFlow & fluidSolver =
    *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseFlow*>());

  solidSolver.ImplicitStepSetup( time_n, dt, domain, blockSystem );
  fluidSolver.ImplicitStepSetup( time_n, dt, domain, blockSystem );
}

void HydrofractureSolver::ImplicitStepComplete( real64 const& time_n,
                                              real64 const& dt,
                                              DomainPartition * const domain)
{
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
    GEOS_ERROR("invalid coupling type option");
  }

}

void HydrofractureSolver::InitializePostInitialConditions_PreSubGroups(ManagedGroup * const problemManager)
{

}

HydrofractureSolver::~HydrofractureSolver()
{
  // TODO Auto-generated destructor stub
}

void HydrofractureSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{

}

real64 HydrofractureSolver::SolverStep( real64 const & time_n,
                                      real64 const & dt,
                                      int const cycleNumber,
                                      DomainPartition * domain )
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

void HydrofractureSolver::UpdateDeformationForCoupling( DomainPartition * const domain )
{
  MeshLevel * const meshLevel = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();

  arrayView1d<R1Tensor> const & u = nodeManager->getReference< array1d<R1Tensor> >( keys::TotalDisplacement );
  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();
  arrayView1d<real64 const> const & faceArea = faceManager->faceArea();
  array1d< array1d<localIndex> > const & facesToNodes = faceManager->nodeList();

  elemManager->forElementRegions<FaceElementRegion>([&]( FaceElementRegion * const faceElemRegion )
  {
    faceElemRegion->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion * const subRegion )
    {
      arrayView1d<real64> const & aperture = subRegion->getElementAperture();
      arrayView1d<real64> const & volume = subRegion->getElementVolume();
      arrayView1d<real64> const & area = subRegion->getElementArea();
      array1d< array1d<localIndex> > const & elemsToNodes = subRegion->nodeList();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();

      for( localIndex kfe=0 ; kfe<subRegion->size() ; ++kfe )
      {
        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const kf1 = elemsToFaces[kfe][1];
        localIndex const numNodesPerFace=facesToNodes[kf0].size();
        localIndex const * const nodelist0 = facesToNodes[kf0];
        localIndex const * const nodelist1 = facesToNodes[kf1];
        R1Tensor temp;
        for( localIndex a=0 ; a<numNodesPerFace ; ++a )
        {
          temp += u[nodelist0[a]];
          temp -= u[nodelist1[a]];
        }
        area[kfe] = faceArea[kfe];
        // TODO this needs a proper contact based strategy for aperture
        aperture[kfe] = -Dot(temp,faceNormal[kf0]) / numNodesPerFace+0.0001;
        volume[kfe] = aperture[kfe] * area[kfe];
        //std::cout<<"kfe, area, aperture, volume = "<<kfe<<", "<<area[kfe]<<", "<<aperture[kfe]<<", "<<volume[kfe]<<std::endl;
      }

    });
  });

}

real64 HydrofractureSolver::SplitOperatorStep( real64 const& time_n,
                                             real64 const& dt,
                                             integer const cycleNumber,
                                             DomainPartition * const domain)
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary = dtReturn;

  SolverBase &
  solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolverBase*>());

  SinglePhaseFlow &
  fluidSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseFlow*>());

  fluidSolver.ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
  solidSolver.ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
  this->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );

  this->UpdateDeformationForCoupling(domain);

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
                                                          getLinearSystemRepository() );

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
                                                          getLinearSystemRepository() );
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

real64 HydrofractureSolver::ExplicitStep( real64 const& time_n,
                                          real64 const& dt,
                                          const int cycleNumber,
                                          DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  SolverBase & solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolverBase*>());
  SinglePhaseFlow & fluidSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseFlow*>());

  solidSolver.ExplicitStep( time_n, dt, cycleNumber, domain );
  fluidSolver.SolverStep( time_n, dt, cycleNumber, domain );

  return dt;
}


REGISTER_CATALOG_ENTRY( SolverBase, HydrofractureSolver, std::string const &, ManagedGroup * const )

} /* namespace geosx */
