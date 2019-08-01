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
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "../solidMechanics/SolidMechanicsLagrangianFEM.hpp"

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
  MeshLevel * const meshLevel = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementRegions<FaceElementRegion>([&]( FaceElementRegion * const faceElemRegion )
  {
    faceElemRegion->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion * const subRegion )
    {
      arrayView1d<real64> const & volume = subRegion->getElementVolume();
      arrayView1d<real64> const & deltaVolume = subRegion->getReference<array1d<real64> >(SinglePhaseFlow::viewKeyStruct::deltaVolumeString);

      for( localIndex kfe=0 ; kfe<subRegion->size() ; ++kfe )
      {
        volume[kfe] += deltaVolume[kfe];
      }

    });
  });

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
      arrayView1d<real64> const & deltaFluidPressure = subRegion->getReference<array1d<real64> >("deltaPressure");
      arrayView1d<real64> const & aperture = subRegion->getElementAperture();
      arrayView1d<real64> const & volume = subRegion->getElementVolume();
      arrayView1d<real64> const & deltaVolume = subRegion->getReference<array1d<real64> >(SinglePhaseFlow::viewKeyStruct::deltaVolumeString);
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

        real64 const oldAperture = aperture[kfe];


        // TODO this needs a proper contact based strategy for aperture
        aperture[kfe] = -Dot(temp,faceNormal[kf0]) / numNodesPerFace + 1.0e-4;

        real64 const K = 2.0e9;
        real64 const dP_dAper = -K / oldAperture;//aperture[kfe];

        std::cout<<"kfe, aperture, oldAperture, deltaP = "<<kfe<<", "<<aperture[kfe]<<", "<<oldAperture<<", "<<deltaFluidPressure[kfe]<<std::endl;
        std::cout<<"P correction option 1:  K * ( 1 + aperture[kfe] / oldAperture ) = "<<dP_dAper * ( aperture[kfe] - oldAperture )<<std::endl;
        std::cout<<"P correction option 2:  K * (  oldAperture / aperture[kfe] - 1 ) = "<<K * ( oldAperture / aperture[kfe] - 1 )<<std::endl;
//        deltaFluidPressure[kfe] += dP_dAper * ( aperture[kfe] - oldAperture );
        std::cout<<"    new deltaP = "<<deltaFluidPressure[kfe]<<std::endl;

        deltaVolume[kfe] = aperture[kfe] * area[kfe] - volume[kfe];
      }

    });
  });

}


void HydrofractureSolver::ApplyFractureFluidCoupling( DomainPartition * const domain,
                                                      systemSolverInterface::EpetraBlockSystem & blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  FaceManager const * const faceManager = mesh->getFaceManager();
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  arrayView1d<R1Tensor> const & u = nodeManager->getReference< array1d<R1Tensor> >( keys::TotalDisplacement );

  arrayView1d<real64 const>   const & faceArea   = faceManager->faceArea();
  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();
  array1d<localIndex_array> const & facesToNodes = faceManager->nodeList();
  arrayView1d<R1Tensor> const & fext = nodeManager->getReference< array1d<R1Tensor> >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternal );
  fext = {0,0,0};



  arrayView1d<globalIndex> const &
  blockLocalDofNumber =  nodeManager->getReference<globalIndex_array>(SolidMechanicsLagrangianFEM::viewKeyStruct::globalDofNumberString);
  Epetra_FEVector * const rhs = blockSystem.GetResidualVector( systemSolverInterface::BlockIDs::displacementBlock );
  Epetra_FECrsMatrix * const matrix = blockSystem.GetMatrix( systemSolverInterface::BlockIDs::displacementBlock,
                                                              systemSolverInterface::BlockIDs::displacementBlock );

  elemManager->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasView("pressure") )
    {

      arrayView1d<real64> const & aperture = subRegion->getElementAperture();

      arrayView1d<real64 const> const & fluidPressure = subRegion->getReference<array1d<real64> >("pressure");
      arrayView1d<real64 const> const & deltaFluidPressure = subRegion->getReference<array1d<real64> >("deltaPressure");

      arrayView1d<real64> const & volume = subRegion->getElementVolume();
      arrayView1d<real64> const & deltaVolume = subRegion->getReference<array1d<real64> >(SinglePhaseFlow::viewKeyStruct::deltaVolumeString);
      arrayView1d<real64> const & area = subRegion->getElementArea();
      array1d< array1d<localIndex> > const & elemsToNodes = subRegion->nodeList();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();


      forall_in_range<serialPolicy>( 0,
                                   subRegion->size(),
                                   GEOSX_LAMBDA ( localIndex const kfe )
      {

        R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
        Nbar -= faceNormal[elemsToFaces[kfe][1]];
        Nbar.Normalize();

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
        area[kfe] = 0.5 * ( faceArea[elemsToFaces[kfe][0]] + faceArea[elemsToFaces[kfe][1]] );
        // TODO this needs a proper contact based strategy for aperture
        aperture[kfe] = -Dot(temp,Nbar) / numNodesPerFace + 1.0e-4;
        deltaVolume[kfe] = aperture[kfe] * area[kfe] - volume[kfe];
        std::cout<<"kfe, area, aperture, volume, dVolume = "<<kfe<<", "<<area[kfe]<<", "<<aperture[kfe]<<", "<<volume[kfe]<<", "<<deltaVolume[kfe]<<std::endl;




        globalIndex rowDOF[24];
        globalIndex colDOF[24];
        real64 nodeRHS[24];
        stackArray2d<real64, 12*12> dRdU(numNodesPerFace*3, numNodesPerFace*3);

        real64 nodalForceMag = ( fluidPressure[kfe]+deltaFluidPressure[kfe] ) * faceArea[elemsToFaces[kfe][0]] / numNodesPerFace;
        std::cout<<"fluidPressure[kfe]+deltaFluidPressure[kfe] = "<<fluidPressure[kfe]<<" + "<<deltaFluidPressure[kfe]<<std::endl;
        R1Tensor nodalForce(Nbar);
        nodalForce *= nodalForceMag;
        for( localIndex kf=0 ; kf<2 ; ++kf )
        {
          localIndex const faceIndex = elemsToFaces[kfe][kf];

          for( localIndex a=0 ; a<numNodesPerFace ; ++a )
          {
            for( int component=0 ; component<3 ; ++component )
            {
              rowDOF[3*a+component] = 3*blockLocalDofNumber[facesToNodes[faceIndex][a]]+component;
              nodeRHS[3*a+component] = - nodalForce[component] * pow(-1,kf);
              fext[facesToNodes[faceIndex][a]][component] += - nodalForce[component] * pow(-1,kf);
            }
          }

          rhs->SumIntoGlobalValues( integer_conversion<int>(numNodesPerFace*3), rowDOF, nodeRHS );
        }


        real64 dAper_dU[2][4][3];
        for( localIndex kf=0 ; kf<2 ; ++kf )
        {
          for( localIndex a=0 ; a<numNodesPerFace ; ++a )
          {
            for( int component=0 ; component<3 ; ++component )
            {
              dAper_dU[kf][a][component] = - pow(-1,kf) * Nbar[component] / numNodesPerFace;
            }
          }
        }

        real64 const K = 2e9;
        real64 const dP_dAper = -K / aperture[kfe];

        for( localIndex kfa=0 ; kfa<2 ; ++kfa )
        {
          localIndex const faceIndexA = elemsToFaces[kfe][kfa];
          for( localIndex kfb=0 ; kfb<2 ; ++kfb )
          {
            localIndex const faceIndexb = elemsToFaces[kfe][kfb];
            for( int a=0 ; a<numNodesPerFace ; ++a )
            {
              for( int i=0 ; i<3 ; ++i )
              {
                rowDOF[3*a+i] = 3*blockLocalDofNumber[facesToNodes[faceIndexA][a]]+i;
                real64 const integrationFactor = - faceArea[elemsToNodes[kfe][0]] / numNodesPerFace * Nbar[i] * pow(-1,kfa);

                for( int b=0 ; b<numNodesPerFace ; ++b )
                {
                  for( int j=0 ; j<3 ; ++j )
                  {
                    colDOF[3*b+j] = 3*blockLocalDofNumber[facesToNodes[faceIndexb][b]]+j;
                    dRdU(3*a+i,3*b+j) = integrationFactor * dP_dAper * dAper_dU[kfb][b][j];
                  }
                }
              }
            }
            matrix->SumIntoGlobalValues( integer_conversion<int>( numNodesPerFace * 3),
                                         rowDOF,
                                         integer_conversion<int>( numNodesPerFace * 3),
                                         colDOF,
                                         dRdU.data() );

          }
        }


      });
    }
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


  SystemSolverParameters const * const solverParams = getSystemSolverParameters();

  this->UpdateDeformationForCoupling(domain);

  int iter = 0;
  while (iter < solverParams->maxIterNewton() )
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
//    dtReturnTemporary = fluidSolver.NonlinearImplicitStep( time_n,
//                                                          dtReturn,
//                                                          cycleNumber,
//                                                          domain,
//                                                          getLinearSystemRepository() );

    // call assemble to fill the matrix and the rhs
    fluidSolver.AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );

    // apply boundary conditions to system
    fluidSolver.ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );

    // call the default linear solver on the system
    fluidSolver.SolveSystem( getLinearSystemRepository(),
                 getSystemSolverParameters() );

    // apply the system solution to the fields/variables
    fluidSolver.ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );


    if (dtReturnTemporary < dtReturn)
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

//    if (fluidSolver.getSystemSolverParameters()->numNewtonIterations() == 0 && iter > 0 && this->verboseLevel() >= 1)
//    {
//      GEOS_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
//      break;
//    }

    if (this->verboseLevel() >= 1)
    {
      GEOS_LOG_RANK_0( "\tIteration: " << iter+1  << ", MechanicsSolver: " );
    }
//    dtReturnTemporary = solidSolver.NonlinearImplicitStep( time_n,
//                                                          dtReturn,
//                                                          cycleNumber,
//                                                          domain,
//                                                          getLinearSystemRepository() );

    // call assemble to fill the matrix and the rhs
    solidSolver.AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );


    ApplyFractureFluidCoupling( domain, *getLinearSystemRepository() );

    // apply boundary conditions to system
    solidSolver.ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );

    // call the default linear solver on the system
    solidSolver.SolveSystem( getLinearSystemRepository(),
                 getSystemSolverParameters() );

    // apply the system solution to the fields/variables
    solidSolver.ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );

    if( fluidSolver.CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() &&
        solidSolver.CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() )
    {
      GEOS_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
      break;
    }

    if (dtReturnTemporary < dtReturn)
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }
//    if (solidSolver.getSystemSolverParameters()->numNewtonIterations() > 0)
    {
      this->UpdateDeformationForCoupling(domain);
//      fluidSolver.UpdateState(domain);
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
