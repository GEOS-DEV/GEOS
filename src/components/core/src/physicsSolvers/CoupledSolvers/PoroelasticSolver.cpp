/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
//#include "../../../../PhysicsSolverPackage1/src/SolidMechanicsLagrangianFEM.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

PoroelasticSolver::PoroelasticSolver( const std::string& name,
                                      ManagedGroup * const parent ):
  SolverBase(name,parent)
{
  this->RegisterViewWrapper(viewKeyStruct::solidSolverNameString, &m_solidSolverName, 0);
  this->RegisterViewWrapper(viewKeyStruct::fluidSolverNameString, &m_flowSolverName, 0);
}


void PoroelasticSolver::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  SolverBase::FillDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("An example single phase flow solver");


  docNode->AllocateChildNode( viewKeyStruct::fluidSolverNameString,
                              viewKeyStruct::fluidSolverNameString,
                              -1,
                              "string",
                              "string",
                              "name of fluid solver",
                              "name of fluid solver",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::solidSolverNameString,
                              viewKeyStruct::solidSolverNameString,
                              -1,
                              "string",
                              "string",
                              "name of solid solver",
                              "name of solid solver",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::couplingTypeOptionString,
                              viewKeyStruct::couplingTypeOptionString,
                              -1,
                              "string",
                              "string",
                              "option for default coupling method",
                              "option for default coupling method",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::biotCoefficienttring,
                              viewKeyStruct::biotCoefficienttring,
                              -1,
                              "real64",
                              "real64",
                              "Biot's Coefficient",
                              "Biot's Coefficient",
                              "",
                              "",
                              0,
                              1,
                              0 );
}

void PoroelasticSolver::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup )
{

  SolverBase::FillOtherDocumentationNodes( rootGroup );
  DomainPartition * domain  = rootGroup->GetGroup<DomainPartition>(keys::domain);

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody*>(mesh.second)->getMeshLevel(0);

    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks( [&]( CellBlockSubRegion * const cellBlock ) -> void
      {
        cxx_utilities::DocumentationNode * const docNode = cellBlock->getDocumentationNode();

        docNode->AllocateChildNode( viewKeyStruct::effectiveStressString,
                                    viewKeyStruct::effectiveStressString,
                                    -1,
                                    "r2Sym_array",
                                    "r2Sym_array",
                                    "Effective Stress",
                                    "Effective Stress",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeyStruct::deltaEffectiveStressString,
                                    viewKeyStruct::deltaEffectiveStressString,
                                    -1,
                                    "r2Sym_array",
                                    "r2Sym_array",
                                    "Change in Effective Stress",
                                    "Change in Effective Stress",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeyStruct::deltaVolumetricStrainString,
                                    viewKeyStruct::deltaVolumetricStrainString,
                                    -1,
                                    "r2Sym_array",
                                    "r2Sym_array",
                                    "Change in Volumetric Strain",
                                    "Change in Volumetric Strain",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );


        docNode->AllocateChildNode( viewKeyStruct::deltaPorosityString,
                                    viewKeyStruct::deltaPorosityString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Change in Porosity",
                                    "Change in Porosity",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeyStruct::dPorosity_dPressureString,
                                    viewKeyStruct::dPorosity_dPressureString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Derivative of Porosity wrt Pressure",
                                    "Derivative of Porosity wrt Pressure",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );

      });
  }
}

void PoroelasticSolver::ImplicitStepSetup( real64 const& time_n,
                                           real64 const& dt,
                                           DomainPartition * const domain,
                                           systemSolverInterface::EpetraBlockSystem * const blockSystem)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto dEffStresss = elemManager->ConstructViewAccessor<r2Sym_array>(viewKeyStruct::deltaEffectiveStressString);
  auto dVolumetricStrain = elemManager->ConstructViewAccessor<r2Sym_array>(viewKeyStruct::deltaVolumetricStrainString);

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    dEffStresss[er][esr][k] = 0.0;
    dVolumetricStrain[er][esr][k] = 0.0;
  });
}

void PoroelasticSolver::ImplicitStepComplete( real64 const& time_n,
                                              real64 const& dt,
                                              DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto effectiveStress = elemManager->ConstructViewAccessor<r2Sym_array>(viewKeyStruct::effectiveStressString);
  auto dEffStresss = elemManager->ConstructViewAccessor<r2Sym_array>(viewKeyStruct::deltaEffectiveStressString);

  //***** loop over all elements and update the derivative arrays *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const k)->void
  {
    effectiveStress[er][esr][k] += dEffStresss[er][esr][k];
  });
}

void PoroelasticSolver::ReadXML_PostProcess()
{
  string ctOption = this->getReference<string>(viewKeyStruct::couplingTypeOptionString);

  if( ctOption == "FixedStress" )
  {
    this->m_couplingTypeOption = couplingTypeOption::FixedStress;
  }
  else if( ctOption == "FullyImplicit" )
  {
    this->m_couplingTypeOption = couplingTypeOption::FullyImplicit;
  }
  else
  {
    GEOS_ERROR("invalid coupling type option");
  }
}

void PoroelasticSolver::FinalInitialization( ManagedGroup * const problemManager )
{
  m_biotCoef = this->getReference<integer>(viewKeyStruct::biotCoefficienttring);
}

PoroelasticSolver::~PoroelasticSolver()
{
  // TODO Auto-generated destructor stub
}

real64 PoroelasticSolver::SolverStep( real64 const & time_n,
                                      real64 const & dt,
                                      int const cycleNumber,
                                      DomainPartition * domain )
{
  real64 dtReturn = dt;
  if( m_couplingTypeOption == couplingTypeOption::FixedStress )
  {
    dtReturn = SplitOperatorStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>() );
  }
  else if( m_couplingTypeOption == couplingTypeOption::FullyImplicit )
  {
    dtReturn = SplitOperatorStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>() );
  }
  return dtReturn;
}

void PoroelasticSolver::UpdateDeformationForCoupling( DomainPartition * const domain )
{

  SolverBase &
  solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolverBase*>());

  SinglePhaseFlow &
  fluidSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseFlow*>());

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  NumericalMethodsManager const * const
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteElementSpaceManager const * const
  feSpaceManager = numericalMethodManager->GetGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);


  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  view_rtype<r1_array> X     = nodeManager->getData<r1_array>(nodeManager->viewKeys.referencePosition);
  view_rtype<r1_array> u     = nodeManager->getData<r1_array>(keys::TotalDisplacement);
  view_rtype<r1_array> uhat  = nodeManager->getData<r1_array>(keys::IncrementalDisplacement);

  auto elemsToNodes = elemManager->
                       ConstructViewAccessor<FixedOneToManyRelation>( CellBlockSubRegion::viewKeyStruct::nodeListString );

  auto effectiveStress = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::effectiveStressString);
  auto dEffStresss = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaEffectiveStressString);
  auto dVolumetricStrain = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaVolumetricStrainString);

  auto dPres = elemManager->ConstructViewAccessor<real64_array>( SinglePhaseFlow::
                                                                 viewKeyStruct::deltaPressureString);


  auto poro = elemManager->ConstructViewAccessor<real64_array>( SinglePhaseFlow::
                                                                viewKeyStruct::
                                                                porosityString);

  auto dPoro = elemManager->ConstructViewAccessor<real64_array>( viewKeyStruct::
                                                                 deltaPorosityString);

  auto volume    = elemManager->ConstructViewAccessor<real64_array>( SinglePhaseFlow::
                                                                     viewKeyStruct::
                                                                     volumeString);

  auto dVol      = elemManager->ConstructViewAccessor<real64_array>( SinglePhaseFlow::
                                                                     viewKeyStruct::
                                                                     deltaVolumeString);

  auto dPoro_dPres = elemManager->ConstructViewAccessor<real64_array>( viewKeyStruct::
                                                                       dPorosity_dPressureString);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  meanStress = elemManager->
               ConstructMaterialViewAccessor< array2d<real64> >( "meanStress",
                                                                 constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  bulkModulus = elemManager->
                ConstructMaterialViewAccessor< array2d<real64> >( "BulkModulus",
                                                                  constitutiveManager);

  localIndex const solidIndex = fluidSolver.solidIndex();

  // TODO
  //   dVolumetricStrain[er][esr][k] += dUhatdX.trace() / numQuadraturePoints;
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);

    auto const & numMethodName = elemRegion->getData<string>(keys::numericalMethod);
    FiniteElementSpace const * feSpace = feSpaceManager->GetGroup<FiniteElementSpace>(numMethodName);

    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const cellBlockSubRegion = elemRegion->GetSubRegion(esr);

      multidimensionalArray::ManagedArray<R1Tensor, 3> const &
      dNdX = cellBlockSubRegion->getReference< multidimensionalArray::ManagedArray<R1Tensor, 3> >(keys::dNdX);

      localIndex const numNodesPerElement = elemsToNodes[er][esr].get().size(1);
      r1_array uhat_local( numNodesPerElement );
      auto const &
      constitutiveMap = cellBlockSubRegion->getReference< std::pair< array2d<localIndex>,array2d<localIndex> > >(cellBlockSubRegion->viewKeys().constitutiveMap);

      for( localIndex ei=0 ; ei<cellBlockSubRegion->size() ; ++ei )
      {

        dVolumetricStrain[er][esr][ei] = 0.0;
        dEffStresss[er][esr][ei] = 0.0;
        localIndex const numQuadraturePoints = feSpace->m_finiteElement->n_quadrature_points() ;
        for( localIndex q=0 ; q<numQuadraturePoints; ++q )
        {
          R2Tensor dUhatdX;
          CalculateGradient( dUhatdX, uhat_local, dNdX[ei][q] );

          dVolumetricStrain[er][esr][ei] += dUhatdX.Trace();
          dEffStresss[er][esr][ei] += ( meanStress[er][esr][solidIndex][ei][0]
                                      - effectiveStress[er][esr][ei] ) ;
        }
        dVolumetricStrain[er][esr][ei] /= numQuadraturePoints;
        dEffStresss[er][esr][ei] /= numQuadraturePoints;

        dPoro[er][esr][ei] = (m_biotCoef - poro[er][esr][ei])
                          / bulkModulus[er][esr][solidIndex][ei][0]
                          * (dEffStresss[er][esr][ei] + (1 - m_biotCoef) * dPres[er][esr][ei]);

        dPoro_dPres[er][esr][ei] = (m_biotCoef - poro[er][esr][ei]) / bulkModulus[er][esr][solidIndex][ei][0];


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
  SolverBase &
  solidSolver = *(this->getParent()->GetGroup(m_solidSolverName)->group_cast<SolverBase*>());

  SinglePhaseFlow &
  fluidSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<SinglePhaseFlow*>());

  fluidSolver.ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
  solidSolver.ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
  this->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );

  int iter = 0;
  while (iter < (*(this->getSystemSolverParameters())).maxIterNewton() )
  {
    dtReturn = fluidSolver.NonlinearImplicitStep( time_n,
                                                  dtReturn,
                                                  cycleNumber,
                                                  domain,
                                                  getLinearSystemRepository() );

    if (fluidSolver.getSystemSolverParameters()->numNewtonIterations() == 0 && iter > 0)
    {
      break;
    }

    dtReturn = solidSolver.NonlinearImplicitStep( time_n,
                                                  dtReturn,
                                                  cycleNumber,
                                                  domain,
                                                  getLinearSystemRepository() );
    this->UpdateDeformationForCoupling(domain);
    ++iter;
  }

  fluidSolver.ImplicitStepComplete( time_n, dt, domain );
  solidSolver.ImplicitStepComplete( time_n, dt, domain );
  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}


REGISTER_CATALOG_ENTRY( SolverBase, PoroelasticSolver, std::string const &, ManagedGroup * const )

} /* namespace geosx */
