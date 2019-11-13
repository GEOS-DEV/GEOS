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
 * @file SolidMechanicsLagrangianSSLE.hpp
 */

#include "SolidMechanicsLagrangianSSLE.hpp"

#include "codingUtilities/Utilities.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/NumericalMethodsManager.hpp"



namespace geosx
{

using namespace constitutive;

SolidMechanicsLagrangianSSLE::SolidMechanicsLagrangianSSLE( string const & name,
                                                            Group * const parent ):
  SolidMechanicsLagrangianFEM( name, parent )
{
  this->m_strainTheory = 0;
}

SolidMechanicsLagrangianSSLE::~SolidMechanicsLagrangianSSLE()
{}

void SolidMechanicsLagrangianSSLE::ApplySystemSolution( DofManager const & dofManager,
                                                        ParallelVector const & solution,
                                                        real64 const scalingFactor,
                                                        DomainPartition * const domain  )
{

  SolidMechanicsLagrangianFEM::ApplySystemSolution( dofManager, solution, scalingFactor, domain );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  Group * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(dataRepository::keys::ConstitutiveManager);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(dataRepository::keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(dataRepository::keys::finiteElementDiscretizations);

  ElementRegionManager::MaterialViewAccessor<real64> const
  biotCoefficient = elemManager->ConstructFullMaterialViewAccessor<real64>( "BiotCoefficient", constitutiveManager);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const
  fluidPres = elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("pressure");

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const
  dPres = elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("deltaPressure");

  r1_array const& disp = nodeManager->getReference<r1_array>(dataRepository::keys::TotalDisplacement);


  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  // begin region loop
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegionBase * const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const *
    feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr,
                                                                        CellElementSubRegion const * const elementSubRegion )
    {
      array3d<R1Tensor> const &
      dNdX = elementSubRegion->getReference< array3d<R1Tensor> >(dataRepository::keys::dNdX);

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(dataRepository::keys::detJ);

      arrayView2d< localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM > const & elemsToNodes = elementSubRegion->nodeList();
      localIndex const numNodesPerElement = elemsToNodes.size(1);

      std::unique_ptr<FiniteElementBase>
      fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

      // space for element matrix and rhs

      using Kernels = SolidMechanicsLagrangianSSLEKernels::StressCalculationKernel;
      return SolidMechanicsLagrangianFEMKernels::
             ElementKernelLaunchSelector<Kernels>( numNodesPerElement,
                                                   fe->n_quadrature_points(),
                                                   constitutiveRelations[er][esr][m_solidMaterialFullIndex],
                                                   elementSubRegion->size(),
                                                   elemsToNodes,
                                                   dNdX,
                                                   detJ,
                                                   disp );

    });
  }


}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianSSLE, string const &, dataRepository::Group * const )
} /* namespace geosx */

