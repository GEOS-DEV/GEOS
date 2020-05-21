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
                                                        DomainPartition * const domain )
{
  SolidMechanicsLagrangianFEM::ApplySystemSolution( dofManager, solution, scalingFactor, domain );
}


void
SolidMechanicsLagrangianSSLE::updateStress( DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodeManager = *mesh.getNodeManager();

  NumericalMethodsManager const & numericalMethodManager = domain->getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const &
  feDiscretization = *feDiscretizationManager.GetGroup< FiniteElementDiscretization >( m_discretizationName );

  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & incDisp = nodeManager.incrementalDisplacement();

  // begin region loop
  forTargetSubRegions< CellElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                          CellElementSubRegion & elementSubRegion )
  {
    arrayView4d< real64 const > const & dNdX = elementSubRegion.dNdX();

    arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();

    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    localIndex const numNodesPerElement = elemsToNodes.size( 1 );

    std::unique_ptr< FiniteElementBase >
    fe = feDiscretization.getFiniteElement( elementSubRegion.GetElementTypeString() );

    SolidBase & constitutiveRelation = GetConstitutiveModel< SolidBase >( elementSubRegion, m_solidMaterialNames[targetIndex] );

    // space for element matrix and rhs

    using Kernels = SolidMechanicsLagrangianSSLEKernels::StressCalculationKernel;
    return SolidMechanicsLagrangianFEMKernels::
      ElementKernelLaunchSelector< Kernels >( numNodesPerElement,
                                              fe->n_quadrature_points(),
                                              &constitutiveRelation,
                                              elementSubRegion.size(),
                                              elemsToNodes,
                                              dNdX,
                                              detJ,
                                              incDisp );
  } );
}


REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianSSLE, string const &, dataRepository::Group * const )
} /* namespace geosx */
