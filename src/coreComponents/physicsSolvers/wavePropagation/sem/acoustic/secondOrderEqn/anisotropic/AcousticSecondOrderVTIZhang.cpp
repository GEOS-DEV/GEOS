/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A / TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "AcousticSecondOrderVTIZhang.hpp"

// For the kernel factory used in applyStiffnessKernels
#include "AcousticVTIZhangWaveEquationSEMKernel.hpp"
#include "constitutive/NullModel.hpp"

// For initializeMatrices implementation
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticMatricesSEMKernel.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"

// Field includes
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"
#include "AcousticVTIFields.hpp"

namespace geos
{

using namespace fields;

AcousticSecondOrderVTIZhang::AcousticSecondOrderVTIZhang( const std::string & name, Group * const parent )
  : AcousticSecondOrderAnisotropicVTISEM( name, parent )
{}

void AcousticSecondOrderVTIZhang::registerDataOnMesh( Group & meshBodies )
{
  // First call parent implementation to register common fields
  AcousticSecondOrderAnisotropicVTISEM::registerDataOnMesh( meshBodies );

  // Zhang VTI doesn't need to register any additional fields beyond what the parent registers
  // (no equivalent to AcousticSigma like Fletcher has)
}

void AcousticSecondOrderVTIZhang::applyStiffnessKernels( const real64 & dt, MeshLevel & mesh, const string_array & regionNames )
{
  auto kernelFactory = acousticVTIZhangWaveEquationSEMKernels::ExplicitAcousticVTIZhangSEMFactory( dt );

  finiteElement::regionBasedKernelApplication< EXEC_POLICY,
                                               constitutive::NullModel,
                                               CellElementSubRegion >( mesh,
                                                                       regionNames,
                                                                       getDiscretizationName(),
                                                                       "",
                                                                       kernelFactory );
}

void AcousticSecondOrderVTIZhang::initializeMatrices( MeshLevel & mesh, string_array const & regionNames )
{
  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
  arrayView1d< integer const > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
  arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

  /// get face to node map
  ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

  // mass matrix to be computed in this function
  arrayView1d< real32 > const mass = nodeManager.getField< acousticfields::AcousticMassVector >();
  mass.zero();

  /// damping matrices to be computed for each dof in the boundary of the mesh
  arrayView1d< real32 > const damping_pp = nodeManager.getField< acousticvtifields::DampingVector_pp >();
  arrayView1d< real32 > const damping_pq = nodeManager.getField< acousticvtifields::DampingVector_pq >();
  arrayView1d< real32 > const damping_qp = nodeManager.getField< acousticvtifields::DampingVector_qp >();
  arrayView1d< real32 > const damping_qq = nodeManager.getField< acousticvtifields::DampingVector_qq >();
  damping_pp.zero();
  damping_pq.zero();
  damping_qp.zero();
  damping_qq.zero();

  arrayView1d< real32 > const dofEpsilon = nodeManager.getField< acousticvtifields::AcousticDofEpsilon >();
  arrayView1d< real32 > const dofDelta   = nodeManager.getField< acousticvtifields::AcousticDofDelta >();
  arrayView1d< real32 > const dofOrder   = nodeManager.getField< acousticvtifields::AcousticDofOrder >();
  dofEpsilon.zero();
  dofDelta.zero();
  dofOrder.zero();   // number of Hexa countaining a dof


  /// get array of indicators: 1 if face is on the free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< acousticfields::AcousticFreeSurfaceFaceIndicator >();
  arrayView1d< localIndex const > const lateralSurfaceFaceIndicator = faceManager.getField< acousticvtifields::AcousticLateralSurfaceFaceIndicator >();
  arrayView1d< localIndex const > const bottomSurfaceFaceIndicator = faceManager.getField< acousticvtifields::AcousticBottomSurfaceFaceIndicator >();


  elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                              CellElementSubRegion & elementSubRegion )
  {
    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

    computeTargetNodeSet( elemsToNodes, elementSubRegion.size(), fe.getNumQuadraturePoints() );

    arrayView1d< real32 const > const velocity = elementSubRegion.getField< acousticfields::AcousticVelocity >();
    arrayView1d< real32 const > const density = elementSubRegion.getField< acousticfields::AcousticDensity >();
    arrayView1d< real32 const > const vti_epsilon = elementSubRegion.getField< acousticvtifields::AcousticEpsilon >();
    arrayView1d< real32 const > const vti_delta = elementSubRegion.getField< acousticvtifields::AcousticDelta >();
    arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
    /// Partial gradient if gradient as to be computed
    arrayView1d< real32 > grad = elementSubRegion.getField< acousticfields::PartialGradient >();

    grad.zero();

    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );


      AcousticMatricesSEM::DofArrays< FE_TYPE > kernelDebug( finiteElement );
      kernelDebug.template computeDofArraysVTI< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                              nodeCoords,
                                                                              elemsToNodes,
                                                                              vti_epsilon,
                                                                              vti_delta,
                                                                              dofEpsilon,
                                                                              dofDelta,
                                                                              dofOrder,
                                                                              sourceCoordinates,
                                                                              m_radiusIsoAroundSource );


      AcousticMatricesSEM::MassMatrix< FE_TYPE > kernelM( finiteElement );
      kernelM.template computeMassMatrix< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                        nodeCoords,
                                                                        elemsToNodes,
                                                                        velocity,
                                                                        density,
                                                                        mass );

      AcousticMatricesSEM::DampingMatrix< FE_TYPE > kernelD( finiteElement );
      kernelD.template computeVTIZhangDampingMatrices< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                                     nodeCoords,
                                                                                     elemsToFaces,
                                                                                     facesToNodes,
                                                                                     facesDomainBoundaryIndicator,
                                                                                     freeSurfaceFaceIndicator,
                                                                                     lateralSurfaceFaceIndicator,
                                                                                     bottomSurfaceFaceIndicator,
                                                                                     velocity,
                                                                                     density,
                                                                                     dofEpsilon,
                                                                                     dofDelta,
                                                                                     damping_pp,
                                                                                     damping_pq,
                                                                                     damping_qp,
                                                                                     damping_qq );

    } );
  } );
}

// Registration
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, AcousticSecondOrderVTIZhang, string const &, dataRepository::Group * const )
} // namespace geos
