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

#include "AcousticSecondOrderTTIFletcher.hpp"

// For the kernel factory used in applyStiffnessKernels
#include "AcousticTTIFletcherWaveEquationSEMKernel.hpp"
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
#include "AcousticTTIFields.hpp"

namespace geos
{

using namespace fields;

AcousticSecondOrderTTIFletcher::AcousticSecondOrderTTIFletcher( const std::string & name, Group * const parent )
  : AcousticSecondOrderAnisotropicTTISEM( name, parent )
{}

void AcousticSecondOrderTTIFletcher::registerDataOnMesh( Group & meshBodies )
{
  // First, call the TTI parent to register all TTI, VTI, and common base fields.
  AcousticSecondOrderAnisotropicTTISEM::registerDataOnMesh( meshBodies );

  // Now, register the single field that is specific to the Fletcher approximation.
  forDiscretizationOnMeshTargets( meshBodies, [&]( std::string const &, MeshLevel & mesh, string_array const & )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< acousticvtifields::AcousticSigma >( getName());
    } );
  } );
}


void AcousticSecondOrderTTIFletcher::applyStiffnessKernels( const real64 & dt, MeshLevel & mesh, const string_array & regionNames )
{
  //ok
  auto kernelFactory = acousticTTIFletcherWaveEquationSEMKernels::ExplicitAcousticTTIFletcherSEMFactory( dt );

  finiteElement::regionBasedKernelApplication< EXEC_POLICY,
                                               constitutive::NullModel,
                                               CellElementSubRegion >( mesh,
                                                                       regionNames,
                                                                       getDiscretizationName(),
                                                                       "",
                                                                       kernelFactory );
}

void AcousticSecondOrderTTIFletcher::initializeMatrices( MeshLevel & mesh, string_array const & regionNames )
{
  //ok
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

  // DOF arrays (these will already be zeroed by parent classes)
  arrayView1d< real32 > const dofEpsilon = nodeManager.getField< acousticvtifields::AcousticDofEpsilon >();
  arrayView1d< real32 > const dofDelta = nodeManager.getField< acousticvtifields::AcousticDofDelta >();
  arrayView1d< real32 > const dofOrder = nodeManager.getField< acousticvtifields::AcousticDofOrder >();
  arrayView1d< real32 > const dofTilt = nodeManager.getField< acousticttifields::AcousticDofTilt >();
  arrayView1d< real32 > const dofAzimuth = nodeManager.getField< acousticttifields::AcousticDofAzimuth >();

  /// get array of indicators: 1 if face is on the free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< acousticfields::AcousticFreeSurfaceFaceIndicator >();
  arrayView1d< localIndex const > const lateralSurfaceFaceIndicator = faceManager.getField< acousticvtifields::AcousticLateralSurfaceFaceIndicator >();
  arrayView1d< localIndex const > const bottomSurfaceFaceIndicator = faceManager.getField< acousticvtifields::AcousticBottomSurfaceFaceIndicator >();

  elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                              CellElementSubRegion & elementSubRegion )
  {
    finiteElement::FiniteElementBase const & fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName());

    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

    computeTargetNodeSet( elemsToNodes, elementSubRegion.size(), fe.getNumQuadraturePoints());

    arrayView1d< real32 const > const velocity = elementSubRegion.getField< acousticfields::AcousticVelocity >();
    arrayView1d< real32 const > const density = elementSubRegion.getField< acousticfields::AcousticDensity >();
    arrayView1d< real32 const > const vti_epsilon = elementSubRegion.getField< acousticvtifields::AcousticEpsilon >();
    arrayView1d< real32 const > const vti_delta = elementSubRegion.getField< acousticvtifields::AcousticDelta >();
    arrayView1d< real32 const > const tti_dipx = elementSubRegion.getField< acousticttifields::AcousticDipX >();
    arrayView1d< real32 const > const tti_dipy = elementSubRegion.getField< acousticttifields::AcousticDipY >();
    arrayView1d< real32 const > const vti_sigma = elementSubRegion.getField< acousticvtifields::AcousticSigma >();
    arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();

    /// Partial gradient if gradient has to be computed
    arrayView1d< real32 > grad = elementSubRegion.getField< acousticfields::PartialGradient >();
    grad.zero();

    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&]( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      // DOF arrays computation
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

      kernelDebug.template computeDofArraysTTI< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                              nodeCoords,
                                                                              elemsToNodes,
                                                                              tti_dipx,
                                                                              tti_dipy,
                                                                              dofTilt,
                                                                              dofAzimuth,
                                                                              dofOrder,
                                                                              sourceCoordinates,
                                                                              m_radiusIsoAroundSource );

      // Mass matrix computation
      AcousticMatricesSEM::MassMatrix< FE_TYPE > kernelM( finiteElement );
      kernelM.template computeMassMatrix< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                        nodeCoords,
                                                                        elemsToNodes,
                                                                        velocity,
                                                                        density,
                                                                        mass );

      // Fletcher-specific damping matrix computation
      AcousticMatricesSEM::DampingMatrix< FE_TYPE > kernelD( finiteElement );
      kernelD.template computeVTIFletcherDampingMatrices< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                                        nodeCoords,
                                                                                        elemsToFaces,
                                                                                        facesToNodes,
                                                                                        facesDomainBoundaryIndicator,
                                                                                        freeSurfaceFaceIndicator,
                                                                                        lateralSurfaceFaceIndicator,
                                                                                        bottomSurfaceFaceIndicator,
                                                                                        velocity,
                                                                                        density,
                                                                                        vti_epsilon,
                                                                                        vti_delta,
                                                                                        vti_sigma,
                                                                                        damping_pp,
                                                                                        damping_pq,
                                                                                        damping_qp,
                                                                                        damping_qq );
    } );
  } );
}


REGISTER_CATALOG_ENTRY( PhysicsSolverBase, AcousticSecondOrderTTIFletcher, string const &, dataRepository::Group * const )

} // namespace geos
