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

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_DAMPINGMATRIXCOMPUTERS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_DAMPINGMATRIXCOMPUTERS_HPP_

#include "common/DataTypes.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticMatricesSEMKernel.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#include "AcousticVTIFields.hpp"

namespace geos
{
using namespace fields;
using wsCoordType = WaveSolverUtils::wsCoordType;

/**
 * @brief Base class for damping matrix computers
 */
class DampingMatrixComputerBase
{
public:
  virtual ~DampingMatrixComputerBase() = default;
};

/**
 * @brief Zhang-style damping computation (works for both VTI and TTI)
 */
class ZhangDampingComputer : public DampingMatrixComputerBase
{
public:
  template< typename FE_TYPE, typename EXEC_POLICY, typename ATOMIC_POLICY >
  void computeDampingMatrices(
    FE_TYPE const & finiteElement,
    CellElementSubRegion & elementSubRegion,
    arrayView2d< localIndex const > const & GEOS_UNUSED_PARAM( elemsToNodes ),
    arrayView2d< localIndex const > const & elemsToFaces,
    ArrayOfArraysView< localIndex const > const & facesToNodes,
    arrayView1d< integer const > const & facesDomainBoundaryIndicator,
    arrayView1d< localIndex const > const & freeSurfaceFaceIndicator,
    arrayView1d< localIndex const > const & lateralSurfaceFaceIndicator,
    arrayView1d< localIndex const > const & bottomSurfaceFaceIndicator,
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const & nodeCoords,
    arrayView1d< real32 const > const & velocity,
    arrayView1d< real32 const > const & density,
    arrayView1d< real32 const > const & GEOS_UNUSED_PARAM( vti_epsilon ),
    arrayView1d< real32 const > const & GEOS_UNUSED_PARAM( vti_delta ),
    arrayView1d< real32 > const & damping_pp,
    arrayView1d< real32 > const & damping_pq,
    arrayView1d< real32 > const & damping_qp,
    arrayView1d< real32 > const & damping_qq ) const
  {
    // Get DOF arrays needed for Zhang method
    NodeManager & nodeManager = elementSubRegion.getParent().getParent().getGroup< MeshLevel >( 0 ).getNodeManager();
    arrayView1d< real32 const > const dofEpsilon = nodeManager.getField< acousticvtifields::AcousticDofEpsilon >();
    arrayView1d< real32 const > const dofDelta = nodeManager.getField< acousticvtifields::AcousticDofDelta >();

    AcousticMatricesSEM::DampingMatrix< FE_TYPE > kernelD( finiteElement );

    kernelD.template computeVTIZhangDampingMatrices< EXEC_POLICY, ATOMIC_POLICY >(
      elementSubRegion.size(),
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
      damping_qq
      );
  }
};

/**
 * @brief Fletcher-style damping computation (works for both VTI and TTI)
 */
class FletcherDampingComputer : public DampingMatrixComputerBase
{
public:
  template< typename FE_TYPE, typename EXEC_POLICY, typename ATOMIC_POLICY >
  void computeDampingMatrices(
    FE_TYPE const & finiteElement,
    CellElementSubRegion & elementSubRegion,
    arrayView2d< localIndex const > const & GEOS_UNUSED_PARAM( elemsToNodes ),
    arrayView2d< localIndex const > const & elemsToFaces,
    ArrayOfArraysView< localIndex const > const & facesToNodes,
    arrayView1d< integer const > const & facesDomainBoundaryIndicator,
    arrayView1d< localIndex const > const & freeSurfaceFaceIndicator,
    arrayView1d< localIndex const > const & lateralSurfaceFaceIndicator,
    arrayView1d< localIndex const > const & bottomSurfaceFaceIndicator,
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const & nodeCoords,
    arrayView1d< real32 const > const & velocity,
    arrayView1d< real32 const > const & density,
    arrayView1d< real32 const > const & vti_epsilon,
    arrayView1d< real32 const > const & vti_delta,
    arrayView1d< real32 > const & damping_pp,
    arrayView1d< real32 > const & damping_pq,
    arrayView1d< real32 > const & damping_qp,
    arrayView1d< real32 > const & damping_qq ) const
  {
    // Fletcher always needs sigma field
    arrayView1d< real32 const > const vti_sigma = elementSubRegion.getField< acousticvtifields::AcousticSigma >();

    AcousticMatricesSEM::DampingMatrix< FE_TYPE > kernelD( finiteElement );

    kernelD.template computeVTIFletcherDampingMatrices< EXEC_POLICY, ATOMIC_POLICY >(
      elementSubRegion.size(),
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
      damping_qq
      );
  }
};

/**
 * @brief Template deduction for computeDampingMatrices
 *
 * This wrapper function explicitly specifies all template parameters,
 * eliminating template argument deduction issues that can occur when
 * the compiler has trouble matching argument types.
 */
template< typename DampingComputer, typename FE_TYPE, typename EXEC_POLICY, typename ATOMIC_POLICY >
void callComputeDampingMatrices(
  DampingComputer & dampingComputer,
  FE_TYPE const & finiteElement,
  CellElementSubRegion & elementSubRegion,
  arrayView2d< localIndex const > const & elemsToNodes,
  arrayView2d< localIndex const > const & elemsToFaces,
  ArrayOfArraysView< localIndex const > const & facesToNodes,
  arrayView1d< integer const > const & facesDomainBoundaryIndicator,
  arrayView1d< localIndex const > const & freeSurfaceFaceIndicator,
  arrayView1d< localIndex const > const & lateralSurfaceFaceIndicator,
  arrayView1d< localIndex const > const & bottomSurfaceFaceIndicator,
  arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const & nodeCoords,
  arrayView1d< real32 const > const & velocity,
  arrayView1d< real32 const > const & density,
  arrayView1d< real32 const > const & vti_epsilon,
  arrayView1d< real32 const > const & vti_delta,
  arrayView1d< real32 > const & damping_pp,
  arrayView1d< real32 > const & damping_pq,
  arrayView1d< real32 > const & damping_qp,
  arrayView1d< real32 > const & damping_qq )
{
  // Forward the call with explicitly specified template parameters
  dampingComputer.template computeDampingMatrices< FE_TYPE, EXEC_POLICY, ATOMIC_POLICY >(
    finiteElement,
    elementSubRegion,
    elemsToNodes,
    elemsToFaces,
    facesToNodes,
    facesDomainBoundaryIndicator,
    freeSurfaceFaceIndicator,
    lateralSurfaceFaceIndicator,
    bottomSurfaceFaceIndicator,
    nodeCoords,
    velocity,
    density,
    vti_epsilon,
    vti_delta,
    damping_pp,
    damping_pq,
    damping_qp,
    damping_qq
    );
}

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_DAMPINGMATRIXCOMPUTERS_HPP_
