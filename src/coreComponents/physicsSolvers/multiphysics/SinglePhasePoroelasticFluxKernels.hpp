/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoroelasticFluxKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICFLUXKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICFLUXKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

namespace SinglePhasePoroelasticFluxKernels
{

/******************************** PermeabilityKernel ********************************/

struct PermeabilityKernel
{
  template< typename POLICY,
            typename PERM_WRAPPER,
            typename FE_TYPE >
  static void
  launch( localIndex const size,
          FE_TYPE finiteElementSpace,
          PERM_WRAPPER const & permWrapper,
          arrayView1d< real64 const > const & pressure,
          arrayView1d< real64 const > const & deltaPressure,
          arrayView2d< real64 const > const & porosity,
          arrayView2d< real64 const > const & dPorosity_dVolStrain,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodeLocation,
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const displacement,
          arrayView2d< localIndex const > const & elemsToNodes,
          arrayView3d< real64 > const & dPerm_dDisplacement )
  {

    GEOSX_UNUSED_VAR( dPorosity_dVolStrain );
    static constexpr int numNodesPerElem = FE_TYPE::numNodes;

    static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      real64 displacementLocal[numNodesPerElem][3];
      real64 xLocal[numNodesPerElem][3];

      for( localIndex a=0; a<numNodesPerElem; ++a )
      {
        localIndex const localNodeIndex = elemsToNodes( k, a );

        for( int i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeLocation[localNodeIndex][i];
          displacementLocal[a][i] = displacement[localNodeIndex][i];
        }
      }

      for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
      {
        real64 N[numNodesPerElem];
        real64 dNdX[numNodesPerElem][3];
        FE_TYPE::calcN( q, N );
        finiteElementSpace.template getGradN< FE_TYPE >( k, q, xLocal, dNdX );

        real64 strainIncrement[6] = {0};
        real64 dPerm_dVolStrain[3] = {0};

        FE_TYPE::symmetricGradient( dNdX, displacementLocal, strainIncrement );

        real64 const volStrain =   strainIncrement[0] + strainIncrement[1] + strainIncrement[2];

        permWrapper.updatePorosity( k, q, porosity[k][q] );

        permWrapper.updatePressureStrain( k, q, pressure[k] + deltaPressure[k], volStrain, dPerm_dVolStrain );

        // TODO: chain rule all dependencies to get to dPerm_dDisplacement
      }
    } );
  }
};


/******************************** EmbeddedSurfaceFluxKernel ********************************/

struct EmbeddedSurfaceFluxKernel
{
  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam SurfaceElementStencilWrapper The type of the stencil that is being used.
   * @param[in] stencil The stencil object.
   * @param[in] dt The timestep for the integration step.
   * @param[in] dofNumber The dofNumbers for each element
   * @param[in] pres The pressures in each element
   * @param[in] dPres The change in pressure for each element
   * @param[in] gravCoef The factor for gravity calculations (g*H)
   * @param[in] dens The material density in each element
   * @param[in] dDens_dPres The change in material density for each element
   * @param[in] mob The fluid mobility in each element
   * @param[in] dMob_dPres The derivative of mobility wrt pressure in each element
   * @param[in] permeability
   * @param[in] dPerm_dPres The derivative of permeability wrt pressure in each element
   * @param[in] transTMultiplier
   * @param[in] gravityVector
   * @param[out] localMatrix The linear system matrix
   * @param[out] localRhs The linear system residual
   */
  template< typename STENCIL_WRAPPER_TYPE >
  static void
  launch( STENCIL_WRAPPER_TYPE const & stencilWrapper,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
          ElementViewConst< arrayView1d< globalIndex const > > const & jumpDofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
          ElementViewConst< arrayView2d< real64 const > > const & transTMultiplier,
          R1Tensor const & gravityVector,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

  /**
   * @brief Compute flux and its derivatives for a given tpfa connector.
   *
   *
   */
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const numFluxElems,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           real64 const (&transmissibility)[2],
           real64 const (&dTrans_dPres)[2],
           real64 const (&dTrans_dAper)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & dens,
           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
           ElementViewConst< arrayView1d< real64 const > > const & mob,
           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian );
};


/******************************** FaceElementFluxKernel ********************************/

struct FaceElementFluxKernel
{
  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam SurfaceElementStencilWrapper The type of the stencil that is being used.
   * @param[in] stencil The stencil object.
   * @param[in] dt The timestep for the integration step.
   * @param[in] dofNumber The dofNumbers for each element
   * @param[in] pres The pressures in each element
   * @param[in] dPres The change in pressure for each element
   * @param[in] gravCoef The factor for gravity calculations (g*H)
   * @param[in] dens The material density in each element
   * @param[in] dDens_dPres The change in material density for each element
   * @param[in] mob The fluid mobility in each element
   * @param[in] dMob_dPres The derivative of mobility wrt pressure in each element
   * @param[in] permeability
   * @param[in] dPerm_dPres The derivative of permeability wrt pressure in each element
   * @param[in] transTMultiplier
   * @param[in] gravityVector
   * @param[out] localMatrix The linear system matrix
   * @param[out] localRhs The linear system residual
   */
  template< typename STENCIL_WRAPPER_TYPE >
  static void
  launch( STENCIL_WRAPPER_TYPE const & stencilWrapper,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
          ElementViewConst< arrayView2d< real64 const > > const & transTMultiplier,
          R1Tensor const & gravityVector,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs,
          CRSMatrixView< real64, localIndex const > const & dR_dAper );

  /**
   * @brief Compute flux and its derivatives for a given tpfa connector.
   *
   *
   */
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const numFluxElems,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           real64 const (&transmissibility)[2],
           real64 const (&dTrans_dPres)[2],
           real64 const (&dTrans_dAper)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & dens,
           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
           ElementViewConst< arrayView1d< real64 const > > const & mob,
           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian,
           arraySlice2d< real64 > const & dFlux_dAperture );
};


} // namespace PoroelasticFluxKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICFLUXKERNELS_HPP
