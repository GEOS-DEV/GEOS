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

template< typename REGIONTYPE >
struct PermeabilityKernel
{};

template<>
struct PermeabilityKernel< CellElementSubRegion >
{
  template< typename POLICY, typename PERM_WRAPPER, typename FE_TYPE >
  static void
  launch( localIndex const size,
          FE_TYPE finiteElementSpace,
          PERM_WRAPPER const & permWrapper,
          arrayView1d< real64 const > const & pressure,
          arrayView2d< real64 const > const & displacement,
          arrayView3d< real64 > const & dPerm_dDisplacement )
  {
    // TODO: pass also the porosity and chain rule all dependencies.
    static constexpr int numNodesPerElem = FE_TYPE::numNodes;

    static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      real64 displacementLocal = LV_ARRAY_INIT( displacement[k] );
      for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
      {
        real64 N[numNodesPerElem];
        real64 dNdX[numNodesPerElem][3];
        FE_TYPE::calcN( q, N );
        real64 const detJxW = finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

        real64 strainIncrement[6] = {0};

        FE_TYPE::symmetricGradient( dNdX, displacementLocal, strainIncrement );

        permWrapper.updatePressureStrain( k, q, pressure, volStrain, dPerm_dVolStrain );
      }
    } );
  }

  template< typename POLICY, typename PERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          PERM_WRAPPER const & permWrapper,
          arrayView2d< real64 const > const & displacement )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < permWrapper.numGauss(); ++q )
      {
        permWrapper.updatePressureStrain( k, q, pressure, volStrain, dPerm_dVolStrain );
      }
    } );
  }
};

// TODO is there a way to use the FlowSolverBaseKernel here ?

template<>
struct PermeabilityKernel< SurfaceElementSubRegion >
{
  template< typename POLICY, typename PERM_WRAPPER >
  static void
  launch( localIndex const size,
          PERM_WRAPPER const & permWrapper,
          arrayView1d< real64 const > const & effectiveAperture )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < permWrapper.numGauss(); ++q )
      {
        permWrapper.update( k, q, effectiveAperture[k] );
      }
    } );
  }

  template< typename POLICY, typename PERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          PERM_WRAPPER const & permWrapper,
          arrayView1d< real64 const > const & effectiveAperture )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < permWrapper.numGauss(); ++q )
      {
        permWrapper.update( k, q, effectiveAperture[k] );
      }
    } );
  }
};


/******************************** FluxKernel ********************************/

struct FluxKernel
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
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
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
  template< typename STENCILWRAPPER_TYPE >
  static void
  launch( STENCILWRAPPER_TYPE const & stencilWrapper,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
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
          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
          R1Tensor const & GEOSX_UNUSED_PARAM ( gravityVector ),
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    constexpr localIndex maxNumFluxElems = STENCILWRAPPER_TYPE::NUM_POINT_IN_FLUX;
    constexpr localIndex maxStencilSize = STENCILWRAPPER_TYPE::MAX_STENCIL_SIZE;

    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();


    forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
    {
      localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
      localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );

      // working arrays
      stackArray1d< globalIndex, maxNumFluxElems > dofColIndices( stencilSize );
      stackArray1d< real64, maxNumFluxElems > localFlux( numFluxElems );
      stackArray2d< real64, maxNumFluxElems * maxStencilSize > localFluxJacobian( numFluxElems, stencilSize );

      // compute transmissibility
      real64 transmissiblity[2], dTrans_dPres[2];
      stencilWrapper.computeTransmissibility( iconn, permeability, transmissiblity );
      stencilWrapper.dTrans_dPressure( iconn, dPerm_dPres, dTrans_dPres );

      compute( stencilSize,
               seri[iconn],
               sesri[iconn],
               sei[iconn],
               transmissiblity,
               dTrans_dPres,
               pres,
               dPres,
               gravCoef,
               dens,
               dDens_dPres,
               mob,
               dMob_dPres,
               dt,
               localFlux,
               localFluxJacobian );


      // extract DOF numbers
      for( localIndex i = 0; i < stencilSize; ++i )
      {
        dofColIndices[i] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
      }

      for( localIndex i = 0; i < numFluxElems; ++i )
      {

        if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
        {
          globalIndex const globalRow = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GT( localMatrix.numRows(), localRow );

          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux[i] );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                            dofColIndices.data(),
                                                                            localFluxJacobian[i].dataIfContiguous(),
                                                                            stencilSize );

        }
      }

    } );
  }

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


} // namespace PoroelasticFluxKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICFLUXKERNELS_HPP
