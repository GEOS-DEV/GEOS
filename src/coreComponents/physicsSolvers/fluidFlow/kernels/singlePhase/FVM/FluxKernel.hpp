/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FluxKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP

//#include "common/DataLayouts.hpp"
//#include "common/DataTypes.hpp"
//#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidBase.hpp"
//#include "constitutive/fluid/singlefluid/SlurryFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
//#include "fieldSpecification/AquiferBoundaryCondition.hpp"
//#include "finiteVolume/BoundaryStencil.hpp"
//#include "finiteVolume/FluxApproximationBase.hpp"
//#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/fields/FlowSolverBaseFields.hpp"
#include "FluxKernelsHelper.hpp"
#include "physicsSolvers/fluidFlow/fields/SinglePhaseBaseFields.hpp"
//#include "physicsSolvers/fluidFlow/kernels/singlePhase/HydrostaticPressureKernel.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
//#include "physicsSolvers/fluidFlow/kernels/singlePhase/FVM/MobilityKernel.hpp"
//#include "ElementBasedAssemblyKernel.hpp"
//#include "physicsSolvers/fluidFlow/kernels/singlePhase/FVM/FaceBasedAssemblyKernelBase.hpp"
//#include "physicsSolvers/fluidFlow/kernels/singlePhase/FVM/FaceBasedAssemblyKernel.hpp"

namespace geos
{

namespace singlePhaseFVMKernels
{

using namespace constitutive;

using namespace fluxKernelsHelper;

/******************************** FluxKernel ********************************/

// To delete it after verifying FaceBasedAssemblyKernel
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

  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::pressure,
                      fields::flow::pressure_n,
                      fields::flow::gravityCoefficient,
                      fields::flow::mobility,
                      fields::flow::dMobility_dPressure >;

  using SinglePhaseFluidAccessors =
    StencilMaterialAccessors< SingleFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity_dPressure >;

  using SlurryFluidAccessors =
    StencilMaterialAccessors< SlurryFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity_dPressure >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  using ProppantPermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure,
                              fields::permeability::dPerm_dDispJump,
                              fields::permeability::permeabilityMultiplier >;


  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
   * @param[in] stencil The stencil object.
   * @param[in] dt The timestep for the integration step.
   * @param[in] dofNumber The dofNumbers for each element
   * @param[in] pres The pressures in each element
   * @param[in] gravCoef The factor for gravity calculations (g*H)
   * @param[in] dens The material density in each element
   * @param[in] dDens_dPres The change in material density for each element
   * @param[in] mob The fluid mobility in each element
   * @param[in] dMob_dPres The derivative of mobility wrt pressure in each element
   * @param[in] permeability
   * @param[in] dPerm_dPres The derivative of permeability wrt pressure in each element
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
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

    constexpr localIndex maxNumElems = STENCILWRAPPER_TYPE::maxNumPointsInFlux;
    constexpr localIndex maxStencilSize = STENCILWRAPPER_TYPE::maxStencilSize;

    forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [stencilWrapper, dt, rankOffset, dofNumber, ghostRank,
                                                              pres, gravCoef, dens, dDens_dPres, mob,
                                                              dMob_dPres, permeability, dPerm_dPres,
                                                              seri, sesri, sei, localMatrix, localRhs] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
      localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );

      // working arrays
      stackArray1d< globalIndex, maxNumElems > dofColIndices( stencilSize );
      stackArray1d< real64, maxNumElems > localFlux( numFluxElems );
      stackArray2d< real64, maxNumElems * maxStencilSize > localFluxJacobian( numFluxElems, stencilSize );


      // compute transmissibility
      real64 transmissibility[STENCILWRAPPER_TYPE::maxNumConnections][2];
      real64 dTrans_dPres[STENCILWRAPPER_TYPE::maxNumConnections][2];

      stencilWrapper.computeWeights( iconn,
                                     permeability,
                                     dPerm_dPres,
                                     transmissibility,
                                     dTrans_dPres );

      compute( numFluxElems,
               seri[iconn],
               sesri[iconn],
               sei[iconn],
               transmissibility,
               dTrans_dPres,
               pres,
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
          GEOS_ASSERT_GE( localRow, 0 );
          GEOS_ASSERT_GT( localMatrix.numRows(), localRow );

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
  template< localIndex maxNumConnections >
  GEOS_HOST_DEVICE
  static void
  compute( localIndex const numFluxElems,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           real64 const (&transmissibility)[maxNumConnections][2],
           real64 const (&dTrans_dPres)[maxNumConnections][2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & dens,
           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
           ElementViewConst< arrayView1d< real64 const > > const & mob,
           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian )
  {

    localIndex k[2];
    localIndex connectionIndex = 0;;
    for( k[0]=0; k[0]<numFluxElems; ++k[0] )
    {
      for( k[1]=k[0]+1; k[1]<numFluxElems; ++k[1] )
      {
        real64 fluxVal = 0.0;
        real64 dFlux_dTrans = 0.0;
        real64 const trans[2] = {transmissibility[connectionIndex][0], transmissibility[connectionIndex][1]};
        real64 const dTrans[2] = { dTrans_dPres[connectionIndex][0], dTrans_dPres[connectionIndex][1] };
        real64 dFlux_dP[2] = {0.0, 0.0};
        localIndex const regionIndex[2]    = {seri[k[0]], seri[k[1]]};
        localIndex const subRegionIndex[2] = {sesri[k[0]], sesri[k[1]]};
        localIndex const elementIndex[2]   = {sei[k[0]], sei[k[1]]};


        computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
                                trans,
                                dTrans,
                                pres,
                                gravCoef,
                                dens,
                                dDens_dPres,
                                mob,
                                dMob_dPres,
                                fluxVal,
                                dFlux_dP,
                                dFlux_dTrans );

        // populate local flux vector and derivatives
        flux[k[0]] +=  dt * fluxVal;
        flux[k[1]] -=  dt * fluxVal;

        fluxJacobian[k[0]][k[0]] += dt * dFlux_dP[0];
        fluxJacobian[k[0]][k[1]] += dt * dFlux_dP[1];
        fluxJacobian[k[1]][k[0]] -= dt * dFlux_dP[0];
        fluxJacobian[k[1]][k[1]] -= dt * dFlux_dP[1];

        connectionIndex++;
      }
    }
  }
};

} // namespace singlePhaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
