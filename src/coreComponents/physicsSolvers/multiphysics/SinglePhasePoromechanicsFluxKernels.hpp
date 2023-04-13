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
 * @file SinglePhasePoromechanicsFluxKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSFLUXKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSFLUXKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"


namespace geosx
{

namespace singlePhasePoromechanicsFluxKernels
{


template< integer NUM_DOF, typename STENCILWRAPPER >
class EmbeddedAssemblyKernel : public singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_DOF, STENCILWRAPPER >
{
 public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using AbstractBase = singlePhaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_mob;
  using AbstractBase::m_dens;

  using Base = singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_DOF, STENCILWRAPPER >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  EmbeddedAssemblyKernel( globalIndex const rankOffset,
                          STENCILWRAPPER const & stencilWrapper,
                          DofNumberAccessor const & dofNumberAccessor,
                          SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                          ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                          SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                          ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors,
                          PermeabilityAccessors const & permeabilityAccessors,
                          ThermalConductivityAccessors const & thermalConductivityAccessors,
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
    : Base( rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_dispJumpDofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_dPerm_dDispJump( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) )
  {} 

 private:
    
    ElementViewConst< arrayView1d< globalIndex const > > const m_dispJumpDofNumber;

    ElementViewConst< arrayView4d< real64 const > > const m_dPerm_dDispJump;
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
  template< typename STENCIL_WRAPPER_TYPE >
  static void
  launch( STENCIL_WRAPPER_TYPE const & stencilWrapper,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
          ElementViewConst< arrayView1d< globalIndex const > > const & jumpDofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPerm_dDispJump,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

  /**
   * @brief Compute flux and its derivatives for a given tpfa connector.
   *
   *
   */
  template< localIndex MAX_NUM_CONNECTIONS >
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const numFluxElems,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           real64 const (&transmissibility)[MAX_NUM_CONNECTIONS][2],
           real64 const (&dTrans_dPres)[MAX_NUM_CONNECTIONS][2],
           real64 const (&dTrans_dDispJump)[MAX_NUM_CONNECTIONS][2][3],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
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
 * @param[in] stencilWrapper The stencil wrapper object.
 * @param[in] dt The timestep for the integration step.
 * @param[in] rankOffset The rank offset
 * @param[in] pressureDofNumber The pressure dof number for each element
 * @param[in] ghostRank The ghost rank
 * @param[in] pres The pressures in each element
 * @param[in] gravCoef The factor for gravity calculations (g*H)
 * @param[in] dens The material density in each element
 * @param[in] dDens_dPres The change in material density for each element
 * @param[in] mob The fluid mobility in each element
 * @param[in] dMob_dPres The derivative of mobility wrt pressure in each element
 * @param[in] permeability The permeability in each element
 * @param[in] dPerm_dPres The derivative of permeability wrt pressure in each element
 * @param[in] dPerm_dDispJump The derivative of permeability wrt aperture in each element
 * @param[in] permeabilityMultiplier The permeability multiplier
 * @param[in] gravityVector The gravity vector
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
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPerm_dDispJump,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs,
          CRSMatrixView< real64, localIndex const > const & dR_dAper );

  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam SurfaceElementStencilWrapper The type of the stencil that is being used.
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
   * @param[in] permeabilityMultiplier
   * @param[in] gravityVector
   * @param[out] localMatrix The linear system matrix
   * @param[out] localRhs The linear system residual
   */
  static void
  launch( SurfaceElementStencilWrapper const & stencilWrapper,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPerm_dDispJump,
          ElementViewConst< arrayView3d< real64 const > > const & permeabilityMultiplier,
          R1Tensor const & gravityVector,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );



  /**
   * @brief Compute flux and its derivatives for a given tpfa connector.
   *
   *
   */
  template< localIndex MAX_NUM_CONNECTIONS >
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const numFluxElems,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           real64 const (&transmissibility)[MAX_NUM_CONNECTIONS][2],
           real64 const (&dTrans_dPres)[MAX_NUM_CONNECTIONS][2],
           real64 const (&dTrans_dDispJump)[MAX_NUM_CONNECTIONS][2][3],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
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


} // namespace singlePhasePoromechanicsFluxKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSFLUXKERNELS_HPP
