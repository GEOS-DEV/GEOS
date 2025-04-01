/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALSINGLEPHASEWELLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALSINGLEPHASEWELLKERNELS_HPP

#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"

namespace geos
{

namespace thermalSinglePhaseWellKernels
{



/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam IS_THERMAL thermal switch
 * @brief Define the interface for the assembly kernel in charge of accumulation and energy balance
 */
template< integer IS_THERMAL >
class ElementBasedAssemblyKernel : public singlePhaseWellKernels::ElementBasedAssemblyKernel< IS_THERMAL >
{
public:

  // Well jacobian column and row indicies
  using FLUID_PROP_COFFSET = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;
  using WJ_COFFSET = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
  using WJ_ROFFSET = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;

  using Base = singlePhaseWellKernels::ElementBasedAssemblyKernel< IS_THERMAL >;
  using Base::m_rankOffset;
  using Base::m_dofNumber;
  using Base::m_elemGhostRank;
  using Base::m_wellElemVolume;
  using Base::m_wellElemDensity;
  using Base::m_wellElemDensity_n;
  using Base::m_dWellElemDensity;
  using Base::m_localMatrix;
  using Base::m_localRhs;

  /**
   * @brief Constructor
   * @param[in] isProducer well type
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( integer const isProducer,
                              globalIndex const rankOffset,
                              string const dofKey,
                              WellElementSubRegion const & subRegion,
                              constitutive::SingleFluidBase const & fluid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :    Base( isProducer, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs ),
    m_internalEnergy( fluid.internalEnergy() ),
    m_internalEnergy_n( fluid.internalEnergy_n() ),
    m_dInternalEnergy( fluid.dInternalEnergy() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const iwelem ) const
  { return m_elemGhostRank( iwelem ); }


  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] iwelem the element index
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const iwelem ) const
  {
    Base::computeAccumulation( iwelem, [&]( globalIndex const presDofColIndex, globalIndex const tempDofColIndex )
    {
      // assemble the accumulation term of the energy equation
      localIndex const eqnRowIndex = m_dofNumber[iwelem] + WJ_ROFFSET::ENERGYBAL - m_rankOffset;

      real64 const localEnergyAccum = m_wellElemVolume[iwelem] * ( m_wellElemDensity[iwelem][0]*m_internalEnergy[iwelem][0] - m_wellElemDensity_n[iwelem][0]*m_internalEnergy_n[iwelem][0]);
      real64 const localEnergyAccumDP = m_wellElemVolume[iwelem] *(m_internalEnergy[iwelem][0] *m_dWellElemDensity[iwelem][0][FLUID_PROP_COFFSET::dP] +
                                                                   m_wellElemDensity[iwelem][0]*m_dInternalEnergy[iwelem][0][FLUID_PROP_COFFSET::dP]);
      real64 const localEnergyAccumDT = m_wellElemVolume[iwelem] *(m_internalEnergy[iwelem][0] *m_dWellElemDensity[iwelem][0][FLUID_PROP_COFFSET::dT] +
                                                                   m_wellElemDensity[iwelem][0]*m_dInternalEnergy[iwelem][0][FLUID_PROP_COFFSET::dT]);

      m_localRhs[eqnRowIndex] += localEnergyAccum;
      m_localMatrix.template addToRow< serialAtomic >( eqnRowIndex, &presDofColIndex, &localEnergyAccumDP, 1 );

      m_localMatrix.template addToRow< serialAtomic >( eqnRowIndex, &tempDofColIndex, &localEnergyAccumDT, 1 );

    } );
  }



  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
    {
      if( kernelComponent.elemGhostRank( iwelem ) >= 0 )
      {
        return;
      }
      kernelComponent.computeAccumulation( iwelem );
    } );
  }

protected:

  /// Views on fluid internal energy
  arrayView2d< real64 const > const m_internalEnergy;
  arrayView2d< real64 const > const m_internalEnergy_n;
  arrayView3d< real64 const > const m_dInternalEnergy;

};


/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:
  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] isProducer well type
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const isProducer,
                   globalIndex const rankOffset,
                   string const dofKey,
                   WellElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr isThermal = 1;
    ElementBasedAssemblyKernel< isThermal >
    kernel( isProducer, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs );
    ElementBasedAssemblyKernel< isThermal >::template
    launch< POLICY, ElementBasedAssemblyKernel< isThermal > >( subRegion.size(), kernel );

  }
};

/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @class FaceBasedAssemblyKernel
 * @tparam IS_THERMAL flag to include temperature dependencies and energy balance
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer IS_THERMAL >
class FaceBasedAssemblyKernel : public singlePhaseWellKernels::FaceBasedAssemblyKernel< IS_THERMAL >
{
public:

  using Base  = singlePhaseWellKernels::FaceBasedAssemblyKernel< IS_THERMAL >;

  // Well jacobian column and row indicies
  using COFFSET = singlePhaseWellKernels::ColOffset;
  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using TAG = singlePhaseWellKernels::ElemTag;

  using FLUID_PROP_COFFSET = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;
  using WJ_COFFSET = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
  using WJ_ROFFSET = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;

  using CP_Deriv = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;


  using Base::m_isProducer;
  using Base::m_dt;
  using Base::m_localRhs;
  using Base::m_localMatrix;
  using Base::m_rankOffset;
  using Base::maxNumElems;
  using Base::maxStencilSize;


  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = WJ_COFFSET::nDer;

/// Compile time value for the number of equations except volume and momentum
  static constexpr integer numEqn = WJ_ROFFSET::nEqn - 2;

  /**
   * @brief Constructor for the kernel interface
   *
   * @param[in] dt time step size
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor
   * @param[in] wellControls well information
   * @param[in] subRegion  region containing well
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( real64 const dt,
                           globalIndex const rankOffset,
                           string const wellDofKey,
                           WellControls const & wellControls,
                           WellElementSubRegion const & subRegion,
                           constitutive::SingleFluidBase const & fluid,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    : Base( dt
            , rankOffset
            , wellDofKey
            , wellControls
            , subRegion
            , localMatrix
            , localRhs ),

    m_globalWellElementIndex( subRegion.getGlobalWellElementIndex() ),
    m_enthalpy( fluid.enthalpy()),
    m_dEnthalpy( fluid.dEnthalpy())
  { }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] ie the element index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */

  GEOS_HOST_DEVICE
  inline
  void computeFlux( localIndex const iwelem ) const
  {
    Base::computeFlux ( iwelem, [&] ( localIndex const & iwelemNext
                                      , localIndex const & iwelemUp
                                      , real64 const & currentConnRate )
    {} );

  }


  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElements the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack
   * variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElements,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const ie )
    {
      kernelComponent.computeFlux( ie );

    } );
  }

protected:

  /// Global index of local element
  arrayView1d< globalIndex const > m_globalWellElementIndex;

  /// Views on enthalpy
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > m_enthalpy;
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > m_dEnthalpy;


};

/**
 * @class FaceBasedAssemblyKernelFactory
 */
class FaceBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] dt time step size
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] wellControls object holding well control/constraint information
   * @param[in] subregion well subregion
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch(
    real64 const dt,
    globalIndex const rankOffset,
    string const dofKey,
    WellControls const & wellControls,
    WellElementSubRegion const & subRegion,
    constitutive::SingleFluidBase const & fluid,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs )
  {
    integer constexpr isThermal=1;
    using kernelType = FaceBasedAssemblyKernel< isThermal >;
    kernelType kernel( dt, rankOffset, dofKey, wellControls, subRegion, fluid, localMatrix, localRhs );
    kernelType::template launch< POLICY >( subRegion.size(), kernel );
  }
};


} // end namespace singlePhaseWellKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP
