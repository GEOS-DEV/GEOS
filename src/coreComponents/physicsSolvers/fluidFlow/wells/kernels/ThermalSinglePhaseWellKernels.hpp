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
    :    Base( rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs ),
    m_isProducer( isProducer ),
    m_iwelemControl( subRegion.getTopWellElementIndex() ),
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
      real64 localEnergyAccum;
      real64 localEnergyAccumDP;
      real64 localEnergyAccumDT;
      if( iwelem == m_iwelemControl && !m_isProducer )
      {
        integer const numRows = 1+ IS_THERMAL;
        // For top segment energy balance eqn replaced with  T(n+1) - T = 0
        // No other energy balance derivatives
        // Assumption is global index == 0 is top segment with fixed temp BC

        localEnergyAccum = 0.0;
        localEnergyAccumDT = 1.0;
        localEnergyAccumDP = 0.0;
      }
      else
      {
        localEnergyAccum = m_wellElemVolume[iwelem] * ( m_wellElemDensity[iwelem][0]*m_internalEnergy[iwelem][0] - m_wellElemDensity_n[iwelem][0]*m_internalEnergy_n[iwelem][0]);
        localEnergyAccumDP = m_wellElemVolume[iwelem] *(m_internalEnergy[iwelem][0] *m_dWellElemDensity[iwelem][0][FLUID_PROP_COFFSET::dP] +
                                                        m_wellElemDensity[iwelem][0]*m_dInternalEnergy[iwelem][0][FLUID_PROP_COFFSET::dP]);
        localEnergyAccumDT = m_wellElemVolume[iwelem] *(m_internalEnergy[iwelem][0] *m_dWellElemDensity[iwelem][0][FLUID_PROP_COFFSET::dT] +
                                                        m_wellElemDensity[iwelem][0]*m_dInternalEnergy[iwelem][0][FLUID_PROP_COFFSET::dT]);
      }
      m_localRhs[eqnRowIndex]  = localEnergyAccum;
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
  /// Well type
  integer const m_isProducer;
  /// Index of the element where the control is enforced
  localIndex const m_iwelemControl;
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

  using FLUID_PROP_COFFSET = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;
  using WJ_COFFSET = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
  using WJ_ROFFSET = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;
  using TAG = singlePhaseWellKernels::ElemTag;


  using Base::m_isProducer;
  using Base::m_dt;
  using Base::m_localRhs;
  using Base::m_localMatrix;
  using Base::m_rankOffset;
  using Base::m_wellElemDofNumber;
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
    {

      if( iwelemNext < 0 )
      {
        globalIndex dofRowIndices  =  m_wellElemDofNumber[iwelemUp] + WJ_ROFFSET::ENERGYBAL - m_rankOffset;

        if( dofRowIndices  >= 0 && dofRowIndices < m_localMatrix.numRows() )
        {
          globalIndex dofColIndices_dRate =   m_wellElemDofNumber[iwelemUp] + WJ_COFFSET::dQ;
          globalIndex dofColIndices[FLUID_PROP_COFFSET::nDer]{};
          dofColIndices[0]= m_wellElemDofNumber[iwelemUp] + WJ_COFFSET::dP;
          dofColIndices[1]= m_wellElemDofNumber[iwelemUp] + WJ_COFFSET::dT;

          real64 localEnergyFlux[1]{};
          real64 localEnergyFluxJacobian_dQ[1]{};
          real64 localEnergyFluxJacobian[numDof]{};

          real64 eflux = -m_dt * m_enthalpy[iwelemUp][0];
          localEnergyFlux[0]= eflux * currentConnRate;
          localEnergyFluxJacobian_dQ[0] =   -m_dt *m_dEnthalpy[iwelemUp][0];
          localEnergyFluxJacobian[FLUID_PROP_COFFSET::dP] =  eflux * currentConnRate * m_dEnthalpy[iwelemUp][0][FLUID_PROP_COFFSET::dP];
          localEnergyFluxJacobian[FLUID_PROP_COFFSET::dT] =  eflux * currentConnRate * m_dEnthalpy[iwelemUp][0][FLUID_PROP_COFFSET::dT];
          if( !m_isProducer && m_globalWellElementIndex[iwelem] == 0 )
          {
            eflux= 0.0;
            localEnergyFluxJacobian_dQ[0]  = 0.0;
            localEnergyFluxJacobian[FLUID_PROP_COFFSET::dP] = 0.0;
            localEnergyFluxJacobian[FLUID_PROP_COFFSET::dT] = 0.0;
          }

          m_localMatrix.template addToRow< parallelDeviceAtomic >( dofRowIndices,
                                                                   &dofColIndices_dRate,
                                                                   localEnergyFluxJacobian_dQ[0],
                                                                   1 );
          m_localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dofRowIndices,
                                                                                       dofColIndices,
                                                                                       localEnergyFluxJacobian[0],
                                                                                       FLUID_PROP_COFFSET::nDer );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[dofRowIndices], localEnergyFlux[0] );
        }
      }
      else
      {
        globalIndex row_current = m_wellElemDofNumber[iwelemUp] + WJ_ROFFSET::ENERGYBAL   - m_rankOffset;
        globalIndex row_next = m_wellElemDofNumber[iwelemNext] + WJ_ROFFSET::ENERGYBAL   - m_rankOffset;

        // Setup Jacobian global row indicies
        // equations for  ENERGY balances
        globalIndex eqnRowIndices[2]{};
        eqnRowIndices[TAG::CURRENT ] = row_current;
        eqnRowIndices[TAG::NEXT ]    = row_next;

        // Setup Jacobian global col indicies  ( Mapping from local jac order to well jac order)
        globalIndex dofColIndices[FLUID_PROP_COFFSET::nDer]{};
        dofColIndices[0]= m_wellElemDofNumber[iwelemUp] + WJ_COFFSET::dP;
        dofColIndices[1]= m_wellElemDofNumber[iwelemUp] + WJ_COFFSET::dT;

        globalIndex dofColIndices_dRate =  m_wellElemDofNumber[iwelemUp]   + WJ_COFFSET::dQ;

        real64 localEnergyFlux[2]{};
        real64 localEnergyFluxJacobian_dQ[2]{};
        real64 localEnergyFluxJacobian[2*numDof]{};

        real64 eflux    =  m_dt * m_enthalpy[iwelemUp][0];
        real64 eflux_dq =  m_dt * m_dEnthalpy[iwelemUp][0];
        real64 dprop_dp =  eflux * currentConnRate  *m_dEnthalpy[iwelemUp][0][FLUID_PROP_COFFSET::dP];
        real64 dprop_dt =  eflux * currentConnRate * m_dEnthalpy[iwelemUp][0][FLUID_PROP_COFFSET::dT];

        localEnergyFlux[TAG::NEXT   ]   =   eflux * currentConnRate;
        localEnergyFlux[TAG::CURRENT  ] = -eflux * currentConnRate;

        localEnergyFluxJacobian_dQ [TAG::NEXT   ]  =   eflux_dq;
        localEnergyFluxJacobian_dQ [TAG::CURRENT]  =  -eflux_dq;

        localEnergyFluxJacobian[TAG::NEXT ][FLUID_PROP_COFFSET::dP] = dprop_dp;
        localEnergyFluxJacobian[TAG::NEXT][FLUID_PROP_COFFSET::dT]  = dprop_dt;

        localEnergyFluxJacobian[TAG::CURRENT][FLUID_PROP_COFFSET::dP] = -dprop_dp;
        localEnergyFluxJacobian[TAG::CURRENT][FLUID_PROP_COFFSET::dT] = -dprop_dt;

        // Note this updates diag and offdiag
        for( integer i = 0; i < 2; ++i )
        {
          if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < m_localMatrix.numRows() )
          {
            m_localMatrix.template addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                     &dofColIndices_dRate,
                                                                     localEnergyFluxJacobian_dQ[i],
                                                                     1 );
            m_localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                                         dofColIndices,
                                                                                         localEnergyFluxJacobian[i],
                                                                                         FLUID_PROP_COFFSET::nDer );
            RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[eqnRowIndices[i]], localEnergyFlux[i] );
          }
        }
      }

    } );

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

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */

class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 2 >
{
public:


  using WJ_ROFFSET = singlePhaseWellKernels::RowOffset_WellJac< 1 >;

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 2 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      WellElementSubRegion const & subRegion,
                      constitutive::SingleFluidBase const & fluid,
                      WellControls const & wellControls,
                      real64 const time,
                      real64 const dt,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_dt( dt ),
    m_isLocallyOwned( subRegion.isLocallyOwned() ),
    m_iwelemControl( subRegion.getTopWellElementIndex() ),
    m_isProducer( wellControls.isProducer() ),
    m_currentControl( wellControls.getControl() ),
    m_targetBHP( wellControls.getTargetBHP( time ) ),
    m_targetRate( wellControls.getTargetTotalRate( time ) ),
    m_targetMassRate( wellControls.getTargetMassRate( time ) ),
    m_volume( subRegion.getElementVolume() ),
    m_density_n( fluid.density_n() ),

    m_internalEnergy_n( fluid.internalEnergy_n() )
  {}


  GEOS_HOST_DEVICE
  void computeMassEnergyNormalizers( localIndex const iwelem,
                                     real64 & massNormalizer,
                                     real64 & energyNormalizer ) const
  {
    massNormalizer = LvArray::math::max( m_minNormalizer, m_density_n[iwelem][0] *  m_volume[iwelem] );
    energyNormalizer = m_internalEnergy_n[iwelem][0] * m_density_n[iwelem][0]  *    m_volume[iwelem];

    // warning: internal energy can be negative
    energyNormalizer = LvArray::math::max( m_minNormalizer, LvArray::math::abs( energyNormalizer ) );
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const iwelem,
                            LinfStackVariables & stack ) const override
  {
    real64 normalizer = 0.0;
    for( integer idof = 0; idof < WJ_ROFFSET::nEqn; ++idof )
    {
      // Step 1: compute a normalizer for the control or pressure equation

      // for the control equation, we distinguish two cases
      if( idof == WJ_ROFFSET::CONTROL )
      {

        // for the top well element, normalize using the current control
        if( m_isLocallyOwned && iwelem == m_iwelemControl )
        {
          if( m_currentControl == WellControls::Control::BHP )
          {
            // the residual entry is in pressure units
            normalizer = m_targetBHP;
          }
          else if( m_currentControl == WellControls::Control::TOTALVOLRATE )
          {
            // the residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetRate ), m_minNormalizer );
          }
          else if( m_currentControl == WellControls::Control::MASSRATE )
          {
            // the residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetMassRate ), m_minNormalizer );
          }
        }
        // for the pressure difference equation, always normalize by the BHP
        else
        {
          normalizer = m_targetBHP;
        }
      }
      // Step 2: compute a normalizer for the mass balance equation
      else if( idof == WJ_ROFFSET::MASSBAL )
      {

        // this residual entry is in mass units
        normalizer = m_dt * LvArray::math::abs( m_targetRate ) * m_density_n[iwelem][0];

        // to make sure that everything still works well if the rate is zero, we add this check
        normalizer = LvArray::math::max( normalizer, m_volume[iwelem] * m_density_n[iwelem][0] );


        // we have the normalizer now, we can compute a dimensionless Linfty norm contribution
        real64 const val = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / normalizer;
        if( val > stack.localValue[0] )
        {
          stack.localValue[0] = val;
        }

      }
      // step 3: energy residual
      else if( idof == WJ_ROFFSET::ENERGYBAL )
      {
        real64 massNormalizer = 0.0, energyNormalizer = 0.0;
        computeMassEnergyNormalizers( iwelem, massNormalizer, energyNormalizer );
        real64 const valEnergy = LvArray::math::abs( m_localResidual[stack.localRow + WJ_ROFFSET::ENERGYBAL] ) / energyNormalizer;
        if( valEnergy > stack.localValue[1] )
        {
          stack.localValue[1] = valEnergy;
        }
      }
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const iwelem,
                          L2StackVariables & stack ) const override
  {
    GEOS_UNUSED_VAR( iwelem, stack );
    GEOS_ERROR( "The L2 norm is not implemented for SinglePhaseWell" );
  }


protected:

  /// Time step size
  real64 const m_dt;

  /// Flag indicating whether the well is locally owned or not
  bool const m_isLocallyOwned;

  /// Index of the element where the control is enforced
  localIndex const m_iwelemControl;

  /// Flag indicating whether the well is a producer or an injector
  bool const m_isProducer;

  /// Controls
  WellControls::Control const m_currentControl;
  real64 const m_targetBHP;
  real64 const m_targetRate;
  real64 const m_targetMassRate;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on phase/total density at the previous converged time step
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_density_n;
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_internalEnergy_n;

};

/*
 *@class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp number of fluid components
   * @param[in] numDof number of dofs per well element
   * @param[in] targetPhaseIndex the index of the target phase (for phase volume control)
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the well element subregion
   * @param[in] fluid the fluid model
   * @param[in] wellControls the controls
   * @param[in] time the time
   * @param[in] dt the time step size
   * @param[out] residualNorm the residual norm on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   WellElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   WellControls const & wellControls,
                   real64 const time,
                   real64 const dt,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[2] )
  {
    using kernelType = ResidualNormKernel;
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    kernelType kernel( rankOffset, localResidual, dofNumber, ghostRank,
                       subRegion, fluid, wellControls, time, dt, minNormalizer );
    kernelType::template launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );

  }

};

} // end namespace singlePhaseWellKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP
