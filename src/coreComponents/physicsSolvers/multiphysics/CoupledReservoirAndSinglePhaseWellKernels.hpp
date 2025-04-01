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
 * @file ThermalSinglePhaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALSINGLEPHASEWELLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALSINGLEPHASEWELLKERNELS_HPP


#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellTags.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"
namespace geos
{

namespace coupledReservoirAndSinglePhaseWellKernels
{
/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer IS_THERMAL >
class IsothermalSinglePhaseFluxKernel
{
public:

  /// Compile time value for the number of components
  static constexpr integer resNumDOF  = 1+IS_THERMAL;

  // Well jacobian column and row indicies
  using WJ_COFFSET = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
  using WJ_ROFFSET = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  using CP_Deriv = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

  using TAG = singlePhaseWellKernels::SubRegionTag;



  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = WJ_COFFSET::nDer;

  /// Compile time value for the number of equations except volume and momentum
  static constexpr integer numEqn = WJ_ROFFSET::nEqn - 2;  // tjb revisit

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  IsothermalSinglePhaseFluxKernel( real64 const dt,
                                   globalIndex const rankOffset,
                                   string const wellDofKey,
                                   ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber,
                                   WellElementSubRegion const & subRegion,
                                   PerforationData const * const perforationData,
                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                   arrayView1d< real64 > const & localRhs )
    :
    m_dt( dt ),
    m_rankOffset( rankOffset ),
    m_perfRate( perforationData->getField< fields::well::perforationRate >() ),
    m_dPerfRate( perforationData->getField< fields::well::dPerforationRate >() ),
    m_perfWellElemIndex( perforationData->getField< fields::perforation::wellElementIndex >() ),
    m_wellElemDofNumber( subRegion.getReference< array1d< globalIndex > >( wellDofKey ) ),
    m_resElemDofNumber( resDofNumber ),
    m_resElementRegion( perforationData->getField< fields::perforation::reservoirElementRegion >() ),
    m_resElementSubRegion( perforationData->getField< fields::perforation::reservoirElementSubRegion >() ),
    m_resElementIndex( perforationData->getField< fields::perforation::reservoirElementIndex >() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )

  { }


  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] ie the element index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */

  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void computeFlux( localIndex const iperf,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {

    // local working variables and arrays
    localIndex eqnRowIndices[ 2 ] = { -1 };
    globalIndex dofColIndices[ 2 ] = { -1 };


    real64 localPerf[ 2 ]{};
    real64 localPerfJacobian[ 2 ][ 2 ]{};

    // get the reservoir (sub)region and element indices
    localIndex const er = m_resElementRegion[iperf];
    localIndex const esr = m_resElementSubRegion[iperf];
    localIndex const ei = m_resElementIndex[iperf];

    // get the well element index for this perforation
    localIndex const iwelem = m_perfWellElemIndex[iperf];
    globalIndex const elemOffset = m_wellElemDofNumber[iwelem];

    // row index on reservoir side
    eqnRowIndices[TAG::RES] = m_resElemDofNumber[er][esr][ei] - m_rankOffset;
    // column index on reservoir side
    dofColIndices[TAG::RES] = m_resElemDofNumber[er][esr][ei];

    // row index on well side
    eqnRowIndices[TAG::WELL] = LvArray::integerConversion< localIndex >( elemOffset - m_rankOffset ) + ROFFSET::MASSBAL;
    // column index on well side
    dofColIndices[TAG::WELL] = elemOffset + COFFSET::DPRES;

    // populate local flux vector and derivatives
    localPerf[TAG::RES] = m_dt * m_perfRate[iperf];
    localPerf[TAG::WELL] = -localPerf[TAG::RES];

    for( integer ke = 0; ke < 2; ++ke )
    {
      localPerfJacobian[TAG::RES][ke] = m_dt * m_dPerfRate[iperf][ke][0];     // tjb tag
      localPerfJacobian[TAG::WELL][ke] = -localPerfJacobian[TAG::RES][ke];
    }

    for( integer i = 0; i < 2; ++i )
    {
      if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < m_localMatrix.numRows() )
      {
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                            &dofColIndices[0],
                                                                            &localPerfJacobian[0][0] + 2 * i,
                                                                            2 );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[eqnRowIndices[i]], localPerf[i] );
      }
    }
    //compFluxKernelOp( resOffset, wellElemOffset, dofColIndices, iwelem );

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

  /// Time step size
  real64 const m_dt;

  globalIndex const m_rankOffset;
  // Perfoation variables
  arrayView1d< real64 const > const m_perfRate;
  arrayView3d< real64 const > const m_dPerfRate;
  arrayView1d< localIndex const > const m_perfWellElemIndex;

  // Element region, subregion, index
  arrayView1d< globalIndex const > const m_wellElemDofNumber;
  ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const m_resElemDofNumber;
  arrayView1d< localIndex const > const m_resElementRegion;
  arrayView1d< localIndex const > const m_resElementSubRegion;
  arrayView1d< localIndex const > const m_resElementIndex;

  // RHS and Jacobian
  CRSMatrixView< real64, globalIndex const >  m_localMatrix;
  arrayView1d< real64 > const m_localRhs;

};

/**
 * @class FaceBasedAssemblyKernelFactory
 */
class IsothermalSinglePhaseFluxKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComps the number of fluid components
   * @param[in] dt time step size
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] useTotalMassEquation flag specifying whether to replace one component bal eqn with total mass eqn
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] wellControls object holding well control/constraint information
   * @param[in] subregion well subregion
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( real64 const dt,
                   globalIndex const rankOffset,
                   string const wellDofKey,

                   ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber,
                   WellElementSubRegion const & subRegion,
                   PerforationData const * const perforationData,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs
                   )
  {
    using kernelType = IsothermalSinglePhaseFluxKernel< 0 >;
    kernelType kernel( dt, rankOffset, wellDofKey, resDofNumber, subRegion, perforationData, localMatrix,
                       localRhs );
    kernelType::template launch< POLICY >( perforationData->size(), kernel );
  }
};


/**
 * @class FaceBasedAssemblyKernel
 * @tparam IS_THERMAL flag to include temperature derivatives and energy balance
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer IS_THERMAL >
class ThermalSinglePhaseFluxKernel : public IsothermalSinglePhaseFluxKernel< IS_THERMAL >
{
public:
  using Base = IsothermalSinglePhaseFluxKernel< IS_THERMAL >;
  // Reservoir degrees of freedom
  static constexpr integer resNumDOF  =  1+IS_THERMAL;

  // Well jacobian column and row indicies
  using WJ_COFFSET = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
  using WJ_ROFFSET = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;

  //using ROFFSET = singlePhaseWellKernels::RowOffset;
  //using COFFSET = singlePhaseWellKernels::ColOffset;

  using CP_Deriv = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

  //using TAG = singlePhaseWellKernels::SubRegionTag;

  using Base::m_dt;
  using Base::m_localRhs;
  using Base::m_localMatrix;
  using Base::m_rankOffset;



  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = WJ_COFFSET::nDer;

  /// Compile time value for the number of equations except volume and momentum
  static constexpr integer numEqn = WJ_ROFFSET::nEqn - 2;  // tjb revisity

  /**
   * @brief Constructor for the kernel interface   // tjb clean up
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ThermalSinglePhaseFluxKernel( integer const isProducer,
                                real64 const dt,

                                globalIndex const rankOffset,
                                string const wellDofKey,

                                ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber,
                                WellElementSubRegion const & subRegion,
                                PerforationData const * const perforationData,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs )
    : Base( dt,
            rankOffset,
            wellDofKey,
            resDofNumber,
            subRegion,
            perforationData,
            localMatrix,
            localRhs ),
    m_isProducer( isProducer ),
    m_globalWellElementIndex( subRegion.getGlobalWellElementIndex() ),
    m_energyPerfFlux( perforationData->getField< fields::well::energyPerforationFlux >()),
    m_dEnergyPerfFlux( perforationData->getField< fields::well::dEnergyPerforationFlux >())

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
  void computeFlux( localIndex const iperf ) const
  {
    Base::computeFlux( iperf, [&] ( globalIndex const & resOffset,
                                    globalIndex const & wellElemOffset,
                                    stackArray1d< globalIndex, 2*resNumDOF > & dofColIndices,
                                    localIndex const iwelem )
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

  /// Well type
  integer const m_isProducer;

  /// Global index of local element
  arrayView1d< globalIndex const >  m_globalWellElementIndex;

  /// Views on energy flux
  arrayView1d< real64 const > const m_energyPerfFlux;
  arrayView3d< real64 const > const m_dEnergyPerfFlux;
};

/**
 * @class ThermalSinglePhaseFluxKernelFactory
 */
class ThermalSinglePhaseFluxKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComps the number of fluid components
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
  createAndLaunch( integer const isProducer,
                   real64 const dt,
                   globalIndex const rankOffset,
                   string const wellDofKey,
                   ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber,
                   WellElementSubRegion const & subRegion,
                   PerforationData const * const perforationData,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    using kernelType = ThermalSinglePhaseFluxKernel< 1 >;
    kernelType kernel( isProducer, dt, rankOffset, wellDofKey, resDofNumber, subRegion, perforationData, localMatrix,
                       localRhs );
    kernelType::template launch< POLICY >( perforationData->size(), kernel );
  }
};

} // end namespace coupledReservoirAndWellKernels

} // end namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALSINGLEPHASEWELLKERNELS_HPP
