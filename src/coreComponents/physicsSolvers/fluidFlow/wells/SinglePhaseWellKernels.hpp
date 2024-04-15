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
 * @file SinglePhaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP

#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/SolverBaseKernels.hpp"

namespace geos
{

namespace singlePhaseWellKernels
{


template< integer IS_THERMAL >
struct ColOffset_WellJac {};


template<>
struct ColOffset_WellJac< 0 >
{
  static constexpr integer dP = 0;
  static constexpr integer dQ = 1;
};
template<>
struct ColOffset_WellJac< 1 >
{
  static constexpr integer dP = 0;
  static constexpr integer dQ = 1;
  static constexpr integer dT = 2;
};


// tag to access well and reservoir elements in perforation rates computation
struct SubRegionTag
{
  static constexpr integer RES  = 0;
  static constexpr integer WELL = 1;
};

// tag to access the next and current well elements of a connection
struct ElemTag
{
  static constexpr integer CURRENT = 0;
  static constexpr integer NEXT    = 1;
};

// define the column offset of the derivatives
struct ColOffset
{
  static constexpr integer DPRES = 0;
  static constexpr integer DRATE = 1;
};

// define the row offset of the residual equations
struct RowOffset
{
  static constexpr integer CONTROL = 0;
  static constexpr integer MASSBAL = 1;
};


/******************************** ControlEquationHelper ********************************/

struct ControlEquationHelper
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  // add an epsilon to the checks to avoid control changes due to tiny pressure/rate updates
  static constexpr real64 EPS = 1e-15;

  GEOS_HOST_DEVICE
  inline
  static
  void
  switchControl( bool const isProducer,
                 WellControls::Control const & currentControl,
                 real64 const & targetBHP,
                 real64 const & targetRate,
                 real64 const & currentBHP,
                 real64 const & currentVolRate,
                 WellControls::Control & newControl );

  GEOS_HOST_DEVICE
  inline
  static
  void
  compute( globalIndex const rankOffset,
           WellControls::Control const currentControl,
           real64 const & targetBHP,
           real64 const & targetRate,
           real64 const & currentBHP,
           real64 const & dCurrentBHP_dPres,
           real64 const & currentVolRate,
           real64 const & dCurrentVolRate_dPres,
           real64 const & dCurrentVolRate_dRate,
           globalIndex const dofNumber,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs );

};


/******************************** FluxKernel ********************************/

struct FluxKernel
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;
  using TAG = singlePhaseWellKernels::ElemTag;

  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};


/******************************** PressureRelationKernel ********************************/

struct PressureRelationKernel
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;
  using TAG = singlePhaseWellKernels::ElemTag;

  static localIndex
  launch( localIndex const size,
          globalIndex const rankOffset,
          bool const isLocallyOwned,
          localIndex const iwelemControl,
          WellControls const & wellControls,
          real64 const & timeAtEndOfStep,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView2d< real64 const > const & wellElemDensity,
          arrayView2d< real64 const > const & dWellElemDensity_dPres,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};


/******************************** PerforationKernel ********************************/

struct PerforationKernel
{

  using TAG = singlePhaseWellKernels::SubRegionTag;

  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::pressure >;

  using SingleFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity_dPressure,
                              fields::singlefluid::viscosity,
                              fields::singlefluid::dViscosity_dPressure >;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  GEOS_HOST_DEVICE
  inline
  static
  void
  compute( real64 const & resPressure,
           real64 const & resDensity,
           real64 const & dResDensity_dPres,
           real64 const & resViscosity,
           real64 const & dResViscosity_dPres,
           real64 const & wellElemGravCoef,
           real64 const & wellElemPressure,
           real64 const & wellElemDensity,
           real64 const & dWellElemDensity_dPres,
           real64 const & wellElemViscosity,
           real64 const & dWellElemViscosity_dPres,
           real64 const & perfGravCoef,
           real64 const & trans,
           real64 & perfRate,
           arraySlice1d< real64 > const & dPerfRate_dPres );

  static void
  launch( localIndex const size,
          ElementViewConst< arrayView1d< real64 const > > const & resPressure,
          ElementViewConst< arrayView2d< real64 const > > const & resDensity,
          ElementViewConst< arrayView2d< real64 const > > const & dResDensity_dPres,
          ElementViewConst< arrayView2d< real64 const > > const & resViscosity,
          ElementViewConst< arrayView2d< real64 const > > const & dResViscosity_dPres,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView2d< real64 const > const & wellElemDensity,
          arrayView2d< real64 const > const & dWellElemDensity_dPres,
          arrayView2d< real64 const > const & wellElemViscosity,
          arrayView2d< real64 const > const & dWellElemViscosity_dPres,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< localIndex const > const & perfWellElemIndex,
          arrayView1d< real64 const > const & perfTransmissibility,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 > const & perfRate,
          arrayView2d< real64 > const & dPerfRate_dPres );

};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemVolume,
          arrayView2d< real64 const > const & wellElemDensity,
          arrayView2d< real64 const > const & dWellElemDensity_dPres,
          arrayView2d< real64 const > const & wellElemDensity_n,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};

/******************************** PressureInitializationKernel ********************************/

struct PresInitializationKernel
{

  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::pressure >;

  using SingleFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density >;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  static void
  launch( localIndex const perforationSize,
          localIndex const subRegionSize,
          localIndex const numPerforations,
          WellControls const & wellControls,
          real64 const & currentTime,
          ElementViewConst< arrayView1d< real64 const > > const & resPressure,
          ElementViewConst< arrayView2d< real64 const > > const & resDensity,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 > const & wellElemPressure );

};

/******************************** RateInitializationKernel ********************************/

struct RateInitializationKernel
{

  static void
  launch( localIndex const subRegionSize,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView2d< real64 const > const & wellElemDens,
          arrayView1d< real64 > const & connRate );

};


/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public solverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = solverBaseKernels::ResidualNormKernelBase< 1 >;
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
                      real64 const timeAtEndOfStep,
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
    m_currentControl( wellControls.getControl() ),
    m_targetBHP( wellControls.getTargetBHP( timeAtEndOfStep ) ),
    m_targetRate( wellControls.getTargetTotalRate( timeAtEndOfStep ) ),
    m_volume( subRegion.getElementVolume() ),
    m_density_n( fluid.density_n() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const iwelem,
                            LinfStackVariables & stack ) const override
  {
    for( localIndex idof = 0; idof < 2; ++idof )
    {
      real64 normalizer = 0.0;
      if( idof == singlePhaseWellKernels::RowOffset::CONTROL )
      {
        // for the top well element, normalize using the current control
        if( m_isLocallyOwned && iwelem == m_iwelemControl )
        {
          if( m_currentControl == WellControls::Control::BHP )
          {
            // this residual entry is in pressure units
            normalizer = m_targetBHP;
          }
          else if( m_currentControl == WellControls::Control::TOTALVOLRATE )
          {
            // this residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetRate ), m_minNormalizer );
          }
        }
        // for the pressure difference equation, always normalize by the BHP
        else
        {
          // this residual is in pressure units
          normalizer = m_targetBHP;
        }
      }
      else // SinglePhaseWell::RowOffset::MASSBAL
      {
        // this residual entry is in mass units
        normalizer = m_dt * LvArray::math::abs( m_targetRate ) * m_density_n[iwelem][0];

        // to make sure that everything still works well if the rate is zero, we add this check
        normalizer = LvArray::math::max( normalizer, m_volume[iwelem] * m_density_n[iwelem][0] );
      }

      // we have the normalizer now, we can compute a dimensionless Linfty norm contribution
      real64 const val = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / normalizer;
      if( val > stack.localValue[0] )
      {
        stack.localValue[0] = val;
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

  /// Controls
  WellControls::Control const m_currentControl;
  real64 const m_targetBHP;
  real64 const m_targetRate;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on total density at the previous converged time step
  arrayView2d< real64 const > const m_density_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the well element subregion
   * @param[in] fluid the fluid model
   * @param[in] wellControls the controls
   * @param[in] timeAtEndOfStep the time at the end of the step (time_n + dt)
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
                   real64 const timeAtEndOfStep,
                   real64 const dt,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, subRegion, fluid, wellControls, timeAtEndOfStep, dt, minNormalizer );
    ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
  }

};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & presDofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & pres,
          real64 const scalingFactor )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, localIndex > minVal( 1 );

    forAll< POLICY >( presDofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 && presDofNumber[ei] >= 0 )
      {
        localIndex const lid = presDofNumber[ei] - rankOffset;
        real64 const newPres = pres[ei] + scalingFactor * localSolution[lid];

        if( newPres < 0.0 )
        {
          minVal.min( 0 );
        }
      }
    } );
    return minVal.get();
  }
};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of accumulation and volume balance
 */
template< integer NUM_DOF >
class ElementBasedAssemblyKernel
{
public:
  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_DOF;

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              constitutive::SingleFluidBase const & fluid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :
    m_rankOffset( rankOffset ),
    m_wellElemDofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_wellElemVolume( subRegion.getElementVolume() ),
    m_wellElemDensity( fluid.density() ),
    m_wellElemDensity_n( fluid.density_n() ),

    m_dWellElemDensity_dPressure( fluid.dDensity_dPressure()  ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    //  volume information (used by both accumulation and volume balance)
    real64 volume = 0.0;
    real64 density = 0.0;
    real64 density_n = 0.0;
    real64 dDensity_dPres = 0.0;


    // Residual information

    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Indices of the matrix rows/columns corresponding to the dofs in this element
    globalIndex dofIndices[numDof]{};
    globalIndex eqnRowIndices[numDof]{};
    globalIndex dofColIndices[numDof]{};

    /// C-array storage for the element local residual vector (all equations )
    real64 localResidual[numDof]{};

    /// C-array storage for the element local Jacobian matrix (all equations  , all dofs)
    real64 localJacobian[numDof][numDof]{};

  };
  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    // initialize the volume
    stack.volume = m_wellElemVolume[ei];
    stack.density = m_wellElemDensity[ei][0];
    stack.density_n = m_wellElemDensity_n[ei][0];
    stack.dDensity_dPres = m_dWellElemDensity_dPressure[ei][0];

    // set row index and degrees of freedom indices for this element (mass + vol bal)
    for( integer ic = 0; ic < numDof; ++ic )
    {
      stack.eqnRowIndices[ic] = m_wellElemDofNumber[ei] +  ic - m_rankOffset;
    }

    // set DOF col indices for this block ( mass + vol bal)
    for( integer idof = 0; idof < numDof; ++idof )
    {
      stack.dofColIndices[idof] = m_wellElemDofNumber[ei]   + idof;
    }

    for( integer jc = 0; jc < numDof; ++jc )
    {
      stack.localResidual[jc] = 0.0;
      for( integer ic = 0; ic < numDof; ++ic )
      {
        stack.localJacobian[jc][ic] = 0.0;
      }

    }

  }


  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] phaseAmountKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const iwelem,
                            StackVariables & stack,
                            FUNC && KernelOp = NoOpFunc{} ) const
  {

    localIndex const eqnRowIndex = m_wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - m_rankOffset;
    globalIndex const presDofColIndex = m_wellElemDofNumber[iwelem] + COFFSET::DPRES;

    stack.localResidual[0] = stack.volume * ( stack.density - stack.density_n );
    stack.localJacobian[0][1] = stack.volume * stack.dDensity_dPres;

    KernelOp();
    // check zero diagonal (works only in debug)
    /*
       for( integer ic = 0; ic < numComp; ++ic )
       {
       GEOS_ASSERT_MSG ( LvArray::math::abs( stack.localJacobian[ic][ic] ) > minDensForDivision,
                        GEOS_FMT( "Zero diagonal in Jacobian: equation {}, value = {}", ic, stack.localJacobian[ic][ic] ) );
       }
     */
  }


  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const GEOS_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {

    // add contribution to residual and jacobian into:
    // - mass bal
    // - pressure eqn
    // - energy if thermal
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels

    for( integer i = 0; i < NUM_DOF; ++i )
    {
      m_localRhs[stack.eqnRowIndices[i]]  += stack.localResidual[i];
      m_localMatrix.template addToRow< serialAtomic >( stack.eqnRowIndices[i],
                                                       stack.dofColIndices,
                                                       stack.localJacobian[i],
                                                       numDof );
    }

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
      typename KERNEL_TYPE::StackVariables stack;
      kernelComponent.setup( iwelem, stack );
      kernelComponent.computeAccumulation( iwelem, stack );
      kernelComponent.complete( iwelem, stack );

    } );
  }

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_wellElemDofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_wellElemVolume;

  /// Views on the densities
  arrayView2d< real64 const > const m_wellElemDensity;
  arrayView2d< real64 const > const m_wellElemDensity_n;
  arrayView2d< real64 const > const m_dWellElemDensity_dPressure;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;


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
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 2;
    ElementBasedAssemblyKernel< NUM_DOF >
    kernel( rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs );
    ElementBasedAssemblyKernel< NUM_DOF >::template
    launch< POLICY, ElementBasedAssemblyKernel< NUM_DOF > >( subRegion.size(), kernel );

  }
};
} // end namespace singlePhaseWellKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP
