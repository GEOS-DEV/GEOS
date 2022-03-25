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
 * @file IsothermalCompositionalMultiphaseBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "functions/TableFunction.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"

namespace geosx
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

using namespace constitutive;

static constexpr real64 minDensForDivision = 1e-10;

/******************************** PropertyKernelBase ********************************/

/**
 * @brief Internal struct to provide no-op defaults used in the inclusion
 *   of lambda functions into kernel component functions.
 * @struct NoOpFunc
 */
struct NoOpFunc
{
  template< typename ... Ts >
  GEOSX_HOST_DEVICE
  constexpr void
  operator()( Ts && ... ) const {}
};

/**
 * @class PropertyKernelBase
 * @tparam NUM_COMP number of fluid components
 * @brief Define the base interface for the property update kernels
 */
template< integer NUM_COMP >
class PropertyKernelBase
{
public:

  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      kernelComponent.compute( ei );
    } );
  }

  /**
   * @brief Performs the kernel launch on a sorted array
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] targetSet the indices of the elements in which we compute the property
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex const ei = targetSet[ i ];
      kernelComponent.compute( ei );
    } );
  }

};

namespace internal
{

template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: type should be integral" );

  switch( value )
  {
    case 1:
    { lambda( std::integral_constant< T, 1 >() ); return; }
    case 2:
    { lambda( std::integral_constant< T, 2 >() ); return; }
    case 3:
    { lambda( std::integral_constant< T, 3 >() ); return; }
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return; }
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return; }
    default:
    { GEOSX_ERROR( "Unsupported number of components: " << value ); }
  }
}

} // namespace internal


/******************************** ComponentFractionKernel ********************************/

/**
 * @class ComponentFractionKernel
 * @tparam NUM_COMP number of fluid components
 * @brief Define the interface for the update kernel in charge of computing the phase volume fractions
 */
template< integer NUM_COMP >
class ComponentFractionKernel : public PropertyKernelBase< NUM_COMP >
{
public:

  using Base = PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  ComponentFractionKernel( ObjectManagerBase & subRegion )
    : Base(),
    m_compDens( subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >() ),
    m_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >() ),
    m_compFrac( subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >() ),
    m_dCompFrac_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity >() )
  {}

  /**
   * @brief Compute the phase volume fractions in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseVolFractionKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && compFractionKernelOp = NoOpFunc{} ) const
  {
    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compDens = m_compDens[ei];
    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const dCompDens = m_dCompDens[ei];
    arraySlice1d< real64, compflow::USD_COMP - 1 > const compFrac = m_compFrac[ei];
    arraySlice2d< real64, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];

    real64 totalDensity = 0.0;

    for( integer ic = 0; ic < numComp; ++ic )
    {
      totalDensity += compDens[ic] + dCompDens[ic];
    }

    real64 const totalDensityInv = 1.0 / totalDensity;

    for( integer ic = 0; ic < numComp; ++ic )
    {
      compFrac[ic] = (compDens[ic] + dCompDens[ic]) * totalDensityInv;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dCompFrac_dCompDens[ic][jc] = -compFrac[ic] * totalDensityInv;
      }
      dCompFrac_dCompDens[ic][ic] += totalDensityInv;
    }

    compFractionKernelOp( compFrac, dCompFrac_dCompDens );
  }

protected:

  // inputs

  // Views on component densities
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;
  arrayView2d< real64 const, compflow::USD_COMP > m_dCompDens;

  // outputs

  // Views on component fraction
  arrayView2d< real64, compflow::USD_COMP > m_compFrac;
  arrayView3d< real64, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

};

/**
 * @class ComponentFractionKernelFactory
 */
class ComponentFractionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   ObjectManagerBase & subRegion )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      ComponentFractionKernel< NUM_COMP > kernel( subRegion );
      ComponentFractionKernel< NUM_COMP >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }

};

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @class PhaseVolumeFractionKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase volume fractions
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseVolumeFractionKernel : public PropertyKernelBase< NUM_COMP >
{
public:

  using Base = PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  PhaseVolumeFractionKernel( ObjectManagerBase & subRegion,
                             MultiFluidBase const & fluid )
    : Base(),
    m_phaseVolFrac( subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac_dPres( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure >() ),
    m_dPhaseVolFrac_dComp( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity >() ),
    m_compDens( subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >() ),
    m_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >() ),
    m_dCompFrac_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseFrac( fluid.phaseFraction() ),
    m_dPhaseFrac( fluid.dPhaseFraction() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() )
  {}

  /**
   * @brief Compute the phase volume fractions in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseVolFractionKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && phaseVolFractionKernelOp = NoOpFunc{} ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compDens = m_compDens[ei];
    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const dCompDens = m_dCompDens[ei];
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const phaseFrac = m_phaseFrac[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const dPhaseFrac = m_dPhaseFrac[ei][0];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres[ei];
    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac_dComp = m_dPhaseVolFrac_dComp[ei];

    real64 work[numComp]{};

    // compute total density from component partial densities
    real64 totalDensity = 0.0;
    real64 const dTotalDens_dCompDens = 1.0;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      totalDensity += compDens[ic] + dCompDens[ic];
    }

    for( integer ip = 0; ip < numPhase; ++ip )
    {

      // set the saturation to zero if the phase is absent
      bool const phaseExists = (phaseFrac[ip] > 0);
      if( !phaseExists )
      {
        phaseVolFrac[ip] = 0.;
        dPhaseVolFrac_dPres[ip] = 0.;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseVolFrac_dComp[ip][jc] = 0.;
        }
        continue;
      }

      // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
      real64 const phaseDensInv = 1.0 / phaseDens[ip];

      // compute saturation and derivatives except multiplying by the total density
      phaseVolFrac[ip] = phaseFrac[ip] * phaseDensInv;

      dPhaseVolFrac_dPres[ip] =
        (dPhaseFrac[ip][Deriv::dP] - phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dP]) * phaseDensInv;

      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseVolFrac_dComp[ip][jc] =
          (dPhaseFrac[ip][Deriv::dC+jc] - phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dC+jc]) * phaseDensInv;
      }

      // apply chain rule to convert derivatives from global component fractions to densities
      applyChainRuleInPlace( numComp, dCompFrac_dCompDens, dPhaseVolFrac_dComp[ip], work );

      // call the lambda in the phase loop to allow the reuse of the phaseVolFrac and totalDensity
      // possible use: assemble the derivatives wrt temperature
      phaseVolFractionKernelOp( ip, phaseVolFrac[ip], phaseDensInv, totalDensity );

      // now finalize the computation by multiplying by total density
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseVolFrac_dComp[ip][jc] *= totalDensity;
        dPhaseVolFrac_dComp[ip][jc] += phaseVolFrac[ip] * dTotalDens_dCompDens;
      }

      phaseVolFrac[ip] *= totalDensity;
      dPhaseVolFrac_dPres[ip] *= totalDensity;
    }
  }

protected:

  // outputs

  /// Views on phase volume fractions
  arrayView2d< real64, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView2d< real64, compflow::USD_PHASE > m_dPhaseVolFrac_dPres;
  arrayView3d< real64, compflow::USD_PHASE_DC > m_dPhaseVolFrac_dComp;

  // inputs

  /// Views on component densities
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;
  arrayView2d< real64 const, compflow::USD_COMP > m_dCompDens;
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

  /// Views on phase fractions
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseFrac;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseFrac;

  /// Views on phase densities
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseDens;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseDens;

};

/**
 * @class PhaseVolumeFractionKernelFactory
 */
class PhaseVolumeFractionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   MultiFluidBase const & fluid )
  {
    if( numPhase == 2 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        PhaseVolumeFractionKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        PhaseVolumeFractionKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};


/******************************** RelativePermeabilityUpdateKernel ********************************/

struct RelativePermeabilityUpdateKernel
{
  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( localIndex const size,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** CapillaryPressureUpdateKernel ********************************/

struct CapillaryPressureUpdateKernel
{
  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( localIndex const size,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of accumulation and volume balance
 */
template< integer NUM_COMP, integer NUM_DOF >
class ElementBasedAssemblyKernel
{
public:

  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;

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
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( localIndex const numPhases,
                              globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              MultiFluidBase const & fluid,
                              CoupledSolidBase const & solid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    : m_numPhases( numPhases ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_porosityOld( solid.getOldPorosity() ),
    m_porosityNew( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    m_dCompFrac_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseVolFracOld( subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFractionOld >() ),
    m_phaseVolFrac( subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac_dPres( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure >() ),
    m_dPhaseVolFrac_dCompDens( subRegion.getExtrinsicData< extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity >() ),
    m_phaseDensOld( fluid.phaseDensityOld() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseCompFracOld( fluid.phaseCompFractionOld() ),
    m_phaseCompFrac( fluid.phaseCompFraction() ),
    m_dPhaseCompFrac( fluid.dPhaseCompFraction() ),
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

    // Pore volume information (used by both accumulation and volume balance)

    /// Pore volume at time n+1
    real64 poreVolumeNew = 0.0;

    /// Pore volume at the previous converged time step
    real64 poreVolumeOld = 0.0;

    /// Derivative of pore volume with respect to pressure
    real64 dPoreVolume_dPres = 0.0;

    // Residual information

    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Indices of the matrix rows/columns corresponding to the dofs in this element
    globalIndex dofIndices[numDof]{};

    /// C-array storage for the element local residual vector (all equations except volume balance)
    real64 localResidual[numEqn]{};

    /// C-array storage for the element local Jacobian matrix (all equations except volume balance, all dofs)
    real64 localJacobian[numEqn][numDof]{};

  };

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOSX_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    // initialize the pore volume
    stack.poreVolumeNew = m_volume[ei] * m_porosityNew[ei][0];
    stack.poreVolumeOld = m_volume[ei] * m_porosityOld[ei][0];
    stack.dPoreVolume_dPres = m_volume[ei] * m_dPoro_dPres[ei][0];

    // set row index and degrees of freedom indices for this element
    stack.localRow = m_dofNumber[ei] - m_rankOffset;
    for( integer idof = 0; idof < numDof; ++idof )
    {
      stack.dofIndices[idof] = m_dofNumber[ei] + idof;
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
  GEOSX_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack,
                            FUNC && phaseAmountKernelOp = NoOpFunc{} ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    // construct the slices for variables accessed multiple times
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFracOld = m_phaseVolFracOld[ei];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[ei];

    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDensOld = m_phaseDensOld[ei][0];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens = m_dPhaseDens[ei][0];

    arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFracOld = m_phaseCompFracOld[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac = m_phaseCompFrac[ei][0];
    arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac = m_dPhaseCompFrac[ei][0];

    // temporary work arrays
    real64 dPhaseAmount_dC[numComp]{};
    real64 dPhaseCompFrac_dC[numComp]{};

    // sum contributions to component accumulation from each phase
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      real64 const phaseAmountNew = stack.poreVolumeNew * phaseVolFrac[ip] * phaseDens[ip];
      real64 const phaseAmountOld = stack.poreVolumeOld * phaseVolFracOld[ip] * phaseDensOld[ip];

      real64 const dPhaseAmount_dP = stack.dPoreVolume_dPres * phaseVolFrac[ip] * phaseDens[ip]
                                     + stack.poreVolumeNew * ( dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                                               + phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dP] );

      // assemble density dependence
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseDens[ip], dPhaseAmount_dC, Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                              + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
        dPhaseAmount_dC[jc] *= stack.poreVolumeNew;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for( integer ic = 0; ic < numComp; ++ic )
      {
        real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ip][ic];
        real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ip][ic];

        real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                           + phaseAmountNew * dPhaseCompFrac[ip][ic][Deriv::dP];

        stack.localResidual[ic] += phaseCompAmountNew - phaseCompAmountOld;
        stack.localJacobian[ic][0] += dPhaseCompAmount_dP;

        // jc - index of component w.r.t. whose compositional var the derivative is being taken
        // (i.e. col number in local matrix)

        // assemble phase composition dependence
        applyChainRule( numComp, dCompFrac_dCompDens, dPhaseCompFrac[ip][ic], dPhaseCompFrac_dC, Deriv::dC );
        for( integer jc = 0; jc < numComp; ++jc )
        {
          real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                             + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
          stack.localJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
        }

      }

      // call the lambda in the phase loop to allow the reuse of the phase amounts and their derivatives
      // possible use: assemble the derivatives wrt temperature, and the accumulation term of the energy equation for this phase
      phaseAmountKernelOp( ip, phaseAmountNew, phaseAmountOld, dPhaseAmount_dP, dPhaseAmount_dC );

    }
  }

  /**
   * @brief Compute the local volume balance contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] phaseVolFractionSumKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeVolumeBalance( localIndex const ei,
                             StackVariables & stack,
                             FUNC && phaseVolFractionSumKernelOp = NoOpFunc{} ) const
  {
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[ei];

    real64 oneMinusPhaseVolFracSum = 1.0;

    // sum contributions to component accumulation from each phase
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      oneMinusPhaseVolFracSum -= phaseVolFrac[ip];
      stack.localJacobian[numComp][0] -= dPhaseVolFrac_dPres[ip];

      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localJacobian[numComp][jc+1] -= dPhaseVolFrac_dCompDens[ip][jc];
      }
    }

    // call the lambda in the phase loop to allow the reuse of the phase amounts and their derivatives
    // possible use: assemble the derivatives wrt temperature, and use oneMinusPhaseVolFracSum if poreVolumeNew depends on temperature
    phaseVolFractionSumKernelOp( oneMinusPhaseVolFracSum );

    // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
    stack.localResidual[numComp] = stack.poreVolumeNew * oneMinusPhaseVolFracSum;
    for( integer idof = 0; idof < numDof; ++idof )
    {
      stack.localJacobian[numComp][idof] *= stack.poreVolumeNew;
    }
    stack.localJacobian[numComp][0] += stack.dPoreVolume_dPres * oneMinusPhaseVolFracSum;
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void complete( localIndex const GEOSX_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {
    using namespace compositionalMultiphaseUtilities;

    // apply equation/variable change transformation to the component mass balance equations
    real64 work[numDof]{};
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numDof, stack.localJacobian, work );
    shiftElementsAheadByOneAndReplaceFirstElementWithSum( numComp, stack.localResidual );

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numComp-1)
    // - the volume balance equations (i = numComp)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    for( integer i = 0; i < numComp+1; ++i )
    {
      m_localRhs[stack.localRow + i] += stack.localResidual[i];
      m_localMatrix.addToRow< serialAtomic >( stack.localRow + i,
                                              stack.dofIndices,
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
    GEOSX_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.elemGhostRank( ei ) >= 0 )
      {
        return;
      }

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.computeAccumulation( ei, stack );
      kernelComponent.computeVolumeBalance( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosityOld;
  arrayView2d< real64 const > const m_porosityNew;
  arrayView2d< real64 const > const m_dPoro_dPres;

  /// Views on the derivatives of comp fractions wrt component density
  arrayView3d< real64 const, compflow::USD_COMP_DC > const m_dCompFrac_dCompDens;

  /// Views on the phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFracOld;
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFrac;
  arrayView2d< real64 const, compflow::USD_PHASE > const m_dPhaseVolFrac_dPres;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > const m_dPhaseVolFrac_dCompDens;

  /// Views on the phase densities
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseDensOld;
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseDens;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const m_dPhaseDens;

  /// Views on the phase component fraction
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const m_phaseCompFracOld;
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const m_phaseCompFrac;
  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const m_dPhaseCompFrac;

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
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   MultiFluidBase const & fluid,
                   CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC()+1;
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
      kernel( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }

};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{

  template< typename POLICY, typename REDUCE_POLICY >
  static void launch( arrayView1d< real64 const > const & localResidual,
                      globalIndex const rankOffset,
                      integer const numComponents,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< integer const > const & ghostRank,
                      arrayView1d< real64 const > const & refPoro,
                      arrayView1d< real64 const > const & volume,
                      arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDensOld,
                      real64 & localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > localSum( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        real64 const normalizer = totalDensOld[ei][0] * refPoro[ei] * volume[ei];

        for( integer idof = 0; idof < numComponents + 1; ++idof )
        {
          real64 const val = localResidual[localRow + idof] / normalizer;
          localSum += val * val;
        }
      }
    } );
    localResidualNorm += localSum.get();
  }

};


/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY, typename REDUCE_POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          localIndex const numComponents,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          integer const allowCompDensChopping,
          real64 const scalingFactor )
  {
    real64 constexpr eps = minDensForDivision;

    RAJA::ReduceMin< REDUCE_POLICY, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        {
          real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[localRow];
          check.min( newPres >= 0.0 );
        }

        // if component density chopping is not allowed, the time step fails if a component density is negative
        // otherwise, we just check that the total density is positive, and negative component densities
        // will be chopped (i.e., set to zero) in ApplySystemSolution)
        if( !allowCompDensChopping )
        {
          for( integer ic = 0; ic < numComponents; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            check.min( newDens >= 0.0 );
          }
        }
        else
        {
          real64 totalDens = 0.0;
          for( integer ic = 0; ic < numComponents; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            totalDens += (newDens > 0.0) ? newDens : 0.0;
          }
          check.min( totalDens >= eps );
        }
      }
    } );
    return check.get();
  }

};

/******************************** HydrostaticPressureKernel ********************************/

struct HydrostaticPressureKernel
{

  // TODO: this type of constants should be centralized somewhere or provided by fluid model
  static real64 constexpr MIN_FOR_PHASE_PRESENCE = 1e-12;

  enum class ReturnType : integer
  {
    FAILED_TO_CONVERGE = 0,
    DETECTED_MULTIPHASE_FLOW = 1,
    SUCCESS = 2
  };

  template< typename FLUID_WRAPPER >
  static ReturnType
  computeHydrostaticPressure( integer const numComps,
                              integer const numPhases,
                              integer const ipInit,
                              integer const maxNumEquilIterations,
                              real64 const & equilTolerance,
                              real64 const (&gravVector)[ 3 ],
                              FLUID_WRAPPER fluidWrapper,
                              arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                              TableFunction::KernelWrapper tempTableWrapper,
                              real64 const & refElevation,
                              real64 const & refPres,
                              arraySlice1d< real64 const > const & refPhaseMassDens,
                              real64 const & newElevation,
                              real64 & newPres,
                              arraySlice1d< real64 > const & newPhaseMassDens )
  {
    // fluid properties at this elevation
    StackArray< real64, 2, constitutive::MultiFluidBase::MAX_NUM_COMPONENTS, compflow::LAYOUT_COMP > compFrac( 1, numComps );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseFrac( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseDens( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseMassDens( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseVisc( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseEnthalpy( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, multifluid::LAYOUT_PHASE > phaseInternalEnergy( 1, 1, numPhases );
    StackArray< real64, 4, constitutive::MultiFluidBase::MAX_NUM_PHASES *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS,
                multifluid::LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhases, numComps );
    real64 totalDens = 0.0;

    bool isSinglePhaseFlow = true;

    // Step 1: compute the hydrostatic pressure at the current elevation

    real64 const gravCoef = gravVector[2] * ( refElevation - newElevation );
    real64 const temp = tempTableWrapper.compute( &newElevation );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      compFrac[0][ic] = compFracTableWrappers[ic].compute( &newElevation );
    }

    // Step 2: guess the pressure with the refPhaseMassDensity

    real64 pres0 = refPres - refPhaseMassDens[ipInit] * gravCoef;
    real64 pres1 = 0.0;

    // Step 3: compute the mass density at this elevation using the guess, and update pressure

    fluidWrapper.compute( pres0,
                          temp,
                          compFrac[0],
                          phaseFrac[0][0],
                          phaseDens[0][0],
                          phaseMassDens[0][0],
                          phaseVisc[0][0],
                          phaseEnthalpy[0][0],
                          phaseInternalEnergy[0][0],
                          phaseCompFrac[0][0],
                          totalDens );
    pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + phaseMassDens[0][0][ipInit] ) * gravCoef;

    // Step 4: fixed-point iteration until convergence

    bool equilHasConverged = false;
    for( integer eqIter = 0; eqIter < maxNumEquilIterations; ++eqIter )
    {

      // check convergence
      equilHasConverged = ( LvArray::math::abs( pres0 - pres1 ) < equilTolerance );
      pres0 = pres1;

      // if converged, check number of phases and move on
      if( equilHasConverged )
      {
        // make sure that the fluid is single-phase, other we have to issue a warning (for now)
        // if only one phase is mobile, we are in good shape (unfortunately it is hard to access relperm from here)
        localIndex numberOfPhases = 0;
        for( integer ip = 0; ip < numPhases; ++ip )
        {
          if( phaseFrac[0][0][ip] > MIN_FOR_PHASE_PRESENCE )
          {
            numberOfPhases++;
          }
        }
        if( numberOfPhases > 1 )
        {
          isSinglePhaseFlow = false;
        }

        break;
      }

      // compute the mass density at this elevation using the previous pressure, and compute the new pressure
      fluidWrapper.compute( pres0,
                            temp,
                            compFrac[0],
                            phaseFrac[0][0],
                            phaseDens[0][0],
                            phaseMassDens[0][0],
                            phaseVisc[0][0],
                            phaseEnthalpy[0][0],
                            phaseInternalEnergy[0][0],
                            phaseCompFrac[0][0],
                            totalDens );
      pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + phaseMassDens[0][0][ipInit] ) * gravCoef;
    }

    // Step 5: save the hydrostatic pressure and the corresponding density

    newPres = pres1;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      newPhaseMassDens[ip] = phaseMassDens[0][0][ip];
    }

    if( !equilHasConverged )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( !isSinglePhaseFlow )
    {
      return ReturnType::DETECTED_MULTIPHASE_FLOW;
    }
    else
    {
      return ReturnType::SUCCESS;
    }
  }

  template< typename FLUID_WRAPPER >
  static ReturnType
  launch( localIndex const size,
          integer const numComps,
          integer const numPhases,
          integer const ipInit,
          integer const maxNumEquilIterations,
          real64 const equilTolerance,
          real64 const (&gravVector)[ 3 ],
          real64 const & minElevation,
          real64 const & elevationIncrement,
          real64 const & datumElevation,
          real64 const & datumPres,
          FLUID_WRAPPER fluidWrapper,
          arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
          TableFunction::KernelWrapper tempTableWrapper,
          arrayView1d< arrayView1d< real64 > const > elevationValues,
          arrayView1d< real64 > pressureValues )
  {

    ReturnType returnVal = ReturnType::SUCCESS;

    // Step 1: compute the phase mass densities at datum

    // datum fluid properties
    array2d< real64, compflow::LAYOUT_COMP > datumCompFrac( 1, numComps );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseFrac( 1, 1, numPhases );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseDens( 1, 1, numPhases );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseMassDens( 1, 1, numPhases );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseVisc( 1, 1, numPhases );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseEnthalpy( 1, 1, numPhases );
    array3d< real64, multifluid::LAYOUT_PHASE > datumPhaseInternalEnergy( 1, 1, numPhases );
    array4d< real64, multifluid::LAYOUT_PHASE_COMP > datumPhaseCompFrac( 1, 1, numPhases, numComps );
    real64 datumTotalDens = 0.0;

    real64 const datumTemp = tempTableWrapper.compute( &datumElevation );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      datumCompFrac[0][ic] = compFracTableWrappers[ic].compute( &datumElevation );
    }
    fluidWrapper.compute( datumPres,
                          datumTemp,
                          datumCompFrac[0],
                          datumPhaseFrac[0][0],
                          datumPhaseDens[0][0],
                          datumPhaseMassDens[0][0],
                          datumPhaseVisc[0][0],
                          datumPhaseEnthalpy[0][0],
                          datumPhaseInternalEnergy[0][0],
                          datumPhaseCompFrac[0][0],
                          datumTotalDens );

    // Step 2: find the closest elevation to datumElevation

    forAll< parallelHostPolicy >( size, [=] ( localIndex const i )
    {
      real64 const elevation = minElevation + i * elevationIncrement;
      elevationValues[0][i] = elevation;
    } );
    integer const iRef = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                 elevationValues[0].size(),
                                                                 datumElevation );

    // Step 3: compute the mass density and pressure at the reference elevation

    array2d< real64 > phaseMassDens( pressureValues.size(), numPhases );
    // temporary array without permutation to compile on Lassen
    array1d< real64 > datumPhaseMassDensTmp( numPhases );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      datumPhaseMassDensTmp[ip] = datumPhaseMassDens[0][0][ip];
    }

    ReturnType const refReturnVal =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ipInit,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  datumElevation,
                                  datumPres,
                                  datumPhaseMassDensTmp,
                                  elevationValues[0][iRef],
                                  pressureValues[iRef],
                                  phaseMassDens[iRef] );
    if( refReturnVal == ReturnType::FAILED_TO_CONVERGE )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( refReturnVal == ReturnType::DETECTED_MULTIPHASE_FLOW )
    {
      returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
    }

    // Step 4: for each elevation above the reference elevation, compute the pressure

    localIndex const numEntriesAboveRef = size - iRef - 1;
    forAll< serialPolicy >( numEntriesAboveRef, [=, &returnVal] ( localIndex const i )
    {
      ReturnType const returnValAboveRef =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ipInit,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][iRef+i],
                                    pressureValues[iRef+i],
                                    phaseMassDens[iRef+i],
                                    elevationValues[0][iRef+i+1],
                                    pressureValues[iRef+i+1],
                                    phaseMassDens[iRef+i+1] );
      if( returnValAboveRef == ReturnType::FAILED_TO_CONVERGE )
      {
        returnVal = ReturnType::FAILED_TO_CONVERGE;
      }
      else if( ( returnValAboveRef == ReturnType::DETECTED_MULTIPHASE_FLOW ) &&
               ( returnVal != ReturnType::FAILED_TO_CONVERGE ) )
      {
        returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
      }

    } );

    // Step 5: for each elevation below the reference elevation, compute the pressure

    localIndex const numEntriesBelowRef = iRef;
    forAll< serialPolicy >( numEntriesBelowRef, [=, &returnVal] ( localIndex const i )
    {
      ReturnType const returnValBelowRef =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ipInit,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][iRef-i],
                                    pressureValues[iRef-i],
                                    phaseMassDens[iRef-i],
                                    elevationValues[0][iRef-i-1],
                                    pressureValues[iRef-i-1],
                                    phaseMassDens[iRef-i-1] );
      if( returnValBelowRef == ReturnType::FAILED_TO_CONVERGE )
      {
        returnVal = ReturnType::FAILED_TO_CONVERGE;
      }
      else if( ( returnValBelowRef == ReturnType::DETECTED_MULTIPHASE_FLOW ) &&
               ( returnVal != ReturnType::FAILED_TO_CONVERGE ) )
      {
        returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
      }

    } );

    return returnVal;
  }

};


/******************************** Kernel launch machinery ********************************/

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector1( integer const numComp, ARGS && ... args )
{
  internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    KERNELWRAPPER::template launch< NC() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector2( integer const numComp, integer const numPhase, ARGS && ... args )
{
  // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
  if( numPhase == 2 )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 2 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else if( numPhase == 3 )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 3 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else
  {
    GEOSX_ERROR( "Unsupported number of phases: " << numPhase );
  }
}

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_ISOTHERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP
