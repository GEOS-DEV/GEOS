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
 * @file DARTSSuperEngineKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_DARTSSUPERENGINEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_DARTSSUPERENGINEKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "functions/MultivariableTableFunctionKernels.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/DARTSSuperEngineExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"



namespace geosx
{

namespace DARTSSuperEngineKernels
{

using namespace constitutive;

static constexpr real64 minDensForDivision = 1e-10;


// The number of operators in use depends on:
// 1. Number of phases
// 2. Number of components
// 3. Features required in simulation (now its only one - ENergy balance)
// This number needs to be used in solver and in kernels (as a template parameter)
// IMHO, this number is too big ( order of 10-100) to be treated via a kernelLaunchSelectorSwitch construct
// Hence, a way to define it once and for all is needed.
// Could be constexpr member of solver, but passing constexpr lambdas require -std=c++17

constexpr integer COMPUTE_NUM_OPS ( integer const NP, integer const NC, bool ENERGY )
{
  auto DOF = NC + ENERGY;
  return DOF /*accumulation*/ +
         DOF * NP /*flux*/ +
         NP /*up_constant*/ +
         DOF * NP /*gradient*/ +
         DOF /*kinetic rate*/ +
         2 /*rock internal energy and conduction*/ +
         2 * NP /*gravity and capillarity*/ +
         1 /*rock porosity*/ +
         1;
}



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
  static constexpr integer numComps = NUM_COMP;

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

// template< integer NUM_COMP, typename LAMBDA >
// void kernelLaunchSelectorOpSwitch( T value, LAMBDA && lambda )

template< bool ENABLE_ENERGY, integer NUM_PHASES, typename T, typename LAMBDA >
void kernelLaunchSelectorCompSwitch( T numComps, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: type should be integral" );

  switch( numComps )
  {
    case 1:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 1 >(), std::integral_constant< bool, ENABLE_ENERGY >() ); return; }
    case 2:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 2 >(), std::integral_constant< bool, ENABLE_ENERGY >() ); return; }
    case 3:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 3 >(), std::integral_constant< bool, ENABLE_ENERGY >() ); return; }
    case 4:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 4 >(), std::integral_constant< bool, ENABLE_ENERGY >() ); return; }
    case 5:
    { lambda( std::integral_constant< T, NUM_PHASES >(), std::integral_constant< T, 5 >(), std::integral_constant< bool, ENABLE_ENERGY >()); return; }
    default:
    { GEOSX_ERROR( "Unsupported number of components: " << numComps ); }
  }
}

template< bool ENABLE_ENERGY, typename T, typename LAMBDA >
void kernelLaunchSelectorPhaseSwitch( T numPhases, T numComps, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorPhaseSwitch: type should be integral" );

  switch( numPhases )
  {
    case 1:
    { kernelLaunchSelectorCompSwitch< ENABLE_ENERGY, 1 >( numComps, lambda ); return; }
    case 2:
    { kernelLaunchSelectorCompSwitch< ENABLE_ENERGY, 2 >( numComps, lambda ); return; }
    case 3:
    { kernelLaunchSelectorCompSwitch< ENABLE_ENERGY, 3 >( numComps, lambda ); return; }
    default:
    { GEOSX_ERROR( "Unsupported number of phases: " << numPhases ); }
  }
}

template< typename T, typename LAMBDA >
void kernelLaunchSelectorEnergySwitch( T numPhases, T numComps, bool enableEnergyBalance, LAMBDA && lambda )
{
  if( enableEnergyBalance )
  {
    kernelLaunchSelectorPhaseSwitch< true >( numPhases, numComps, lambda );
  }
  else
  {
    kernelLaunchSelectorPhaseSwitch< false >( numPhases, numComps, lambda );
  }
}

} // namespace internal


/******************************** OBLOperatorsKernel ********************************/

/**
 * @class OBLOperatorsKernel
 * @tparam NUM_DIMS number of degrees of freedom
 * @tparam NUM_OPS number of degrees of freedom
 * @brief Compute OBL Operators and derivatives for specified region
 */
template< integer NUM_PHASES, integer NUM_COMPS, bool ENABLE_ENERGY >
class OBLOperatorsKernel
{
public:
  /// Compile time value for the energy balance switch
  static constexpr integer enableEnergyBalance = ENABLE_ENERGY;
  /// Compile time value for the number of components
  static constexpr integer numComps = NUM_COMPS;

  /// Compile time value for the number of dimensions
  static constexpr integer numDofs = numComps + enableEnergyBalance;
  // /// Compile time value for the number of operators
  static constexpr integer numOps = COMPUTE_NUM_OPS( NUM_PHASES, NUM_COMPS, ENABLE_ENERGY );

  static constexpr real64 barToPascalMult = 1e5;
  static constexpr real64 pascalToBarMult = 1 / 1e5;

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
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  OBLOperatorsKernel( ObjectManagerBase & subRegion,
                      MultivariableTableFunctionStaticKernel< numDofs, numOps > OBLOperatorsTable )
    :
    m_OBLOperatorsTable( OBLOperatorsTable ),
    m_pressure( subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >() ),
    m_compFrac( subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >() ),
    m_temperature( subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature >() ),
    m_dPressure( subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >() ),
    m_dCompFrac( subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompFraction >() ),
    m_dTemperature( subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaTemperature >() ),
    m_OBLOperatorValues ( subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValues >()),
    m_OBLOperatorDerivatives ( subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorDerivatives >())
  {}

  /**
   * @brief Compute the phase volume fractions in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseVolFractionKernelOp the function used to customize the kernel
   */
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei ) const
  {
    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compFrac = m_compFrac[ei];
    arraySlice1d< real64 const, compflow::USD_COMP - 1 > const dCompFrac = m_dCompFrac[ei];
    arraySlice1d< real64, compflow::USD_OBL_VAL - 1 > const & OBLVals = m_OBLOperatorValues[ei];
    arraySlice2d< real64, compflow::USD_OBL_DER - 1 > const & OBLDers = m_OBLOperatorDerivatives[ei];
    real64 state[numDofs];

    // we need to convert pressure from Pa (internal unit in GEOSX) to bar (internal unit in DARTS)
    state[0] = (m_pressure[ei] + m_dPressure[ei]) * pascalToBarMult;

    for( integer i = 1; i < numComps - 1; i++ )
    {
      state[i] = compFrac[i - 1] + dCompFrac[i - 1];
      //printf ( " %lf,", state[i] );
    }

    if( enableEnergyBalance )
    {
      state[numDofs - 1] = (m_temperature[ei] + m_dTemperature[ei]);
    }

    m_OBLOperatorsTable.compute( state, OBLVals, OBLDers );

    // perform pressure unit conversion back
    for( integer i = 0; i < numOps; i++ )
    {
      OBLDers[i][0] *= barToPascalMult;
    }

    if( ei == 0 )
    {
      printf ( "\n Ops after:\n" );
      for( integer i = 0; i < numOps; i++ )
      {
        printf ( "%lf [", OBLVals[i] );
        for( integer j = 0; j < numDofs; j++ )
        {
          printf ( " %lf,", OBLDers[i][j] );
        }
        printf ( "]\n" );
      }
    }
  }

protected:

  // inputs
  MultivariableTableFunctionStaticKernel< numDofs, numOps > m_OBLOperatorsTable;

  // Views on primary variables and their updates
  arrayView1d< real64 const > m_pressure;
  arrayView2d< real64 const, compflow::USD_COMP > m_compFrac;
  arrayView1d< real64 const > m_temperature;

  arrayView1d< real64 const > m_dPressure;
  arrayView2d< real64 const, compflow::USD_COMP > m_dCompFrac;
  arrayView1d< real64 const > m_dTemperature;

  // outputs

  // Views on component fraction
  arrayView2d< real64, compflow::USD_OBL_VAL > m_OBLOperatorValues;
  arrayView3d< real64, compflow::USD_OBL_DER > m_OBLOperatorDerivatives;


};

/**
 * @class ComponentFractionKernelFactory
 */
class OBLOperatorsKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhases the number of phases
   * @param[in] numDofs the number of degrees of freedom
   * @param[in] subRegion the element subregion
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numPhases,
                   integer const numComponents,
                   bool const enableEnergyBalance,
                   ObjectManagerBase & subRegion,
                   MultivariableTableFunction const & function )
  {
    internal::kernelLaunchSelectorEnergySwitch( numPhases, numComponents, enableEnergyBalance, [&] ( auto NP, auto NC, auto E )
    {
      integer constexpr ENABLE_ENERGY = E();
      integer constexpr NUM_PHASES = NP();
      integer constexpr NUM_COMPS = NC();
      integer constexpr NUM_DIMS = ENABLE_ENERGY + NUM_COMPS;
      integer constexpr NUM_OPS  = COMPUTE_NUM_OPS( NUM_PHASES, NUM_COMPS, ENABLE_ENERGY );

      OBLOperatorsKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >
      kernel( subRegion,
              MultivariableTableFunctionStaticKernel< NUM_DIMS, NUM_OPS >( function.getAxisMinimums(),
                                                                           function.getAxisMaximums(),
                                                                           function.getAxisPoints(),
                                                                           function.getAxisSteps(),
                                                                           function.getAxisStepInvs(),
                                                                           function.getAxisHypercubeMults(),
                                                                           function.getHypercubeData()
                                                                           ) );
      OBLOperatorsKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >::template launch< POLICY >( subRegion.size(), kernel );
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
template< integer NUM_PHASES, integer NUM_COMPS, bool ENABLE_ENERGY >
class ElementBasedAssemblyKernel
{
public:

  /// Compile time value for the number of phases
  static constexpr integer numPhases = NUM_PHASES;
  /// Compile time value for the number of components
  static constexpr integer numComps = NUM_COMPS;
  /// Compile time value for the energy balance switch
  static constexpr integer enableEnergyBalance = ENABLE_ENERGY;

  /// Compile time value for the number of dimensions
  static constexpr integer numDofs = numComps + enableEnergyBalance;
  // /// Compile time value for the number of operators
  static constexpr integer numOps = COMPUTE_NUM_OPS( NUM_PHASES, NUM_COMPS, ENABLE_ENERGY );

  // order of operators:
  static constexpr integer ACC_OP = 0;
  static constexpr integer FLUX_OP = numDofs;
  // diffusion
  static constexpr integer UPSAT_OP = numDofs + numDofs * numPhases;
  static constexpr integer GRAD_OP = numDofs + numDofs * numPhases + numPhases;
  // kinumDofstic reaction
  static constexpr integer KIN_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases;


  static constexpr integer RE_INTER_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs;
  static constexpr integer RE_TEMP_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 1;
  static constexpr integer ROCK_COND = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 2;
  static constexpr integer GRAV_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 3;
  static constexpr integer PC_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 3 + numPhases;
  static constexpr integer PORO_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 3 + 2 * numPhases;

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
  ElementBasedAssemblyKernel( real64 const dt,
                              globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              CoupledSolidBase const & solid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :
    m_dt( dt ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_pressure ( subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >() ),
    m_compFrac ( subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >() ),
    m_temperature ( subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature >() ),
    m_dPressure( subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >() ),
    m_dCompFrac( subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompFraction >() ),
    m_dTemperature( subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaTemperature >() ),
    m_referencePoreVolume( subRegion.getExtrinsicData< extrinsicMeshData::flow::referencePoreVolume >() ),
    m_referenceRockVolume( subRegion.getExtrinsicData< extrinsicMeshData::flow::referenceRockVolume >() ),
    m_rockVolumetricHeatCapacity( subRegion.getExtrinsicData< extrinsicMeshData::flow::rockVolumetricHeatCapacity >() ),
    m_rockThermalConductivity( subRegion.getExtrinsicData< extrinsicMeshData::flow::rockThermalConductivity >() ),
    m_rockKineticRateFactor( subRegion.getExtrinsicData< extrinsicMeshData::flow::rockKineticRateFactor >() ),
    m_OBLOperatorValues ( subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValues >()),
    m_OBLOperatorValuesOld ( subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValuesOld >()),
    m_OBLOperatorDerivatives ( subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorDerivatives >()),
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

    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Indices of the matrix rows/columns corresponding to the dofs in this element
    globalIndex dofIndices[numDofs]{};

    /// C-array storage for the element local residual vector (all equations except volume balance)
    real64 localResidual[numDofs]{};

    /// C-array storage for the element local Jacobian matrix (all equations except volume balance, all dofs)
    real64 localJacobian[numDofs][numDofs]{};

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

    // set row index and degrees of freedom indices for this element
    stack.localRow = m_dofNumber[ei] - m_rankOffset;

    for( integer idof = 0; idof < numDofs; ++idof )
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
  GEOSX_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLVals = m_OBLOperatorValues[ei];
    arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLValsOld = m_OBLOperatorValuesOld[ei];
    arraySlice2d< real64 const, compflow::USD_OBL_DER - 1 > const & OBLDers = m_OBLOperatorDerivatives[ei];

    // [1] fill diagonal part for both mass (and energy equations if needed, only fluid energy is involved here)
    for( integer c = 0; c < numDofs; c++ )
    {
      stack.localResidual[c] = m_referencePoreVolume[ei] * (OBLVals[ACC_OP + c] - OBLValsOld[ACC_OP + c]);   // acc operators
      // only

      // Add reaction term to diagonal of reservoir cells (here the volume is pore volume or block volume):
      stack.localResidual[c] += (m_referencePoreVolume[ei] + m_referenceRockVolume[ei]) * m_dt * OBLVals[KIN_OP + c] * m_rockKineticRateFactor[ei]; // kinetics

      for( integer v = 0; v < numDofs; v++ )
      {
        stack.localJacobian[c][v] = m_referencePoreVolume[ei] * OBLDers[ACC_OP + c][v];   // der of accumulation term

        // Include derivatives for reaction term if part of reservoir cells:
        stack.localJacobian[c][v] += (m_referencePoreVolume[ei] + m_referenceRockVolume[ei]) * m_dt * OBLDers[KIN_OP + c][v] * m_rockKineticRateFactor[ei];     // derivative
      }
    }

    // + rock energy
    if( enableEnergyBalance )
    {
      stack.localResidual[numDofs-1] += m_referenceRockVolume[ei] * (OBLVals[RE_INTER_OP] - OBLValsOld[RE_INTER_OP]) * m_rockVolumetricHeatCapacity[ei];

      for( integer v = 0; v < numDofs; v++ )
      {
        stack.localJacobian[numDofs-1][v] += m_referenceRockVolume[ei] * OBLDers[RE_INTER_OP][v] * m_rockVolumetricHeatCapacity[ei];
      }   // end of fill offdiagonal part + contribute to diagonal
    }
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
    using namespace CompositionalMultiphaseUtilities;

    // add contribution to residual and jacobian into:
    // - the component mass balance equations
    // - the volume balance equations
    for( integer i = 0; i < numDofs; ++i )
    {
      m_localRhs[stack.localRow + i] += stack.localResidual[i];
      m_localMatrix.addToRow< serialAtomic >( stack.localRow + i,
                                              stack.dofIndices,
                                              stack.localJacobian[i],
                                              numDofs );
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
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

// Time step size
  real64 m_dt;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;


  // Views on primary variables and their updates
  arrayView1d< real64 const > const m_pressure;
  arrayView2d< real64 const, compflow::USD_COMP > const m_compFrac;
  arrayView1d< real64 const > const m_temperature;

  arrayView1d< real64 const > const m_dPressure;
  arrayView2d< real64 const, compflow::USD_COMP > const m_dCompFrac;
  arrayView1d< real64 const > const m_dTemperature;

  // views on solid properties
  arrayView1d< real64 const > const m_referencePoreVolume;
  arrayView1d< real64 const > const m_referenceRockVolume;
  arrayView1d< real64 const > const m_rockVolumetricHeatCapacity;
  arrayView1d< real64 const > const m_rockThermalConductivity;
  arrayView1d< real64 const > const m_rockKineticRateFactor;


  // Views on OBL operators and their derivatives
  arrayView2d< real64 const, compflow::USD_OBL_VAL > const m_OBLOperatorValues;
  arrayView2d< real64 const, compflow::USD_OBL_VAL > const m_OBLOperatorValuesOld;
  arrayView3d< real64 const, compflow::USD_OBL_DER > const m_OBLOperatorDerivatives;

  // outputs

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
  createAndLaunch(
    integer const numPhases,
    integer const numComps,
    bool const enableEnergyBalance,
    real64 const dt,
    globalIndex const rankOffset,
    string const dofKey,
    ElementSubRegionBase const & subRegion,
    CoupledSolidBase const & solid,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs )
  {
    internal::kernelLaunchSelectorEnergySwitch( numPhases, numComps, enableEnergyBalance, [&] ( auto NP, auto NC, auto E )
    {
      integer constexpr ENABLE_ENERGY = E();
      integer constexpr NUM_PHASES = NP();
      integer constexpr NUM_COMPS = NC();

      ElementBasedAssemblyKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >
      kernel( dt, rankOffset, dofKey, subRegion, solid, localMatrix, localRhs );
      ElementBasedAssemblyKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY >::template launch< POLICY >( subRegion.size(), kernel );
    } );

  }

};

/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @brief Base class for FaceBasedAssemblyKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of components/dofs).
 */
class FaceBasedAssemblyKernelBase
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

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using CompFlowAccessors =
    StencilAccessors< extrinsicMeshData::ghostRank,
                      extrinsicMeshData::flow::gravityCoefficient,
                      extrinsicMeshData::flow::pressure,
                      extrinsicMeshData::flow::deltaPressure,
                      extrinsicMeshData::flow::globalCompFraction,
                      extrinsicMeshData::flow::deltaGlobalCompFraction,
                      extrinsicMeshData::flow::temperature,
                      extrinsicMeshData::flow::deltaTemperature,
                      extrinsicMeshData::flow::referencePorosity,
                      extrinsicMeshData::flow::rockThermalConductivity,
                      extrinsicMeshData::flow::OBLOperatorValues,
                      extrinsicMeshData::flow::OBLOperatorDerivatives >;

  using PermeabilityAccessors =
    StencilAccessors< extrinsicMeshData::permeability::permeability,
                      extrinsicMeshData::permeability::dPerm_dPressure >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] capPressureFlag flag specifying whether capillary pressure is used or not
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernelBase( globalIndex const rankOffset,
                               DofNumberAccessor const & dofNumberAccessor,
                               CompFlowAccessors const & compFlowAccessors,
                               PermeabilityAccessors const & permeabilityAccessors,
                               real64 const & dt,
                               real64 const & transMultExp,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
    :     m_rankOffset( rankOffset ),
    m_dt( dt ),
    m_transMultExp (),
    m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_permeability( permeabilityAccessors.get( extrinsicMeshData::permeability::permeability {} ) ),
    m_dPerm_dPres( permeabilityAccessors.get( extrinsicMeshData::permeability::dPerm_dPressure {} ) ),
    m_referencePorosity( compFlowAccessors.get( extrinsicMeshData::flow::referencePorosity {} ) ),
    m_rockThermalConductivity( compFlowAccessors.get( extrinsicMeshData::flow::rockThermalConductivity {} ) ),
    m_ghostRank( compFlowAccessors.get( extrinsicMeshData::ghostRank {} ) ),
    m_gravCoef( compFlowAccessors.get( extrinsicMeshData::flow::gravityCoefficient {} ) ),
    m_pres( compFlowAccessors.get( extrinsicMeshData::flow::pressure {} ) ),
    m_dPres( compFlowAccessors.get( extrinsicMeshData::flow::deltaPressure {} ) ),
    m_compFrac( compFlowAccessors.get( extrinsicMeshData::flow::globalCompFraction {} ) ),
    m_dCompFrac( compFlowAccessors.get( extrinsicMeshData::flow::deltaGlobalCompFraction {} ) ),
    m_temp( compFlowAccessors.get( extrinsicMeshData::flow::temperature {} ) ),
    m_dTemp( compFlowAccessors.get( extrinsicMeshData::flow::deltaTemperature {} ) ),
    m_OBLOperatorValues ( compFlowAccessors.get( extrinsicMeshData::flow::OBLOperatorValues {} ) ),
    m_OBLOperatorDerivatives ( compFlowAccessors.get( extrinsicMeshData::flow::OBLOperatorDerivatives {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Time step size
  real64 const m_dt;

  /// trans multiplier
  real64 m_transMultExp;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > m_permeability;
  ElementViewConst< arrayView3d< real64 const > > m_dPerm_dPres;
  ElementViewConst< arrayView1d< real64 const > > m_referencePorosity;
  ElementViewConst< arrayView1d< real64 const > > m_rockThermalConductivity;


  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  // Primary variables

  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres;
  ElementViewConst< arrayView1d< real64 const > > const m_dPres;

  /// Views on global fraction
  ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const m_compFrac;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const m_dCompFrac;

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;
  ElementViewConst< arrayView1d< real64 const > > const m_dTemp;


  // Views on OBL operators and their derivatives
  ElementViewConst< arrayView2d< real64 const, compflow::USD_OBL_VAL > > const m_OBLOperatorValues;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_OBL_DER > > const m_OBLOperatorDerivatives;

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;
};

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */


template< integer NUM_PHASES, integer NUM_COMPS, bool ENABLE_ENERGY, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public FaceBasedAssemblyKernelBase
{
public:

  /// Compile time value for the number of phases
  static constexpr integer numPhases = NUM_PHASES;
  /// Compile time value for the number of components
  static constexpr integer numComps = NUM_COMPS;
  /// Compile time value for the energy balance switch
  static constexpr integer enableEnergyBalance = ENABLE_ENERGY;

  /// Compile time value for the number of primary variables
  static constexpr integer numDofs = numComps + enableEnergyBalance;

  /// Compile time value for the number of equations
  static constexpr integer numEqns = numDofs;

  // /// Compile time value for the number of operators
  static constexpr integer numOps = COMPUTE_NUM_OPS( NUM_PHASES, NUM_COMPS, ENABLE_ENERGY );

  // order of operators:
  static constexpr integer ACC_OP = 0;
  static constexpr integer FLUX_OP = numDofs;
  // diffusion
  static constexpr integer UPSAT_OP = numDofs + numDofs * numPhases;
  static constexpr integer GRAD_OP = numDofs + numDofs * numPhases + numPhases;
  // kinumDofstic reaction
  static constexpr integer KIN_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases;


  static constexpr integer RE_INTER_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs;
  static constexpr integer RE_TEMP_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 1;
  static constexpr integer ROCK_COND = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 2;
  static constexpr integer GRAV_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 3;
  static constexpr integer PC_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 3 + numPhases;
  static constexpr integer PORO_OP = numDofs + numDofs * numPhases + numPhases + numDofs * numPhases + numDofs + 3 + 2 * numPhases;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::NUM_POINT_IN_FLUX;

  /// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::MAX_NUM_OF_CONNECTIONS;

  /// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::MAX_STENCIL_SIZE;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] capPressureFlag flag specifying whether capillary pressure is used or not
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
  FaceBasedAssemblyKernel( globalIndex const rankOffset,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           CompFlowAccessors const & compFlowAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           real64 const & dt,
                           real64 const & transMultExp,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    : FaceBasedAssemblyKernelBase( rankOffset,
                                   dofNumberAccessor,
                                   compFlowAccessors,
                                   permeabilityAccessors,
                                   dt,
                                   transMultExp,
                                   localMatrix,
                                   localRhs ),
    m_stencilWrapper( stencilWrapper ),
    m_seri( stencilWrapper.getElementRegionIndices() ),
    m_sesri( stencilWrapper.getElementSubRegionIndices() ),
    m_sei( stencilWrapper.getElementIndices() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOSX_HOST_DEVICE
    StackVariables()
      : dofColIndices( maxNumElems * numDofs ),
      localFlux( numEqns ),
      localFluxJacobian( numEqns, maxNumElems * numDofs )
    {}

    /// Transmissibility
    real64 transmissibility[maxNumConns][2]{};
    /// Thermal Transmissibility
    real64 diffusiveTransmissibility[maxNumConns][2]{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dPres[maxNumConns][2]{};

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    stackArray1d< globalIndex, maxNumElems * numDofs > dofColIndices;

    /// Storage for the face local residual vector
    stackArray1d< real64, numEqns > localFlux;
    /// Storage for the face local Jacobian matrix
    stackArray2d< real64, maxNumElems * numEqns *  numDofs > localFluxJacobian;

  };



  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // The kernel is only designed for TPFA
    GEOSX_ASSERT_EQ( maxNumElems, 2 );
    GEOSX_ASSERT_EQ( maxNumConns, 1 );
    GEOSX_ASSERT_EQ( maxStencilSize, 2 );

    // set degrees of freedom indices for this face
    for( integer i = 0; i < maxStencilSize; ++i )
    {
      globalIndex const offset = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];

      for( integer jdof = 0; jdof < numDofs; ++jdof )
      {
        stack.dofColIndices[i * numDofs + jdof] = offset + jdof;
      }
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC1 the type of the function that can be used to customize the computation of the phase fluxes
   * @tparam FUNC2 the type of the function that can be used to customize the assembly into the local Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] phaseFluxKernelOp the function used to customize the computation of the phase fluxes
   * @param[in] localFluxJacobianKernelOp the function used to customize the computation of the assembly into the local Jacobian
   */
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {



    // first, compute the transmissibilities at this face
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    m_stencilWrapper.computeWeights( iconn,
                                     stack.diffusiveTransmissibility,
                                     stack.dTrans_dPres );


    // As a first iteration, do everything in TPFA style: i is the left (minus) element, j is the right(plus) element.
    // Consider inflow and outflow to and from element i.

    localIndex const erI  = m_seri( iconn, 0 );
    localIndex const esrI = m_sesri( iconn, 0 );
    localIndex const eiI  = m_sei( iconn, 0 );

    localIndex const erJ  = m_seri( iconn, 1 );
    localIndex const esrJ = m_sesri( iconn, 1 );
    localIndex const eiJ  = m_sei( iconn, 1 );

    arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLValsI = m_OBLOperatorValues[erI][esrI][eiI];
    arraySlice2d< real64 const, compflow::USD_OBL_DER - 1 > const & OBLDersI = m_OBLOperatorDerivatives[erI][esrI][eiI];

    arraySlice1d< real64 const, compflow::USD_OBL_VAL - 1 > const & OBLValsJ = m_OBLOperatorValues[erJ][esrJ][eiJ];
    arraySlice2d< real64 const, compflow::USD_OBL_DER - 1 > const & OBLDersJ = m_OBLOperatorDerivatives[erJ][esrJ][eiJ];

    // single gravity coefficient and transimissibility values correspond to the current connection
    real64 const gravCoef = m_gravCoef[erI][esrI][eiI] - m_gravCoef[erJ][esrJ][eiJ];
    real64 const trans = stack.transmissibility[0][0];
    real64 const transD = stack.diffusiveTransmissibility[0][0];


    real64 transMult = 1;
    real64 transMultDerI[numDofs];
    real64 transMultDerJ[numDofs];
    if( m_transMultExp > 0 )
    {
      // Calculate transmissibility multiplier:

      // Take average interface porosity:
      real64 const poroAverage = (OBLValsI[PORO_OP] + OBLValsJ[PORO_OP]) * 0.5;
      real64 const transMultD = m_transMultExp * pow( poroAverage, m_transMultExp - 1 ) * 0.5;

      for( integer v = 0; v < numDofs; v++ )
      {
        transMultDerI[v] = transMultD * OBLDersI[PORO_OP][v];
        transMultDerJ[v] = transMultD * OBLDersJ[PORO_OP][v];
      }
      transMult = pow( poroAverage, m_transMultExp );
    }
    else
    {
      for( integer v = 0; v < numDofs; v++ )
      {
        transMultDerI[v] = 0;
        transMultDerJ[v] = 0;
      }
    }

    real64 const pDiff = (m_pres[erJ][esrJ][eiJ] + m_dPres[erJ][esrJ][eiJ]) - (m_pres[erI][esrI][eiI] + m_dPres[erI][esrI][eiI]);


    // [2] fill offdiagonal part + contribute to diagonal, only fluid part is considered in energy equation
    for( integer p = 0; p < numPhases; p++ )
    {     // loop over number of phases for convective operator

      // calculate gravity term for phase p
      real64 const avg_density = (OBLValsI[GRAV_OP + p] + OBLValsJ[GRAV_OP + p]) / 2;

      // p = 1 means oil phase, it's reference phase. pw=po-pcow, pg=po-(-pcog).
      real64 const phase_p_diff = pDiff + avg_density * gravCoef - OBLValsJ[PC_OP + p] + OBLValsI[PC_OP + p];

      // calculate partial derivatives for gravity and capillary terms
      real64 gravPcDerI[numDofs];
      real64 gravPcDerJ[numDofs];
      for( integer v = 0; v < numDofs; v++ )
      {
        gravPcDerI[v] = -(OBLDersI[GRAV_OP + p][v]) * gravCoef / 2 - OBLDersI[PC_OP + p][v];
        gravPcDerJ[v] = -(OBLDersJ[GRAV_OP + p][v]) * gravCoef / 2 + OBLDersJ[PC_OP + p][v];
      }

      real64 phase_gamma_p_diff = transMult * trans * m_dt * phase_p_diff;

      if( phase_p_diff < 0 )
      {
        // mass and energy outflow with effect of gravity and capillarity
        for( integer c = 0; c < numDofs; c++ )
        {
          real64 c_flux = transMult * trans * m_dt * OBLValsI[FLUX_OP + p * numDofs + c];

          stack.localFlux[c] -= phase_p_diff * c_flux;     // flux operators only
          for( integer v = 0; v < numDofs; v++ )
          {
            stack.localFluxJacobian[c][v] -= (phase_gamma_p_diff * OBLDersI[FLUX_OP + p * numDofs + c][v] +
                                              trans * m_dt * phase_p_diff * transMultDerI[v] * OBLValsI[FLUX_OP + p * numDofs + c]);
            stack.localFluxJacobian[c][v] += c_flux * gravPcDerI[v];
            stack.localFluxJacobian[c][numDofs + v] += c_flux * gravPcDerJ[v];

            if( v == 0 )
            {
              stack.localFluxJacobian[c][numDofs + v] -= c_flux;
              stack.localFluxJacobian[c][v] += c_flux;
            }
          }
        }
      }
      else
      {
        // mass and energy inflow with effect of gravity and capillarity
        for( integer c = 0; c < numDofs; c++ )
        {
          real64 c_flux = transMult * trans * m_dt * OBLValsJ[FLUX_OP + p * numDofs + c];

          stack.localFlux[c] -= phase_p_diff * c_flux;     // flux operators only
          for( integer v = 0; v < numDofs; v++ )
          {
            stack.localFluxJacobian[c][numDofs + v] -= (phase_gamma_p_diff * OBLDersJ[FLUX_OP + p * numDofs + c][v] +
                                                        trans * m_dt * phase_p_diff * transMultDerJ[v] * OBLValsJ[FLUX_OP + p * numDofs + c]);
            stack.localFluxJacobian[c][v] += c_flux * gravPcDerI[v];
            stack.localFluxJacobian[c][numDofs + v] += c_flux * gravPcDerJ[v];
            if( v == 0 )
            {
              stack.localFluxJacobian[c][v] += c_flux;     //-= Jac[jac_idx + c * numDofs];
              stack.localFluxJacobian[c][numDofs + v] -= c_flux;      // -trans * m_dt * op_vals[NC + c];
            }
          }
        }
      }

    }     // end of loop over number of phases for convective operator with gravity and capillarity

    // [3] Additional diffusion code here:   (phi_p * S_p) * (rho_p * D_cp * Delta_x_cp)  or (phi_p * S_p) * (kappa_p * Delta_T)

    real64 const poroAverage = (m_referencePorosity[erI][esrI][eiI] + m_referencePorosity[erJ][esrJ][eiJ]) * 0.5;     // diffusion term
                                                                                                                      // depends on total
    // porosity!


    // Add diffusion term to the residual:
    for( integer c = 0; c < numDofs; c++ )
    {
      for( integer p = 0; p < numPhases; p++ )
      {
        real64 grad_con = OBLValsJ[GRAD_OP + c * numPhases + p] - OBLValsI[GRAD_OP + c * numPhases + p];

        if( grad_con < 0 )
        {
          // Diffusion flows from cell i to j (high to low), use upstream quantity from cell i for compressibility and saturation (mass or
          // energy):
          real64 diff_mob_ups_m = m_dt * transD * poroAverage * OBLValsI[UPSAT_OP + p];

          stack.localFlux[c] -= diff_mob_ups_m * grad_con;       // diffusion term

          // Add diffusion terms to Jacobian:
          for( integer v = 0; v < numDofs; v++ )
          {
            stack.localFluxJacobian[c][v] += diff_mob_ups_m * OBLDersI[GRAD_OP + c * numPhases + p][v];
            stack.localFluxJacobian[c][numDofs + v] -= diff_mob_ups_m * OBLDersJ[GRAD_OP + c * numPhases + p][v];

            stack.localFluxJacobian[c][v] -= grad_con * m_dt * transD * poroAverage * OBLDersI[UPSAT_OP + p][v];
          }
        }
        else
        {
          // Diffusion flows from cell j to i (high to low), use upstream quantity from cell j for density and saturation:
          real64 diff_mob_ups_m = m_dt * transD * poroAverage * OBLValsJ[UPSAT_OP + p];

          stack.localFlux[c] -= diff_mob_ups_m * grad_con;       // diffusion term

          // Add diffusion terms to Jacobian:
          for( integer v = 0; v < numDofs; v++ )
          {
            stack.localFluxJacobian[c][v] += diff_mob_ups_m * OBLDersI[GRAD_OP + c * numPhases + p][v];
            stack.localFluxJacobian[c][numDofs + v] -= diff_mob_ups_m * OBLDersJ[GRAD_OP + c * numPhases + p][v];

            stack.localFluxJacobian[c][numDofs + v] -= grad_con * m_dt * transD * poroAverage * OBLDersJ[UPSAT_OP + p][v];
          }
        }
      }
    }

    // [4] add rock conduction
    if( enableEnergyBalance )
    {
      real64 const tDiff = OBLValsJ[RE_TEMP_OP] - OBLValsI[RE_TEMP_OP];
      real64 const gammaTDiff = transD * m_dt * tDiff;

      if( tDiff < 0 )
      {
        real64 const poroI = m_referencePorosity[erI][esrI][eiI];
        real64 const rockCondI = m_rockThermalConductivity[erI][esrI][eiI];

        // rock heat transfers flows from cell i to j
        stack.localFlux[numComps] -= gammaTDiff * OBLValsI[ROCK_COND] * (1 - poroI) * rockCondI;
        for( integer v = 0; v < numDofs; v++ )
        {
          stack.localFluxJacobian[numComps][v] -= gammaTDiff * OBLDersI[ROCK_COND][v] * (1 - poroI) * rockCondI;
          // the last variable - T
          if( v == numDofs - 1 )
          {
            stack.localFluxJacobian[numComps][numDofs + v] -= transD * m_dt * OBLValsI[ROCK_COND] * (1 - poroI) * rockCondI;
            stack.localFluxJacobian[numComps][v] += transD * m_dt * OBLValsI[ROCK_COND] * (1 - poroI) * rockCondI;
          }
        }
      }
      else
      {
        real64 const poroJ = m_referencePorosity[erJ][esrJ][eiJ];
        real64 const rockCondJ = m_rockThermalConductivity[erJ][esrJ][eiJ];

        // rock heat transfers flows from cell j to i
        stack.localFlux[numComps] -= gammaTDiff * OBLValsJ[ROCK_COND] * (1 - poroJ) * rockCondJ;     // energy
        // cond
        // operator
        for( integer v = 0; v < numDofs; v++ )
        {
          stack.localFluxJacobian[numComps][numDofs + v] -= gammaTDiff * OBLDersJ[ROCK_COND][v] * (1 - poroJ) * rockCondJ;
          // the last variable - T
          if( v == numDofs - 1 )
          {
            stack.localFluxJacobian[numComps][v] += transD * m_dt * OBLValsJ[ROCK_COND] * (1 - poroJ) * rockCondJ;
            stack.localFluxJacobian[numComps][numDofs + v] -= transD * m_dt * OBLValsJ[ROCK_COND] * (1 - poroJ) * rockCondJ;
          }
        }
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    // Add to residual/jacobian
    for( integer i = 0; i < maxNumElems; ++i )
    {
      // during first loop iteration, we fill the equations for element i as is.
      if( i == 1 )
      {
        // during second (and the last) loop iteration, we fill now the equations for element j.
        // here, we need to multiply all values by -1
        for( integer ic = 0; ic < numComps; ++ic )
        {
          stack.localFlux[ic] *= -1;
          for( integer v = 0; v < 2 * numDofs; ++v )
          {
            stack.localFluxJacobian[ic][v] *= -1;
          }
        }
      }

      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GE( m_localMatrix.numRows(), localRow + numEqns );

        for( integer ic = 0; ic < numComps; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + ic], stack.localFlux[ic] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
            ( localRow + ic,
            stack.dofColIndices.data(),
            stack.localFluxJacobian[ic].dataIfContiguous(),
            maxNumElems * numDofs );
        }
      }
    }
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numConnections the number of connections
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numConnections,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    if( numConnections )
    {
      forAll< POLICY >( numConnections, [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
      {
        typename KERNEL_TYPE::StackVariables stack;

        kernelComponent.setup( iconn, stack );
        kernelComponent.computeFlux( iconn, stack );
        kernelComponent.complete( iconn, stack );
      } );
    }
  }

private:

  // Stencil information

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
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
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] capPressureFlag flag specifying whether capillary pressure is used or not
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] targetRegionNames names of the target regions
   * @param[in] fluidModelNames names of the fluid models
   * @param[in] capPresModelNames names of the capillary pressure models
   * @param[in] permeabilityModelNames names of the permeability models
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numPhases,
                   integer const numComps,
                   bool const enableEnergyBalance,
                   real64 const & transMultExp,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   arrayView1d< string const > const & targetRegionNames,
                   arrayView1d< string const > const & permeabilityModelNames,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    internal::kernelLaunchSelectorEnergySwitch( numPhases, numComps, enableEnergyBalance, [&] ( auto NP, auto NC, auto E )
    {
      integer constexpr ENABLE_ENERGY = E();
      integer constexpr NUM_PHASES = NP();
      integer constexpr NUM_COMPS = NC();

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      using KERNEL_TYPE = FaceBasedAssemblyKernel< NUM_PHASES, NUM_COMPS, ENABLE_ENERGY, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName, targetRegionNames, permeabilityModelNames );

      KERNEL_TYPE kernel( rankOffset, stencilWrapper, dofNumberAccessor, compFlowAccessors, permeabilityAccessors,
                          dt, transMultExp, localMatrix, localRhs );
      KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
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
                      arrayView1d< real64 const > const & totalDensOld,
                      real64 & localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > localSum( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        real64 const normalizer = totalDensOld[ei] * refPoro[ei] * volume[ei];

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



/******************************** Kernel launch machinery ********************************/

// template< typename KERNELWRAPPER, typename ... ARGS >
// void KernelLaunchSelector1( integer const numComps, ARGS && ... args )
// {
//   internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
//   {
//     KERNELWRAPPER::template launch< NC() >( std::forward< ARGS >( args )... );
//   } );
// }

// template< typename KERNELWRAPPER, typename ... ARGS >
// void KernelLaunchSelector2( integer const numComps, integer const numPhase, ARGS && ... args )
// {
//   // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
//   if( numPhase == 2 )
//   {
//     internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
//     {
//       KERNELWRAPPER::template launch< NC(), 2 >( std::forward< ARGS >( args ) ... );
//     } );
//   }
//   else if( numPhase == 3 )
//   {
//     internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
//     {
//       KERNELWRAPPER::template launch< NC(), 3 >( std::forward< ARGS >( args ) ... );
//     } );
//   }
//   else
//   {
//     GEOSX_ERROR( "Unsupported number of phases: " << numPhase );
//   }
// }

} // namespace CompositionalMultiphaseBaseKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_DARTSSUPERENGINEKERNELS_HPP
