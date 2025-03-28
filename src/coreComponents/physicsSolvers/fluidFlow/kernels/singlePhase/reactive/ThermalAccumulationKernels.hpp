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
 * @file ThermalAccumulationKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALACCUMULATIONKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALACCUMULATIONKERNELS_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/AccumulationKernels.hpp"

namespace geos
{

namespace thermalSinglePhaseReactiveBaseKernels
{

/******************************** AccumulationKernel ********************************/

/**
 * @class AccumulationKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */
template< typename SUBREGION_TYPE, integer NUM_DOF, integer NUM_SPECIES >
class AccumulationKernel : public singlePhaseReactiveBaseKernels::AccumulationKernel< SUBREGION_TYPE, NUM_DOF, NUM_SPECIES >
{

public:

  using Base = singlePhaseReactiveBaseKernels::AccumulationKernel< SUBREGION_TYPE, NUM_DOF, NUM_SPECIES >;
  using Base::numDof;
  using Base::numEqn;
  using Base::numSpecies;
  using Base::m_rankOffset;
  using Base::m_dofNumber;
  using Base::m_elemGhostRank;
  using Base::m_localMatrix;
  using Base::m_localRhs;
  using Base::m_dMass;
  using Base::m_volume;
  using Base::m_deltaVolume;
  using Base::m_primarySpeciesAggregateConcentration;

  /// Note: Derivative lineup only supports dP & dT, not component terms
  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
  /**
   * @brief Constructor
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  AccumulationKernel( globalIndex const rankOffset,
                      string const dofKey,
                      SUBREGION_TYPE const & subRegion,
                      constitutive::ReactiveSingleFluid const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      real64 const & dt,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, fluid, solid, dt, localMatrix, localRhs ),
    m_energy( subRegion.template getField< fields::flow::energy >() ),
    m_energy_n( subRegion.template getField< fields::flow::energy_n >() ),
    m_dEnergy( subRegion.template getField< fields::flow::dEnergy >() ),
    m_dPoro_dTemp( solid.getDporosity_dTemperature() )
    // m_dPrimarySpeciesAggregateConcentration_dTemp( fluid.dPrimarySpeciesAggregateConcentration_dTemp() ),
    // m_dPrimarySpeciesTotalKineticRate_dTemp( fluid.dPrimarySpeciesTotalKineticRate_dTemp() ),
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables()
      : Base::StackVariables()
    {}

    using Base::StackVariables::localRow;
    using Base::StackVariables::dofIndices;
    using Base::StackVariables::localResidual;
    using Base::StackVariables::localJacobian;
    using Base::StackVariables::poreVolume;
    using Base::StackVariables::dPoreVolume_dLogPrimaryConc;

    // Pore volume information

    /// Derivative of pore volume with respect to temperature
    real64 dPoreVolume_dTemp = 0.0;
  };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

    stack.dPoreVolume_dTemp = ( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dTemp[ei][0];
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    Base::computeAccumulation( ei, stack );

    // Step 1: assemble the derivatives of the mass balance equation w.r.t temperature
    stack.localJacobian[0][numDof-numSpecies-1] = m_dMass[ei][DerivOffset::dT];

    // Step 2: assemble the accumulation term of the energy equation
    // Step 2.1: assemble the residual and derivatives wrt pressure and temperature
    stack.localResidual[numEqn-numSpecies-1] = m_energy[ei] - m_energy_n[ei];
    stack.localJacobian[numEqn-numSpecies-1][0] += m_dEnergy[ei][DerivOffset::dP];
    stack.localJacobian[numEqn-numSpecies-1][numDof-numSpecies-1] += m_dEnergy[ei][DerivOffset::dT];

    // Step 2.2: assemble the derivatives of the energy equation w.r.t log primary species concentration
    // for( integer is = 0; is < numSpecies; ++is )
    // {
    //   stack.localJacobian[numEqn-numSpecies-1][is+numDof-numSpecies] += stack.dPoreVolume_dLogPrimaryConc[is] * m_density[ei][0] *
    // m_fluidInternalEnergy[ei][0]
    //                                                                     - stack.dPoreVolume_dLogPrimaryConc[is] *
    // m_rockInternalEnergy[ei][0]
    //                                                                     + stack.poreVolume * m_dDensity_dLogPrimaryConc[ei][is] *
    // m_fluidInternalEnergy[ei][0]
    //                                                                     + stack.poreVolume * m_density[ei][0] *
    // m_dFluidInternalEnergy_dLogPrimaryConc[ei][is];
    // }

    // Step 3: assemble the derivatives of the species amount balance equation w.r.t temperature
    for( integer is = 0; is < numSpecies; ++is )
    {
      // Drivative of primary species amount in pore volume wrt temperature
      stack.localJacobian[is+numEqn-numSpecies][numDof-numSpecies-1] += stack.dPoreVolume_dTemp * m_primarySpeciesAggregateConcentration[ei][is]
                                                                        /* + stack.poreVolume *
                                                                           m_dPrimarySpeciesAggregateConcentration_dTemp[ei][is] */;
      // // Derivative of reaction term wrt temperature
      // stack.localJacobian[is+numEqn-numSpecies][numDof-numSpecies-1] -= m_dt * ( m_volume[ei] + m_deltaVolume[ei] ) *
      // m_dPrimarySpeciesTotalKineticRate_dTemp[is];
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // Step 1: assemble the total mass balance equation (i = 0)
    // and species amount balance equation (i = numEqn-numSpecies to i = numEqn-1)
    Base::complete( ei, stack );

    // Step 2: assemble the energy equation (i = numEqn-numSpecies-1)
    m_localRhs[stack.localRow + numEqn-numSpecies-1] += stack.localResidual[numEqn-numSpecies-1];
    m_localMatrix.template addToRow< serialAtomic >( stack.localRow + numEqn-numSpecies-1,
                                                     stack.dofIndices,
                                                     stack.localJacobian[numEqn-numSpecies-1],
                                                     numDof );
  }

protected:

  /// View on energy
  arrayView1d< real64 const > const m_energy;
  arrayView1d< real64 const > const m_energy_n;
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_dEnergy;

  /// Views on the porosity derivative
  arrayView2d< real64 const > const m_dPoro_dTemp;

  // // View on the derivatives of aggregate concentration for the primary species wrt temperature
  // arrayView2d< real64 const, compflow::USD_COMP > m_dPrimarySpeciesAggregateConcentration_dTemp;

  // // View on the derivatives of total kinetic rate of primary species wrt temperature
  // arrayView2d< real64 const, compflow::USD_COMP > m_dPrimarySpeciesTotalKineticRate_dTemp;

};

/**
 * @class AccumulationKernelFactory
 */
class AccumulationKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numSpecies the number of primary species
   * @param[in] dt time step
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename SUBREGION_TYPE >
  static void
  createAndLaunch( integer const numSpecies,
                   real64 const dt,
                   globalIndex const rankOffset,
                   string const dofKey,
                   SUBREGION_TYPE const & subRegion,
                   constitutive::ReactiveSingleFluid const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    singlePhaseReactiveBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numSpecies, [&] ( auto NS )
    {
      integer constexpr NUM_SPECIES = NS();
      integer constexpr NUM_DOF = 2+NS();
      AccumulationKernel< SUBREGION_TYPE, NUM_DOF, NUM_SPECIES > kernel( rankOffset, dofKey, subRegion, fluid, solid, dt, localMatrix, localRhs );
      AccumulationKernel< SUBREGION_TYPE, NUM_DOF, NUM_SPECIES >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }
};

} // namespace thermalSinglePhaseReactiveBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_THERMALACCUMULATIONKERNELS_HPP
