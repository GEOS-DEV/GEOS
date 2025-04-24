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
 * @file AccumulationKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_ACCUMULATIONKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_ACCUMULATIONKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/reactivefluid/ReactiveSinglePhaseFluid.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/AccumulationKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/KernelLaunchSelectors.hpp"

namespace geos
{

namespace singlePhaseReactiveBaseKernels
{

/******************************** AccumulationKernel ********************************/

/**
 * @class AccumulationKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */
template< typename SUBREGION_TYPE, integer NUM_DOF, integer NUM_SPECIES >
class AccumulationKernel : public singlePhaseBaseKernels::AccumulationKernel< SUBREGION_TYPE, NUM_DOF >
{

public:

  using Base = singlePhaseBaseKernels::AccumulationKernel< SUBREGION_TYPE, NUM_DOF >;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_rankOffset;
  using Base::m_dofNumber;
  using Base::m_elemGhostRank;
  using Base::m_localMatrix;
  using Base::m_localRhs;
  using Base::m_dMass;

  /// Note: Derivative lineup only supports dP & dT, not component terms
  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 0 >;

  /// Compile time value for the number of primary species
  static constexpr integer numSpecies = NUM_SPECIES;

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
                      constitutive::reactivefluid::ReactiveSinglePhaseFluid< constitutive::CompressibleSinglePhaseFluid > const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      real64 const & dt,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, localMatrix, localRhs ),
    m_dt( dt ),
    m_volume( subRegion.getElementVolume() ),
    m_deltaVolume( subRegion.template getField< fields::flow::deltaVolume >() ),
    m_porosity( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    // m_dDensity_dLogPrimaryConc( fluid.dDensity_dLogPrimaryConc() ),
    // m_dPoro_dLogPrimaryConc( solid.getDporosity_dLogPrimaryConc() ),
    m_primarySpeciesAggregateConcentration( fluid.primarySpeciesAggregateConcentration() ),
    // m_dPrimarySpeciesAggregateConcentration_dPres( fluid.dPrimarySpeciesAggregateConcentration_dPres() ),
    m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc( fluid.dPrimarySpeciesAggregateConcentration_dLogPrimaryConc() ),
    // m_primarySpeciesAggregateKineticRate( fluid.kineticReactionRates() ),
    // m_dPrimarySpeciesAggregateKineticRate_dPres( fluid.dPrimarySpeciesAggregateKineticRate_dPres() ),
    // m_dPrimarySpeciesAggregateKineticRate_dLogPrimaryConc( fluid.dPrimarySpeciesAggregateKineticRate_dLogPrimaryConc() ),
    m_primarySpeciesAggregateMole_n( subRegion.template getField< fields::flow::primarySpeciesAggregateMole_n >() )
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

    // Pore volume information

    /// Pore volume at time n+1
    real64 poreVolume = 0.0;

    /// Derivative of pore volume with respect to pressure
    real64 dPoreVolume_dPres = 0.0;

    /// Derivative of pore volume with respect to each primary species concentration
    real64 dPoreVolume_dLogPrimaryConc[numSpecies]{};

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
    // initialize the pore volume
    stack.poreVolume = ( m_volume[ei] + m_deltaVolume[ei] ) * m_porosity[ei][0];
    stack.dPoreVolume_dPres = ( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dPres[ei][0];

    Base::setup( ei, stack );

    // // is - index of the primary species
    // for( integer is = 0; is < numSpecies; ++is )
    // {
    //   stack.dPoreVolume_dLogPrimaryConc[is] = ( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dLogPrimaryConc[ei][is]
    // }
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
    // Residual[is] += (primarySpeciesAggregateConcentration[is] * stack.poreVolume - primarySpeciesAggregateMole_n[is])
    //                 - dt * m_volume * primarySpeciesKineticRate[is] // To Check: what's the unit of the kinetic rate

    Base::computeAccumulation( ei, stack );

    // Step 1: assemble the derivatives of the total mass equation w.r.t log of primary species concentration
    //   for( integer is = 0; is < numSpecies; ++is )
    //   {
    //     stack.localJacobian[0][is+numDof-numSpecies] = stack.poreVolume * m_dDensity_dLogPrimaryConc[ei][is] +
    // stack.dPoreVolume_dLogPrimaryConc[is] * m_density[ei][0];
    //   }

    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dPrimarySpeciesAggregateConcentration_dLogPrimaryConc = m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc[ei];

    for( integer is = 0; is < numSpecies; ++is )
    {
      // Step 2: assemble the accumulation term of the species mass balance equation
      // Step 2.1: residual
      // Primary species mole amount in pore volume
      stack.localResidual[is+numEqn-numSpecies] -= m_primarySpeciesAggregateMole_n[ei][is];
      stack.localResidual[is+numEqn-numSpecies] += m_primarySpeciesAggregateConcentration[ei][is] * stack.poreVolume;

      // // Reaction term
      // stack.localResidual[is+numEqn-numSpecies] -= m_dt * ( m_volume[ei] + m_deltaVolume[ei] ) *
      // m_primarySpeciesAggregateKineticRate[is];

      // Step 2.1: jacobian
      // Drivative of primary species amount in pore volume wrt pressure
      stack.localJacobian[is+numEqn-numSpecies][0] += stack.dPoreVolume_dPres * m_primarySpeciesAggregateConcentration[ei][is]
                                                      /* + stack.poreVolume * m_dTotalPrimarySpeciesConcentration_dPres[ei][is] */;
      // // Derivative of reaction term wrt pressure
      // stack.localJacobian[is+numEqn-numSpecies][0] -= m_dt * ( m_volume[ei] + m_deltaVolume[ei] ) *
      // m_dPrimarySpeciesTotalKineticRate_dPres[is];

      // Derivative wrt log of primary species concentration
      // arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dPrimarySpeciesTotalKineticRate_dLogPrimaryConc =
      // m_dPrimarySpeciesTotalKineticRate_dLogPrimaryConc[ei];

      for( integer js = 0; js < numSpecies; ++js )
      {
        stack.localJacobian[is+numEqn-numSpecies][js+numDof-numSpecies] += /* stack.dPoreVolume_dLogPrimaryConc[js] *
                                                                              m_primarySpeciesAggregateConcentration[ei][is]
                                                                            + */stack.poreVolume * dPrimarySpeciesAggregateConcentration_dLogPrimaryConc[is][js]; // To
                                                                                                                                                                  // check
                                                                                                                                                                  // if
                                                                                                                                                                  // the
                                                                                                                                                                  // permutation
                                                                                                                                                                  // is
                                                                                                                                                                  // consistent

        // stack.localJacobian[is+numEqn-numSpecies][js+numDof-numSpecies] -= m_dt * ( m_volume[ei] + m_deltaVolume[ei] ) *
        // dPrimarySpeciesTotalKineticRate_dLogPrimaryConc[is][js];
      }
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
    // Step 1: assemble the total mass balance equation
    // - the total mass balance equations (i = 0)
    Base::complete( ei, stack );

    // Step 2: assemble the primary species mole amount balance equation
    // - the species mole amount balance equations (i = numEqn-numSpecies to i = numEqn-1)
    integer const beginRowSpecies = numEqn-numSpecies;
    for( integer i = 0; i < numSpecies; ++i )
    {
      m_localRhs[stack.localRow + beginRowSpecies + i] += stack.localResidual[beginRowSpecies+i];
      m_localMatrix.template addToRow< serialAtomic >( stack.localRow + beginRowSpecies + i,
                                                       stack.dofIndices,
                                                       stack.localJacobian[beginRowSpecies + i],
                                                       numDof );
    }
  }

protected:

  /// Time step size
  real64 const m_dt;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;
  arrayView1d< real64 const > const m_deltaVolume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosity;
  arrayView2d< real64 const > const m_dPoro_dPres;

  // // View on the derivatives of fluid density wrt log of primary species concentration
  // arrayView2d< real64 const, compflow::USD_COMP > m_dDensity_dLogPrimaryConc;

  // // View on the derivatives of porosity wrt log of primary species concentration
  // arrayView2d< real64 const, compflow::USD_COMP > m_dPoro_dLogPrimaryConc;

  // View on the total concentration of ions that contain the primary species
  arrayView2d< real64 const, compflow::USD_COMP > m_primarySpeciesAggregateConcentration;

  // // View on the derivatives of aggregate concentration for the primary species wrt pressure
  // arrayView2d< real64 const, compflow::USD_COMP > m_dPrimarySpeciesAggregateConcentration_dPres;

  // View on the derivatives of total ion concentration for the primary species wrt log of primary species concentration
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dPrimarySpeciesAggregateConcentration_dLogPrimaryConc;

  // // View on the aggregate kinetic rate of primary species from all reactions
  // arrayView2d< real64 const, compflow::USD_COMP > m_primarySpeciesAggregateKineticRate;

  // // View on the derivatives of aggregate kinetic rate of primary species wrt pressure
  // arrayView2d< real64 const, compflow::USD_COMP > m_dPrimarySpeciesAggregateKineticRate_dPres;

  // // View on the derivatives of aggregate kinetic rate of primary species wrt log of primary species concentration
  // arrayView3d< real64 const, compflow::USD_COMP_DC > m_dPrimarySpeciesAggregateKineticRate_dLogPrimaryConc;

  // View on primary species mole amount from previous time step
  arrayView2d< real64 const, compflow::USD_COMP > m_primarySpeciesAggregateMole_n;
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
                   constitutive::reactivefluid::ReactiveSinglePhaseFluid< constitutive::CompressibleSinglePhaseFluid > const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    internal::kernelLaunchSelectorCompSwitch( numSpecies, [&] ( auto NS )
    {
      integer constexpr NUM_SPECIES = NS();
      integer constexpr NUM_DOF = 1+NS();
      AccumulationKernel< SUBREGION_TYPE, NUM_DOF, NUM_SPECIES > kernel( rankOffset, dofKey, subRegion, fluid, solid, dt, localMatrix, localRhs );
      AccumulationKernel< SUBREGION_TYPE, NUM_DOF, NUM_SPECIES >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }
};

} // namespace singlePhaseReactiveBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVE_ACCUMULATIONKERNELS_HPP
