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
 * @file ElementBasedAssemblyKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMAL_ELEMENTBASEDASSEMBLYKERNEL_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMAL_ELEMENTBASEDASSEMBLYKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/FVM/ElementBasedAssemblyKernel.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

namespace thermalSinglePhaseBaseKernels
{

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */
template< typename SUBREGION_TYPE, integer NUM_DOF >
class ElementBasedAssemblyKernel : public singlePhaseBaseKernels::ElementBasedAssemblyKernel< SUBREGION_TYPE, NUM_DOF >
{

public:

  using Base = singlePhaseBaseKernels::ElementBasedAssemblyKernel< SUBREGION_TYPE, NUM_DOF >;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_rankOffset;
  using Base::m_dofNumber;
  using Base::m_elemGhostRank;
  using Base::m_volume;
  using Base::m_deltaVolume;
  using Base::m_porosity_n;
  using Base::m_porosityNew;
  using Base::m_dPoro_dPres;
  using Base::m_density_n;
  using Base::m_density;
  using Base::m_dDensity_dPres;
  using Base::m_localMatrix;
  using Base::m_localRhs;

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
  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              string const dofKey,
                              SUBREGION_TYPE const & subRegion,
                              constitutive::SingleFluidBase const & fluid,
                              constitutive::CoupledSolidBase const & solid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs ),
    m_dDensity_dTemp( fluid.dDensity_dTemperature() ),
    m_dPoro_dTemp( solid.getDporosity_dTemperature() ),
    m_internalEnergy_n( fluid.internalEnergy_n() ),
    m_internalEnergy( fluid.internalEnergy() ),
    m_dInternalEnergy_dPres( fluid.dInternalEnergy_dPressure() ),
    m_dInternalEnergy_dTemp( fluid.dInternalEnergy_dTemperature() ),
    m_rockInternalEnergy_n( solid.getInternalEnergy_n() ),
    m_rockInternalEnergy( solid.getInternalEnergy() ),
    m_dRockInternalEnergy_dTemp( solid.getDinternalEnergy_dTemperature() )
  { }

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
    { }

    using Base::StackVariables::poreVolume;
    using Base::StackVariables::poreVolume_n;
    using Base::StackVariables::dPoreVolume_dPres;
    using Base::StackVariables::localRow;
    using Base::StackVariables::dofIndices;
    using Base::StackVariables::localResidual;
    using Base::StackVariables::localJacobian;

    /// Derivative of pore volume with respect to temperature
    real64 dPoreVolume_dTemp = 0.0;

    // Solid energy

    /// Solid energy at time n+1
    real64 solidEnergy = 0.0;

    /// Solid energy at the previous converged time step
    real64 solidEnergy_n = 0.0;

    /// Derivative of solid internal energy with respect to pressure
    real64 dSolidEnergy_dPres = 0.0;

    /// Derivative of solid internal energy with respect to temperature
    real64 dSolidEnergy_dTemp = 0.0;
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
    ElementBasedAssemblyKernel< SUBREGION_TYPE, NUM_DOF >::setup( ei, stack );

    stack.dPoreVolume_dTemp = ( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dTemp[ei][0];

    // initialize the solid volume
    real64 const solidVolume = ( m_volume[ei] + m_deltaVolume[ei] ) * ( 1.0 - m_porosityNew[ei][0] );
    real64 const solidVolume_n = m_volume[ei] * ( 1.0 - m_porosity_n[ei][0] );
    real64 const dSolidVolume_dPres = -( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dPres[ei][0];
    real64 const dSolidVolume_dTemp = -( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dTemp[ei][0];

    // initialize the solid internal energy
    stack.solidEnergy = solidVolume * m_rockInternalEnergy[ei][0];
    stack.solidEnergy_n = solidVolume_n * m_rockInternalEnergy_n[ei][0];
    stack.dSolidEnergy_dPres = dSolidVolume_dPres * m_rockInternalEnergy[ei][0];
    stack.dSolidEnergy_dTemp = solidVolume * m_dRockInternalEnergy_dTemp[ei][0] + dSolidVolume_dTemp * m_rockInternalEnergy[ei][0];
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    Base::computeAccumulation( ei, stack, [&]()
    {
      // Step 1: assemble the derivatives of the mass balance equation w.r.t temperature
      stack.localJacobian[0][numDof - 1] = stack.poreVolume * m_dDensity_dTemp[ei][0] + stack.dPoreVolume_dTemp * m_density[ei][0];

      // Step 2: assemble the fluid part of the accumulation term of the energy equation
      real64 const fluidEnergy = stack.poreVolume * m_density[ei][0] * m_internalEnergy[ei][0];
      real64 const fluidEnergy_n = stack.poreVolume_n * m_density_n[ei][0] * m_internalEnergy_n[ei][0];

      real64 const dFluidEnergy_dP = stack.dPoreVolume_dPres * m_density[ei][0] * m_internalEnergy[ei][0]
                                     + stack.poreVolume * m_dDensity_dPres[ei][0] * m_internalEnergy[ei][0]
                                     + stack.poreVolume * m_density[ei][0] * m_dInternalEnergy_dPres[ei][0];

      real64 const dFluidEnergy_dT = stack.poreVolume * m_dDensity_dTemp[ei][0] * m_internalEnergy[ei][0]
                                     + stack.poreVolume * m_density[ei][0] * m_dInternalEnergy_dTemp[ei][0]
                                     + stack.dPoreVolume_dTemp * m_density[ei][0] * m_internalEnergy[ei][0];

      // local accumulation
      stack.localResidual[numEqn - 1] = fluidEnergy - fluidEnergy_n;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[numEqn - 1][0] = dFluidEnergy_dP;
      stack.localJacobian[numEqn - 1][numDof - 1] = dFluidEnergy_dT;
    } );

    // Step 3: assemble the solid part of the accumulation term of the energy equation
    stack.localResidual[numEqn - 1] += stack.solidEnergy - stack.solidEnergy_n;
    stack.localJacobian[numEqn - 1][0] += stack.dSolidEnergy_dPres;
    stack.localJacobian[numEqn - 1][numDof - 1] += stack.dSolidEnergy_dTemp;
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
    // Step 1: assemble the mass balance equation
    ElementBasedAssemblyKernel< SUBREGION_TYPE, NUM_DOF >::complete( ei, stack );

    // Step 2: assemble the energy equation
    m_localRhs[stack.localRow + numEqn - 1] += stack.localResidual[numEqn - 1];
    m_localMatrix.template addToRow< serialAtomic >( stack.localRow + numEqn - 1,
                                                     stack.dofIndices,
                                                     stack.localJacobian[numEqn - 1],
                                                     numDof );


  }

protected:

  /// View on derivative of fluid density w.r.t temperature
  arrayView2d< real64 const > const m_dDensity_dTemp;

  /// View on derivative of porosity w.r.t temperature
  arrayView2d< real64 const > const m_dPoro_dTemp;

  /// Views on fluid internal energy
  arrayView2d< real64 const > const m_internalEnergy_n;
  arrayView2d< real64 const > const m_internalEnergy;
  arrayView2d< real64 const > const m_dInternalEnergy_dPres;
  arrayView2d< real64 const > const m_dInternalEnergy_dTemp;

  /// Views on rock internal energy
  arrayView2d< real64 const > m_rockInternalEnergy_n;
  arrayView2d< real64 const > m_rockInternalEnergy;
  arrayView2d< real64 const > m_dRockInternalEnergy_dTemp;

};

/**
 * @class SurfaceElementBasedAssemblyKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation in SurfaceElementSubRegion
 */
class SurfaceElementBasedAssemblyKernel : public ElementBasedAssemblyKernel< SurfaceElementSubRegion, 2 >
{

public:

  using Base = ElementBasedAssemblyKernel< SurfaceElementSubRegion, 2 >;

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
  SurfaceElementBasedAssemblyKernel( globalIndex const rankOffset,
                                     string const dofKey,
                                     SurfaceElementSubRegion const & subRegion,
                                     constitutive::SingleFluidBase const & fluid,
                                     constitutive::CoupledSolidBase const & solid,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs )
#if ALLOW_CREATION_MASS
    ,
    m_creationMass( subRegion.getReference< array1d< real64 > >( SurfaceElementSubRegion::viewKeyStruct::creationMassString() ) )
#endif
  {
#if !defined(ALLOW_CREATION_MASS)
    static_assert( true, "must have ALLOW_CREATION_MASS defined" );
#endif
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    Base::computeAccumulation( ei, stack );

#if ALLOW_CREATION_MASS
    if( m_volume[ei] * m_density_n[ei][0] > 1.1 * m_creationMass[ei] )
    {
      stack.localResidual[0] += m_creationMass[ei] * 0.25;
    }
#endif
  }

protected:

#if ALLOW_CREATION_MASS
  arrayView1d< real64 const > const m_creationMass;
#endif

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
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   CellElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 2;

    ElementBasedAssemblyKernel< CellElementSubRegion, NUM_DOF >
    kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
    ElementBasedAssemblyKernel< CellElementSubRegion, NUM_DOF >::template
    launch< POLICY, ElementBasedAssemblyKernel< CellElementSubRegion, NUM_DOF > >( subRegion.size(), kernel );
  }

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
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
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   SurfaceElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    SurfaceElementBasedAssemblyKernel
      kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
    SurfaceElementBasedAssemblyKernel::launch< POLICY >( subRegion.size(), kernel );
  }


};

}

}

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMAL_ELEMENTBASEDASSEMBLYKERNEL_HPP
