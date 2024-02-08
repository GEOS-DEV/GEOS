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
 * @file ThermalCompositionalMultiphaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALCOMPOSITIONALMULTIPHASEWELLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALCOMPOSITIONALMULTIPHASEWELLKERNELS_HPP

#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp"

namespace geos
{

  namespace thermalCompositionalMultiphaseWellKernels
  {

    using namespace constitutive;

    /******************************** TotalMassDensityKernel ****************************/

    /**
     * @class TotalMassDensityKernel
     * @tparam NUM_COMP number of fluid components
     * @tparam NUM_PHASE number of fluid phases
     * @brief Define the interface for the property kernel in charge of computing the total mass density
     */
    template <integer NUM_COMP, integer NUM_PHASE>
    class TotalMassDensityKernel : public compositionalMultiphaseWellKernels::TotalMassDensityKernel<NUM_COMP, NUM_PHASE>
    {
    public:
      using Base = compositionalMultiphaseWellKernels::TotalMassDensityKernel<NUM_COMP, NUM_PHASE>;
      using Base::m_dCompFrac_dCompDens;
      using Base::m_dPhaseMassDens;
      using Base::m_dPhaseVolFrac;
      using Base::m_dTotalMassDens_dCompDens;
      using Base::m_dTotalMassDens_dPres;
      using Base::m_phaseMassDens;
      using Base::m_phaseVolFrac;
      using Base::m_totalMassDens;
      using Base::numComp;
      using Base::numPhase;

      /**
       * @brief Constructor
       * @param[in] subRegion the element subregion
       * @param[in] fluid the fluid model
       */
      TotalMassDensityKernel(ObjectManagerBase &subRegion,
                             MultiFluidBase const &fluid)
          : Base(subRegion, fluid),
            m_dTotalMassDens_dTemp(subRegion.getField<fields::well::dTotalMassDensity_dTemperature>())
      {
      }

      /**
       * @brief Compute the total mass density in an element
       * @tparam FUNC the type of the function that can be used to customize the kernel
       * @param[in] ei the element index
       * @param[in] totalMassDensityKernelOp the function used to customize the kernel
       */
      template <typename FUNC = NoOpFunc>
      GEOS_HOST_DEVICE inline void compute(localIndex const ei,
                                           FUNC &&totalMassDensityKernelOp = NoOpFunc{}) const
      {
        using Deriv = multifluid::DerivativeOffset;

        arraySlice1d<real64 const, compflow::USD_PHASE - 1> phaseVolFrac = m_phaseVolFrac[ei];
        arraySlice2d<real64 const, compflow::USD_PHASE_DC - 1> dPhaseVolFrac = m_dPhaseVolFrac[ei];
        arraySlice1d<real64 const, multifluid::USD_PHASE - 2> phaseMassDens = m_phaseMassDens[ei][0];
        arraySlice2d<real64 const, multifluid::USD_PHASE_DC - 2> dPhaseMassDens = m_dPhaseMassDens[ei][0];
        real64 &totalMassDens = m_totalMassDens[ei];
        real64 &dTotalMassDens_dTemp = m_dTotalMassDens_dTemp[ei];

        // Call the base compute the compute the total mass density and derivatives
        return Base::compute(ei, [&](localIndex const ip)
                             { dTotalMassDens_dTemp += dPhaseVolFrac[ip][Deriv::dT] * phaseMassDens[ip] + phaseVolFrac[ip] * dPhaseMassDens[ip][Deriv::dT]; });
      }

    protected:
      // outputs
      arrayView1d<real64> m_dTotalMassDens_dTemp;
    };

    /**
     * @class TotalMassDensityKernelFactory
     */
    class TotalMassDensityKernelFactory
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
      template <typename POLICY>
      static void
      createAndLaunch(integer const numComp,
                      integer const numPhase,
                      ObjectManagerBase &subRegion,
                      MultiFluidBase const &fluid)
      {
        if (numPhase == 2)
        {
          isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch(numComp, [&](auto NC)
                                                                                                 {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel ); });
        }
        else if (numPhase == 3)
        {
          isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch(numComp, [&](auto NC)
                                                                                                 {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel ); });
        }
      }
    };

    /******************************** ElementBasedAssemblyKernel ********************************/

    /**
     * @class ElementBasedAssemblyKernel
     * @tparam NUM_COMP number of fluid components
     * @tparam NUM_DOF number of degrees of freedom
     * @brief Define the interface for the assembly kernel in charge of thermal accumulation and volume balance
     */
    template <localIndex NUM_COMP, localIndex NUM_DOF>
    class ElementBasedAssemblyKernel : public compositionalMultiphaseWellKernels::ElementBasedAssemblyKernel<NUM_COMP, NUM_DOF>
    {
    public:
      using Base = compositionalMultiphaseWellKernels::ElementBasedAssemblyKernel<NUM_COMP, NUM_DOF>;
      using Base::m_dCompFrac_dCompDens;
      using Base::m_dofNumber;
      using Base::m_dPhaseCompFrac;
      using Base::m_dPhaseDens;
      using Base::m_dPhaseVolFrac;
      using Base::m_dPoro_dPres;
      using Base::m_elemGhostRank;
      using Base::m_localMatrix;
      using Base::m_localRhs;
      using Base::m_numPhases;
      using Base::m_phaseCompFrac;
      using Base::m_phaseCompFrac_n;
      using Base::m_phaseDens;
      using Base::m_phaseDens_n;
      using Base::m_phaseVolFrac;
      using Base::m_phaseVolFrac_n;
      using Base::m_porosity;
      using Base::m_porosity_n;
      using Base::m_rankOffset;
      using Base::m_volume;
      using Base::numComp;
      using Base::numDof;
      using Base::numEqn;

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
      ElementBasedAssemblyKernel(localIndex const numPhases,
                                 globalIndex const rankOffset,
                                 string const dofKey,
                                 ElementSubRegionBase const &subRegion,
                                 MultiFluidBase const &fluid,
                                 CRSMatrixView<real64, globalIndex const> const &localMatrix,
                                 arrayView1d<real64> const &localRhs,
                                 BitFlags<isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags> const kernelFlags)
          : Base(numPhases, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs, kernelFlags),
            m_phaseInternalEnergy_n(fluid.phaseInternalEnergy_n()),
            m_phaseInternalEnergy(fluid.phaseInternalEnergy()),
            m_dPhaseInternalEnergy(fluid.dPhaseInternalEnergy())
      {
      }

      struct StackVariables : public Base::StackVariables
      {
      public:
        GEOS_HOST_DEVICE
        StackVariables()
            : Base::StackVariables()
        {
        }

        using Base::StackVariables::dofIndices;
        using Base::StackVariables::localJacobian;
        using Base::StackVariables::localResidual;
        using Base::StackVariables::localRow;
        using Base::StackVariables::volume;
 


      };
      /**
       * @brief Performs the setup phase for the kernel.
       * @param[in] ei the element index
       * @param[in] stack the stack variables
       */
      GEOS_HOST_DEVICE
      void setup(localIndex const ei,
                 StackVariables &stack) const
      {
        Base::setup(ei, stack);

      }

      /**
       * @brief Compute the local accumulation contributions to the residual and Jacobian
       * @tparam FUNC the type of the function that can be used to customize the kernel
       * @param[in] ei the element index
       * @param[inout] stack the stack variables
       */
      GEOS_HOST_DEVICE
      void computeAccumulation(localIndex const ei,
                               StackVariables &stack) const
      {
        using Deriv = multifluid::DerivativeOffset;

        Base::computeAccumulation(ei, stack, [&](integer const ip, real64 const &phaseAmount, real64 const &phaseAmount_n, real64 const &dPhaseAmount_dP, real64 const(&dPhaseAmount_dC)[numComp])
                                  {
      // We are in the loop over phases, ip provides the current phase index.
      // We have to do two things:
      //   1- Assemble the derivatives of the component mass balance equations with respect to temperature
      //   2- Assemble the phase-dependent part of the accumulation term of the energy equation

      real64 dPhaseInternalEnergy_dC[numComp]{};

      // construct the slices
      arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens = m_dPhaseDens[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac = m_phaseCompFrac[ei][0];
      arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac = m_dPhaseCompFrac[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy_n = m_phaseInternalEnergy_n[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy = m_phaseInternalEnergy[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseInternalEnergy = m_dPhaseInternalEnergy[ei][0];

      // Step 1: assemble the derivatives of the component mass balance equations with respect to temperature

      real64 const dPhaseAmount_dT = stack.volume * (dPhaseVolFrac[ip][Deriv::dT] * phaseDens[ip] + phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dT] );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        stack.localJacobian[ic][numDof-1] += dPhaseAmount_dT * phaseCompFrac[ip][ic]
                                             + phaseAmount * dPhaseCompFrac[ip][ic][Deriv::dT];
      }

      // Step 2: assemble the phase-dependent part of the accumulation term of the energy equation

      real64 const phaseEnergy = phaseAmount * phaseInternalEnergy[ip];
      real64 const phaseEnergy_n = phaseAmount_n * phaseInternalEnergy_n[ip];
      real64 const dPhaseEnergy_dP = dPhaseAmount_dP * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dP];
      real64 const dPhaseEnergy_dT = dPhaseAmount_dT * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dT];

      // local accumulation
      stack.localResidual[numEqn-1] += phaseEnergy - phaseEnergy_n;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[numEqn-1][0]        += dPhaseEnergy_dP;
      stack.localJacobian[numEqn-1][numDof-1] += dPhaseEnergy_dT;

      // derivatives w.r.t. component densities
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseInternalEnergy[ip], dPhaseInternalEnergy_dC, Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localJacobian[numEqn-1][jc + 1] += phaseInternalEnergy[ip] * dPhaseAmount_dC[jc]
                                                 + dPhaseInternalEnergy_dC[jc] * phaseAmount;
      } });

 
      }

      /**
       * @brief Compute the local volume balance contributions to the residual and Jacobian
       * @tparam FUNC the type of the function that can be used to customize the kernel
       * @param[in] ei the element index
       * @param[inout] stack the stack variables
       */
      GEOS_HOST_DEVICE
      void computeVolumeBalance(localIndex const ei,
                                StackVariables &stack) const
      {
        using Deriv = multifluid::DerivativeOffset;

        Base::computeVolumeBalance(ei, stack, [&](real64 const &oneMinusPhaseVolFraction)
                                   {
      GEOS_UNUSED_VAR( oneMinusPhaseVolFraction );

      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];

      for( integer ip = 0; ip < m_numPhases; ++ip )
      {
        stack.localJacobian[numEqn-2][numDof-1] -= dPhaseVolFrac[ip][Deriv::dT];
      } });
      }

      GEOS_HOST_DEVICE
      void complete(localIndex const ei,
                    StackVariables &stack) const
      {
        // Step 1: assemble the component mass balance equations and volume balance equations
        Base::complete(ei, stack);

        // Step 2: assemble the energy equation
        m_localRhs[stack.localRow + numEqn - 1] += stack.localResidual[numEqn - 1];
        m_localMatrix.template addToRow<serialAtomic>(stack.localRow + numEqn - 1,
                                                      stack.dofIndices,
                                                      stack.localJacobian[numEqn - 1],
                                                      numDof);
      }

    protected:
      /// View on derivative of porosity w.r.t temperature
      arrayView2d<real64 const> const m_dPoro_dTemp;

      /// Views on phase internal energy
      arrayView3d<real64 const, multifluid::USD_PHASE> m_phaseInternalEnergy_n;
      arrayView3d<real64 const, multifluid::USD_PHASE> m_phaseInternalEnergy;
      arrayView4d<real64 const, multifluid::USD_PHASE_DC> m_dPhaseInternalEnergy;

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
       * @param[inout] localMatrix the local CRS matrix
       * @param[inout] localRhs the local right-hand side vector
       */
      template <typename POLICY>
      static void
      createAndLaunch(localIndex const numComps,
                      localIndex const numPhases,
                      globalIndex const rankOffset,
                      integer const useTotalMassEquation,
                      string const dofKey,
                      ElementSubRegionBase const &subRegion,
                      MultiFluidBase const &fluid,
                      CRSMatrixView<real64, globalIndex const> const &localMatrix,
                      arrayView1d<real64> const &localRhs)
      {
        isothermalCompositionalMultiphaseBaseKernels::
            internal::kernelLaunchSelectorCompSwitch(numComps, [&](auto NC)
                                                     {
      localIndex constexpr NUM_COMP = NC();
      localIndex constexpr NUM_DOF = NC()+2;

      BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags::TotalMassEquation );

      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
      kernel( numPhases, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs, kernelFlags );
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template
      launch< POLICY, ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF > >( subRegion.size(), kernel ); });
      }
    };

  }; // end namespace thermalCompositionalMultiphaseWellKernels

} // end namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALCOMPOSITIONALMULTIPHASEWELLKERNELS_HPP
