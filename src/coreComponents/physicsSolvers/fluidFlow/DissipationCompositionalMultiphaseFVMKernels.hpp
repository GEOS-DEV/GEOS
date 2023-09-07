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
 * @file DissipationCompositionalMultiphaseFVMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_DISSIPATIONCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_DISSIPATIONCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP

#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"

#include "codingUtilities/Utilities.hpp"
namespace geos
{

namespace DissipationCompositionalMultiphaseFVMKernels
{

using namespace constitutive;


/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >
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

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using DissCompFlowAccessors =
    StencilAccessors< fields::flow::pressure_n , 
                      fields::flow::globalCompDensity, 
                      fields::flow::globalCompFraction,
                      fields::elementVolume >;

  using DissMultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseDensity_n,
                              fields::multifluid::phaseMassDensity,
                              fields::multifluid::phaseCompFraction_n,
                              fields::multifluid::phaseViscosity >;  
                                                                                
  using RelPermAccessors =
    StencilMaterialAccessors< RelativePermeabilityBase, fields::relperm::phaseRelPerm_n >;

  using PorosityAccessors = 
    StencilMaterialAccessors< PorosityBase, fields::porosity::porosity_n >;  

  using AbstractBase::m_dt;
  using AbstractBase::m_numPhases;
  using AbstractBase::m_hasCapPressure;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_phaseMob;
  using AbstractBase::m_dPhaseMassDens;
  using AbstractBase::m_phaseCompFrac;
  using AbstractBase::m_dPhaseCompFrac;
  using AbstractBase::m_dCompFrac_dCompDens;
  using AbstractBase::m_dPhaseCapPressure_dPhaseVolFrac;
  using AbstractBase::m_pres;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::numFluxSupportPoints;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  using Deriv = multifluid::DerivativeOffset;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] compFlowAccessor accessor for wrappers registered by the solver
   * @param[in] dissCompFlowAccessor accessor for wrappers registered by the solver needed for dissipation
   * @param[in] multiFluidAccessor accessor for wrappers registered by the multifluid model
   * @param[in] dissMultiFluidAccessor accessor for wrappers registered by the multifluid model needed for dissipation
   * @param[in] capPressureAccessors accessor for wrappers registered by the cap pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] porosityAccessors accessor for wrappers registered by the porosity model
   * @param[in] relPermAccessors accessor for wrappers registered by the relative permeability model needed for dissipation
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( integer const numPhases,
                           globalIndex const rankOffset,
                           integer const hasCapPressure,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           CompFlowAccessors const & compFlowAccessors,
                           DissCompFlowAccessors const & dissCompFlowAccessors,
                           MultiFluidAccessors const & multiFluidAccessors,
                           DissMultiFluidAccessors const & dissMultiFluidAccessors,
                           CapPressureAccessors const & capPressureAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           PorosityAccessors const & porosityAccessors,
                           RelPermAccessors const & relPermAccessors,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs,
                           real64 const omega,
                           integer const curNewton,
                           integer const continuation, 
                           integer const miscible,
                           real64 const kappamin,
                           real64 const contMultiplier )
    : Base( numPhases,
            rankOffset,
            hasCapPressure,
            //0.0,                  // no C1-PPU
            stencilWrapper,
            dofNumberAccessor,
            compFlowAccessors,
            multiFluidAccessors,
            capPressureAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs),
    m_pres_n( dissCompFlowAccessors.get( fields::flow::pressure_n {} ) ),
    m_compDens( dissCompFlowAccessors.get( fields::flow::globalCompDensity {} ) ),
    m_compFrac( dissCompFlowAccessors.get( fields::flow::globalCompFraction {} ) ),
    m_phaseRelPerm_n( relPermAccessors.get( fields::relperm::phaseRelPerm_n {} ) ),
    //m_phaseVisc_n( dissMultiFluidAccessors.get( fields::multifluid::phaseViscosity_n {} ) ),
    m_phaseDens_n( dissMultiFluidAccessors.get( fields::multifluid::phaseDensity_n {} ) ),
    m_phaseMassDens_n( dissMultiFluidAccessors.get( fields::multifluid::phaseMassDensity {} ) ),
    m_phaseCompFrac_n( dissMultiFluidAccessors.get( fields::multifluid::phaseCompFraction_n {} ) ),
    m_permeability( permeabilityAccessors.get( fields::permeability::permeability {} ) ),
    m_dPerm_dPres( permeabilityAccessors.get( fields::permeability::dPerm_dPressure {} ) ),
    m_omegaDBC(omega),
    m_curNewton(curNewton),
    m_continuationDBC(continuation),
    m_miscibleDBC(miscible),
    m_kappaminDBC(kappamin),
    m_contMultiplierDBC(contMultiplier),
    m_volume( dissCompFlowAccessors.get( fields::elementVolume {} ) ),
    m_porosity_n( porosityAccessors.get( fields::porosity::porosity_n {} ) )
  {} 

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numConnectedElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian, including dissipation
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives,
    //
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to compute dissipation terms
    Base::computeFlux( iconn, stack, [&] ( integer const ip,
                                           localIndex const (&k)[2],
                                           localIndex const (&seri)[2],
                                           localIndex const (&sesri)[2],
                                           localIndex const (&sei)[2],
                                           localIndex const connectionIndex,
                                           localIndex const k_up,
                                           localIndex const er_up,
                                           localIndex const esr_up,
                                           localIndex const ei_up,
                                           real64 const & potGrad,
                                           real64 const & phaseFlux,
                                           real64 const (&dPhaseFlux_dP)[2],
                                           real64 const (&dPhaseFlux_dC)[2][numComp] )
    {
      GEOS_UNUSED_VAR(k_up, potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, er_up, esr_up, ei_up );
      bool dz_config = true;

      /// Storing dissipation flux and its derivatives locally
      real64 dissFlux[numComp]{};
      real64 dDissFlux_dP[numFluxSupportPoints][numComp]{};
      real64 dDissFlux_dC[numFluxSupportPoints][numComp][numComp]{};
      real64 fluxPointCoef[numFluxSupportPoints] = {1.0, -1.0}; // for gradients

      real64 viscosityMult[3] = {1.0, 1.0, 1.0}; // for viscosity
      real64 refMassDensity[3] = {700, 1000, 700}; // for reference density

      /// Step 1. Calculate the continuation parameter based on the current Newton iteration
      real64 kappaDBC = 1.0; // default value

      if (m_continuationDBC) // if continuation is enabled
      {
        if (m_curNewton >= 5) 
        {
          kappaDBC = m_kappaminDBC;
        }
        else 
        {
          for (int mp = 0; mp < m_curNewton; mp++) kappaDBC *= m_contMultiplierDBC; 
          kappaDBC = std::max(kappaDBC, m_kappaminDBC);
        }
      }
      /// Step 2. Collect all contributions
      real64 poreVolume_n = 0; // Pore volume contribution
      real64 trans = 0; // Transmissibility contribution
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke ) 
      {
        trans = std::max(stack.transmissibility[connectionIndex][ke], trans);
        // to do: average porevolume
        poreVolume_n += 0.5 *  m_volume[seri[ke]][sesri[ke]][sei[ke]] * m_porosity_n[seri[ke]][sesri[ke]][sei[ke]][0];
        //poreVolume_n = std::max(m_volume[seri[ke]][sesri[ke]][sei[ke]] * m_porosity_n[seri[ke]][sesri[ke]][sei[ke]][0], poreVolume_n);
      }

      // potential gradient contribution
      // pressure
      real64 pressure_gradient = 0;
      real64 densMean = 0;
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke ) 
      {
          pressure_gradient += fluxPointCoef[ke] * m_pres_n[seri[ke]][sesri[ke]][sei[ke]];
          //densMean += 0.5 * m_phaseMassDens_n[seri[ke]][sesri[ke]][sei[ke]][0][ip];  
          densMean += 0.5 * refMassDensity[ip];
          //real64 const dDensMean_dP = 0.5 * m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];
      }
      // gravity
      real64 gravity_gradient = 0;
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke ) 
      {
          gravity_gradient += fluxPointCoef[ke] * densMean * m_gravCoef[seri[ke]][sesri[ke]][sei[ke]];
      }
      if (gravity_gradient < 0) gravity_gradient = 0;

      //real64 potential_gradient = abs(pressure_gradient) + gravity_gradient;

      //real64 potential_gradient = abs(gravity_gradient);
      real64 potential_gradient = abs(pressure_gradient);

      real64 grad_depth = 0;
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke ) 
      {
          grad_depth += fluxPointCoef[ke] * m_gravCoef[seri[ke]][sesri[ke]][sei[ke]];
      }
      
      real64 directional_coef;
      
      if (m_miscibleDBC) 
      {
        directional_coef = 100.0;
        if (abs(grad_depth) != 0) 
        {
          if (1000.f / abs(grad_depth * grad_depth) < 100)
            directional_coef = 1000.f / abs(grad_depth * grad_depth);
        }
      }
      else
      {
        directional_coef = 1.0;
      }
      
      
      //std::cout << "directional_coef = " << directional_coef << std::endl;

      // multiplier with all contributions
      real64 multiplier_n = kappaDBC * m_omegaDBC * m_dt * trans / poreVolume_n * potential_gradient;


      multiplier_n = kappaDBC * m_omegaDBC * trans / poreVolume_n * m_dt * potential_gradient * directional_coef;
      //multiplier_n = kappaDBC * m_omegaDBC * trans / poreVolume_n * m_dt * directional_coef;

      //if (m_dt < 86400) multiplier_n = 0; // don't use DBC until timestep size is at least 1 day

      //std::cout << "kappaDBC " << kappaDBC << std::endl;
      //std::cout << "multiplier_n " << multiplier_n << std::endl;

      /*      
      std::cout << "kappaDBC " << kappaDBC << std::endl;
      std::cout << "m_omegaDBC " << m_omegaDBC << std::endl;
      std::cout << "trans " << trans << std::endl;
      std::cout << "poreVolume_n " << poreVolume_n << std::endl;
      std::cout << "potential_gradient " << potential_gradient << std::endl;
      std::cout << "directional_coef = " << directional_coef << std::endl;
      std::cout << "multiplier_n " << multiplier_n << std::endl;
      std::cout << "m_contMultiplierDBC = " << m_contMultiplierDBC << std::endl;
      std::cout << "kappamin = " << m_kappaminDBC << std::endl;
      std::cout << "miscible = " << m_miscibleDBC << std::endl;
      */
      //std::cout << "gravCoef of neighbor[" << ke << "] = " << m_gravCoef[seri[ke]][sesri[ke]][sei[ke]] << std::endl;
      //std::cout << "gravCoef gradient = " << gravity_gradient_report << std::endl;
      

      

      /// Step 3. Compute the dissipation flux and its derivative
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke ) {
          for( integer ic = 0; ic < numComp; ++ic ) {
            localIndex const er  = seri[ke];
            localIndex const esr = sesri[ke];
            localIndex const ei  = sei[ke];

            // composition gradient contribution to the dissipation flux
            if (dz_config) // using z gradient
              dissFlux[ic] += multiplier_n * viscosityMult[ip] * fluxPointCoef[ke] * m_compFrac[er][esr][ei][ic];
            else // using rho_c gradeint
              dissFlux[ic] += multiplier_n * viscosityMult[ip] * fluxPointCoef[ke] * m_compDens[er][esr][ei][ic];


            dDissFlux_dP[ke][ic] = 0; 

            for( integer jc = 0; jc < numComp; ++jc ) 
            {
              // composition gradient derivative with respect to component density contribution to the dissipation flux
              if (dz_config) // using z gradient
                dDissFlux_dC[ke][ic][jc] += multiplier_n * viscosityMult[ip] * fluxPointCoef[ke] * m_dCompFrac_dCompDens[er][esr][ei][ic][jc];
              else // using rho_c gradeint
                dDissFlux_dC[ke][ic][jc] += multiplier_n * viscosityMult[ip] * fluxPointCoef[ke];
            }
          }
        }

      
      /// Step 4: add the dissipation flux and its derivatives to the residual and Jacobian
      for( integer ic = 0; ic < numComp; ++ic )
      {
        integer const eqIndex0 = k[0] * numEqn + ic;
        integer const eqIndex1 = k[1] * numEqn + ic;

        stack.localFlux[eqIndex0] += m_dt * dissFlux[ic];
        stack.localFlux[eqIndex1] -= m_dt * dissFlux[ic];

        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          localIndex const localDofIndexPres = k[ke] * numDof;
          stack.localFluxJacobian[eqIndex0][localDofIndexPres] +=  m_dt * dDissFlux_dP[ke][ic];
          stack.localFluxJacobian[eqIndex1][localDofIndexPres] -=  m_dt * dDissFlux_dP[ke][ic];

          for( integer jc = 0; jc < numComp; ++jc )
            {
              localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
              stack.localFluxJacobian[eqIndex0][localDofIndexComp] += m_dt * dDissFlux_dC[ke][ic][jc];
              stack.localFluxJacobian[eqIndex1][localDofIndexComp] -= m_dt * dDissFlux_dC[ke][ic][jc];
            }
        }
      }

    } ); // end call to Base::computeFlux

  }

protected:

  /// Views on flow properties at the previous converged time step
  ElementViewConst< arrayView1d< real64 const > > const m_pres_n;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const m_compDens;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const m_compFrac;

  ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const m_phaseRelPerm_n;
  //ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseVisc_n;

  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseDens_n;
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseMassDens_n;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const m_phaseCompFrac_n;

  ElementViewConst< arrayView3d< real64 const > > const m_permeability;
  ElementViewConst< arrayView3d< real64 const > > const m_dPerm_dPres;
  

  // DBC specific parameters
  real64 m_omegaDBC;
  integer m_curNewton;
  integer m_continuationDBC;
  integer m_miscibleDBC;
  real64 m_kappaminDBC;
  real64 m_contMultiplierDBC;


  ElementViewConst< arrayView1d< real64 const > > const m_volume;
  ElementViewConst< arrayView2d< real64 const > > const m_porosity_n;

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
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const hasCapPressure,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   //ElementSubRegionBase const & subRegion,
                   //CoupledSolidBase const & solid,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   real64 const omega,
                   integer const curNewton,
                   integer const continuation, 
                   integer const miscible,
                   real64 const kappamin,
                   real64 const contMultiplier )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );  

      using KERNEL_TYPE = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::DissCompFlowAccessors dissCompFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::DissMultiFluidAccessors dissMultiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PorosityAccessors porosityAccessors( elemManager, solverName );
      typename KERNEL_TYPE::RelPermAccessors relPermAccessors( elemManager, solverName );

      KERNEL_TYPE kernel( numPhases, rankOffset, hasCapPressure, stencilWrapper, dofNumberAccessor,
                          compFlowAccessors, dissCompFlowAccessors, multiFluidAccessors, dissMultiFluidAccessors,
                          capPressureAccessors, permeabilityAccessors, porosityAccessors, relPermAccessors,
                          dt, localMatrix, localRhs, omega, curNewton, continuation, miscible, kappamin, contMultiplier   );
      KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace DissipationCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_DISSIPATIONCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
