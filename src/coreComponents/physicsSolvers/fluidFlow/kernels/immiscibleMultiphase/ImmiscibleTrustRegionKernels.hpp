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
 * @file ImmiscibleTrustRegionKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_TRUSTREGIONKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_TRUSTREGIONKERNELS_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluidFields.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/StencilBase.hpp"
#include "finiteVolume/TwoPointFluxApproximation.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlow.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/RelativePermeabilityUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/CapillaryPressureUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/KernelLaunchSelectors.hpp"

namespace geos
{
namespace immiscibleMultiphaseKernels
{
using namespace constitutive;

/**************************************** Trust Region Kernels ********************************************/

/****************************************** Kink Computation **********************************************/

/**********************************************************************************************************/


template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class KinkFactorKernel
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

  using ImmiscibleMultiphaseFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::pressure,
                      fields::flow::gravityCoefficient,
                      fields::immiscibleMultiphaseFlow::solutionUpdate >;

  using MultiphaseFluidAccessors =
    StencilMaterialAccessors< constitutive::TwoPhaseImmiscibleFluid,
                              fields::twophaseimmisciblefluid::phaseDensity,
                              fields::twophaseimmisciblefluid::dPhaseDensity >;

  using CapPressureAccessors =
    StencilMaterialAccessors< CapillaryPressureBase,
                              fields::cappres::phaseCapPressure,
                              fields::cappres::dPhaseCapPressure_dPhaseVolFraction >;

  using Deriv = immiscibleFlow::DerivativeOffset;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_EQN;

  /// Maximum number of elements at the face (equivalent to a connection with 4 elements, such as a fracture intersection)
  static constexpr localIndex maxNumConn = 6;

  KinkFactorKernel( integer const numPhases,
                    globalIndex const rankOffset,
                    arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localSolution ),
                    arrayView1d< real64 const > const & localResidual,
                    STENCILWRAPPER const & stencilWrapper,
                    DofNumberAccessor const & dofNumberAccessor,
                    ImmiscibleMultiphaseFlowAccessors const & flowAccessors,
                    MultiphaseFluidAccessors const & fluidAccessors,
                    CapPressureAccessors const & capPressureAccessors,                    
                    integer const hasCapPressure,
                    real64 const resNorm,
                    ImmiscibleMultiphaseFlow::TrustRegionParameters const & trustRegionParams )
    : m_numPhases( numPhases ),
      m_rankOffset( rankOffset ),
      m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
      m_ghostRank( flowAccessors.get( fields::ghostRank {} ) ),
      m_gravCoef( flowAccessors.get( fields::flow::gravityCoefficient {} ) ),
      m_pres( flowAccessors.get( fields::flow::pressure {} ) ),
      m_dens( fluidAccessors.get( fields::twophaseimmisciblefluid::phaseDensity {} ) ),
      m_dDens_dPres( fluidAccessors.get( fields::twophaseimmisciblefluid::dPhaseDensity {} ) ),
      m_phaseCapPressure( capPressureAccessors.get( fields::cappres::phaseCapPressure {} ) ),
      m_dPhaseCapPressure_dPhaseVolFrac( capPressureAccessors.get( fields::cappres::dPhaseCapPressure_dPhaseVolFraction {} ) ),
      m_localSolution ( flowAccessors.get( fields::immiscibleMultiphaseFlow::solutionUpdate {} ) ), // localSolution      
      m_localResidual( localResidual ),
      m_resNorm( resNorm ),
      m_hasCapPressure ( hasCapPressure ),
      m_dPhiMin( trustRegionParams.dPhiMin ),
      m_kinkFactorDelta( trustRegionParams.kinkFactorDelta ),
      m_minFactor( trustRegionParams.minKinkFactor ),
      m_relResThres( trustRegionParams.relResThres ),
      m_absResThres( trustRegionParams.absResThres ),
      m_stencilWrapper( stencilWrapper ),      
      m_seri( stencilWrapper.getElementRegionIndices() ),
      m_sesri( stencilWrapper.getElementSubRegionIndices() ),
      m_sei( stencilWrapper.getElementIndices() )
  {}  

  struct StackVariables
  {
    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex numElems )
      : numFluxElems( numElems ),
        nConn( 0 ),
        localEps( numElems * (numElems - 1) / 2 * numEqn ),   // Number of connections * number of equations
        connFactor( 1.0 )
    {}

    /// Number of elements for a given connection
    localIndex const numFluxElems;    

    /// Number of pairwise connections
    localIndex nConn; 

    /// Storage for the face local damping factors
    stackArray1d< real64, maxNumConn * numEqn > localEps;

    /// Storage for the minimum face local damping factor
    real64 connFactor;
  };

  /**
   * @brief Getter for the number of elements at this connection
   * @param[in] iconn the connection index
   * @return the number of elements at this connection
   */
  GEOS_HOST_DEVICE
  localIndex numPointsInFlux( localIndex const iconn ) const
  { return m_stencilWrapper.numPointsInFlux( iconn ); }

  /**
   * @brief Compute the local kink damping factors for the connection
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */  
  GEOS_HOST_DEVICE
  void computeKink( localIndex const iconn,
                    StackVariables & stack ) const
  {
    localIndex k[2];   
    localIndex connectionIndex = 0; 

    // loop over connection elements
    for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
      {
        // clear working arrays
        real64 densMean[numEqn]{};
        
        real64 dPhi_dS[numEqn][2]{};
        real64 dPhi_dP[numEqn][2]{};

        real64 dPhi[numEqn]{};
        real64 dDPhi[numEqn]{};
        real64 eps[numEqn]{};

        real64 dP[2]{};
        real64 dS[2]{};

        localIndex localRow[2]{-1, -1};

        // cell indices
        localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // get pressure and saturation updates for connection elements
        // for ( integer ke = 0; ke < 2; ++ke )
        // {
        //   globalIndex const globalRow = m_dofNumber[seri[ke]][sesri[ke]][sei[ke]];
        //   localRow[ke] = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        //   GEOS_ASSERT_GE( localRow[ke], 0 );

        //   dP[ke] = m_localSolution[localRow[ke]];     // dP = { dP1 , dP2 }
        //   dS[ke] = m_localSolution[localRow[ke] + 1]; // dS = { dS1 , dS2 }          
        // }                                                                        
        for ( integer ke = 0; ke < 2; ++ke )
        {
          dP[ke] = m_localSolution[seri[ke]][sesri[ke]][sei[ke]][0];     // dP = { dP1 , dP2 }
          dS[ke] = m_localSolution[seri[ke]][sesri[ke]][sei[ke]][1];     // dS = { dS1 , dS2 } 
          
          if ( m_ghostRank[seri[ke]][sesri[ke]][sei[ke]] < 0 )
          {
            globalIndex const globalRow = m_dofNumber[seri[ke]][sesri[ke]][sei[ke]];
            localRow[ke] = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
            GEOS_ASSERT_GE( localRow[ke], 0 );
          }
        }

        // loop over phases
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          // adaptive damping based on significant residual values           
          // if( (fabs( m_localResidual[localRow[0] + ip] ) < m_relResThres * m_resNorm || fabs( m_localResidual[localRow[0] + ip] ) < m_absResThres) &&
          //     (fabs( m_localResidual[localRow[1] + ip] ) < m_relResThres * m_resNorm || fabs( m_localResidual[localRow[1] + ip] ) < m_absResThres) )            
          // {
          //   continue;
          // }
          bool skipConnection = true;
          for( integer ke = 0; ke < 2; ++ke )
          {
            if ( localRow[ke] >= 0 )
            {              
              if( fabs( m_localResidual[localRow[ke] + ip] ) > m_relResThres * m_resNorm || 
                  fabs( m_localResidual[localRow[ke] + ip] ) > m_absResThres )
              {
                skipConnection = false;
                break;
              }
            }
          }
          if ( skipConnection )
          {
            continue;
          }

          constexpr int signPotDiff[2] = {1, -1};

          // compute average density and derivatives
          for( integer ke = 0; ke < 2; ++ke )
          {
            real64 const density  = m_dens[seri[ke]][sesri[ke]][sei[ke]][0][ip];                    // r = rho1 || rho2
            real64 const dDens_dP = m_dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0][ip][Deriv::dP];  // dr/dP = drho1/dP1 || drho2/dP2
           
            densMean[ip] += 0.5 * density;                           // rho = (rho1 + rho2) / 2            
            
            dPhi_dP[ip][ke] = -0.5 * dDens_dP;     // dPhi/dP = { -0.5 * drho1/dP1 , -0.5 * drho2/dP2 }
          }

          // compute potential difference and potential difference update    
          for( integer ke = 0; ke < 2; ++ke )
          {
            localIndex const er  = seri[ke];
            localIndex const esr = sesri[ke];
            localIndex const ei  = sei[ke];

            real64 const pressure = m_pres[er][esr][ei];       // P = P1 || P2
            real64 const gravD = m_gravCoef[er][esr][ei];      // D = g z1 || g z2   

            dPhi_dP[ip][ke] = signPotDiff[ke] + dPhi_dP[ip][ke] * gravD;     // dPhi/dP = { 1 - 0.5 * drho1/dP1 g z1 , -1 - 0.5 * drho2/dP2 g z2 }      
            
            real64 pot = pressure - densMean[ip] * gravD;      // Phi = P1 - rho g z1 || P2 - rho g z2           

            if( m_hasCapPressure )
            {              
              pot -= m_phaseCapPressure[er][esr][ei][0][ip];   // Phi = P1 - rho g z1 - Pc1 || P2 - rho g z2 - Pc2

              dPhi_dS[ip][ke] = -signPotDiff[ke] * m_dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][ip];   // dPhi/dS = { -dPc1/dS1 , dPc2/dS2 }
            }

            dPhi[ip] += signPotDiff[ke] * pot;                                  // dPhi = P1 - P2 - rho g (z1 - z2) - Pc1 + Pc2
            dDPhi[ip] += dPhi_dP[ip][ke] * dP[ke] + dPhi_dS[ip][ke] * dS[ke];   // dDPhi = (1 - 0.5 * drho1/dP1 g z1)  dP1 + (-dPc1/dS1) dS1
                                                                                //       + (-1 - 0.5 * drho2/dP2 g z2) dP2 + (dPc2/dS2)  dS2
          }

          // compute kink damping factor
          if ( signbit( dPhi[ip] + dDPhi[ip]) != signbit( dPhi[ip] ) && 
               fabs( dPhi[ip] )                 > m_dPhiMin && 
               fabs( dPhi[ip] + dDPhi[ip] )     > m_dPhiMin )
          {
            eps[ip] = fmax( fmin( -dPhi[ip] / dDPhi[ip] * m_kinkFactorDelta, 1.0 ), m_minFactor );
          }
          else
          {
            eps[ip] = 1.0;
          }

          stack.localEps[connectionIndex * numEqn + ip] = eps[ip];
          stack.nConn++;

        } // loop over phases
        
        connectionIndex++;        
      }
    } // loop over connection elements
    
  }

  /**
   * @brief Return true for negative numbers and false for non-negative ones
   * @param[in] x number to check sign
   */                                 
  GEOS_HOST_DEVICE
  bool signbit( real64 x ) const
  {
    return x < 0 ? 1 : 0;
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */                                 
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  { 
    GEOS_UNUSED_VAR( iconn );
    for( integer ic = 0; ic < stack.nConn; ++ic )
    {
      stack.connFactor = LvArray::math::min( stack.connFactor, stack.localEps[ic] );  
    }   
  }  

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numConnections the number of connections
   * @param[inout] kernelComponent the kernel component providing access to the compute function and stack variables
   * @param[in] kinkFactor the minimum kink damping factor for the current stencil
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void 
  launch( localIndex const numConnections,
          KERNEL_TYPE const & kernelComponent,
          real64 & kinkFactor )
  {
    GEOS_MARK_FUNCTION;

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minKinkFactor( 1.0 );

    forAll< POLICY >( numConnections, [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent.numPointsInFlux( iconn ) );
      
      kernelComponent.computeKink( iconn, stack );
      kernelComponent.complete( iconn, stack );

      minKinkFactor.min( stack.connFactor );
    } );
    
    kinkFactor = minKinkFactor.get();
  }


protected:    

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;
  
  /// View on pressure 
  ElementViewConst< arrayView1d< real64 const > > const m_pres;

  /// Views on fluid density
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_dens;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dDens_dPres;

  /// Views on capillary pressure
  ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const m_phaseCapPressure;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac;  

  /// View on the local solution and residual   
  ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_localSolution;   
  arrayView1d< real64 const > const m_localResidual;

  /// Residual norm for adaptive residual analysis
  real64 const m_resNorm;  

  /// Flags
  integer const m_hasCapPressure;

  /// Trust region parameters
  /// Minimum potential difference to apply a damping factor
  real64 const m_dPhiMin; //1.0, 7000

  /// Stretching factor applied to damping factor to allow for a small crossing of the kink
  real64 const m_kinkFactorDelta; // 1.05

  /// Minimum dampining factor
  real64 const m_minFactor; // 0.1

   /// Minimum relative residual threshold
  real64 const m_relResThres; // 0.0, 0.05, 0.2

  /// Minimum absolute residual threshold
  real64 const m_absResThres; // 1e2, 0.0, 1e3

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
};

/**
 * @class ResidualNormKernelFactory
 */
class KinkFactorKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localSolution the solution vector on my MPI rank
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] hasCapPressure flag for presence of capillary pressure
   * @param[inout] kinkFactor the kink damping factor on the subRegion
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localSolution,
                   arrayView1d< real64 const > const & localResidual,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper, 
                   integer const hasCapPressure,
                   real64 const resNorm,  
                   ImmiscibleMultiphaseFlow::TrustRegionParameters const trustRegionParams,               
                   real64 & kinkFactor )
  {
    integer constexpr NUM_EQN = 2;
    integer constexpr NUM_DOF = 2;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = KinkFactorKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
    typename kernelType::ImmiscibleMultiphaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::MultiphaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );

    kernelType kernel( numPhases, rankOffset, localSolution, localResidual, stencilWrapper, 
                       dofNumberAccessor, flowAccessors, fluidAccessors, capPressureAccessors,
                       hasCapPressure, resNorm, trustRegionParams );
    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel, kinkFactor );    
  }
};

/*****************************************************************************************************************/

/**************************************** Flux Inflection Computation ********************************************/

/*****************************************************************************************************************/


template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER, typename FLUIDWRAPPER, typename RELPERMWRAPPER, typename CAPPRESWRAPPER >
class FluxInflectionFactorKernel
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

  using ImmiscibleMultiphaseFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::pressure,
                      fields::flow::gravityCoefficient,
                      fields::immiscibleMultiphaseFlow::phaseVolumeFraction >;

  using MultiphaseFluidAccessors =
    StencilMaterialAccessors< constitutive::TwoPhaseImmiscibleFluid,
                              fields::twophaseimmisciblefluid::phaseDensity,
                              fields::twophaseimmisciblefluid::dPhaseDensity,
                              fields::twophaseimmisciblefluid::phaseViscosity,
                              fields::twophaseimmisciblefluid::dPhaseViscosity >;

  using RelPermAccessors =
    StencilMaterialAccessors< RelativePermeabilityBase,
                              fields::relperm::phaseRelPerm,
                              fields::relperm::dPhaseRelPerm_dPhaseVolFraction >;                  

  using CapPressureAccessors =
    StencilMaterialAccessors< CapillaryPressureBase,
                              fields::cappres::phaseCapPressure,
                              fields::cappres::dPhaseCapPressure_dPhaseVolFraction >;

  using Deriv = immiscibleFlow::DerivativeOffset;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_EQN;

  /// Maximum number of elements at the face (equivalent to a connection with 4 elements)
  static constexpr localIndex maxNumConn = 6;

  /// Maximum number of nonlinear iterations
  static constexpr integer m_maxIter = 5;

  /// Minimum potential difference to apply a damping factor
  static constexpr real64 m_d2FMin = 1.0; // 7000

  /// Minimum dampining factor
   static constexpr real64 m_minFactor = 0.1;


  FluxInflectionFactorKernel( integer const numPhases,
                              globalIndex const rankOffset,
                              arrayView1d< real64 const > const & localSolution,
                              real64 const globalKinkFactor,
                              STENCILWRAPPER const & stencilWrapper,
                              FLUIDWRAPPER const & fluidWrapper,
                              RELPERMWRAPPER const & relPermWrapper,
                              CAPPRESWRAPPER const * capPressureWrapper,
                              DofNumberAccessor const & dofNumberAccessor,
                              ImmiscibleMultiphaseFlowAccessors const & flowAccessors,
                              MultiphaseFluidAccessors const & fluidAccessors,
                              RelPermAccessors const & relPermAccessor,
                              CapPressureAccessors const & capPressureAccessors,                    
                              integer const hasCapPressure,
                              ImmiscibleMultiphaseFlow::TrustRegionParameters const & trustRegionParams )
    : m_numPhases( numPhases ),
      m_rankOffset( rankOffset ),
      m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
      m_ghostRank( flowAccessors.get( fields::ghostRank {} ) ),
      m_gravCoef( flowAccessors.get( fields::flow::gravityCoefficient {} ) ),
      m_pres( flowAccessors.get( fields::flow::pressure {} ) ),
      m_phaseVolFrac( flowAccessors.get( fields::immiscibleMultiphaseFlow::phaseVolumeFraction {} ) ),
      m_dens( fluidAccessors.get( fields::twophaseimmisciblefluid::phaseDensity {} ) ),
      m_dDens_dPres( fluidAccessors.get( fields::twophaseimmisciblefluid::dPhaseDensity {} ) ),
      m_phaseVisc( fluidAccessors.get( fields::twophaseimmisciblefluid::phaseViscosity {} ) ),
      m_dPhaseVisc( fluidAccessors.get( fields::twophaseimmisciblefluid::dPhaseViscosity {} ) ),
      m_phaseRelPerm( relPermAccessor.get( fields::relperm::phaseRelPerm {} ) ),
      m_dPhaseRelPerm_dPhaseVolFrac( relPermAccessor.get( fields::relperm::dPhaseRelPerm_dPhaseVolFraction {} ) ),
      m_phaseCapPressure( capPressureAccessors.get( fields::cappres::phaseCapPressure {} ) ),
      m_dPhaseCapPressure_dPhaseVolFrac( capPressureAccessors.get( fields::cappres::dPhaseCapPressure_dPhaseVolFraction {} ) ),
      m_localSolution( localSolution ),
      m_globalKinkFactor( globalKinkFactor ),  
      m_hasCapPressure( hasCapPressure ),
      m_trustRegionParams( trustRegionParams ),
      m_stencilWrapper( stencilWrapper ),
      m_fluidWrapper( fluidWrapper ),
      m_relPermWrapper( relPermWrapper ),   
      m_capPressureWrapper( capPressureWrapper ),
      m_seri( stencilWrapper.getElementRegionIndices() ),
      m_sesri( stencilWrapper.getElementSubRegionIndices() ),
      m_sei( stencilWrapper.getElementIndices() )
  {}  

  struct StackVariables
  {
    /**
     * @brief Constructor for the stack variables
     * @param[in] numElems number of elements in this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex numElems )
      : numFluxElems( numElems ),
        nConn ( numElems * (numElems - 1) / 2 ),
        localEps( numElems * (numElems - 1) / 2 * numEqn ),
        connFactor( 1.0 ),
        dPhi( numElems * (numElems - 1) / 2 * numEqn )
    {}   

    /// Number of elements for a given connection
    localIndex const numFluxElems;    

    /// Number of pairwise connections
    localIndex const nConn; 

    /// Storage for the face local damping factors
    stackArray1d< real64, maxNumConn * numEqn > localEps;

    /// Storage for the minimum face local damping factor
    real64 connFactor;

    /// Storage for the face potential differences
    stackArray1d< real64, maxNumConn * numEqn > dPhi;
  };

  /**
   * @brief Getter for the number of elements at this connection
   * @param[in] iconn the connection index
   * @return the number of elements at this connection
   */
  GEOS_HOST_DEVICE
  localIndex numPointsInFlux( localIndex const iconn ) const
  { return m_stencilWrapper.numPointsInFlux( iconn ); }

  /**
   * @brief Compute the local inflection damping factors for the connection
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */  
  GEOS_HOST_DEVICE
  void computeInflection( localIndex const iconn,
                          StackVariables & stack ) const
  {
    localIndex k[2];   
    localIndex connectionIndex = 0; 

    // loop over connection elements
    for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
      {        
        real64 densMean[numEqn]{}; 
        real64 compressibility[numEqn][2]{};
        real64 viscosibility[numEqn][2]{};

        real64 dPhi[numEqn]{};        
        real64 d2F[numEqn][2]{};
        real64 eps[numEqn]{};       

        // cell indices
        localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // loop over phases
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {          
          real64 phaseEps[2] {0.0, 1.0}; 
          real64 newD2F {0.0}; 
          eps[ip] = 1.0;    

          constexpr int signPotDiff[2] = {1, -1};  

          // compute average density, compressibility and viscosibility
          for( integer ke = 0; ke < 2; ++ke )
          {
            real64 const density  = m_dens[seri[ke]][sesri[ke]][sei[ke]][0][ip];                     // r = rho1 || rho2
            real64 const dDens_dP = m_dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0][ip][Deriv::dP];   // dr/dP = dr1/dP1 || dr2/dP

            real64 const viscosity = m_phaseVisc[seri[ke]][sesri[ke]][sei[ke]][0][ip];               // mu = mu1 || mu2
            real64 const dVisc_dP  = m_dPhaseVisc[seri[ke]][sesri[ke]][sei[ke]][0][ip][Deriv::dP];   // dmu/dP = dmu1/dP1 || dmu2/dP

            densMean[ip] += 0.5 * density;                                                           // rho = (rho1 + rho2) / 2  
            compressibility[ip][ke] = dDens_dP / density;                                            // cf = drho1 / rho1 || drho2 / rho2     
            viscosibility[ip][ke] = dVisc_dP / viscosity;                                            // cmu = dmu1 / mu1  || dmu2 / mu2
          }          

          // compute potential difference before update  
          for( integer ke = 0; ke < 2; ++ke )
          {
            real64 const pressure = m_pres[seri[ke]][sesri[ke]][sei[ke]];       // P = P1 || P2
            real64 const gravD = m_gravCoef[seri[ke]][sesri[ke]][sei[ke]];      // D = g z1 || g z2 

            real64 pot = pressure - densMean[ip] * gravD;      // Phi = P1 - rho g z1 || P2 - rho g z2

            if( m_hasCapPressure )
            {              
              pot -= m_phaseCapPressure[seri[ke]][sesri[ke]][sei[ke]][0][ip];   // Phi = P1 - rho g z1 - Pc1 || P2 - rho g z2 - Pc2
            }

            dPhi[ip] += signPotDiff[ke] * pot;                 // dPhi = P1 - P2 - rho g (z1 - z2) - Pc1 + Pc2
          }

          // compute directional derivatives before and after update
          computeDirectionalDerivative(iconn, k, ip, compressibility[ip], viscosibility[ip], dPhi[ip], phaseEps[0], d2F[ip][0]);
          computeDirectionalDerivative(iconn, k, ip, compressibility[ip], viscosibility[ip], dPhi[ip], phaseEps[1], d2F[ip][1]);

          // search for root of second derivative if there is an inflection surface in the Newton path
          for ( integer iter = 0; iter < m_maxIter &&
                                  signbit( d2F[ip][0] ) != signbit( d2F[ip][1] ) &&
                                  fabs( d2F[ip][0] ) > m_d2FMin &&
                                  fabs( d2F[ip][1] ) > m_d2FMin; ++iter )
          {
            real64 slope = (d2F[ip][1] - d2F[ip][0]) / (phaseEps[1] - phaseEps[0]);
            eps[ip] = phaseEps[1] - d2F[ip][1] / slope;

            computeDirectionalDerivative(iconn, k, ip, compressibility[ip], viscosibility[ip], dPhi[ip], eps[ip], newD2F);

            if ( signbit( d2F[ip][0] ) == signbit( newD2F ) )
            {
              d2F[ip][0] = newD2F;
              phaseEps[0] = eps[ip];
            } 
            else
            {
              d2F[ip][1] = newD2F;
              phaseEps[1] = eps[ip];
            }
          }

          stack.localEps[connectionIndex * numEqn + ip] = eps[ip];
          
        }  // loop over phases

        connectionIndex++;
      }
    } // loop over connection elements
    
  }


  /**
   * @brief Compute the directional second derivative of the phase flux function
   *        at a specific location in the solution space, given by the location 
   *        at the previous iteration plus a restricted update
   * @param[in] iconn the connection index
   * @param[in] ip the phase index
   * @param[in] dPhi phase potential difference before update   * 
   * @param[in] eps restriction factor that determines where in the update direction to evaluate
   *                the directional second derivative 
   * @param[out] d2F directional second derivative at given location in the solution space
   */                                 
  GEOS_HOST_DEVICE
  void computeDirectionalDerivative( localIndex const iconn,
                                     localIndex const * k,
                                     localIndex const ip,
                                     real64 const * comp,
                                     real64 const * visc,
                                     real64 const dPhiv,                                     
                                     real64 const eps,
                                     real64 & d2F ) const
  {
    localIndex up = dPhiv >= 0.0 ? 0 : 1;
    localIndex down = 1 - up;
    
    // clear working variables
    real64 dens[2]{};   
    real64 gravD[2]{};  
    real64 capPres[2]{};  
    real64 dcapPres[2]{};
    real64 d2capPres[2]{};  

    real64 dPhi{};
    real64 D{};
    real64 perm{}; 
    real64 dperm{}; 
    real64 d2perm{};

    real64 pressure[2]{};
    real64 phaseVolFrac[2]{};

    real64 dP[2]{};
    real64 dS[2]{};    

    // cell indices
    localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
    localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
    localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

    // get pressure and saturation updates for connection elements
    for ( integer ke = 0; ke < 2; ++ke )
    {
      globalIndex const globalRow = m_dofNumber[seri[ke]][sesri[ke]][sei[ke]];
      localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
      GEOS_ASSERT_GE( localRow, 0 );      

      dP[ke] = m_localSolution[localRow] * m_globalKinkFactor;     // dP = { dP1 , dP2 }
      dS[ke] = m_localSolution[localRow + 1] * m_globalKinkFactor; // dS = { dS1 , dS2 }    // saturation change for nonspecified phase should be negative of dS
    }

    // compute saturation of updated state
    for ( integer ke = 0; ke < 2; ++ke )
    {  
      phaseVolFrac[ke] = ip == 0 ? m_phaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][ip] + eps * dS[ke] :  // S = S + e.dS
                                   m_phaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][ip] - eps * dS[ke] ;
             
    }

    // get grav coef, rel perm and capillary pressure + derivatives      
    for ( integer ke = 0; ke < 2; ++ke )
    {      
      gravD[ke] = m_gravCoef[seri[ke]][sesri[ke]][sei[ke]];  // D = { g z1 , g z2 }      

      if ( m_hasCapPressure )
      {
        (*m_capPressureWrapper).compute(phaseVolFrac[ke], ip, capPres[ke], dcapPres[ke], d2capPres[ke]);
      }      
    }
    m_relPermWrapper.compute(phaseVolFrac[up], ip, perm, dperm, d2perm);
    D = gravD[up] - gravD[down];   // D = g (zup - zdown)


    // get density and potential difference
    if ( eps < 0.001 )
    {
      // get density     
      for ( integer ke = 0; ke < 2; ++ke )
      {
        dens[ke] = m_dens[seri[ke]][sesri[ke]][sei[ke]][0][ip];       // r = { r1 , r2 }            
      }

      // get potential differece
      dPhi = fabs( dPhiv );
    }
    else // (0 < eps < 1)
    { 
      // compute pressure and density of updated state
      for ( integer ke = 0; ke < 2; ++ke )
      {
        pressure[ke] = m_pres[seri[ke]][sesri[ke]][sei[ke]] + eps * dP[ke];  // P = P1 || P2

        real64 dDens;       
        m_fluidWrapper.compute(pressure[ke], ip, dens[ke], dDens);  // r = { r1 , r2 } , dr = { dr1 , dr2 }             
      }    

      // compute potential difference
      constexpr int signPotDiff[2] = {1, -1}; 

      real64 densMean = (dens[0] + dens[1]) / 2;          // rho = (r1 + r2) / 2
      
      for( integer ke = 0; ke < 2; ++ke )
      {    
        real64 pot = pressure[ke] - densMean * gravD[ke]; // Phi = P1 - rho g z1 || P2 - rho g z2

        if( m_hasCapPressure )
        {              
          pot -= capPres[ke];                             // Phi = P1 - rho g z1 - Pc1 || P2 - rho g z2 - Pc2
        }

        dPhi += signPotDiff[ke] * pot;                    // Phi = P1 - P2 - rho g (z1 - z2) - Pc1 + Pc2
      }
      dPhi = fabs( dPhi );    
    }


    // compute partial second derivatives
    real64 d2Pi2 = perm * pow( comp[up] - visc[up], 2.0) * dPhi + 2 * perm * (comp[up] - visc[up]) + perm * comp[up] * (visc[up] - 1.5 * comp[up]) * dens[up] * D;
    real64 d2Pj2 = -0.5 * perm * pow( comp[down], 2.0 ) * dens[down] * D;
    real64 d2Si2 = d2perm * dPhi;
    real64 d2Sj2 = 0;
    real64 d2PiPj = perm * (comp[up] - visc[up]) * (-1 - 0.5 * comp[down] * dens[down] * D);
    real64 d2PiSi = dperm * (dPhi * (comp[up] - visc[up]) + 1 - 0.5 * comp[up] * dens[up] * D);
    real64 d2PiSj = 0;
    real64 d2PjSi = dperm * (-1 - 0.5 * comp[down] * dens[down] * D);
    real64 d2SiSj = 0;

    if ( m_hasCapPressure )
    {
      d2Si2 -= (2 * dperm * dcapPres[up] + perm * d2capPres[up]);
      d2Sj2 += perm * d2capPres[down];
      d2PiSi -= perm * (comp[up] - visc[up]) * dcapPres[up];
      d2PiSj += perm * (comp[up] - visc[up]) * dcapPres[down];
      d2SiSj += dperm * dcapPres[down];
    }

    // compute directional second derivative
    d2F = d2Pi2 * pow( dP[up], 2.0 ) + d2Pj2 * pow( dP[down], 2.0 ) + d2Si2 * pow( dS[up], 2.0 ) + d2Sj2 * pow( dS[down], 2.0 )
        + d2PiPj * 2 * dP[up] * dP[down] + d2PiSi * 2  * dP[up] * dS[up] + d2PiSj * 2  * dP[up] * dS[down] + d2PjSi * 2  * dP[down] * dS[up] + d2SiSj * 2  * dS[up] * dS[down];
  }

  /**
   * @brief Return true for negative numbers and false for non-negative ones
   * @param[in] x number to check sign
   */                                 
  GEOS_HOST_DEVICE
  bool signbit( real64 x ) const
  {
    return x < 0 ? 1 : 0;
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */                                 
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( iconn );
    for( integer ic = 0; ic < stack.nConn; ++ic )
    {
      for (integer ip = 0; ip < m_numPhases; ++ip )
      {
        stack.connFactor = LvArray::math::min( stack.connFactor, stack.localEps[ic * numEqn + ip] );
      }  
    }
    stack.connFactor = LvArray::math::max( stack.connFactor, m_minFactor ); 
  }  

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numConnections the number of connections
   * @param[inout] kernelComponent the kernel component providing access to the compute function and stack variables
   * @param[in] inflectionFactor the minimum inflection damping factor for the current stencil
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void 
  launch( localIndex const numConnections,
          KERNEL_TYPE const & kernelComponent,
          real64 & inflectionFactor )
  {
    GEOS_MARK_FUNCTION;

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minInflectionFactor( 1.0 );

    forAll< POLICY >( numConnections, [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent.numPointsInFlux( iconn ) );
      
      kernelComponent.computeInflection( iconn, stack );
      kernelComponent.complete( iconn, stack );

      minInflectionFactor.min( stack.connFactor );
    } );

    inflectionFactor = minInflectionFactor.get();
  }


protected:    

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;
  
  /// View on pressure and phase volume fraction
  ElementViewConst< arrayView1d< real64 const > > const m_pres;
  ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_phaseVolFrac;

  /// Views on fluid density
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_dens;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dDens_dPres;

  /// Views on the phase viscosities
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_phaseVisc;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dPhaseVisc;

  /// Views on the phase relative permeabilities
  ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const m_phaseRelPerm;
  ElementViewConst< arrayView4d< real64 const, relperm::USD_RELPERM_DS > > const m_dPhaseRelPerm_dPhaseVolFrac;
  ElementViewConst< arrayView4d< real64 const, relperm::USD_RELPERM_DS > > const m_d2PhaseRelPerm_d2PhaseVolFrac;

  /// Views on capillary pressure
  ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const m_phaseCapPressure;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_d2PhaseCapPressure_d2PhaseVolFrac;  

  /// View on the local solution
  arrayView1d< real64 const > const m_localSolution;

  /// Kink factor computed previously
  real64 const m_globalKinkFactor;

  /// Flags
  integer const m_hasCapPressure;

  /// Trust region parameters
  ImmiscibleMultiphaseFlow::TrustRegionParameters const m_trustRegionParams;

  /// Reference to the wrappers
  STENCILWRAPPER const m_stencilWrapper;
  FLUIDWRAPPER const m_fluidWrapper;
  RELPERMWRAPPER const m_relPermWrapper;
  CAPPRESWRAPPER const * m_capPressureWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
};

/**
 * @class FluxInflectionFactorKernelFactory
 */
class FluxInflectionFactorKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localSolution the residual vector on my MPI rank
   * @param[in] globalKinkFactor kink restriction factor already computed
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] relPermWrapper reference to the relative permeability wrapper
   * @param[in] capPressureWrapper pointer to the capillary pressure wrapper
   * @param[in] hasCapPressure flag for presence of capillary pressure
   * @param[inout] inflectionFactor the inflection damping factor on the subRegion
   */
  template< typename POLICY, typename STENCILWRAPPER, typename FLUIDWRAPPER, typename RELPERMWRAPPER, typename CAPPRESWRAPPER >
  static void
  createAndLaunch( integer const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localSolution,
                   real64 const globalKinkFactor,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   FLUIDWRAPPER const & fluidWrapper,
                   RELPERMWRAPPER const & relPermWrapper,
                   CAPPRESWRAPPER const * capPressureWrapper,
                   integer const hasCapPressure,
                   ImmiscibleMultiphaseFlow::TrustRegionParameters const trustRegionParams,                 
                   real64 & inflectionFactor )
  {
    integer constexpr NUM_EQN = 2;
    integer constexpr NUM_DOF = 2;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = FluxInflectionFactorKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER, FLUIDWRAPPER, RELPERMWRAPPER, CAPPRESWRAPPER >;
    typename kernelType::ImmiscibleMultiphaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::MultiphaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::RelPermAccessors relPermAccessor( elemManager, solverName );
    typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );

    kernelType kernel( numPhases, rankOffset, localSolution, globalKinkFactor,
                       stencilWrapper, fluidWrapper, relPermWrapper, capPressureWrapper,
                       dofNumberAccessor, flowAccessors, fluidAccessors, relPermAccessor,
                       capPressureAccessors, hasCapPressure, trustRegionParams );
    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel, inflectionFactor );    
  }
};


/*****************************************************************************************************************/

/************************************** Residual Inflection Computation ******************************************/

/*****************************************************************************************************************/

template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER, typename FLUIDWRAPPER, typename RELPERMWRAPPER, typename CAPPRESWRAPPER >
class ResidualInflectionFactorKernel
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

  using ImmiscibleMultiphaseFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::pressure,
                      fields::flow::gravityCoefficient,
                      fields::immiscibleMultiphaseFlow::phaseVolumeFraction,
                      fields::immiscibleMultiphaseFlow::solutionUpdate >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  using MultiphaseFluidAccessors =
    StencilMaterialAccessors< TwoPhaseImmiscibleFluid,
                              fields::twophaseimmisciblefluid::phaseDensity,
                              fields::twophaseimmisciblefluid::dPhaseDensity,
                              fields::twophaseimmisciblefluid::phaseViscosity,
                              fields::twophaseimmisciblefluid::dPhaseViscosity >;

  using RelPermAccessors =
    StencilMaterialAccessors< RelativePermeabilityBase,
                              fields::relperm::phaseRelPerm,
                              fields::relperm::dPhaseRelPerm_dPhaseVolFraction >;    
                              
  using PorosityAccessors = 
    StencilMaterialAccessors< PorosityBase,
                              fields::porosity::porosity,
                              fields::porosity::dPorosity_dPressure >;                     

  using CapPressureAccessors =
    StencilMaterialAccessors< CapillaryPressureBase,
                              fields::cappres::phaseCapPressure,
                              fields::cappres::dPhaseCapPressure_dPhaseVolFraction >;

  using Deriv = immiscibleFlow::DerivativeOffset;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_EQN;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumConn = 6;

  ResidualInflectionFactorKernel( integer const numPhases,
                                  globalIndex const rankOffset,
                                  string const dofKey,
                                  arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localSolution ),
                                  arrayView1d< real64 const > const & localResidual,
                                  real64 const resNorm,
                                  real64 const globalKinkFactor,
                                  CellElementSubRegion const & subRegion,
                                  constitutive::CoupledSolidBase const & solid,
                                  STENCILWRAPPER const & stencilWrapper,
                                  FLUIDWRAPPER const & fluidWrapper,
                                  RELPERMWRAPPER const & relPermWrapper,
                                  CAPPRESWRAPPER const * capPressureWrapper,
                                  DofNumberAccessor const & dofNumberAccessor,
                                  PermeabilityAccessors const & permeabilityAccessors,
                                  ImmiscibleMultiphaseFlowAccessors const & flowAccessors,
                                  MultiphaseFluidAccessors const & fluidAccessors,
                                  RelPermAccessors const & relPermAccessors,
                                  PorosityAccessors const & poroAccessors,
                                  CapPressureAccessors const & capPressureAccessors,                    
                                  integer const hasCapPressure,
                                  ImmiscibleMultiphaseFlow::TrustRegionParameters const trustRegionParams )
    : m_numPhases( numPhases ),
      m_rankOffset( rankOffset ),
      m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
      m_dofNumberElem( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
      m_ghostRank( subRegion.ghostRank() ), //  flowAccessors.get( fields::ghostRank {} )
      m_connMap( subRegion.getConnectionMap() ),
      m_gravCoef( flowAccessors.get( fields::flow::gravityCoefficient {} ) ),
      m_pres( flowAccessors.get( fields::flow::pressure {} ) ),
      m_phaseVolFrac( flowAccessors.get( fields::immiscibleMultiphaseFlow::phaseVolumeFraction {} ) ),
      m_permeability( permeabilityAccessors.get( fields::permeability::permeability {} ) ),
      m_dPerm_dPres( permeabilityAccessors.get( fields::permeability::dPerm_dPressure {} ) ),
      m_dens( fluidAccessors.get( fields::twophaseimmisciblefluid::phaseDensity {} ) ),
      m_dDens_dPres( fluidAccessors.get( fields::twophaseimmisciblefluid::dPhaseDensity {} ) ),
      m_phaseVisc( fluidAccessors.get( fields::twophaseimmisciblefluid::phaseViscosity {} ) ),
      m_dPhaseVisc( fluidAccessors.get( fields::twophaseimmisciblefluid::dPhaseViscosity {} ) ),
      m_phaseRelPerm( relPermAccessors.get( fields::relperm::phaseRelPerm {} ) ),
      m_dPhaseRelPerm_dPhaseVolFrac( relPermAccessors.get( fields::relperm::dPhaseRelPerm_dPhaseVolFraction {} ) ),     
      m_phaseCapPressure( capPressureAccessors.get( fields::cappres::phaseCapPressure {} ) ),
      m_dPhaseCapPressure_dPhaseVolFrac( capPressureAccessors.get( fields::cappres::dPhaseCapPressure_dPhaseVolFraction {} ) ),
      m_volume( subRegion.getElementVolume() ),
      m_porosity( solid.getPorosity() ), // poroAccessors.get( fields::porosity::porosity {} )
      m_dPoro_dPres( solid.getDporosity_dPressure() ), // poroAccessors.get( fields::porosity::dPorosity_dPressure {} )
      m_localSolution( flowAccessors.get( fields::immiscibleMultiphaseFlow::solutionUpdate {} ) ), // localSolution
      m_localResidual( localResidual ),
      m_resNorm( resNorm ),
      m_globalKinkFactor( globalKinkFactor ),  
      m_hasCapPressure( hasCapPressure ),
      m_maxIter( trustRegionParams.maxIter ),
      m_d2RMin( trustRegionParams.d2RMin ),
      m_minFactor( trustRegionParams.minInfFactor ),
      m_relResThres( trustRegionParams.relResThres ),
      m_absResThres( trustRegionParams.absResThres ),
      m_stencilWrapper( stencilWrapper ),
      m_fluidWrapper( fluidWrapper ),
      m_relPermWrapper( relPermWrapper ),   
      m_capPressureWrapper( capPressureWrapper ),      
      m_seri( stencilWrapper.getElementRegionIndices() ),
      m_sesri( stencilWrapper.getElementSubRegionIndices() ),
      m_sei( stencilWrapper.getElementIndices() )
  { GEOS_UNUSED_VAR( poroAccessors );}  

  struct StackVariables
  {
  public:
  
    /// Storage for the phase local damping factors    
    real64 localEps[numEqn]{ 1.0, 1.0 };

    /// Storage for the minimum cell local damping factor
    real64 elemFactor = 1.0;

    /// Storage for the face potential differences
    real64 dPhi[maxNumConn][numEqn]{};
  };  

  /**
   * @brief Compute the local inflection damping factors for the connection
   * @param[in] ei the cell index
   * @param[inout] stack the stack variables
   */  
  GEOS_HOST_DEVICE
  void computeInflection( localIndex const ei,
                          StackVariables & stack ) const
  {    
    // get list of with the connection number associated to the faces of the current cell
    arraySlice1d< localIndex const > const connList = m_connMap[ei];
    integer numConns = m_connMap.sizeOfArray( ei );

    localIndex seri[maxNumConn][2]{};
    localIndex sesri[maxNumConn][2]{};
    localIndex sei[maxNumConn][2]{};

    real64 transmissibility[maxNumConn]{};    

    // build cell maps (assuming only 2 cells per connection) and get geometric transmissibilities
    for ( integer conn = 0; conn < numConns; ++conn )
    {
      localIndex const iconn = connList[conn];

      // cell indices
      seri[conn][0] = m_seri( iconn, 0 );   seri[conn][1] = m_seri( iconn, 1 );
      sesri[conn][0] = m_sesri( iconn, 0 ); sesri[conn][1] = m_sesri( iconn, 1 );
      sei[conn][0] = m_sei( iconn, 0 );     sei[conn][1] = m_sei( iconn, 1 );      

      // transmissibilities      
      real64 trans[STENCILWRAPPER::maxNumConnections][2]{};
      real64 dTrans_dPres[STENCILWRAPPER::maxNumConnections][2]{};

      m_stencilWrapper.computeWeights( iconn,
                                       m_permeability,
                                       m_dPerm_dPres,                                    
                                       trans,
                                       dTrans_dPres );
      transmissibility[conn] = trans[0][0];
    }
    
    constexpr int signPotDiff[2] = {1, -1};
    real64 eps[2]{1.0, 1.0};
    
    // analyze residual inflections for each phase
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {        
      // Adaptive residual analysis
      globalIndex const globalRow = m_dofNumberElem[ei];
      localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
      GEOS_ASSERT_GE( localRow, 0 );  
      if( fabs( m_localResidual[localRow + ip] ) < m_relResThres * m_resNorm ||
          fabs( m_localResidual[localRow + ip] ) < m_absResThres )        
      {           
        continue; // skip analysis if phase residual is below minimum threshold
      }      
      
      real64 dPhi[maxNumConn]{};      
      real64 compressibility[maxNumConn][2]{};
      real64 viscosibility[maxNumConn][2]{};              

      real64 phaseEps[2]{0.0, 1.0}; 
      real64 d2F[maxNumConn][2]{};
      real64 d2R[2]{};      
      real64 newD2F[maxNumConn]{};            

      // loop through connections to compute directional derivative at eps = 0 and eps = 1
      for ( integer conn = 0; conn < numConns; ++conn )
      {           
        localIndex const iconn = connList[conn];   

        // compute average density, compressibility and viscosibility
        real64 densMean{};

        for( integer ke = 0; ke < 2; ++ke )
        {
          real64 const density  = m_dens[seri[conn][ke]][sesri[conn][ke]][sei[conn][ke]][0][ip];                              // r = rho1 || rho2
          real64 const dDens_dP = m_dDens_dPres[seri[conn][ke]][sesri[conn][ke]][sei[conn][ke]][0][ip][Deriv::dP];     // dr/dP = dr1/dP1 || dr2/dP

          real64 const viscosity = m_phaseVisc[seri[conn][ke]][sesri[conn][ke]][sei[conn][ke]][0][ip];                        // mu = mu1 || mu2
          real64 const dVisc_dP  = m_dPhaseVisc[seri[conn][ke]][sesri[conn][ke]][sei[conn][ke]][0][ip][Deriv::dP];   // dmu/dP = dmu1/dP1 || dmu2/dP

          densMean += 0.5 * density;                                                                 // rho = (rho1 + rho2) / 2  
          compressibility[conn][ke] = dDens_dP / density;                                            // cf = drho1 / rho1 || drho2 / rho2     
          viscosibility[conn][ke] = dVisc_dP / viscosity;                                            // cmu = dmu1 / mu1  || dmu2 / mu2
        }          

        // compute potential difference before update  
        for( integer ke = 0; ke < 2; ++ke )
        {
          real64 const pressure = m_pres[seri[conn][ke]][sesri[conn][ke]][sei[conn][ke]];       // P = P1 || P2
          real64 const gravD = m_gravCoef[seri[conn][ke]][sesri[conn][ke]][sei[conn][ke]];      // D = g z1 || g z2 

          real64 pot = pressure - densMean * gravD;      // Phi = P1 - rho g z1 || P2 - rho g z2

          if( m_hasCapPressure )
          {              
            pot -= m_phaseCapPressure[seri[conn][ke]][sesri[conn][ke]][sei[conn][ke]][0][ip];   // Phi = P1 - rho g z1 - Pc1 || P2 - rho g z2 - Pc2
          }

          dPhi[conn] += signPotDiff[ke] * pot;           // dPhi = P1 - P2 - rho g (z1 - z2) - Pc1 + Pc2          
        }

        // compute directional derivatives before and after update
        computeDirectionalDerivative(iconn, ip, transmissibility[conn], compressibility[conn], viscosibility[conn], dPhi[conn], phaseEps[0], d2F[conn][0]);
        computeDirectionalDerivative(iconn, ip, transmissibility[conn], compressibility[conn], viscosibility[conn], dPhi[conn], phaseEps[1], d2F[conn][1]);

        d2R[0] += d2F[conn][0];
        d2R[1] += d2F[conn][1];

      } // connection loop for initial derivative computations

      // check for presence of inflection
      for ( integer iter = 0; iter < m_maxIter &&
                              signbit( d2R[0] ) != signbit( d2R[1] ) &&
                              fabs( d2R[0] ) > m_d2RMin &&
                              fabs( d2R[1] ) > m_d2RMin; ++iter )
      {
        real64 newD2R{};
        real64 slope = (d2R[1] - d2R[0]) / (phaseEps[1] - phaseEps[0]);
        eps[ip] = phaseEps[1] - d2R[1] / slope;        

        // loop through connections to compute directional derivative at potential root
        for ( integer conn = 0; conn < numConns; ++conn )
        {        
          localIndex const iconn = connList[conn];
          computeDirectionalDerivative(iconn, ip, transmissibility[conn], compressibility[conn], viscosibility[conn], dPhi[conn], eps[ip], newD2F[conn]);

          newD2R += newD2F[conn];
        } // connection loop

        if ( signbit( d2R[0] ) == signbit( newD2R ) )
        {
          d2R[0] = newD2R;
          phaseEps[0] = eps[ip];
        } 
        else
        {
          d2R[1] = newD2R;
          phaseEps[1] = eps[ip];
        }

      } // iterative root search

      stack.localEps[ip] = eps[ip];
      
    } // loop through phases
    
  }

  /**
   * @brief Compute the directional second derivative of the phase flux function
   *        at a specific location in the solution space, given by the location 
   *        at the previous iteration plus a restricted update
   * @param[in] iconn the connection index
   * @param[in] k the cell index in the connection
   * @param[in] ip the phase index
   * @param[in] trans the connection transmissibility
   * @param[in] comp cell compressibilities
   * @param[in] visc cell viscosibility
   * @param[in] dPhiv phase potential difference before update  
   * @param[in] eps restriction factor that determines where in the update direction to evaluate
   *                the directional second derivative 
   * @param[out] d2F directional second derivative at given location in the solution space
   */                                 
  GEOS_HOST_DEVICE
  void computeDirectionalDerivative( localIndex const iconn,
                                     localIndex const ip,
                                     real64 const trans,
                                     real64 const * comp,
                                     real64 const * visc,
                                     real64 const dPhiv,                                     
                                     real64 const eps,
                                     real64 & d2F ) const
  {
    localIndex up = dPhiv >= 0.0 ? 0 : 1;
    localIndex down = 1 - up;
    
    // clear working variables
    real64 gravD[2]{};
    real64 dens[2]{};
    real64 viscosity[2]{};     
    real64 capPres[2]{};  
    real64 dcapPres[2]{};
    real64 d2capPres[2]{};  

    real64 dPhi{};
    real64 D{};
    real64 perm{}; 
    real64 dperm{}; 
    real64 d2perm{};

    real64 pressure[2]{};
    real64 phaseVolFrac[2]{};

    real64 dP[2]{};
    real64 dS[2]{};    

    // cell indices
    localIndex const seri[2]  = {m_seri( iconn, 0 ), m_seri( iconn, 1 )};
    localIndex const sesri[2] = {m_sesri( iconn, 0 ), m_sesri( iconn, 1 )};
    localIndex const sei[2]   = {m_sei( iconn, 0 ), m_sei( iconn, 1 )};

    // get pressure and saturation updates for connection elements
    // for ( integer ke = 0; ke < 2; ++ke )
    // {
    //   globalIndex const globalRow = m_dofNumber[seri[ke]][sesri[ke]][sei[ke]];
    //   localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
    //   GEOS_ASSERT_GE( localRow, 0 );      

    //   dP[ke] = m_localSolution[localRow] * m_globalKinkFactor;     // dP = { dP1 , dP2 }
    //   dS[ke] = m_localSolution[localRow + 1] * m_globalKinkFactor; // dS = { dS1 , dS2 }
    // }
    for ( integer ke = 0; ke < 2; ++ke )
    {
      dP[ke] = m_localSolution[seri[ke]][sesri[ke]][sei[ke]][0];     // dP = { dP1 , dP2 }
      dS[ke] = m_localSolution[seri[ke]][sesri[ke]][sei[ke]][1];     // dS = { dS1 , dS2 }
    }

    // compute saturation of updated state
    for ( integer ke = 0; ke < 2; ++ke )
    {  
      phaseVolFrac[ke] = ip == 0 ? m_phaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][ip] + eps * dS[ke] :  // S = S +- e.dS
                                   m_phaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][ip] - eps * dS[ke] ;
      
      // check physical bounds for saturation
      if ( phaseVolFrac[ke] < 0.0 )
      {
        phaseVolFrac[ke] = 0.0;
        dS[ke] = -m_phaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][ip]; // dS' = S' - S0
      }
      else if ( phaseVolFrac[ke] > 1.0 )
      {
        phaseVolFrac[ke] = 1.0;
        dS[ke] = 1.0 - m_phaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][ip]; // dS' = S' - S0
      }
    }

    // get grav coef, rel perm and capillary pressure + derivatives      
    for ( integer ke = 0; ke < 2; ++ke )
    {      
      gravD[ke] = m_gravCoef[seri[ke]][sesri[ke]][sei[ke]];              // D = { g z1 , g z2 }      

      if ( m_hasCapPressure )
      {
        (*m_capPressureWrapper).compute(phaseVolFrac[ke], ip, capPres[ke], dcapPres[ke], d2capPres[ke]);    // Pc   = { Pc1   || Pc2   }
                                                                                                            // Pc'  = { Pc1'  || Pc2'  }
                                                                                                            // Pc'' = { Pc1'' || Pc2'' }
      }      
    }
    m_relPermWrapper.compute(phaseVolFrac[up], ip, perm, dperm, d2perm);                                    // kr   = { krup   }
                                                                                                            // kr'  = { krup'  }
                                                                                                            // kr'' = { krup'' }
    D = gravD[up] - gravD[down];                                         // D = g (zup - zdown)

    // get density, viscosity and potential difference
    if ( eps < 0.001 )
    {
      // get density and viscosity   
      for ( integer ke = 0; ke < 2; ++ke )
      {
        dens[ke] = m_dens[seri[ke]][sesri[ke]][sei[ke]][0][ip];                 // r  = { r1  , r2  }
        viscosity[ke] = m_phaseVisc[seri[ke]][sesri[ke]][sei[ke]][0][ip];       // mu = { mu1 , mu2 }            
      }

      // get potential differece
      dPhi = fabs( dPhiv );
    }
    else // (0 < eps < 1)
    { 
      // compute pressure and density of updated state
      for ( integer ke = 0; ke < 2; ++ke )
      {
        pressure[ke] = m_pres[seri[ke]][sesri[ke]][sei[ke]] + eps * dP[ke];  // P = P1 || P2

        // check physical bounds for pressure
        if ( pressure[ke] < 0.0 )
        {
          pressure[ke] = 0.0;
          dP[ke] = -m_pres[seri[ke]][sesri[ke]][sei[ke]]; // dP' = P' - P0
        }        

        real64 dDens;   
        real64 dVisc;    
        m_fluidWrapper.compute(pressure[ke], ip, dens[ke], dDens, viscosity[ke], dVisc);  // r  = { r1  , r2  }
                                                                                          // mu = { mu1 , mu2 }             
      }    

      // compute potential difference
      constexpr int signPotDiff[2] = {1, -1}; 

      real64 densMean = (dens[0] + dens[1]) / 2;          // rho = (r1 + r2) / 2
      
      for( integer ke = 0; ke < 2; ++ke )
      {    
        real64 pot = pressure[ke] - densMean * gravD[ke]; // Phi = P1 - rho g z1 || P2 - rho g z2

        if( m_hasCapPressure )
        {              
          pot -= capPres[ke];                             // Phi = P1 - rho g z1 - Pc1 || P2 - rho g z2 - Pc2
        }

        dPhi += signPotDiff[ke] * pot;                    // Phi = P1 - P2 - rho g (z1 - z2) - Pc1 + Pc2
      }
      dPhi = fabs( dPhi );    
    }

    // compute partial second derivatives
    real64 d2Pi2 = perm * pow( comp[up] - visc[up], 2.0) * dPhi + 2 * perm * (comp[up] - visc[up]) + perm * comp[up] * (visc[up] - 1.5 * comp[up]) * dens[up] * D;
    real64 d2Pj2 = -0.5 * perm * pow( comp[down], 2.0 ) * dens[down] * D;
    real64 d2Si2 = d2perm * dPhi;
    real64 d2Sj2 = 0;
    real64 d2PiPj = perm * (comp[up] - visc[up]) * (-1 - 0.5 * comp[down] * dens[down] * D);
    real64 d2PiSi = dperm * (dPhi * (comp[up] - visc[up]) + 1 - 0.5 * comp[up] * dens[up] * D);
    real64 d2PiSj = 0;
    real64 d2PjSi = dperm * (-1 - 0.5 * comp[down] * dens[down] * D);
    real64 d2SiSj = 0;

    if ( m_hasCapPressure )
    {
      d2Si2 -= (2 * dperm * dcapPres[up] + perm * d2capPres[up]);
      d2Sj2 += perm * d2capPres[down];
      d2PiSi -= perm * (comp[up] - visc[up]) * dcapPres[up];
      d2PiSj += perm * (comp[up] - visc[up]) * dcapPres[down];
      d2SiSj += dperm * dcapPres[down];
    }

    // compute directional second derivative
    d2F = d2Pi2 * pow( dP[up], 2.0 ) + d2Pj2 * pow( dP[down], 2.0 ) + d2Si2 * pow( dS[up], 2.0 ) + d2Sj2 * pow( dS[down], 2.0 )
        + d2PiPj * 2 * dP[up] * dP[down] + d2PiSi * 2 * dP[up] * dS[up] + d2PiSj * 2 * dP[up] * dS[down]
        + d2PjSi * 2 * dP[down] * dS[up] + d2SiSj * 2 * dS[up] * dS[down];
    d2F *= trans * dens[up] / viscosity[up];
  }

  /**
   * @brief Return true for negative numbers and false for non-negative ones
   * @param[in] x number to check sign
   */                                 
  GEOS_HOST_DEVICE
  bool signbit( real64 x ) const
  {
    return x < 0 ? 1 : 0;
  }

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { 
    return m_ghostRank( ei );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the cell index
   * @param[inout] stack the stack variables
   */                                 
  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( ei );
    for ( integer ip = 0; ip < m_numPhases; ++ip )
    {
      stack.elemFactor = LvArray::math::min( stack.elemFactor, stack.localEps[ip] );
    }   
    stack.elemFactor = LvArray::math::max( stack.elemFactor, m_minFactor );
  }  

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[in] kernelComponent the kernel component providing access to the compute function and stack variables
   * @param[inout] inflectionFactor the minimum inflection damping factor for the current stencil
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void 
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent,
          real64 & inflectionFactor )
  {
    GEOS_MARK_FUNCTION;

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minInflectionFactor( 1.0 );

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.elemGhostRank( ei ) >= 0 )
      {
        return;
      }
      
      typename KERNEL_TYPE::StackVariables stack;
      
      kernelComponent.computeInflection( ei, stack );
      kernelComponent.complete( ei, stack );

      minInflectionFactor.min( stack.elemFactor );
    } );

    inflectionFactor = minInflectionFactor.get();
  }


protected:    

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;  

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;
  arrayView1d< globalIndex const > const m_dofNumberElem;

  /// Views on ghost rank numbers and gravity coefficients
  arrayView1d< integer const > const m_ghostRank;

  /// Views on the connection map
  ArrayOfArraysView< localIndex const > const m_connMap;

  /// ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;
  
  /// View on pressure and phase volume fraction
  ElementViewConst< arrayView1d< real64 const > > const m_pres;
  ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_phaseVolFrac;  

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > m_permeability;
  ElementViewConst< arrayView3d< real64 const > > m_dPerm_dPres;

  /// Views on fluid density
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_dens;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dDens_dPres;

  /// Views on the phase viscosities
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_phaseVisc;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dPhaseVisc;

  /// Views on the phase relative permeabilities
  ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const m_phaseRelPerm;
  ElementViewConst< arrayView4d< real64 const, relperm::USD_RELPERM_DS > > const m_dPhaseRelPerm_dPhaseVolFrac;
  ElementViewConst< arrayView4d< real64 const, relperm::USD_RELPERM_DS > > const m_d2PhaseRelPerm_d2PhaseVolFrac;  

  /// Views on capillary pressure
  ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const m_phaseCapPressure;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_d2PhaseCapPressure_d2PhaseVolFrac;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosity;
  arrayView2d< real64 const > const m_dPoro_dPres;

  /// View on the local solution and residual
  //arrayView1d< real64 const > const m_localSolution; 
  ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_localSolution;
  arrayView1d< real64 const > const m_localResidual;

  /// Residual norm for adaptive residual analysis
  real64 const m_resNorm;

  /// Kink factor computed previously
  real64 const m_globalKinkFactor;

  /// Flags
  integer const m_hasCapPressure;

  /// Trust region parameters
  /// Maximum number of nonlinear iterations
  integer const m_maxIter; // 5

  /// Minimum second derivative to apply a damping factor
  real64 const m_d2RMin; // 1.0e-14

  /// Minimum dampining factor
  real64 const m_minFactor; // 0.1

   /// Minimum relative residual threshold
  real64 const m_relResThres; // 0.6 

  /// Minimum absolute residual threshold
  real64 const m_absResThres; // 0.0  

  /// Reference to the wrappers
  STENCILWRAPPER const m_stencilWrapper;
  FLUIDWRAPPER const m_fluidWrapper;
  RELPERMWRAPPER const m_relPermWrapper;
  CAPPRESWRAPPER const * m_capPressureWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
};

/**
 * @class ResidualInflectionFactorKernelFactory
 */
class ResidualInflectionFactorKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localSolution the residual vector on my MPI rank
   * @param[in] globalKinkFactor kink restriction factor already computed
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] subRegion reference to the subregion
   * @param[in] solid the solid model
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] relPermWrapper reference to the relative permeability wrapper
   * @param[in] capPressureWrapper pointer to the capillary pressure wrapper
   * @param[in] hasCapPressure flag for presence of capillary pressure
   * @param[inout] inflectionFactor the inflection damping factor on the subRegion
   */
  template< typename POLICY, typename STENCILWRAPPER, typename FLUIDWRAPPER, typename RELPERMWRAPPER, typename CAPPRESWRAPPER >
  static void
  createAndLaunch( integer const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localSolution,
                   arrayView1d< real64 const > const & localResidual,
                   real64 const globalKinkFactor,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   CellElementSubRegion const & subRegion,
                   constitutive::CoupledSolidBase const & solid,
                   STENCILWRAPPER const & stencilWrapper,
                   FLUIDWRAPPER const & fluidWrapper,
                   RELPERMWRAPPER const & relPermWrapper,
                   CAPPRESWRAPPER const * capPressureWrapper,
                   integer const hasCapPressure,
                   real64 const resNorm,
                   ImmiscibleMultiphaseFlow::TrustRegionParameters const trustRegionParams,             
                   real64 & inflectionFactor )
  {
    if ( subRegion.getConnectionMap().size() == 0 )
    {
      inflectionFactor = 1.0; // no connections, no inflection damping
      return;
    }
    
    integer constexpr NUM_EQN = 2;
    integer constexpr NUM_DOF = 2;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = ResidualInflectionFactorKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER, FLUIDWRAPPER, RELPERMWRAPPER, CAPPRESWRAPPER >;
    typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
    typename kernelType::ImmiscibleMultiphaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::MultiphaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::RelPermAccessors relPermAccessors( elemManager, solverName );
    typename kernelType::PorosityAccessors poroAccessors( elemManager, solverName );
    typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );

    kernelType kernel( numPhases, rankOffset, dofKey, localSolution, localResidual, resNorm, globalKinkFactor, subRegion,
                       solid, stencilWrapper, fluidWrapper, relPermWrapper, capPressureWrapper,
                       dofNumberAccessor, permeabilityAccessors, flowAccessors, fluidAccessors, 
                       relPermAccessors, poroAccessors, capPressureAccessors, hasCapPressure, trustRegionParams );
    kernelType::template launch< POLICY >( subRegion.size(), kernel, inflectionFactor );    
  }
};

} // namespace immiscible multiphasekernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_MULTIPHASEKERNELS_HPP
