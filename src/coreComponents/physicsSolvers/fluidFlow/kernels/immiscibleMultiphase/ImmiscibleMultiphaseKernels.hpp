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
 * @file ImmiscibleMultiphaseKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_MULTIPHASEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_MULTIPHASEKERNELS_HPP

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
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/RelativePermeabilityUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/CapillaryPressureUpdateKernel.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/KernelLaunchSelectors.hpp"

namespace geos
{
namespace immiscibleMultiphaseKernels
{
using namespace constitutive;
using namespace immiscibleFlowUtilities;


/******************************** FluxComputeKernelBase ********************************/

/**
 * @brief Base class for FluxComputeKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of dofs).
 */
class FluxComputeKernelBase
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
                      fields::immiscibleMultiphaseFlow::phaseMobility,
                      fields::immiscibleMultiphaseFlow::dPhaseMobility >;

  using MultiphaseFluidAccessors =
    StencilMaterialAccessors< constitutive::TwoPhaseImmiscibleFluid,
                              fields::twophaseimmisciblefluid::phaseDensity,
                              fields::twophaseimmisciblefluid::dPhaseDensity >;

  using CapPressureAccessors =
    StencilMaterialAccessors< CapillaryPressureBase,
                              fields::cappres::phaseCapPressure,
                              fields::cappres::dPhaseCapPressure_dPhaseVolFraction >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  using Deriv = immiscibleFlow::DerivativeOffset;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor accessor for the dof numbers
   * @param[in] multiPhaseFlowAccessors accessor for wrappers registered by the solver
   * @param[in] fluidAccessors accessor for wrappers registered by the fluid model
   * @param[in] capPressureAccessors accessor for wrappers registered by the capillary pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[inout] hasCapPressure flag to indicate whether problem includes capillarity
   * @param[inout] useTotalMassEquation flag to indicate whether to use the total mass formulation
   */
  FluxComputeKernelBase( integer const numPhases,
                         globalIndex const rankOffset,
                         DofNumberAccessor const & dofNumberAccessor,
                         ImmiscibleMultiphaseFlowAccessors const & multiPhaseFlowAccessors,
                         MultiphaseFluidAccessors const & fluidAccessors,
                         CapPressureAccessors const & capPressureAccessors,
                         PermeabilityAccessors const & permeabilityAccessors,
                         real64 const & dt,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         integer const hasCapPressure,
                         integer const useTotalMassEquation,
                         integer const checkPhasePresenceInGravity )
    : m_numPhases ( numPhases ),
    m_rankOffset( rankOffset ),
    m_dt( dt ),
    m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_permeability( permeabilityAccessors.get( fields::permeability::permeability {} ) ),
    m_dPerm_dPres( permeabilityAccessors.get( fields::permeability::dPerm_dPressure {} ) ),
    m_ghostRank( multiPhaseFlowAccessors.get( fields::ghostRank {} ) ),
    m_gravCoef( multiPhaseFlowAccessors.get( fields::flow::gravityCoefficient {} ) ),
    m_pres( multiPhaseFlowAccessors.get( fields::flow::pressure {} ) ),
    m_phaseVolFrac( multiPhaseFlowAccessors.get( fields::immiscibleMultiphaseFlow::phaseVolumeFraction {} ) ),
    m_mob( multiPhaseFlowAccessors.get( fields::immiscibleMultiphaseFlow::phaseMobility {} ) ),
    m_dMob( multiPhaseFlowAccessors.get( fields::immiscibleMultiphaseFlow::dPhaseMobility {} ) ),
    m_dens( fluidAccessors.get( fields::twophaseimmisciblefluid::phaseDensity {} ) ),
    m_dDens_dPres( fluidAccessors.get( fields::twophaseimmisciblefluid::dPhaseDensity {} ) ),
    m_phaseCapPressure( capPressureAccessors.get( fields::cappres::phaseCapPressure {} ) ),
    m_dPhaseCapPressure_dPhaseVolFrac( capPressureAccessors.get( fields::cappres::dPhaseCapPressure_dPhaseVolFraction {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs ),
    m_hasCapPressure ( hasCapPressure ),
    m_useTotalMassEquation ( useTotalMassEquation ),
    m_checkPhasePresenceInGravity ( checkPhasePresenceInGravity )
  {}

protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Time step size
  real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > m_permeability;
  ElementViewConst< arrayView3d< real64 const > > m_dPerm_dPres;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  // Primary and secondary variables
  /// Views on pressure and phase volume fraction
  ElementViewConst< arrayView1d< real64 const > > const m_pres;
  ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_phaseVolFrac;

  /// Views on fluid mobility
  ElementViewConst< arrayView2d< real64 const, immiscibleFlow::USD_PHASE > > const m_mob;
  ElementViewConst< arrayView3d< real64 const, immiscibleFlow::USD_PHASE_DS > > const m_dMob;

  /// Views on fluid density
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_dens;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dDens_dPres;

  /// Views on capillary pressure
  ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const m_phaseCapPressure;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac;

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  // Flags
  integer const m_hasCapPressure;
  integer const m_useTotalMassEquation;
  integer const m_checkPhasePresenceInGravity;
};

/***************************************** */

/**
 * @class FluxComputeKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class FluxComputeKernel : public FluxComputeKernelBase
{
public:

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_EQN;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::maxNumPointsInFlux;

  /// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::maxNumConnections;

  /// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::maxStencilSize;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] multiPhaseFlowAccessors
   * @param[in] fluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] hasCapPressure flags for capillary pressure
   * @param[in] useTotalMassEquation flags for using total velocity formulation
   */
  FluxComputeKernel( integer const numPhases,
                     globalIndex const rankOffset,
                     STENCILWRAPPER const & stencilWrapper,
                     DofNumberAccessor const & dofNumberAccessor,
                     ImmiscibleMultiphaseFlowAccessors const & multiPhaseFlowAccessors,
                     MultiphaseFluidAccessors const & fluidAccessors,
                     CapPressureAccessors const & capPressureAccessors,
                     PermeabilityAccessors const & permeabilityAccessors,
                     real64 const & dt,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs,
                     integer const hasCapPressure,
                     integer const useTotalMassEquation,
                     integer const checkPhasePresenceInGravity )
    : FluxComputeKernelBase( numPhases,
                             rankOffset,
                             dofNumberAccessor,
                             multiPhaseFlowAccessors,
                             fluidAccessors,
                             capPressureAccessors,
                             permeabilityAccessors,
                             dt,
                             localMatrix,
                             localRhs,
                             hasCapPressure,
                             useTotalMassEquation,
                             checkPhasePresenceInGravity ),
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
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : stencilSize( size ),
      numFluxElems( numElems ),
      dofColIndices( size * numDof ),
      localFlux( numElems * numEqn ),
      localFluxJacobian( numElems * numEqn, size * numDof )
    {}

    // Stencil information

    /// Stencil size for a given connection
    localIndex const stencilSize;

    /// Number of elements for a given connection
    localIndex const numFluxElems;

    // Transmissibility and derivatives

    /// Transmissibility
    real64 transmissibility[maxNumConns][2]{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dPres[maxNumConns][2]{};

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    stackArray1d< globalIndex, maxNumElems * numDof > dofColIndices;

    /// Storage for the face local residual vector (all equations except volume balance)
    stackArray1d< real64, maxNumElems * numEqn > localFlux;
    /// Storage for the face local Jacobian matrix
    stackArray2d< real64, maxNumElems * numEqn * maxStencilSize * numDof > localFluxJacobian;
  };

  /**
   * @brief Getter for the stencil size at this connection
   * @param[in] iconn the connection index
   * @return the size of the stencil at this connection
   */
  GEOS_HOST_DEVICE
  localIndex stencilSize( localIndex const iconn ) const
  { return m_sei[iconn].size(); }

  /**
   * @brief Getter for the number of elements at this connection
   * @param[in] iconn the connection index
   * @return the number of elements at this connection
   */
  GEOS_HOST_DEVICE
  localIndex numPointsInFlux( localIndex const iconn ) const
  { return m_stencilWrapper.numPointsInFlux( iconn ); }

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // set degrees of freedom indices for this face
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      globalIndex const offset = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];

      for( integer jdof = 0; jdof < numDof; ++jdof )
      {
        stack.dofColIndices[i * numDof + jdof] = offset + jdof;
      }
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] NoOpFunc the function used to customize the computation of the flux
   */
  template< typename FUNC = NoOpFunc >   // should change to multiphase
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = NoOpFunc{} ) const
  {
    // first, compute the transmissibilities at this face                                             // get k and dk/dP from global arrays
    // and place in stack
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    localIndex k[2];
    localIndex connectionIndex = 0;

    for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
      {
        // clear working arrays
        real64 densMean[numEqn]{};
        real64 dDensMean_dP[numEqn][2]{};

        real64 presGrad[numEqn]{};
        real64 dPresGrad_dP[numEqn][2]{};

        real64 gravHead[numEqn]{};
        real64 dGravHead_dP[numEqn][2]{};

        real64 capGrad[numEqn]{};
        real64 dCapGrad_dP[numEqn][2]{};
        real64 dCapGrad_dS[numEqn][2]{};

        real64 fluxVal[numEqn]{};
        real64 dFlux_dP[numEqn][2]{};
        real64 dFlux_dS[numEqn][2]{};

        real64 mobility[numEqn]{};
        real64 dMob_dP[numEqn][2]{};
        real64 dMob_dS[numEqn][2]{};

        real64 const trans[2] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };
        real64 const dTrans_dP[2] = { stack.dTrans_dPres[connectionIndex][0], stack.dTrans_dPres[connectionIndex][1] };

        // cell indices
        localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // loop over phases
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          // calculate quantities on primary connected cells
          integer denom = 0;
          for( integer ke = 0; ke < 2; ++ke )
          {
            // density
            bool const phaseExists = (m_phaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][ip] > 0);
            if( m_checkPhasePresenceInGravity && !phaseExists )
            {
              continue;
            }

            real64 const density  = m_dens[seri[ke]][sesri[ke]][sei[ke]][0][ip];                    // r = rho1 || rho2
            real64 const dDens_dP = m_dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0][ip][Deriv::dP];  // dr/dP = dr1/dP1 || dr2/dP

            // average density and derivatives
            densMean[ip] += density;          // rho = (rho1 + rho2)
            dDensMean_dP[ip][ke] = dDens_dP;  // drho/dP = { (dr1/dP1) , (dr2/dP2) }

            denom++;
          }

          if( denom > 1 )
          {
            densMean[ip] /= denom; // rho = (rho1 + rho2) / denom
            for( integer ke = 0; ke < 2; ++ke )
            {
              dDensMean_dP[ip][ke] /= denom; // drho/dP = { (dr1/dP1) / denom , (dr2/dP2) / denom }
            }
          }

          //***** calculation of flux *****

          // compute potential difference
          real64 potScale = 0.0;
          real64 dPresGrad_dTrans = 0.0;
          real64 dGravHead_dTrans = 0.0;
          real64 dCapGrad_dTrans = 0.0;
          constexpr int signPotDiff[2] = {1, -1};

          for( integer ke = 0; ke < 2; ++ke )
          {
            localIndex const er  = seri[ke];
            localIndex const esr = sesri[ke];
            localIndex const ei  = sei[ke];

            real64 const pressure = m_pres[er][esr][ei];      // P = P1 || P2
            presGrad[ip] += trans[ke] * pressure;             // DPv = T (P1 - P2)
            dPresGrad_dTrans += signPotDiff[ke] * pressure;   // dDPv/dT = (P1 - P2)
            dPresGrad_dP[ip][ke] = trans[ke];                 // dDPv/dP = { T , -T }

            real64 const gravD = trans[ke] * m_gravCoef[er][esr][ei];       // D = T g z1 || -T g z2
            real64 pot = trans[ke] * pressure - densMean[ip] * gravD; // Phi = T P1 - rho T g z1 || -T P2 + rho T g z2

            gravHead[ip] += densMean[ip] * gravD;                                         // DPg = rho (T g z1 - T g z2) = T rho g (z1 - z2)
            dGravHead_dTrans += signPotDiff[ke] * densMean[ip] * m_gravCoef[er][esr][ei]; // dDPg/dT = rho g z1 - rho g z2 = rho g (z1 - z2)

            for( integer i = 0; i < 2; ++i )
            {
              dGravHead_dP[ip][i] += dDensMean_dP[ip][i] * gravD; // dDPg/dP = {drho/dP1 * T g (z1 - z2) , drho/dP2 * T g (z1 - z2)}
            }

            if( m_hasCapPressure )  // check sign convention
            {
              real64 const capPres = m_phaseCapPressure[er][esr][ei][0][ip];  // Pc = Pc1 || Pc2
              dCapGrad_dTrans -= signPotDiff[ke] * capPres;                   // dDPc/dT = (-Pc1 + Pc2)
              pot -= trans[ke] * capPres;                                     // Phi = T P1 - rho T g z1 - T Pc1 || -T P2 + rho T g z2 + T
                                                                              // Pc2
              capGrad[ip] -= trans[ke] * capPres;                             // DPc = T (-Pc1 + Pc2)
            }

            potScale = fmax( potScale, fabs( pot ) ); // maxPhi = Phi1 > Phi2 ? Phi1 : Phi2
          }

          for( integer ke = 0; ke < 2; ++ke )
          {
            dPresGrad_dP[ip][ke] += dTrans_dP[ke] * dPresGrad_dTrans;   // dDPv/dP = { T + dT/dP1 * (P1 - P2) , -T + dT/dP2 * (P1 - P2)}
            dGravHead_dP[ip][ke] += dTrans_dP[ke] * dGravHead_dTrans;   // dDPg/dP = { drho/dP1 * T g (z1 - z2) + dT/dP1 * rho g (z1 - z2) ,
                                                                        //             drho/dP2 * T g (z1 - z2) + dT/dP2 * rho g (z1 - z2) }
            if( m_hasCapPressure )
            {
              real64 const dCapPres_dS = m_dPhaseCapPressure_dPhaseVolFrac[seri[ke]][sesri[ke]][sei[ke]][0][ip][ip]; // dPc/dS = dPc1/dS1 ||
                                                                                                                     // dPc2/dS2
              dCapGrad_dP[ip][ke] += dTrans_dP[ke] * dCapGrad_dTrans;                                                // dDPc/dP = { dT/dP1 *
                                                                                                                     // (-Pc1 + Pc2) ,
                                                                                                                     //             dT/dP2 *
                                                                                                                     // (-Pc1 + Pc2) }
              dCapGrad_dS[ip][ke] -= trans[ke] * dCapPres_dS;                                                        // dDPc/dS = { -T *
                                                                                                                     // dPc1/dS1 , T *
                                                                                                                     // dPc2/dS2 }
            }
          }

          // *** upwinding ***

          // compute potential gradient
          real64 potGrad = presGrad[ip] - gravHead[ip]; // DPhi = T (P1 - P2) - T rho g (z1 - z2)
          if( m_hasCapPressure )
          {
            potGrad += capGrad[ip]; // DPhi = T (P1 - P2) - T rho g (z1 - z2) + T (-Pc1 + Pc2)
          }

          // compute upwinding tolerance
          real64 constexpr upwRelTol = 1e-8;
          real64 const upwAbsTol = fmax( potScale * upwRelTol, LvArray::NumericLimits< real64 >::epsilon ); // abstol = maxPhi * tol > eps ?
                                                                                                            // maxPhi * tol : eps

          // decide mobility coefficients - smooth variation in [-upwAbsTol; upwAbsTol]
          real64 const alpha = ( potGrad + upwAbsTol ) / ( 2 * upwAbsTol );   // alpha = (DPhi + abstol) / abstol / 2

          // choose upstream cell
          if( alpha <= 0.0 || alpha >= 1.0 )  // no smoothing needed
          {
            localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) ); // 1 upwind -> k_up = 0 || 2 upwind -> k_up = 1

            mobility[ip] = m_mob[seri[k_up]][sesri[k_up]][sei[k_up]][ip];                     // M = Mupstream
            dMob_dP[ip][k_up] = m_dMob[seri[k_up]][sesri[k_up]][sei[k_up]][ip][Deriv::dP];    // dM/dP = {dM/dP1 , 0} OR {0 , dM/dP2}
            dMob_dS[ip][k_up] = m_dMob[seri[k_up]][sesri[k_up]][sei[k_up]][ip][Deriv::dS];    // dM/dS = {dM/dS1 , 0} OR {0 , dM/dS2}
          }
          else  // perform smoothing
          {
            real64 const mobWeights[2] = { alpha, 1.0 - alpha };
            for( integer ke = 0; ke < 2; ++ke )
            {
              mobility[ip] += mobWeights[ke] * m_mob[seri[ke]][sesri[ke]][sei[ke]][ip];               // M = alpha * M1 + (1 - alpha) * M2
              dMob_dP[ip][ke] = mobWeights[ke] * m_dMob[seri[ke]][sesri[ke]][sei[ke]][ip][Deriv::dP]; // dM/dP = {alpha * dM1/dP1 , (1 -
                                                                                                      // alpha) * dM2/dP2}
              dMob_dS[ip][ke] = mobWeights[ke] * m_dMob[seri[ke]][sesri[ke]][sei[ke]][ip][Deriv::dS]; // dM/dP = {alpha * dM1/dS1 , (1 -
                                                                                                      // alpha) * dM2/dS2}
            }
          }

          // pressure gradient depends on all points in the stencil
          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] += dPresGrad_dP[ip][ke]; // dF/dP = { T + dT/dP1 * (P1 - P2) ,
                                                      //          -T + dT/dP2 * (P1 - P2) }
          }

          // gravitational head depends only on the two cells connected (same as mean density)
          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] -= dGravHead_dP[ip][ke]; // dF/dP = { T + dT/dP1 * (P1 - P2) - drho/dP1 * T g (z1 - z2) - dT/dP1 * rho g (z1 -
                                                      // z2) ,
                                                      //          -T + dT/dP2 * (P1 - P2) - drho/dP2 * T g (z1 - z2) - dT/dP2 * rho g (z1 -
                                                      // z2) }
          }

          // capillary pressure contribution
          if( m_hasCapPressure )
          {
            for( integer ke = 0; ke < 2; ++ke )
            {
              dFlux_dP[ip][ke] += dCapGrad_dP[ip][ke];  // dF/dP = { T + dT/dP1 * (P1 - P2) - drho/dP1 * T g (z1 - z2) - dT/dP1 * rho g (z1
                                                        // - z2) + dT/dP1 * (-Pc1 + Pc2) ,
                                                        //          -T + dT/dP2 * (P1 - P2) - drho/dP2 * T g (z1 - z2) - dT/dP2 * rho g (z1
                                                        // - z2) + dT/dP2 * (-Pc1 + Pc2) }

              dFlux_dS[ip][ke] += dCapGrad_dS[ip][ke];  // dF/dS = { T * -dPc/dS1 , T * dPc/dS2 }
            }
          }

          // compute the flux and derivatives using upstream cell mobility
          fluxVal[ip] = mobility[ip] * potGrad; // F = M * DPhi

          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] *= mobility[ip];   // dF/dP = { M [ T + dT/dP1 * (P1 - P2) - drho/dP1 * T g (z1 - z2) - dT/dP1 * rho g (z1 -
                                                // z2) + dT/dP1 * (-Pc1 + Pc2)] ,
                                                //           M [-T + dT/dP2 * (P1 - P2) - drho/dP2 * T g (z1 - z2) - dT/dP2 * rho g (z1 -
                                                // z2) + dT/dP2 * (-Pc1 + Pc2)] }

            dFlux_dS[ip][ke] *= mobility[ip];   // dF/dS = { M [T * -dPc/dS1] , M [T * dPc/dS2] }
          }

          // add contribution from upstream cell mobility derivatives
          for( integer ke = 0; ke < 2; ++ke )
          {
            dFlux_dP[ip][ke] += dMob_dP[ip][ke] * potGrad;  // dF/dP = { M [ T + dT/dP1 * (P1 - P2) - drho1/dP * T g (z1 - z2) - dT/dP1 *
                                                            // rho g (z1 - z2) + dT/dP1 * (-Pc1 + Pc2)] + dM/dP1 * DPhi ,
                                                            //           M [-T + dT/dP2 * (P1 - P2) - drho2/dP * T g (z1 - z2) - dT/dP2 *
                                                            // rho g (z1 - z2) + dT/dP2 * (-Pc1 + Pc2)] + dM/dP2 * DPhi }

            dFlux_dS[ip][ke] += dMob_dS[ip][ke] * potGrad;  // dF/dS = { M [T * -dPc/dS1] + dM/dS1 * DPhi , M [T * dPc/dS2] + dM/dS2 * DPhi
                                                            // }
          }

          // populate local flux vector and derivatives
          stack.localFlux[k[0]*numEqn + ip] += m_dt * fluxVal[ip];
          stack.localFlux[k[1]*numEqn + ip] -= m_dt * fluxVal[ip];

          for( integer ke = 0; ke < 2; ++ke )
          {
            // pressure
            localIndex const localDofIndexPres = k[ke] * numDof;
            stack.localFluxJacobian[k[0]*numEqn + ip][localDofIndexPres] += m_dt * dFlux_dP[ip][ke];
            stack.localFluxJacobian[k[1]*numEqn + ip][localDofIndexPres] -= m_dt * dFlux_dP[ip][ke];

            // saturation (hard-coded for 2-phase currently)
            localIndex const localDofIndexSat = k[ke] * numDof + 1;
            stack.localFluxJacobian[k[0]*numEqn + ip][localDofIndexSat] += m_dt * dFlux_dS[ip][ke];
            stack.localFluxJacobian[k[1]*numEqn + ip][localDofIndexSat] -= m_dt * dFlux_dS[ip][ke];
          }

          // Customize the kernel with this lambda
          kernelOp( k, seri, sesri, sei, connectionIndex, alpha, mobility, potGrad, fluxVal, dFlux_dP, dFlux_dS );  // Not sure what this
                                                                                                                    // does

        } // loop over phases

        connectionIndex++;
      }
    } // loop over connection elements
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = NoOpFunc >                                  // should change to multiphase
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && kernelOp = NoOpFunc{} ) const
  {
    using namespace compositionalMultiphaseUtilities;

    if( m_useTotalMassEquation )
    {
      // Apply equation/variable change transformation(s)
      stackArray1d< real64, maxStencilSize * numDof > work( stack.stencilSize * numDof );
      shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numPhases, numEqn, numDof * stack.stencilSize, stack.numFluxElems,
                                                               stack.localFluxJacobian, work );
      shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( m_numPhases, numEqn, stack.numFluxElems,
                                                                 stack.localFlux );
    }

    // add contribution to residual and jacobian into:
    // - the mass balance equation
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    for( integer i = 0; i < stack.numFluxElems; ++i )
    {
      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOS_ASSERT_GE( localRow, 0 );

        GEOS_ASSERT_GT( m_localMatrix.numRows(), localRow + numEqn - 1 );

        for( integer ic = 0; ic < numEqn; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + ic],
                           stack.localFlux[i * numEqn + ic] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
            ( localRow + ic,
            stack.dofColIndices.data(),
            stack.localFluxJacobian[i * numEqn + ic].dataIfContiguous(),
            stack.stencilSize * numDof );
        }

        // call the lambda to assemble additional terms, such as thermal terms
        kernelOp( i, localRow );
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
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numConnections, [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent.stencilSize( iconn ),
                                                  kernelComponent.numPointsInFlux( iconn ) );

      kernelComponent.setup( iconn, stack );
      kernelComponent.computeFlux( iconn, stack );
      kernelComponent.complete( iconn, stack );
    } );
  }

protected:

  // Stencil information

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
};

/****************************************** */

/**
 * @class FluxComputeKernelFactory
 */
class FluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const hasCapPressure,
                   integer const useTotalMassEquation,
                   integer const checkPhasePresenceInGravity,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_EQN = 2;
    integer constexpr NUM_DOF = 2;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = FluxComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
    typename kernelType::ImmiscibleMultiphaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::MultiphaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );

    kernelType kernel( numPhases, rankOffset, stencilWrapper, dofNumberAccessor,
                       flowAccessors, fluidAccessors, capPressureAccessors, permAccessors,
                       dt, localMatrix, localRhs, hasCapPressure, useTotalMassEquation,
                       checkPhasePresenceInGravity );
    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};


/******************************** AccumulationKernel ********************************/

enum class KernelFlags
{
  TotalMassEquation = 1 << 0, // 1

  /// Add more flags like that if needed:
  // Flag2 = 1 << 1, // 2
  // Flag3 = 1 << 2, // 4
  // Flag4 = 1 << 3, // 8
  // Flag5 = 1 << 4, // 16
  // Flag6 = 1 << 5, // 32
  // Flag7 = 1 << 6, // 64
  // Flag8 = 1 << 7  //128
};
/**
 * @class AccumulationKernel
 * @tparam NUM_EQN number of fluid phases
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */

template< integer NUM_EQN, integer NUM_DOF >
class AccumulationKernel
{
public:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_EQN;

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
   * @param[inout] kernelFlags the kernel options
   */
  AccumulationKernel( localIndex const numPhases,
                      globalIndex const rankOffset,
                      string const dofKey,
                      ElementSubRegionBase const & subRegion,
                      constitutive::TwoPhaseImmiscibleFluid const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs,
                      BitFlags< KernelFlags > const kernelFlags )
    : m_numPhases( numPhases ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    m_phaseVolFrac( subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseMass_n( subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs ),
    m_kernelFlags( kernelFlags )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    // Pore volume information (used by both accumulation)

    /// Pore volume at time n+1
    real64 poreVolume = 0.0;

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
    // initialize the pore volume
    stack.poreVolume = m_volume[ei] * m_porosity[ei][0];
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
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    constexpr int sign[2] = {1, -1};
    // ip - index of phase/component whose conservation equation is assembled
    // (i.e. row number in local matrix)
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      real64 const phaseMass = stack.poreVolume * m_phaseVolFrac[ei][ip] * m_phaseDens[ei][0][ip];
      real64 const phaseMass_n = m_phaseMass_n[ei][ip];

      stack.localResidual[ip] += phaseMass - phaseMass_n;

      real64 const dPhaseMass_dP = stack.dPoreVolume_dPres * m_phaseVolFrac[ei][ip] * m_phaseDens[ei][0][ip]
                                   + stack.poreVolume * m_phaseVolFrac[ei][ip] * m_dPhaseDens[ei][0][ip][0];
      stack.localJacobian[ip][0] += dPhaseMass_dP;

      real64 const dPhaseMass_dS = stack.poreVolume * m_phaseDens[ei][0][ip];
      stack.localJacobian[ip][1] += sign[ip] * dPhaseMass_dS;
    }
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
    using namespace compositionalMultiphaseUtilities;

    if( m_kernelFlags.isSet( KernelFlags::TotalMassEquation ) )
    {
      // apply equation/variable change transformation to the component mass balance equations
      real64 work[numDof]{};
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numPhases, numDof, stack.localJacobian, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( m_numPhases, stack.localResidual );
    }

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numPhase - 1)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    integer const numRows = m_numPhases;
    for( integer i = 0; i < numRows; ++i )
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
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
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

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosity;
  arrayView2d< real64 const > const m_dPoro_dPres;

  /// Views on the phase volume fractions
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const m_phaseVolFrac;

  /// Views on the phase densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const m_dPhaseDens;

  // View on component amount (mass) from previous time step
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_phaseMass_n;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  BitFlags< KernelFlags > const m_kernelFlags;
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
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[inout] useTotalMassEquation option for using total mass equation
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numPhases,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   constitutive::TwoPhaseImmiscibleFluid const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {

    geos::immiscibleMultiphaseKernels::kernelLaunchSelectorPhaseSwitch( numPhases, [&] ( auto NP )
    {
      integer constexpr NUM_EQN = NP();
      integer constexpr NUM_DOF = NP();

      BitFlags< KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( KernelFlags::TotalMassEquation );

      AccumulationKernel< NUM_EQN, NUM_DOF > kernel( numPhases, rankOffset, dofKey, subRegion,
                                                     fluid, solid, localMatrix, localRhs, kernelFlags );
      AccumulationKernel< NUM_EQN, NUM_DOF >::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }

};



/******************************** PhaseMobilityKernel ********************************/

/**
 * @class PhaseMobilityKernel
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase mobilities
 */
template< integer NUM_PHASE >
class PhaseMobilityKernel
{
public:

  //using Base = MultiphaseFluidAccessors::PropertyKernelBase< NUM_COMP >;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  PhaseMobilityKernel( ObjectManagerBase & subRegion,
                       TwoPhaseImmiscibleFluid const & fluid,
                       RelativePermeabilityBase const & relperm )
    :
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseVisc( fluid.phaseViscosity() ),
    m_dPhaseVisc( fluid.dPhaseViscosity() ),
    m_phaseRelPerm( relperm.phaseRelPerm() ),
    m_dPhaseRelPerm_dPhaseVolFrac( relperm.dPhaseRelPerm_dPhaseVolFraction() ),
    m_phaseMob( subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMobility >() ),
    m_dPhaseMob( subRegion.getField< fields::immiscibleMultiphaseFlow::dPhaseMobility >() )
  {}

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
    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      kernelComponent.compute( ei );
    } );
  }

  /**
   * @brief Compute the phase mobilities in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseMobilityKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && phaseMobilityKernelOp = NoOpFunc{} ) const
  {
    using Deriv = immiscibleFlow::DerivativeOffset;

    arraySlice1d< real64, immiscibleFlow::USD_PHASE - 1 > const phaseMob = m_phaseMob[ei];
    arraySlice2d< real64, immiscibleFlow::USD_PHASE_DS - 1 > const dPhaseMob = m_dPhaseMob[ei];
    constexpr int sign[2] = {1, -1};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      real64 const density = m_phaseDens[ei][0][ip];
      real64 const dDens_dP = m_dPhaseDens[ei][0][ip][Deriv::dP];
      real64 const viscosity = m_phaseVisc[ei][0][ip];
      real64 const dVisc_dP = m_dPhaseVisc[ei][0][ip][Deriv::dP];

      real64 const relPerm = m_phaseRelPerm[ei][0][ip];

      real64 const mobility = relPerm * density / viscosity;

      phaseMob[ip] = mobility;
      dPhaseMob[ip][Deriv::dP] = mobility * (dDens_dP / density - dVisc_dP / viscosity);

      // for( integer jp = 0; jp < numPhase-1; ++jp )
      // {
      // Derivative matrix is currently diagonal. Implementation below handles missing off-diagonal entry.
      real64 const dRelPerm_dS = sign[ip] * m_dPhaseRelPerm_dPhaseVolFrac[ei][0][ip][ip];
      dPhaseMob[ip][Deriv::dS] = dRelPerm_dS * density / viscosity;
      // }


      // call the lambda in the phase loop to allow the reuse of the relperm, density, viscosity, and mobility
      // possible use: assemble the derivatives wrt temperature
      phaseMobilityKernelOp( ip, phaseMob[ip], dPhaseMob[ip] );
    }
  }

protected:

  // inputs

  /// Views on the phase densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseDens;

  //arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_phaseDens;
  //arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_dPhaseDens;

  /// Views on the phase viscosities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseVisc;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseVisc;
  //arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_phaseVisc;
  //arrayView2d< real64 const, immiscibleFlow::USD_PHASE > m_dPhaseVisc;

  /// Views on the phase relative permeabilities
  arrayView3d< real64 const, relperm::USD_RELPERM > m_phaseRelPerm;
  arrayView4d< real64 const, relperm::USD_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;

  // outputs

  /// Views on the phase mobilities
  arrayView2d< real64, immiscibleFlow::USD_PHASE > m_phaseMob;
  arrayView3d< real64, immiscibleFlow::USD_PHASE_DS > m_dPhaseMob;
};

/**
 * @class PhaseMobilityKernelFactory
 */
class PhaseMobilityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] relperm the relperm model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numPhase,
                   ObjectManagerBase & subRegion,
                   TwoPhaseImmiscibleFluid const & fluid,
                   RelativePermeabilityBase const & relperm )
  {
    if( numPhase == 2 )
    {
      PhaseMobilityKernel< 2 > kernel( subRegion, fluid, relperm );
      PhaseMobilityKernel< 2 >::template launch< POLICY >( subRegion.size(), kernel );
    }
  }
};


struct FluidUpdateKernel
{
  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k] );
      }
    } );
  }
};

//******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 *
 * @brief Define the interface for the property kernel in charge of computing the residual norm
 */

/// Compile time value for the number of norms to compute
static constexpr integer numNorm = 1;

class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< numNorm >
{
public:


  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< numNorm >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      integer const numPhases,
                      ElementSubRegionBase const & subRegion,
                      constitutive::CoupledSolidBase const & solid,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numPhases( numPhases ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_phaseMass_n( subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 massNormalizer = 0;
    for( integer idof = 0; idof < m_numPhases; ++idof )
    {
      massNormalizer += LvArray::math::max( m_minNormalizer, m_phaseMass_n[ei][idof] );
    }

    for( integer idof = 0; idof < m_numPhases; ++idof )
    {
      real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / massNormalizer;
      if( valMass > stack.localValue[0] )
      {
        stack.localValue[0] = valMass;
      }
    }

    // for( integer idof = 0; idof < m_numPhases; ++idof )
    // {
    //   real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_phaseMass_n[ei][idof] );
    //   real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / massNormalizer;
    //   if( valMass > stack.localValue[0] )
    //   {
    //     stack.localValue[0] = valMass;
    //   }
    // }

  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    for( integer idof = 0; idof < m_numPhases; ++idof )
    {
      real64 const massNormalizer = LvArray::math::max( m_minNormalizer, m_phaseMass_n[ei][idof] );
      stack.localValue[0] += m_localResidual[stack.localRow + idof] * m_localResidual[stack.localRow + idof];
      stack.localNormalizer[0] += massNormalizer;
    }
  }


protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on porosity at the previous converged time step
  arrayView2d< real64 const > const m_porosity_n;
  /// View on mass at the previous converged time step
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const m_phaseMass_n;
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
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[in] solid the solid model
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( physicsSolverBaseKernels::NormType const normType,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   constitutive::CoupledSolidBase const & solid,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1],
                   real64 (& residualNormalizer)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, numPhases, subRegion, solid, minNormalizer );
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

/******************************** SolutionCheckKernel ********************************/

/**
 * @class SolutionCheckKernel
 * @brief Define the kernel for checking the updated solution
 */
class SolutionCheckKernel
{
public:

  /**
 * @brief Create a new kernel instance
 * @param[in] numPhases the number of fluid phases
 * @param[in] subRegion the element subregion
 * @param[in] rankOffset the offset of my MPI rank
 * @param[in] dofKey the string key to retrieve the degress of freedom numbers
 * @param[in] localSolution the solution vector on my MPI rank
 * @param[in] pressure pressure vector
 * @param[in] phaseVolFrac phase volume fraction vector
 * @param[in] pressureScalingFactor scaling factor for pressure
 * @param[in] compDensScalingFactor scaling factor for component density
 * @param[in] allowNegativePressure flag to allow negative pressure
 * @param[in] allowOutOfBoundsSaturation flag to allow out of bounds saturation
 * @param[out] StackVariables information about the solution update inside the subRegion
 */
  SolutionCheckKernel( integer const numPhases,
                       ElementSubRegionBase const & subRegion,
                       globalIndex const rankOffset,
                       string const dofKey,
                       arrayView1d< real64 const > const & localSolution,                       
                       arrayView1d< real64 const > const pressure,
                       arrayView2d< real64 const > const phaseVolFrac,
                       ScalingType const scalingType,
                       real64 const scalingFactor,
                       arrayView1d< real64 > pressureScalingFactor,
                       arrayView1d< real64 > saturationScalingFactor,
                       integer const allowNegativePressure,
                       integer const allowOutOfBoundsSaturation )
    : m_numPhase( numPhases ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_ghostRank( subRegion.ghostRank() ),
    m_localSolution( localSolution ),
    m_pressure( pressure ),
    m_phaseVolFrac( phaseVolFrac ),
    m_scalingType( scalingType ),
    m_scalingFactor( scalingFactor ),
    m_pressureScalingFactor( pressureScalingFactor ),
    m_saturationScalingFactor( saturationScalingFactor ),
    m_allowNegativePressure( allowNegativePressure ),
    m_allowOutOfBoundsSaturation( allowOutOfBoundsSaturation )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
    { } 

    StackVariables( real64 _localMinVal,
                    real64 _localMinPres,
                    real64 _localMinSat,
                    real64 _localMaxSat,
                    integer _localNumNegPressures,
                    integer _localNumNegSat,
                    integer _localNumUnitySat )
      :
      localMinVal( _localMinVal ),
      localRow( -1 ),
      localMinPres( _localMinPres ),
      localMinSat( _localMinSat ),
      localMaxSat( _localMaxSat ),
      localNumNegPressures( _localNumNegPressures ),
      localNumNegSat( _localNumNegSat ),
      localNumUnitySat( _localNumUnitySat )
    { }

    integer localMinVal;
    localIndex localRow;

    real64 localMinPres;
    real64 localMinSat;
    real64 localMaxSat;

    integer localNumNegPressures;
    integer localNumNegSat;
    integer localNumUnitySat;

  };

  /**
   * @brief Getter for the ghost rank
   * @param[in] i the looping index of the element/node/face
   * @return the ghost rank of the element/node/face
   */
  GEOS_HOST_DEVICE
  integer ghostRank( localIndex const i ) const
  { return m_ghostRank( i ); }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static StackVariables
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, integer > globalMinVal( 1 );

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minPres( 0.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minSat( 0.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxSat( 1.0 );

    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegPressures( 0 );
    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegSat( 0 );
    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numUnitySat( 0 );

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.ghostRank( ei ) >= 0 )
      {
        return;
      }

      StackVariables stack;
      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );

      globalMinVal.min( stack.localMinVal );

      minPres.min( stack.localMinPres );
      minSat.min( stack.localMinSat );
      maxSat.max( stack.localMaxSat );

      numNegPressures += stack.localNumNegPressures;
      numNegSat += stack.localNumNegSat;
      numUnitySat += stack.localNumUnitySat;
    } );

    return StackVariables( globalMinVal.get(),
                           minPres.get(),
                           minSat.get(),
                           maxSat.get(),
                           numNegPressures.get(),
                           numNegSat.get(),
                           numUnitySat.get() );
  }

  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    stack.localMinVal = 1;

    stack.localRow = m_dofNumber[ei] - m_rankOffset;

    stack.localMinPres = 0.0;
    stack.localMinSat = 0.0;
    stack.localMaxSat = 1.0;

    stack.localNumNegPressures = 0;
    stack.localNumNegSat = 0;
    stack.localNumUnitySat = 0;
  }

  /**
   * @brief Compute the local value
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    bool const localScaling = m_scalingType == ScalingType::Local;

    real64 const newPres = m_pressure[ei] + (localScaling ? m_pressureScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow];
    if( newPres < 0 )
    {
      if( !m_allowNegativePressure )
      {
        stack.localMinVal = 0;
      }
      stack.localNumNegPressures += 1;
      if( newPres < stack.localMinPres )
      {
        stack.localMinPres = newPres;
      }
    }

    for( integer ip = 0; ip < m_numPhase; ++ip )
    {
      real64 const newSat = ip == 0 ? m_phaseVolFrac[ei][ip] + (localScaling ? m_saturationScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow + 1] :  // S = S + dS || S = S - dS
                                      m_phaseVolFrac[ei][ip] - (localScaling ? m_saturationScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow + 1] ;
      
      if( newSat < 0.0 )
      {        
        if( !m_allowOutOfBoundsSaturation )
        {
          stack.localMinVal = 0;
        }
        stack.localNumNegSat += 1;
        if( newSat < stack.localMinSat )
        {
          stack.localMinSat = newSat;
        }
      }
      else if( newSat > 1.0 )
      {
        if( !m_allowOutOfBoundsSaturation )
        {
          stack.localMinVal = 0;
        }
        stack.localNumUnitySat += 1;
        if( newSat > stack.localMaxSat )
        {
          stack.localMaxSat = newSat;
        }
      }      
    }
  }

protected:

  /// Number of components
  integer const m_numPhase;  

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;  

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_ghostRank;

  /// View on the local residual
  arrayView1d< real64 const > const m_localSolution;

  /// View on the primary variables
  arrayView1d< real64 const > const m_pressure;
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const m_phaseVolFrac;

  /// scaling type (global or local)
  ScalingType const m_scalingType;

  /// scaling factor
  real64 const m_scalingFactor;

  /// View on the scaling factors
  arrayView1d< real64 > const m_pressureScalingFactor;
  arrayView1d< real64 > const m_saturationScalingFactor;  

  /// flag to allow negative pressure values
  integer const m_allowNegativePressure;

  /// flag to allow the component density chopping
  integer const m_allowOutOfBoundsSaturation;
};

/**
 * @class SolutionCheckKernelFactory
 */
class SolutionCheckKernelFactory
{
public:

/**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhases the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localSolution the solution vector on my MPI rank
   * @param[in] pressure pressure vector
   * @param[in] phaseVolFrac phase volume fraction vector
   * @param[in] pressureScalingFactor scaling factor for pressure
   * @param[in] compDensScalingFactor scaling factor for component density
   * @param[in] allowNegativePressure flag to allow negative pressure
   * @param[in] allowOutOfBoundsSaturation flag to allow out of bounds saturation
   * @param[in] allowSatChopping flag to allow saturation chopping
   * @param[in] satChangeLimit saturation change limit
   * @param[out] StackVariables information about the solution update inside the subRegion
   */
  template< typename POLICY >
  static SolutionCheckKernel::StackVariables
  createAndLaunch( integer const numPhases,
                   ElementSubRegionBase const & subRegion,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localSolution,                   
                   arrayView1d< real64 const > const pressure,
                   arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac,
                   ScalingType const scalingType,
                   real64 const scalingFactor,
                   arrayView1d< real64 > pressureScalingFactor,
                   arrayView1d< real64 > saturationScalingFactor,
                   integer const allowNegativePressure,
                   integer const allowOutOfBoundsSaturation )
  {                                                
    SolutionCheckKernel kernel( numPhases, subRegion, rankOffset, dofKey, localSolution,
                                pressure, phaseVolFrac, scalingType, scalingFactor,
                                pressureScalingFactor, saturationScalingFactor,
                                allowNegativePressure, allowOutOfBoundsSaturation );
    return SolutionCheckKernel::launch< POLICY >( subRegion.size(), kernel );                                           
  }

};


/******************************** SolutionScalingKernel ********************************/

/**
 * @class SolutionScalingKernel
 * @brief Define the kernel for scaling the Newton update
 */
class SolutionScalingKernel
{
public:

/**
   * @brief Create a new kernel instance
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] maxRelativeCompDensChange the max allowed comp density change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] compDens the component density vector
   * @param[in] pressureScalingFactor the pressure local scaling factor
   * @param[in] compDensScalingFactor the component density local scaling factor
   */
  SolutionScalingKernel( real64 const maxRelativePresChange,
                         real64 const maxAbsolutePresChange,
                         real64 const maxRelativeSatChange,
                         real64 const maxAbsoluteSatChange,
                         globalIndex const rankOffset,
                         integer const numPhases,
                         string const dofKey,
                         ElementSubRegionBase const & subRegion,
                         arrayView1d< real64 const > const localSolution,
                         arrayView1d< real64 const > const pressure,
                         arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac,
                         arrayView1d< real64 > pressureScalingFactor,
                         arrayView1d< real64 > saturationScalingFactor )
    : m_rankOffset( rankOffset ),
    m_numPhases( numPhases ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_ghostRank( subRegion.ghostRank() ),
    m_localSolution( localSolution ),
    m_pressure( pressure ),   
    m_phaseVolFrac( phaseVolFrac ),  
    m_pressureScalingFactor( pressureScalingFactor ),
    m_saturationScalingFactor( saturationScalingFactor ),
    m_maxRelativePresChange( maxRelativePresChange ),
    m_maxAbsolutePresChange( maxAbsolutePresChange ),
    m_maxRelativeSatChange( maxRelativeSatChange ),
    m_maxAbsoluteSatChange( maxAbsoluteSatChange )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
    { }

    StackVariables( real64 _localMinVal,
                    real64 _localMaxDeltaPres,
                    localIndex _localMaxDeltaPresLoc,                    
                    real64 _localMaxDeltaSat,
                    localIndex _localMaxDeltaSatLoc,
                    real64 _localMinPresScalingFactor,
                    real64 _localMinSatScalingFactor )
      :
      localRow( -1 ),
      localMinVal( _localMinVal ),
      localMaxDeltaPres( _localMaxDeltaPres ),
      localMaxDeltaPresLoc( _localMaxDeltaPresLoc ),      
      localMaxDeltaSat( _localMaxDeltaSat ),
      localMaxDeltaSatLoc( _localMaxDeltaSatLoc ),
      localMinPresScalingFactor( _localMinPresScalingFactor ),
      localMinSatScalingFactor( _localMinSatScalingFactor )
    { }

    localIndex localRow;
    real64 localMinVal;    

    real64 localMaxDeltaPres;
    localIndex localMaxDeltaPresLoc;    
    real64 localMaxDeltaSat;
    localIndex localMaxDeltaSatLoc;

    real64 localMinPresScalingFactor;
    real64 localMinSatScalingFactor;

  };

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static StackVariables
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > globalScalingFactor( 1.0 );

    RAJA::ReduceMaxLoc< ReducePolicy< POLICY >, real64 > maxDeltaPres( std::numeric_limits< real64 >::min(), -1 );    
    RAJA::ReduceMaxLoc< ReducePolicy< POLICY >, real64 > maxDeltaSat( std::numeric_limits< real64 >::min(), -1 );

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minPresScalingFactor( 1.0 );    
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minSatScalingFactor( 1.0 );

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.ghostRank( ei ) >= 0 )
      {
        return;
      }

      StackVariables stack;
      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );

      globalScalingFactor.min( stack.localMinVal );

      maxDeltaPres.maxloc( stack.localMaxDeltaPres, ei );      
      maxDeltaSat.maxloc( stack.localMaxDeltaSat, ei );

      minPresScalingFactor.min( stack.localMinPresScalingFactor );      
      minSatScalingFactor.min( stack.localMinSatScalingFactor );
    } );

    return StackVariables( globalScalingFactor.get(),
                           maxDeltaPres.get(),
                           maxDeltaPres.getLoc(),                           
                           maxDeltaSat.get(),
                           maxDeltaSat.getLoc(),
                           minPresScalingFactor.get(),                           
                           minSatScalingFactor.get() );
  }

  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    stack.localRow = m_dofNumber[ei] - m_rankOffset;
    stack.localMinVal = 1.0;

    stack.localMaxDeltaPres = 0.0;
    stack.localMaxDeltaPresLoc = -1;
    stack.localMaxDeltaSat = 0.0;
    stack.localMaxDeltaSatLoc =-1;

    stack.localMinPresScalingFactor = 1.0;
    stack.localMinSatScalingFactor = 1.0;
  }

  /**
   * @brief Compute the local value
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    // compute the change in pressure
    real64 const pres = m_pressure[ei];
    real64 const absPresChange = LvArray::math::abs( m_localSolution[stack.localRow] );
    if( stack.localMaxDeltaPres < absPresChange )
    {
      stack.localMaxDeltaPres = absPresChange;
    }

    // compute pressure scaling factor
    real64 presScalingFactor = 1.0;
    // when enabled, absolute change scaling has a priority over relative change
    if( m_maxAbsolutePresChange > 0.0 ) // maxAbsolutePresChange <= 0.0 means that absolute scaling is disabled
    {
      if( absPresChange > m_maxAbsolutePresChange )
      {
        presScalingFactor = m_maxAbsolutePresChange / absPresChange;
      }
    }
    else if( m_maxRelativePresChange > 0.0 && pres > minDensForDivision )
    {
      real64 const relativePresChange = absPresChange / pres;
      if( relativePresChange > m_maxRelativePresChange )
      {
        presScalingFactor = m_maxRelativePresChange / relativePresChange;
      }
    }
    m_pressureScalingFactor[ei] = presScalingFactor;
    if( stack.localMinVal > presScalingFactor )
    {
      stack.localMinVal = presScalingFactor;
    }
    if( stack.localMinPresScalingFactor > presScalingFactor )
    {
      stack.localMinPresScalingFactor = presScalingFactor;
    }


    // compute the change in saturation
    real64 const phaseVolFrac = m_phaseVolFrac[ei][0];
    real64 const absSatChange = LvArray::math::abs( m_localSolution[stack.localRow + 1] );
    if( stack.localMaxDeltaSat < absSatChange )
    {
      stack.localMaxDeltaSat = absSatChange;
    }
    // compute saturation scaling factor
    real64 satScalingFactor = 1.0;
    // when enabled, absolute change scaling has a priority over relative change
    if( m_maxAbsoluteSatChange > 0.0 ) // maxAbsoluteSatChange <= 0.0 means that absolute scaling is disabled
    {
      if( absSatChange > m_maxAbsoluteSatChange )
      {
        satScalingFactor = m_maxAbsoluteSatChange / absSatChange;
      }
    }
    else if( m_maxRelativeSatChange > 0.0 && phaseVolFrac > minDensForDivision )
    {
      real64 const relativeSatChange = absSatChange / phaseVolFrac;
      if( relativeSatChange > m_maxRelativeSatChange )
      {
        satScalingFactor = m_maxRelativeSatChange / relativeSatChange;
      }
    }
    m_saturationScalingFactor[ei] = satScalingFactor;
    if( stack.localMinVal > satScalingFactor )
    {
      stack.localMinVal = satScalingFactor;
    }
    if( stack.localMinSatScalingFactor > satScalingFactor )
    {
      stack.localMinSatScalingFactor = satScalingFactor;
    }
  }

  /**
   * @brief Getter for the ghost rank
   * @param[in] i the looping index of the element/node/face
   * @return the ghost rank of the element/node/face
   */
  GEOS_HOST_DEVICE
  integer ghostRank( localIndex const i ) const
  { return m_ghostRank( i ); }

protected:

  static constexpr real64 minDensForDivision = 1e-10;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Number of phases
  real64 const m_numPhases;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_ghostRank;

  /// View on the local residual
  arrayView1d< real64 const > const m_localSolution;

  /// View on the primary variables
  arrayView1d< real64 const > const m_pressure;
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const m_phaseVolFrac;

  /// View on the scaling factors
  arrayView1d< real64 > const m_pressureScalingFactor;
  arrayView1d< real64 > const m_saturationScalingFactor;

  /// Max allowed changes in primary variables
  real64 const m_maxRelativePresChange;
  real64 const m_maxAbsolutePresChange;
  real64 const m_maxRelativeSatChange;
  real64 const m_maxAbsoluteSatChange;

};


/**
 * @class SolutionScalingKernelFactory
 */
class SolutionScalingKernelFactory
{
public:

  /**
   * @brief Create and launch the kernel computing the scaling factor
   * @tparam POLICY the kernel policy
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxRelativeSatChange the max allowed relative saturation change
   * @param[in] maxAbsoluteSatChange the max allowed absolute saturation change
   * @param[in] pressure the pressure vector
   * @param[in] phaseVolFrac the phase volume fraction vector
   * @param[in] pressureScalingFactor the scaling factor for pressure
   * @param[in] saturationScalingFactor the scaling factor for saturation
   * @param[in] rankOffset the rank offset
   * @param[in] numPhases the number of phases
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[out] StackVariables information about the scaling factors inside the subRegion
   */
  template< typename POLICY >
  static SolutionScalingKernel::StackVariables
  createAndLaunch( real64 const maxRelativePresChange,
                   real64 const maxAbsolutePresChange,
                   real64 const maxRelativeSatChange,
                   real64 const maxAbsoluteSatChange,
                   arrayView1d< real64 const > const pressure,
                   arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac,
                   arrayView1d< real64 > pressureScalingFactor,
                   arrayView1d< real64 > saturationScalingFactor,
                   globalIndex const rankOffset,
                   integer const numPhases,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    SolutionScalingKernel kernel( maxRelativePresChange, maxAbsolutePresChange, maxRelativeSatChange, maxAbsoluteSatChange, rankOffset,
                                  numPhases, dofKey, subRegion, localSolution, pressure, phaseVolFrac, pressureScalingFactor, saturationScalingFactor );
    return SolutionScalingKernel::launch< POLICY >( subRegion.size(), kernel );
  }
};

} // namespace immiscible multiphasekernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_MULTIPHASEKERNELS_HPP
