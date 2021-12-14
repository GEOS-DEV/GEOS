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
 * @file CompositionalMultiphaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVMKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureExtrinsicData.hpp"
#include "constitutive/fluid/MultiFluidExtrinsicData.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityExtrinsicData.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"


namespace geosx
{

namespace CompositionalMultiphaseFVMKernels
{

using namespace constitutive;

/******************************** PhaseMobilityKernel ********************************/

/**
 * @brief Functions to compute phase mobilities and derivatives from density, viscosity and relperm
 */
struct PhaseMobilityKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  compute( arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dComp,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseVisc,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseVisc_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
           arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseMob,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & dPhaseMob_dPres,
           arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseMob,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseMob,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp );
};


/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< localIndex NUM_COMP, localIndex NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel
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

  using DofNumberAccessor =
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using CompFlowAccessors =
    ElementRegionManager::
      StencilAccessors< extrinsicMeshData::ghostRank,   // 0
                        extrinsicMeshData::flow::gravityCoefficient,   // 1
                        extrinsicMeshData::flow::pressure,   // 2
                        extrinsicMeshData::flow::deltaPressure,   // 3
                        extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity,   // 4
                        extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure,   // 5
                        extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity,   // 6
                        extrinsicMeshData::flow::phaseMobility,   // 7
                        extrinsicMeshData::flow::dPhaseMobility_dPressure,   // 8
                        extrinsicMeshData::flow::dPhaseMobility_dGlobalCompDensity >;   // 9

  using MultiFluidAccessors =
    ElementRegionManager::
      StencilMaterialAccessors< extrinsicMeshData::multifluid::phaseMassDensity,   // 0
                                extrinsicMeshData::multifluid::dPhaseMassDensity_dPressure,   // 1
                                extrinsicMeshData::multifluid::dPhaseMassDensity_dGlobalCompFraction,   // 2
                                extrinsicMeshData::multifluid::phaseCompFraction,   // 3
                                extrinsicMeshData::multifluid::dPhaseCompFraction_dPressure,   // 4
                                extrinsicMeshData::multifluid::dPhaseCompFraction_dGlobalCompFraction >;   // 5

  using CapPressureAccessors =
    ElementRegionManager::
      StencilMaterialAccessors< extrinsicMeshData::cappres::phaseCapPressure,   // 0
                                extrinsicMeshData::cappres::dPhaseCapPressure_dPhaseVolFraction >;   // 1

  using PermeabilityAccessors
    = ElementRegionManager::
        StencilMaterialAccessors< extrinsicMeshData::permeability::permeability, // 0
                                  extrinsicMeshData::permeability::dPerm_dPressure >; // 1


  /// Compile time value for the number of components
  static constexpr localIndex numComp = NUM_COMP;

  /// Compute time value for the number of degrees of freedom
  static constexpr localIndex numDof = NUM_DOF;

  /// Compute time value for the number of equations (all of them, except the volume balance equation)
  static constexpr localIndex numEqn = NUM_DOF-1;

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
   * @param[in] dofKey string to get the element degrees of freedom numbers
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
  FaceBasedAssemblyKernel( localIndex const numPhases,
                           globalIndex const rankOffset,
                           integer const capPressureFlag,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           CompFlowAccessors const & compFlowAccessors,
                           MultiFluidAccessors const & multiFluidAccessors,
                           CapPressureAccessors const & capPressureAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    : m_numPhases( numPhases ),
    m_rankOffset( rankOffset ),
    m_capPressureFlag( capPressureFlag ),
    m_dt( dt ),
    m_stencilWrapper( stencilWrapper ),
    m_seri( stencilWrapper.getElementRegionIndices() ),
    m_sesri( stencilWrapper.getElementSubRegionIndices() ),
    m_sei( stencilWrapper.getElementIndices() ),
    m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_ghostRank( std::get< 0 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_gravCoef( std::get< 1 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_permeability( std::get< 0 >( permeabilityAccessors.accessors ).toNestedViewConst() ),
    m_dPerm_dPres( std::get< 1 >( permeabilityAccessors.accessors ).toNestedViewConst() ),
    m_pres( std::get< 2 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_dPres( std::get< 3 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_dCompFrac_dCompDens( std::get< 4 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_dPhaseVolFrac_dPres( std::get< 5 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_dPhaseVolFrac_dCompDens( std::get< 6 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_phaseMob( std::get< 7 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_dPhaseMob_dPres( std::get< 8 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_dPhaseMob_dCompDens( std::get< 9 >( compFlowAccessors.accessors ).toNestedViewConst() ),
    m_phaseMassDens( std::get< 0 >( multiFluidAccessors.accessors ).toNestedViewConst() ),
    m_dPhaseMassDens_dPres( std::get< 1 >( multiFluidAccessors.accessors ).toNestedViewConst() ),
    m_dPhaseMassDens_dComp( std::get< 2 >( multiFluidAccessors.accessors ).toNestedViewConst() ),
    m_phaseCompFrac( std::get< 3 >( multiFluidAccessors.accessors ).toNestedViewConst() ),
    m_dPhaseCompFrac_dPres( std::get< 4 >( multiFluidAccessors.accessors ).toNestedViewConst() ),
    m_dPhaseCompFrac_dComp( std::get< 5 >( multiFluidAccessors.accessors ).toNestedViewConst() ),
    m_phaseCapPressure( std::get< 0 >( capPressureAccessors.accessors ).toNestedViewConst() ),
    m_dPhaseCapPressure_dPhaseVolFrac( std::get< 1 >( capPressureAccessors.accessors ).toNestedViewConst() ),
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

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOSX_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : stencilSize( size ),
      numFluxElems( numElems ),
      compFlux( numComp ),
      dCompFlux_dP( size, numComp ),
      dCompFlux_dC( size, numComp, numComp ),
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

    // Component fluxes and derivatives

    /// Component fluxes
    stackArray1d< real64, numComp > compFlux;
    /// Derivatives of component fluxes wrt pressure
    stackArray2d< real64, maxStencilSize * numComp > dCompFlux_dP;
    /// Derivatives of component fluxes wrt component densities
    stackArray3d< real64, maxStencilSize * numComp * numComp > dCompFlux_dC;

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
  GEOSX_HOST_DEVICE
  localIndex stencilSize( localIndex const iconn ) const
  { return meshMapUtilities::size1( m_sei, iconn ); }

  /**
   * @brief Getter for the number of elements at this connection
   * @param[in] iconn the connection index
   * @return the number of elements at this connection
   */
  GEOSX_HOST_DEVICE
  localIndex numPointsInFlux( localIndex const iconn ) const
  { return m_stencilWrapper.numPointsInFlux( iconn ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
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
   * @tparam FUNC1 the type of the function that can be used to customize the computation of the phase fluxes
   * @tparam FUNC2 the type of the function that can be used to customize the assembly into the local Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] phaseFluxKernelOp the function used to customize the computation of the phase fluxes
   * @param[in] localFluxJacobianKernelOp the function used to customize the computation of the assembly into the local Jacobian
   */
  template< typename FUNC1 = CompositionalMultiphaseBaseKernels::NoOpFunc,
            typename FUNC2 = CompositionalMultiphaseBaseKernels::NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC1 && phaseFluxKernelOp = CompositionalMultiphaseBaseKernels::NoOpFunc{},
                    FUNC2 && localFluxJacobianKernelOp = CompositionalMultiphaseBaseKernels::NoOpFunc{} ) const
  {
    // first, compute the transmissibilities at this face
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      // clear working arrays
      real64 densMean{};
      stackArray1d< real64, maxNumElems > dDensMean_dP( stack.numFluxElems );
      stackArray2d< real64, maxNumElems * numComp > dDensMean_dC( stack.numFluxElems, numComp );

      // create local work arrays
      real64 phaseFlux{};
      real64 dPhaseFlux_dP[maxStencilSize]{};
      real64 dPhaseFlux_dC[maxStencilSize][numComp]{};

      real64 presGrad{};
      stackArray1d< real64, maxStencilSize > dPresGrad_dP( stack.stencilSize );
      stackArray2d< real64, maxStencilSize *numComp > dPresGrad_dC( stack.stencilSize, numComp );

      real64 gravHead{};
      stackArray1d< real64, maxNumElems > dGravHead_dP( stack.numFluxElems );
      stackArray2d< real64, maxNumElems * numComp > dGravHead_dC( stack.numFluxElems, numComp );

      real64 dCapPressure_dC[numComp]{};

      // Working array
      real64 dProp_dC[numComp]{};

      // calculate quantities on primary connected cells
      for( integer i = 0; i < stack.numFluxElems; ++i )
      {
        localIndex const er  = m_seri( iconn, i );
        localIndex const esr = m_sesri( iconn, i );
        localIndex const ei  = m_sei( iconn, i );

        // density
        real64 const density  = m_phaseMassDens[er][esr][ei][0][ip];
        real64 const dDens_dP = m_dPhaseMassDens_dPres[er][esr][ei][0][ip];

        applyChainRule( numComp,
                        m_dCompFrac_dCompDens[er][esr][ei],
                        m_dPhaseMassDens_dComp[er][esr][ei][0][ip],
                        dProp_dC );

        // average density and derivatives
        densMean += 0.5 * density;
        dDensMean_dP[i] = 0.5 * dDens_dP;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dDensMean_dC[i][jc] = 0.5 * dProp_dC[jc];
        }
      }

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      for( integer i = 0; i < stack.stencilSize; ++i )
      {
        localIndex const er  = m_seri( iconn, i );
        localIndex const esr = m_sesri( iconn, i );
        localIndex const ei  = m_sei( iconn, i );

        // capillary pressure
        real64 capPressure     = 0.0;
        real64 dCapPressure_dP = 0.0;

        for( integer ic = 0; ic < numComp; ++ic )
        {
          dCapPressure_dC[ic] = 0.0;
        }

        if( m_capPressureFlag )
        {
          capPressure = m_phaseCapPressure[er][esr][ei][0][ip];

          for( integer jp = 0; jp < m_numPhases; ++jp )
          {
            real64 const dCapPressure_dS = m_dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
            dCapPressure_dP += dCapPressure_dS * m_dPhaseVolFrac_dPres[er][esr][ei][jp];

            for( integer jc = 0; jc < numComp; ++jc )
            {
              dCapPressure_dC[jc] += dCapPressure_dS * m_dPhaseVolFrac_dCompDens[er][esr][ei][jp][jc];
            }
          }
        }

        presGrad += stack.transmissibility[0][i] * (m_pres[er][esr][ei] + m_dPres[er][esr][ei] - capPressure);
        dPresGrad_dP[i] += stack.transmissibility[0][i] * (1 - dCapPressure_dP)
                           + stack.dTrans_dPres[0][i] * (m_pres[er][esr][ei] + m_dPres[er][esr][ei] - capPressure);
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPresGrad_dC[i][jc] += -stack.transmissibility[0][i] * dCapPressure_dC[jc];
        }

        real64 const gravD     = stack.transmissibility[0][i] * m_gravCoef[er][esr][ei];
        real64 const dGravD_dP = stack.dTrans_dPres[0][i] * m_gravCoef[er][esr][ei];

        // the density used in the potential difference is always a mass density
        // unlike the density used in the phase mobility, which is a mass density
        // if useMass == 1 and a molar density otherwise
        gravHead += densMean * gravD;

        // need to add contributions from both cells the mean density depends on
        for( integer j = 0; j < stack.numFluxElems; ++j )
        {
          dGravHead_dP[j] += dDensMean_dP[j] * gravD + dGravD_dP * densMean;
          for( integer jc = 0; jc < numComp; ++jc )
          {
            dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
          }
        }
      }

      // *** upwinding ***

      // compute phase potential gradient
      real64 const potGrad = presGrad - gravHead;

      // choose upstream cell
      localIndex const k_up = (potGrad >= 0) ? 0 : 1;

      localIndex const er_up  = m_seri( iconn, k_up );
      localIndex const esr_up = m_sesri( iconn, k_up );
      localIndex const ei_up  = m_sei( iconn, k_up );

      real64 const mobility = m_phaseMob[er_up][esr_up][ei_up][ip];

      // skip the phase flux if phase not present or immobile upstream
      if( LvArray::math::abs( mobility ) < 1e-20 ) // TODO better constant
      {
        continue;
      }

      // pressure gradient depends on all points in the stencil
      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc];
        }
      }

      // gravitational head depends only on the two cells connected (same as mean density)
      for( integer ke = 0; ke < stack.numFluxElems; ++ke )
      {
        dPhaseFlux_dP[ke] -= dGravHead_dP[ke];
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] -= dGravHead_dC[ke][jc];
        }
      }

      // compute the phase flux and derivatives using upstream cell mobility
      phaseFlux = mobility * potGrad;
      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        dPhaseFlux_dP[ke] *= mobility;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] *= mobility;
        }
      }

      real64 const dMob_dP  = m_dPhaseMob_dPres[er_up][esr_up][ei_up][ip];
      arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 > dPhaseMob_dCompSub =
        m_dPhaseMob_dCompDens[er_up][esr_up][ei_up][ip];

      // add contribution from upstream cell mobility derivatives
      dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[k_up][jc] += dPhaseMob_dCompSub[jc] * potGrad;
      }

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d< real64 const, multifluid::USD_PHASE_COMP-3 > phaseCompFracSub =
        m_phaseCompFrac[er_up][esr_up][ei_up][0][ip];
      arraySlice1d< real64 const, multifluid::USD_PHASE_COMP-3 > dPhaseCompFrac_dPresSub =
        m_dPhaseCompFrac_dPres[er_up][esr_up][ei_up][0][ip];
      arraySlice2d< real64 const, multifluid::USD_PHASE_COMP_DC-3 > dPhaseCompFrac_dCompSub =
        m_dPhaseCompFrac_dComp[er_up][esr_up][ei_up][0][ip];

      // compute component fluxes and derivatives using upstream cell composition
      for( integer ic = 0; ic < numComp; ++ic )
      {
        real64 const ycp = phaseCompFracSub[ic];
        stack.compFlux[ic] += phaseFlux * ycp;

        // derivatives stemming from phase flux
        for( integer ke = 0; ke < stack.stencilSize; ++ke )
        {
          stack.dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
          }
        }

        // additional derivatives stemming from upstream cell phase composition
        stack.dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFrac_dPresSub[ic];

        // convert derivatives of comp fraction w.r.t. comp fractions to derivatives w.r.t. comp densities
        applyChainRule( numComp,
                        m_dCompFrac_dCompDens[er_up][esr_up][ei_up],
                        dPhaseCompFrac_dCompSub[ic],
                        dProp_dC );
        for( integer jc = 0; jc < numComp; ++jc )
        {
          stack.dCompFlux_dC[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
        }
      }

      // call the lambda in the phase loop to allow the reuse of the phase fluxes and their derivatives
      // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase
      // note: more variables may need to be passed, but it is hard to tell which ones will be needed for now
      phaseFluxKernelOp( ip, er_up, esr_up, ei_up, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

    }

    // *** end of upwinding

    // populate local flux vector and derivatives
    for( integer ic = 0; ic < numComp; ++ic )
    {
      stack.localFlux[ic]           =  m_dt * stack.compFlux[ic];
      stack.localFlux[numComp + ic] = -m_dt * stack.compFlux[ic];

      for( integer ke = 0; ke < stack.stencilSize; ++ke )
      {
        localIndex const localDofIndexPres = ke * numDof;
        stack.localFluxJacobian[ic][localDofIndexPres]           =  m_dt * stack.dCompFlux_dP[ke][ic];
        stack.localFluxJacobian[numComp + ic][localDofIndexPres] = -m_dt * stack.dCompFlux_dP[ke][ic];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
          stack.localFluxJacobian[ic][localDofIndexComp]           =  m_dt * stack.dCompFlux_dC[ke][ic][jc];
          stack.localFluxJacobian[numComp + ic][localDofIndexComp] = -m_dt * stack.dCompFlux_dC[ke][ic][jc];
        }
      }
    }

    // call the lambda to allow assembly into the localFluxJacobian
    // possible uses: - assemble the derivatives of compFlux wrt temperature into localFluxJacobian;
    //                - assemble energyFlux into localFlux and localFluxJacobian
    localFluxJacobianKernelOp();

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
    using namespace CompositionalMultiphaseUtilities;

    // Apply equation/variable change transformation(s)
    stackArray1d< real64, maxStencilSize * numDof > work( stack.stencilSize * numDof );
    shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numDof*stack.stencilSize, stack.numFluxElems,
                                                             stack.localFluxJacobian, work );
    shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( numComp, stack.numFluxElems,
                                                               stack.localFlux );

    // Add to residual/jacobian
    for( integer i = 0; i < stack.numFluxElems; ++i )
    {
      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( m_localMatrix.numRows(), localRow + numComp );

        for( integer ic = 0; ic < numComp; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + ic], stack.localFlux[i * numComp + ic] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
            ( localRow + ic,
            stack.dofColIndices.data(),
            stack.localFluxJacobian[i * numComp + ic].dataIfContiguous(),
            stack.stencilSize * numDof );
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

    forAll< POLICY >( numConnections, [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent.stencilSize( iconn ),
                                                  kernelComponent.numPointsInFlux( iconn ) );

      kernelComponent.setup( iconn, stack );
      kernelComponent.computeFlux( iconn, stack );
      kernelComponent.complete( iconn, stack );
    } );
  }

private:

  /// Number of fluid phases
  localIndex const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Flag to specify whether capillary pressure is used or not
  integer const m_capPressureFlag;

  /// Time step size
  real64 const m_dt;

  // Stencil information

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;

  /// Views on dof numbers and ghost rank numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > m_permeability;
  ElementViewConst< arrayView3d< real64 const > > m_dPerm_dPres;

  // Primary and secondary variables

  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres;
  ElementViewConst< arrayView1d< real64 const > > const m_dPres;

  /// Views on derivatives of phase volume fractions and comp fractions
  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const m_dCompFrac_dCompDens;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_dPhaseVolFrac_dPres;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseVolFrac_dCompDens;

  /// Views on phase mobilities
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseMob;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_dPhaseMob_dPres;
  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const m_dPhaseMob_dCompDens;

  /// Views on phase mass densities
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseMassDens;
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_dPhaseMassDens_dPres;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const m_dPhaseMassDens_dComp;

  /// Views on phase component fractions
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const m_phaseCompFrac;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const m_dPhaseCompFrac_dPres;
  ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const m_dPhaseCompFrac_dComp;

  /// Views on phase capillary pressure
  ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const m_phaseCapPressure;
  ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const m_dPhaseCapPressure_dPhaseVolFrac;

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

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
  createAndLaunch( localIndex const numComps,
                   localIndex const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const capPressureFlag,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   arrayView1d< string const > const & targetRegionNames,
                   arrayView1d< string const > const & fluidModelNames,
                   arrayView1d< string const > const & capPresModelNames,
                   arrayView1d< string const > const & permeabilityModelNames,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    CompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      localIndex constexpr NUM_COMP = NC();
      localIndex constexpr NUM_DOF = NC()+1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      using KERNEL_TYPE = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::MultiFluidAccessors multiFluidAccessors( elemManager, solverName, targetRegionNames, fluidModelNames );
      typename KERNEL_TYPE::CapPressureAccessors capPressureAccessors( elemManager, solverName, targetRegionNames, capPresModelNames );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName, targetRegionNames, permeabilityModelNames );

      KERNEL_TYPE kernel( numPhases, rankOffset, capPressureFlag, stencilWrapper, dofNumberAccessor,
                          compFlowAccessors, multiFluidAccessors, capPressureAccessors, permeabilityAccessors,
                          dt, localMatrix, localRhs );
      KERNEL_TYPE::template launch< POLICY, KERNEL_TYPE >( stencilWrapper.size(), kernel );
    } );
  }
};

/******************************** CFLFluxKernel ********************************/

/**
 * @brief Functions to compute the (outflux) total volumetric flux needed in the calculation of CFL numbers
 */
struct CFLFluxKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename VIEWTYPE >
  using ElementView = ElementRegionManager::ElementView< VIEWTYPE >;

  using CompFlowAccessors =
    ElementRegionManager::
      StencilAccessors< extrinsicMeshData::flow::pressure, // 0
                        extrinsicMeshData::flow::gravityCoefficient, // 1
                        extrinsicMeshData::flow::phaseVolumeFraction, // 2
                        extrinsicMeshData::flow::phaseOutflux, // 3
                        extrinsicMeshData::flow::componentOutflux >; // 4

  using MultiFluidAccessors =
    ElementRegionManager::
      StencilMaterialAccessors< extrinsicMeshData::multifluid::phaseViscosity, // 0
                                extrinsicMeshData::multifluid::phaseDensity, // 1
                                extrinsicMeshData::multifluid::phaseMassDensity, // 2
                                extrinsicMeshData::multifluid::phaseCompFraction >; // 3

  using PermeabilityAccessors =
    ElementRegionManager::
      StencilMaterialAccessors< extrinsicMeshData::permeability::permeability, // 0
                                extrinsicMeshData::permeability::dPerm_dPressure >; // 1


  using RelPermAccessors =
    ElementRegionManager::
      StencilMaterialAccessors< extrinsicMeshData::relperm::phaseRelPerm >; // 0

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL_SIZE >
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const numPhases,
           localIndex const stencilSize,
           real64 const & dt,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           real64 const (&transmissibility)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
           ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & phaseRelPerm,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseVisc,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux,
           ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux );

  template< localIndex NC, typename STENCILWRAPPER_TYPE >
  static void
  launch( localIndex const numPhases,
          real64 const & dt,
          STENCILWRAPPER_TYPE const & stencil,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & phaseRelPerm,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseVisc,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux,
          ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux );
};

/******************************** CFLKernel ********************************/

/**
 * @brief Functions to compute the CFL number using the phase volumetric outflux and the component mass outflux in each cell
 */
struct CFLKernel
{

  static constexpr real64 minPhaseMobility = 1e-12;
  static constexpr real64 minComponentFraction = 1e-12;

  template< localIndex NP >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  computePhaseCFL( real64 const & poreVol,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
                   arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > phaseRelPerm,
                   arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > dPhaseRelPerm_dPhaseVolFrac,
                   arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseVisc,
                   arraySlice1d< real64 const, compflow::USD_PHASE- 1 > phaseOutflux,
                   real64 & phaseCFLNumber );

  template< localIndex NC >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  computeCompCFL( real64 const & poreVol,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compDens,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compFrac,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compOutflux,
                  real64 & compCFLNumber );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
          arrayView1d< real64 const > const & volume,
          arrayView2d< real64 const > const & porosity,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux,
          arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux,
          arrayView1d< real64 > const & phaseCFLNumber,
          arrayView1d< real64 > const & compCFLNumber,
          real64 & maxPhaseCFLNumber,
          real64 & maxCompCFLNumber );

};

/******************************** AquiferBCKernel ********************************/

/**
 * @brief Functions to assemble aquifer boundary condition contributions to residual and Jacobian
 */
struct AquiferBCKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using CompFlowAccessors =
    ElementRegionManager::
      StencilAccessors< extrinsicMeshData::ghostRank, // 0
                        extrinsicMeshData::flow::pressure, // 1
                        extrinsicMeshData::flow::deltaPressure, // 2
                        extrinsicMeshData::flow::gravityCoefficient, // 3
                        extrinsicMeshData::flow::phaseVolumeFraction, // 4
                        extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure, // 5
                        extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity, // 6
                        extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity >; // 7

  using MultiFluidAccessors =
    ElementRegionManager::
      StencilMaterialAccessors< extrinsicMeshData::multifluid::phaseDensity, // 0
                                extrinsicMeshData::multifluid::dPhaseDensity_dPressure, // 1
                                extrinsicMeshData::multifluid::dPhaseDensity_dGlobalCompFraction, // 2
                                extrinsicMeshData::multifluid::phaseCompFraction, // 3
                                extrinsicMeshData::multifluid::dPhaseCompFraction_dPressure, // 4
                                extrinsicMeshData::multifluid::dPhaseCompFraction_dGlobalCompFraction >; // 5

  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
    compute( localIndex const numPhases,
             localIndex const ipWater,
             bool const allowAllPhasesIntoAquifer,
             real64 const & aquiferVolFlux,
             real64 const & dAquiferVolFlux_dPres,
             real64 const & aquiferWaterPhaseDens,
             arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
             arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens,
             arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > dPhaseDens_dPres,
             arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens_dCompFrac,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > dPhaseVolFrac_dPres,
             arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac_dCompDens,
             arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac,
             arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > dPhaseCompFrac_dPres,
             arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac_dCompFrac,
             arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens,
             real64 const & dt,
             real64 ( &localFlux )[NC],
             real64 ( &localFluxJacobian )[NC][NC+1] );

  template< localIndex NC >
  static void
  launch( localIndex const numPhases,
          localIndex const ipWater,
          bool const allowAllPhasesIntoAquifer,
          BoundaryStencil const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,

          AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
          real64 const & aquiferWaterPhaseDens,
          arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseVolFrac_dPres,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac,
          real64 const & timeAtBeginningOfStep,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};

} // namespace CompositionalMultiphaseFVMKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVMKERNELS_HPP
