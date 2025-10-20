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
 * @file DirichletMFDFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_DIRICHLETMFDFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_DIRICHLETMFDFLUXCOMPUTEKERNEL_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/cell/CellElementSubRegion.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/NodeManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"

#include "FluxComputeKernel.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

/******************************** DirichletMFDFluxComputeKernel ********************************/

/**
 * @class DirichletMFDFluxComputeKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @tparam IP_TYPE the mimetic inner product type (e.g., TPFAInnerProduct, BdVLMInnerProduct)
 * @tparam NF number of faces per element in the targeted subregion
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms using an MFD inner product
 */
template< integer NUM_COMP, integer NUM_DOF, typename FLUIDWRAPPER, typename IP_TYPE, integer NF >
class DirichletMFDFluxComputeKernel : public FluxComputeKernel< NUM_COMP,
                                                                NUM_DOF,
                                                                BoundaryStencilWrapper >
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

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_numPhases;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;
  using AbstractBase::m_phaseCompFrac;
  using AbstractBase::m_dPhaseCompFrac;
  using AbstractBase::m_dCompFrac_dCompDens;
  using AbstractBase::m_localMatrix;
  using AbstractBase::m_localRhs;
  using AbstractBase::m_kernelFlags;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernel< NUM_COMP, NUM_DOF, BoundaryStencilWrapper >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_stencilWrapper;
  using Base::m_phaseMob;
  using Base::m_dPhaseMob;
  using Base::m_phaseMassDens;
  using Base::m_dPhaseMassDens;
  using Base::m_permeability;
  using Base::m_dPerm_dPres;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] faceManager the face manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] nodePosition reference positions of nodes
   * @param[in] faceToNodes face-to-nodes connectivity
   * @param[in] elemCenter element centers (in the targeted subregion)
   * @param[in] elemVolume element volumes (in the targeted subregion)
   * @param[in] elemToFaces element-to-faces connectivity (in the targeted subregion)
   * @param[in] transMultiplier transmissibility multipliers per face
   * @param[in] lengthTolerance geometric tolerance for IP computation
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  DirichletMFDFluxComputeKernel( integer const numPhases,
                                 globalIndex const rankOffset,
                                 FaceManager const & faceManager,
                                 BoundaryStencilWrapper const & stencilWrapper,
                                 FLUIDWRAPPER const & fluidWrapper,
                                 DofNumberAccessor const & dofNumberAccessor,
                                 CompFlowAccessors const & compFlowAccessors,
                                 MultiFluidAccessors const & multiFluidAccessors,
                                 CapPressureAccessors const & capPressureAccessors,
                                 PermeabilityAccessors const & permeabilityAccessors,
                                 arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                 ArrayOfArraysView< localIndex const > const & faceToNodes,
                                 ElementViewConst< arrayView2d< real64 const > > const & elemCenter,
                                 ElementViewConst< arrayView1d< real64 const > > const & elemVolume,
                                 ElementViewConst< arrayView2d< localIndex const > > const & elemToFaces,
                                 arrayView1d< real64 const > const & transMultiplier,
                                 real64 const lengthTolerance,
                                 real64 const dt,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs,
                                 BitFlags< KernelFlags > kernelFlags )
    : Base( numPhases,
            rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            compFlowAccessors,
            multiFluidAccessors,
            capPressureAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs,
            kernelFlags ),
    m_facePres( faceManager.getField< fields::flow::facePressure >() ),
    m_faceTemp( faceManager.getField< fields::flow::faceTemperature >() ),
    m_faceCompFrac( faceManager.getField< fields::flow::faceGlobalCompFraction >() ),
    m_faceGravCoef( faceManager.getField< fields::flow::gravityCoefficient >() ),
    m_nodePosition( nodePosition ),
    m_faceToNodes( faceToNodes ),
    m_elemCenter( elemCenter ),
    m_elemVolume( elemVolume ),
    m_elemToFaces( elemToFaces ),
    m_transMultiplier( transMultiplier ),
    m_lengthTolerance( lengthTolerance ),
    m_fluidWrapper( fluidWrapper )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
  public:

    GEOS_HOST_DEVICE
    StackVariables( localIndex const GEOS_UNUSED_PARAM( size ),
                    localIndex GEOS_UNUSED_PARAM( numElems )) {}

    // Transmissibility
    real64 transmissibility = 0.0;

    // Component fluxes and derivatives

    /// Component fluxes
    real64 compFlux[numComp]{};
    /// Derivatives of component fluxes wrt pressure
    real64 dCompFlux_dP[numComp]{};
    /// Derivatives of component fluxes wrt component densities
    real64 dCompFlux_dC[numComp][numComp]{};

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    globalIndex dofColIndices[numDof]{};

    /// Storage for the face local residual vector
    real64 localFlux[numEqn]{};
    /// Storage for the face local Jacobian matrix
    real64 localFluxJacobian[numEqn][numDof]{};
  };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    globalIndex const offset =
      m_dofNumber[m_seri( iconn, BoundaryStencil::Order::ELEM )][m_sesri( iconn, BoundaryStencil::Order::ELEM )][m_sei( iconn, BoundaryStencil::Order::ELEM )];

    for( integer jdof = 0; jdof < numDof; ++jdof )
    {
      stack.dofColIndices[jdof] = offset + jdof;
    }
  }

  /**
   * @brief Compute the local Dirichlet face flux contributions to the residual and Jacobian (MFD version)
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;
    using Order = BoundaryStencil::Order;

    localIndex const er  = m_seri( iconn, Order::ELEM );
    localIndex const esr = m_sesri( iconn, Order::ELEM );
    localIndex const ei  = m_sei( iconn, Order::ELEM );
    localIndex const kf  = m_sei( iconn, Order::FACE );

    // Step 1: compute the one-sided MFD transmissibility at the boundary face
    {
      stackArray2d< real64, NF * NF > transMatrix( NF, NF );

      // principal permeabilities (assumed stored in first comp index)
      real64 const perm[3] = { m_permeability[er][esr][ei][0][0],
                               m_permeability[er][esr][ei][0][1],
                               m_permeability[er][esr][ei][0][2] };

      IP_TYPE::template compute< NF >( m_nodePosition,
                                       m_transMultiplier,
                                       m_faceToNodes,
                                       m_elemToFaces[er][esr][ei],
                                       m_elemCenter[er][esr][ei],
                                       m_elemVolume[er][esr][ei],
                                       perm,
                                       m_lengthTolerance,
                                       transMatrix );

      // find local face index corresponding to global face kf
      integer ifaceLoc = 0;
      for( integer j = 0; j < NF; ++j )
      {
        if( m_elemToFaces[er][esr][ei][j] == kf )
        {
          ifaceLoc = j;
          break;
        }
      }
      stack.transmissibility = transMatrix[ifaceLoc][ifaceLoc];
    }

    // Step 2: compute the fluid properties on the face
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseFrac( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseDens( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseMassDens( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseVisc( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseEnthalpy( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseInternalEnergy( 1, 1, m_numPhases );
    StackArray< real64, 4, constitutive::MultiFluidBase::MAX_NUM_PHASES * NUM_COMP,
                constitutive::multifluid::LAYOUT_PHASE_COMP > facePhaseCompFrac( 1, 1, m_numPhases, NUM_COMP );
    real64 faceTotalDens = 0.0;

    constitutive::MultiFluidBase::KernelWrapper::computeValues( m_fluidWrapper,
                                                                m_facePres[kf],
                                                                m_faceTemp[kf],
                                                                m_faceCompFrac[kf],
                                                                facePhaseFrac[0][0],
                                                                facePhaseDens[0][0],
                                                                facePhaseMassDens[0][0],
                                                                facePhaseVisc[0][0],
                                                                facePhaseEnthalpy[0][0],
                                                                facePhaseInternalEnergy[0][0],
                                                                facePhaseCompFrac[0][0],
                                                                faceTotalDens );

    // Step 3: loop over phases, compute and upwind phase flux and sum contributions to each component's flux
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      // working variables
      real64 dDensMean_dC[numComp]{};
      real64 dF_dC[numComp]{};
      real64 dProp_dC[numComp]{};

      real64 phaseFlux = 0.0; // for the lambda
      real64 dPhaseFlux_dP = 0.0;
      real64 dPhaseFlux_dC[numComp]{};

      // Step 3.1: compute the average phase mass density at the face
      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[er][esr][ei],
                      m_dPhaseMassDens[er][esr][ei][0][ip],
                      dProp_dC,
                      Deriv::dC );

      real64 const densMean = 0.5 * ( m_phaseMassDens[er][esr][ei][0][ip] + facePhaseMassDens[0][0][ip] );
      real64 const dDensMean_dP = 0.5 * m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dC[jc] = 0.5 * dProp_dC[jc];
      }

      // Step 3.2: compute the MFD potential difference at the face (same expression; only T differs)
      real64 const gravTimesDz = m_gravCoef[er][esr][ei] - m_faceGravCoef[kf];
      real64 const potDif = m_pres[er][esr][ei] - m_facePres[kf] - densMean * gravTimesDz;

      // Derivative of T wrt pressure unknown for general IP: set to 0 (approximation)
      real64 const dTrans_dPres = 0.0;
      real64 const f = stack.transmissibility * potDif;
      real64 const dF_dP = stack.transmissibility * ( 1.0 - dDensMean_dP * gravTimesDz ) + dTrans_dPres * potDif;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dF_dC[jc] = -stack.transmissibility * dDensMean_dC[jc] * gravTimesDz;
      }

      // Step 3.3: approximation for face mobility using total element mobility
      real64 const facePhaseMob = ( facePhaseFrac[0][0][ip] > 0.0 )
        ? facePhaseFrac[0][0][ip] * faceTotalDens / facePhaseVisc[0][0][ip]
        : 0.0;

      // Step 3.4: upwinding
      if( potDif >= 0 ) // the element is upstream
      {
        phaseFlux = m_phaseMob[er][esr][ei][ip] * f;
        dPhaseFlux_dP = m_phaseMob[er][esr][ei][ip] * dF_dP + m_dPhaseMob[er][esr][ei][ip][Deriv::dP] * f;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[jc] = m_phaseMob[er][esr][ei][ip] * dF_dC[jc] + m_dPhaseMob[er][esr][ei][ip][Deriv::dC+jc] * f;
        }

        arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE_COMP-3 > phaseCompFracSub =
          m_phaseCompFrac[er][esr][ei][0][ip];
        arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC-3 > dPhaseCompFracSub =
          m_dPhaseCompFrac[er][esr][ei][0][ip];

        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const ycp = phaseCompFracSub[ic];
          stack.compFlux[ic] += phaseFlux * ycp;
          stack.dCompFlux_dP[ic] += dPhaseFlux_dP * ycp + phaseFlux * dPhaseCompFracSub[ic][Deriv::dP];

          applyChainRule( numComp,
                          m_dCompFrac_dCompDens[er][esr][ei],
                          dPhaseCompFracSub[ic],
                          dProp_dC,
                          Deriv::dC );
          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dCompFlux_dC[ic][jc] += dPhaseFlux_dC[jc] * ycp + phaseFlux * dProp_dC[jc];
          }
        }
      }
      else // the face is upstream
      {
        phaseFlux = facePhaseMob * f;
        dPhaseFlux_dP = facePhaseMob * dF_dP;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[jc] = facePhaseMob * dF_dC[jc];
        }

        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const ycp = facePhaseCompFrac[0][0][ip][ic];
          stack.compFlux[ic] += phaseFlux * ycp;
          stack.dCompFlux_dP[ic] += dPhaseFlux_dP * ycp;
          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dCompFlux_dC[ic][jc] += dPhaseFlux_dC[jc] * ycp;
          }
        }
      }

      // reuse hook
      compFluxKernelOp( ip, er, esr, ei, kf, f,
                        facePhaseMob, facePhaseEnthalpy[0][0], facePhaseCompFrac[0][0],
                        phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );
    }

    // Step 4: populate local flux vector and derivatives
    for( integer ic = 0; ic < numComp; ++ic )
    {
      stack.localFlux[ic]            = m_dt * stack.compFlux[ic];
      stack.localFluxJacobian[ic][0] = m_dt * stack.dCompFlux_dP[ic];
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localFluxJacobian[ic][jc+1] = m_dt * stack.dCompFlux_dC[ic][jc];
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && assemblyKernelOp = NoOpFunc{} ) const
  {
    using namespace compositionalMultiphaseUtilities;
    using Order = BoundaryStencil::Order;

    if( AbstractBase::m_kernelFlags.isSet( KernelFlags::TotalMassEquation ) )
    {
      // Apply equation/variable change transformation(s)
      real64 work[numDof]{};
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numDof, stack.localFluxJacobian, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( numComp, stack.localFlux );
    }

    if( m_ghostRank[m_seri( iconn, Order::ELEM )][m_sesri( iconn, Order::ELEM )][m_sei( iconn, Order::ELEM )] < 0 )
    {
      globalIndex const globalRow = m_dofNumber[m_seri( iconn, Order::ELEM )][m_sesri( iconn, Order::ELEM )][m_sei( iconn, Order::ELEM )];
      localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
      GEOS_ASSERT_GE( localRow, 0 );
      GEOS_ASSERT_GT( AbstractBase::m_localMatrix.numRows(), localRow + numComp );

      for( integer ic = 0; ic < numComp; ++ic )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + ic], stack.localFlux[ic] );
        AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
          ( localRow + ic,
            stack.dofColIndices,
            stack.localFluxJacobian[ic],
            numDof );
      }

      assemblyKernelOp( localRow );
    }
  }

protected:

  /// Views on face pressure, temperature, and composition
  arrayView1d< real64 const > const m_facePres;
  arrayView1d< real64 const > const m_faceTemp;
  arrayView2d< real64 const, compflow::USD_COMP > const m_faceCompFrac;

  /// View on the face gravity coefficient
  arrayView1d< real64 const > const m_faceGravCoef;

  // Geometry/connectivity for IP computation
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodePosition;
  ArrayOfArraysView< localIndex const > const m_faceToNodes;
  ElementViewConst< arrayView2d< real64 const > > const m_elemCenter;
  ElementViewConst< arrayView1d< real64 const > > const m_elemVolume;
  ElementViewConst< arrayView2d< localIndex const > > const m_elemToFaces;
  arrayView1d< real64 const > const m_transMultiplier;
  real64 const m_lengthTolerance;

  /// Reference to the fluid wrapper
  FLUIDWRAPPER const m_fluidWrapper;
};


/**
 * @class DirichletMFDFluxComputeKernelFactory
 */
class DirichletMFDFluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new MFD Dirichlet kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam IP_TYPE mimetic inner product type
   * @tparam NF number of faces per element in the targeted subregion
   */
  template< typename POLICY, typename IP_TYPE, integer NF >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   BitFlags< KernelFlags > kernelFlags,
                   string const & dofKey,
                   string const & solverName,
                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   BoundaryStencilWrapper const & stencilWrapper,
                   constitutive::MultiFluidBase & fluidBase,
                   real64 const lengthTolerance,
                   real64 const dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    // Select fluid wrapper based on number of components
    constitutive::constitutiveComponentUpdatePassThru( fluidBase, numComps, [&]( auto & fluid, auto NC )
    {
      using FluidType = TYPEOFREF( fluid );
      typename FluidType::KernelWrapper const fluidWrapper = fluid.createKernelWrapper();

      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 1;

      // Identify a representative subregion from the stencil (same assumption as TPFA Dirichlet)
      localIndex const er = stencilWrapper.getElementRegionIndices()( 0, BoundaryStencil::Order::ELEM );
      localIndex const esr = stencilWrapper.getElementSubRegionIndices()( 0, BoundaryStencil::Order::ELEM );
      ElementSubRegionBase const & subRegionBase = elemManager.getRegion( er ).getSubRegion( esr );

      // Element geometry/connectivity arrays
      ElementRegionManager::ElementViewConst< arrayView2d< real64 const > > elemCenter =
        elemManager.constructArrayViewAccessor< real64, 2 >( CellElementSubRegion::viewKeyStruct::elementCenterString() ).toViewConst();
      ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > elemVolume =
        elemManager.constructArrayViewAccessor< real64, 1 >( CellElementSubRegion::viewKeyStruct::elementVolumeString() ).toViewConst();
      ElementRegionManager::ElementViewConst< arrayView2d< localIndex const > > elemToFaces =
        elemManager.constructArrayViewAccessor< localIndex, 2 >( CellElementSubRegion::viewKeyStruct::faceListString() ).toViewConst();

      // Face connectivity and multipliers
      ArrayOfArraysView< localIndex const > faceToNodes = faceManager.nodeList().toViewConst();
      arrayView1d< real64 const > const & transMultiplier = faceManager.getField< fields::flow::transMultiplier >();

      using kernelType = DirichletMFDFluxComputeKernel< NUM_COMP, NUM_DOF, typename FluidType::KernelWrapper, IP_TYPE, NF >;
      typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

      kernelType kernel( numPhases,
                         rankOffset,
                         faceManager,
                         stencilWrapper,
                         fluidWrapper,
                         elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey ),
                         compFlowAccessors,
                         multiFluidAccessors,
                         capPressureAccessors,
                         permeabilityAccessors,
                         nodePosition,
                         faceToNodes,
                         elemCenter,
                         elemVolume,
                         elemToFaces,
                         transMultiplier,
                         lengthTolerance,
                         dt,
                         localMatrix,
                         localRhs,
                         kernelFlags );
      kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
      GEOS_UNUSED_VAR( subRegionBase );
    } );
  }
};

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos


#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_DIRICHLETMFDFLUXCOMPUTEKERNEL_HPP
