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
 * @file SinglePhaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/fluid/SingleFluidExtrinsicData.hpp"
#include "constitutive/fluid/SlurryFluidBase.hpp"
#include "constitutive/fluid/SlurryFluidExtrinsicData.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geosx
{

namespace singlePhaseFVMKernels
{
using namespace constitutive;

using namespace fluxKernelsHelper;

/******************************** FaceBasedAssemblyKernelBase ********************************/

/**
 * @brief Base class for FaceBasedAssemblyKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of dofs).
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

  using SinglePhaseFlowAccessors =
    StencilAccessors< extrinsicMeshData::ghostRank,
                      extrinsicMeshData::flow::pressure,
                      extrinsicMeshData::flow::pressure_n,
                      extrinsicMeshData::flow::gravityCoefficient,
                      extrinsicMeshData::flow::mobility,
                      extrinsicMeshData::flow::dMobility_dPressure >;

  using SinglePhaseFluidAccessors =
    StencilMaterialAccessors< SingleFluidBase,
                              extrinsicMeshData::singlefluid::density,
                              extrinsicMeshData::singlefluid::dDensity_dPressure >;

  using SlurryFluidAccessors =
    StencilMaterialAccessors< SlurryFluidBase,
                              extrinsicMeshData::singlefluid::density,
                              extrinsicMeshData::singlefluid::dDensity_dPressure >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              extrinsicMeshData::permeability::permeability,
                              extrinsicMeshData::permeability::dPerm_dPressure >;

  using ProppantPermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              extrinsicMeshData::permeability::permeability,
                              extrinsicMeshData::permeability::dPerm_dPressure,
                              extrinsicMeshData::permeability::dPerm_dDispJump,
                              extrinsicMeshData::permeability::permeabilityMultiplier >;
                              
  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor accessor for the dof numbers
   * @param[in] singleFlowAccessors accessor for wrappers registered by the solver
   * @param[in] singlePhaseFluidAccessors accessor for wrappers registered by the singlefluid model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernelBase( globalIndex const rankOffset,
                               DofNumberAccessor const & dofNumberAccessor,
                               SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                               SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                               PermeabilityAccessors const & permeabilityAccessors,
                               real64 const & dt,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
    : m_rankOffset( rankOffset ),
    m_dt( dt ),
    m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_permeability( permeabilityAccessors.get( extrinsicMeshData::permeability::permeability {} ) ),
    m_dPerm_dPres( permeabilityAccessors.get( extrinsicMeshData::permeability::dPerm_dPressure {} ) ),
    m_ghostRank( singlePhaseFlowAccessors.get( extrinsicMeshData::ghostRank {} ) ),
    m_gravCoef( singlePhaseFlowAccessors.get( extrinsicMeshData::flow::gravityCoefficient {} ) ),
    m_pres( singlePhaseFlowAccessors.get( extrinsicMeshData::flow::pressure {} ) ),
    m_mob( singlePhaseFlowAccessors.get( extrinsicMeshData::flow::mobility {} ) ),
    m_dMob_dPres( singlePhaseFlowAccessors.get( extrinsicMeshData::flow::dMobility_dPressure {} ) ),
    m_dens( singlePhaseFluidAccessors.get( extrinsicMeshData::singlefluid::density {} ) ),
    m_dDens_dPres( singlePhaseFluidAccessors.get( extrinsicMeshData::singlefluid::dDensity_dPressure {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}    

protected:             

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
  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres;

  /// Views on fluid mobility 
  ElementViewConst< arrayView1d< real64 const > > const m_mob; 
  ElementViewConst< arrayView1d< real64 const > > const m_dMob_dPres;

  /// Views on fluid density
  ElementViewConst< arrayView2d< real64 const > > const m_dens; 
  ElementViewConst< arrayView2d< real64 const > > const m_dDens_dPres; 

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;
}; 

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public FaceBasedAssemblyKernelBase
{
public:

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_DOF;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::maxNumPointsInFlux;

  /// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::maxNumConnections;

  /// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::maxStencilSize;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] singlePhaseFlowAccessors
   * @param[in] singlePhaseFluidAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( globalIndex const rankOffset,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                           SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    : FaceBasedAssemblyKernelBase( rankOffset,
                                   dofNumberAccessor,
                                   singlePhaseFlowAccessors,
                                   singlePhaseFluidAccessors,
                                   permeabilityAccessors,
                                   dt,
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
  GEOSX_HOST_DEVICE
  localIndex stencilSize( localIndex const iconn ) const
  { return m_sei[iconn].size(); }

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
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] NoOpFunc the function used to customize the computation of the flux
   */
  template< typename FUNC = singlePhaseBaseKernels::NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = singlePhaseBaseKernels::NoOpFunc{} ) const
  {
    // first, compute the transmissibilities at this face
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    // clear working arrays
    real64 densMean = 0.0;
    stackArray1d< real64, maxNumElems > dDensMean_dP( stack.numFluxElems );

    // create local work arrays
    real64 fluxVal = 0.0; 
    real64 dFlux_dP[maxStencilSize]{}; 

    real64 presGrad = 0.0; 
    stackArray1d< real64, maxStencilSize > dPresGrad_dP( stack.stencilSize ); 

    real64 gravHead = 0.0; 
    stackArray1d< real64, maxNumElems > dGravHead_dP( stack.numFluxElems ); 

    // calculate quantities on primary connected cells
    for( integer ke = 0; ke < stack.numFluxElems; ++ke )
    {
      localIndex const er  = m_seri( iconn, ke );
      localIndex const esr = m_sesri( iconn, ke );
      localIndex const ei  = m_sei( iconn, ke );

      // density 
      real64 const density  = m_dens[er][esr][ei][0]; 
      real64 const dDens_dP = m_dDens_dPres[er][esr][ei][0]; 

      // average density and derivatives 
      densMean += 0.5 * density; 
      dDensMean_dP[ke] = 0.5 * dDens_dP; 
    }

    //***** calculation of flux *****

    // compute potential difference
    for ( integer i = 0; i < stack.stencilSize; ++i )
    {
      localIndex const er  = m_seri( iconn, i );
      localIndex const esr = m_sesri( iconn, i );
      localIndex const ei  = m_sei( iconn, i );

      presGrad += stack.transmissibility[0][i] * m_pres[er][esr][ei]; 
      dPresGrad_dP[i] += stack.transmissibility[0][i] + stack.dTrans_dPres[0][i] * m_pres[er][esr][ei]; 

      real64 const gravD     = stack.transmissibility[0][i] * m_gravCoef[er][esr][ei]; 
      real64 const dGravD_dP = stack.dTrans_dPres[0][i] * m_gravCoef[er][esr][ei]; 

      gravHead += densMean * gravD;

      for ( integer ke = 0; ke < stack.numFluxElems; ++ke )
      {
        dGravHead_dP[ke] += dDensMean_dP[ke] * gravD + dGravD_dP * densMean; 
      }
    }

    // *** upwinding ***

    // compute potential gradient
    real64 const potGrad = presGrad - gravHead;

    // choose upstream cell 
    localIndex const k_up = (potGrad >= 0) ? 0 : 1;

    localIndex const er_up  = m_seri( iconn, k_up );
    localIndex const esr_up = m_sesri( iconn, k_up );
    localIndex const ei_up  = m_sei( iconn, k_up ); 

    real64 const mobility = m_mob[er_up][esr_up][ei_up]; 
    real64 const dMob_dP = m_dMob_dPres[er_up][esr_up][ei_up]; 

    // pressure gradient depends on all points in the stencil
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      dFlux_dP[i] += dPresGrad_dP[i];
    }

    // gravitational head depends only on the two cells connected (same as mean density)
    for( integer ke = 0; ke < stack.numFluxElems; ++ke )
    {
      dFlux_dP[ke] -= dGravHead_dP[ke];
    }

    // compute the flux and derivatives using upstream cell mobility
    fluxVal = mobility * potGrad;

    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      dFlux_dP[i] *= mobility;  
    }

    // add contribution from upstream cell mobility derivatives
    dFlux_dP[k_up] += dMob_dP * potGrad;

    // populate local flux vector and derivatives
    stack.localFlux[0]      += m_dt * fluxVal; 
    stack.localFlux[numEqn] -= m_dt * fluxVal; 

    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      localIndex const localDofIndexPres = i * numDof; 
      stack.localFluxJacobian[0][localDofIndexPres]      += m_dt * dFlux_dP[i]; 
      stack.localFluxJacobian[numEqn][localDofIndexPres] -= m_dt * dFlux_dP[i]; 
    }

    // Customize the kernel with this lambda
    kernelOp( k_up, er_up, esr_up, ei_up, potGrad, fluxVal, dFlux_dP );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = singlePhaseBaseKernels::NoOpFunc >
  GEOSX_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && kernelOp = singlePhaseBaseKernels::NoOpFunc{} ) const
  {
    // add contribution to residual and jacobian into:
    // - the mass balance equation
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    for( integer i = 0; i < stack.numFluxElems; ++i )
    {
      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( m_localMatrix.numRows(), localRow );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow], stack.localFlux[i * numEqn] );
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                            stack.dofColIndices.data(),
                                                                            stack.localFluxJacobian[i * numEqn].dataIfContiguous(),
                                                                            stack.stencilSize * numDof );

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


protected: 

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
  createAndLaunch( globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 1;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = FaceBasedAssemblyKernel< NUM_DOF, STENCILWRAPPER >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );

    kernelType kernel( rankOffset, stencilWrapper, dofNumberAccessor,
                       flowAccessors, fluidAccessors, permAccessors,
                       dt, localMatrix, localRhs );
    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

/******************************** FluxKernel ********************************/

// To delete it after verifying FaceBasedAssemblyKernel
struct FluxKernel
{
  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using SinglePhaseFlowAccessors =
    StencilAccessors< extrinsicMeshData::ghostRank,
                      extrinsicMeshData::flow::pressure,
                      extrinsicMeshData::flow::pressure_n,
                      extrinsicMeshData::flow::gravityCoefficient,
                      extrinsicMeshData::flow::mobility,
                      extrinsicMeshData::flow::dMobility_dPressure >;

  using SinglePhaseFluidAccessors =
    StencilMaterialAccessors< SingleFluidBase,
                              extrinsicMeshData::singlefluid::density,
                              extrinsicMeshData::singlefluid::dDensity_dPressure >;

  using SlurryFluidAccessors =
    StencilMaterialAccessors< SlurryFluidBase,
                              extrinsicMeshData::singlefluid::density,
                              extrinsicMeshData::singlefluid::dDensity_dPressure >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              extrinsicMeshData::permeability::permeability,
                              extrinsicMeshData::permeability::dPerm_dPressure >;

  using ProppantPermeabilityAccessors =
    StencilMaterialAccessors< PermeabilityBase,
                              extrinsicMeshData::permeability::permeability,
                              extrinsicMeshData::permeability::dPerm_dPressure,
                              extrinsicMeshData::permeability::dPerm_dDispJump,
                              extrinsicMeshData::permeability::permeabilityMultiplier >;


  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
   * @param[in] stencil The stencil object.
   * @param[in] dt The timestep for the integration step.
   * @param[in] dofNumber The dofNumbers for each element
   * @param[in] pres The pressures in each element
   * @param[in] gravCoef The factor for gravity calculations (g*H)
   * @param[in] dens The material density in each element
   * @param[in] dDens_dPres The change in material density for each element
   * @param[in] mob The fluid mobility in each element
   * @param[in] dMob_dPres The derivative of mobility wrt pressure in each element
   * @param[in] permeability
   * @param[in] dPerm_dPres The derivative of permeability wrt pressure in each element
   * @param[out] localMatrix The linear system matrix
   * @param[out] localRhs The linear system residual
   */
  template< typename STENCILWRAPPER_TYPE >
  static void
  launch( STENCILWRAPPER_TYPE const & stencilWrapper,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

    constexpr localIndex maxNumElems = STENCILWRAPPER_TYPE::maxNumPointsInFlux;
    constexpr localIndex maxStencilSize = STENCILWRAPPER_TYPE::maxStencilSize;

    forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [stencilWrapper, dt, rankOffset, dofNumber, ghostRank,
                                                              pres, gravCoef, dens, dDens_dPres, mob,
                                                              dMob_dPres, permeability, dPerm_dPres,
                                                              seri, sesri, sei, localMatrix, localRhs] GEOSX_HOST_DEVICE ( localIndex const iconn )
    {
      localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
      localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );

      // working arrays
      stackArray1d< globalIndex, maxNumElems > dofColIndices( stencilSize );
      stackArray1d< real64, maxNumElems > localFlux( numFluxElems );
      stackArray2d< real64, maxNumElems * maxStencilSize > localFluxJacobian( numFluxElems, stencilSize );


      // compute transmissibility
      real64 transmissibility[STENCILWRAPPER_TYPE::maxNumConnections][2];
      real64 dTrans_dPres[STENCILWRAPPER_TYPE::maxNumConnections][2];

      stencilWrapper.computeWeights( iconn,
                                     permeability,
                                     dPerm_dPres,
                                     transmissibility,
                                     dTrans_dPres );

      compute( numFluxElems,
               seri[iconn],
               sesri[iconn],
               sei[iconn],
               transmissibility,
               dTrans_dPres,
               pres,
               gravCoef,
               dens,
               dDens_dPres,
               mob,
               dMob_dPres,
               dt,
               localFlux,
               localFluxJacobian );


      // extract DOF numbers
      for( localIndex i = 0; i < stencilSize; ++i )
      {
        dofColIndices[i] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];

      }

      for( localIndex i = 0; i < numFluxElems; ++i )
      {

        if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
        {
          globalIndex const globalRow = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GT( localMatrix.numRows(), localRow );

          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux[i] );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                            dofColIndices.data(),
                                                                            localFluxJacobian[i].dataIfContiguous(),
                                                                            stencilSize );

        }
      }

    } );
  }

  /**
   * @brief Compute flux and its derivatives for a given tpfa connector.
   *
   *
   */
  template< localIndex maxNumConnections >
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const numFluxElems,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           real64 const (&transmissibility)[maxNumConnections][2],
           real64 const (&dTrans_dPres)[maxNumConnections][2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & dens,
           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
           ElementViewConst< arrayView1d< real64 const > > const & mob,
           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian )
  {

    localIndex k[2];
    localIndex connectionIndex = 0;;
    for( k[0]=0; k[0]<numFluxElems; ++k[0] )
    {
      for( k[1]=k[0]+1; k[1]<numFluxElems; ++k[1] )
      {
        real64 fluxVal = 0.0;
        real64 dFlux_dTrans = 0.0;
        real64 const trans[2] = {transmissibility[connectionIndex][0], transmissibility[connectionIndex][1]};
        real64 const dTrans[2] = { dTrans_dPres[connectionIndex][0], dTrans_dPres[connectionIndex][1] };
        real64 dFlux_dP[2] = {0.0, 0.0};
        localIndex const regionIndex[2]    = {seri[k[0]], seri[k[1]]};
        localIndex const subRegionIndex[2] = {sesri[k[0]], sesri[k[1]]};
        localIndex const elementIndex[2]   = {sei[k[0]], sei[k[1]]};


        computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
                                trans,
                                dTrans,
                                pres,
                                gravCoef,
                                dens,
                                dDens_dPres,
                                mob,
                                dMob_dPres,
                                fluxVal,
                                dFlux_dP,
                                dFlux_dTrans );

        // populate local flux vector and derivatives
        flux[k[0]] +=  dt * fluxVal;
        flux[k[1]] -=  dt * fluxVal;

        fluxJacobian[k[0]][k[0]] += dt * dFlux_dP[0];
        fluxJacobian[k[0]][k[1]] += dt * dFlux_dP[1];
        fluxJacobian[k[1]][k[0]] -= dt * dFlux_dP[0];
        fluxJacobian[k[1]][k[1]] -= dt * dFlux_dP[1];

        connectionIndex++;
      }
    }
  }
};

struct FaceDirichletBCKernel
{
  template< typename VIEWTYPE >
  using ElementViewConst = FluxKernel::ElementViewConst< VIEWTYPE >;

  template< typename FLUID_WRAPPER >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void compute( arraySlice1d< localIndex const > const & seri,
                       arraySlice1d< localIndex const > const & sesri,
                       arraySlice1d< localIndex const > const & sefi,
                       real64 const trans,
                       real64 const dTrans_dPres,
                       ElementViewConst< arrayView1d< real64 const > > const & pres,
                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                       ElementViewConst< arrayView2d< real64 const > > const & dens,
                       ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                       ElementViewConst< arrayView1d< real64 const > > const & mob,
                       ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                       arrayView1d< real64 const > const & presFace,
                       arrayView1d< real64 const > const & gravCoefFace,
                       FLUID_WRAPPER const & fluidWrapper,
                       real64 const dt,
                       real64 & flux,
                       real64 & dFlux_dP )
  {
    using Order = BoundaryStencil::Order;
    localIndex constexpr numElems = BoundaryStencil::maxNumPointsInFlux;

    stackArray1d< real64, numElems > mobility( numElems );
    stackArray1d< real64, numElems > dMobility_dP( numElems );

    localIndex const er  = seri[ Order::ELEM ];
    localIndex const esr = sesri[ Order::ELEM ];
    localIndex const ei  = sefi[ Order::ELEM ];
    localIndex const kf  = sefi[ Order::FACE ];

    // Get flow quantities on the elem/face
    real64 faceDens, faceVisc;
    fluidWrapper.compute( presFace[kf], faceDens, faceVisc );

    mobility[Order::ELEM] = mob[er][esr][ei];
    singlePhaseBaseKernels::MobilityKernel::compute( faceDens, faceVisc, mobility[Order::FACE] );

    dMobility_dP[Order::ELEM] = dMob_dPres[er][esr][ei];
    dMobility_dP[Order::FACE] = 0.0;

    // Compute average density
    real64 const densMean = 0.5 * ( dens[er][esr][ei][0] + faceDens );
    real64 const dDens_dP = 0.5 * dDens_dPres[er][esr][ei][0];

    // Evaluate potential difference
    real64 const potDif = (pres[er][esr][ei] - presFace[kf])
                          - densMean * (gravCoef[er][esr][ei] - gravCoefFace[kf]);
    real64 const dPotDif_dP = 1.0 - dDens_dP * gravCoef[er][esr][ei];

    real64 const f = trans * potDif;
    real64 const dF_dP = trans * dPotDif_dP + dTrans_dPres * potDif;

    // Upwind mobility
    localIndex const k_up = ( potDif >= 0 ) ? Order::ELEM : Order::FACE;
    flux = dt * mobility[k_up] * f;
    dFlux_dP = dt * ( mobility[k_up] * dF_dP + dMobility_dP[k_up] * f );
  }

  template< typename FLUID_WRAPPER >
  static void launch( BoundaryStencilWrapper const & stencil,
                      ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                      ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
                      globalIndex const rankOffset,
                      ElementViewConst< arrayView3d< real64 const > > const & perm,
                      ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                      ElementViewConst< arrayView1d< real64 const > > const & pres,
                      ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                      ElementViewConst< arrayView2d< real64 const > > const & dens,
                      ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                      ElementViewConst< arrayView1d< real64 const > > const & mob,
                      ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                      arrayView1d< real64 const > const & presFace,
                      arrayView1d< real64 const > const & gravCoefFace,
                      FLUID_WRAPPER const & fluidWrapper,
                      real64 const dt,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs )
  {
    BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();

    forAll< parallelDevicePolicy<> >( seri.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
    {
      localIndex const er  = seri( iconn, BoundaryStencil::Order::ELEM );
      localIndex const esr = sesri( iconn, BoundaryStencil::Order::ELEM );
      localIndex const ei  = sefi( iconn, BoundaryStencil::Order::ELEM );

      real64 trans;
      real64 dTrans_dPerm[3];
      stencil.computeWeights( iconn, perm, trans, dTrans_dPerm );
      real64 const dTrans_dPres = LvArray::tensorOps::AiBi< 3 >( dTrans_dPerm, dPerm_dPres[er][esr][ei][0] );

      real64 flux, fluxJacobian;

      compute( seri[iconn],
               sesri[iconn],
               sefi[iconn],
               trans,
               dTrans_dPres,
               pres,
               gravCoef,
               dens,
               dDens_dPres,
               mob,
               dMob_dPres,
               presFace,
               gravCoefFace,
               fluidWrapper,
               dt,
               flux,
               fluxJacobian );

      if( ghostRank[er][esr][ei] < 0 )
      {
        // Add to global residual/jacobian
        globalIndex const dofIndex = dofNumber[er][esr][ei];
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofIndex - rankOffset );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], flux );
        localMatrix.addToRow< parallelDeviceAtomic >( localRow, &dofIndex, &fluxJacobian, 1 );
      }
    } );

  }
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

  GEOSX_HOST_DEVICE
  static void
  compute( real64 const & aquiferVolFlux,
           real64 const & dAquiferVolFlux_dPres,
           real64 const & aquiferDens,
           real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & dt,
           real64 & localFlux,
           real64 & localFluxJacobian )
  {
    if( aquiferVolFlux > 0 ) // aquifer is upstream
    {
      localFlux -= dt * aquiferVolFlux * aquiferDens;
      localFluxJacobian -= dt * dAquiferVolFlux_dPres * aquiferDens;
    }
    else // reservoir is upstream
    {
      localFlux -= dt * aquiferVolFlux * dens;
      localFluxJacobian -= dt * (dAquiferVolFlux_dPres * dens + aquiferVolFlux * dDens_dPres);
    }
  }

  static void
  launch( BoundaryStencil const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
          real64 const & aquiferDens,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & pres_n,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          real64 const & timeAtBeginningOfStep,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    using Order = BoundaryStencil::Order;

    BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
    BoundaryStencil::WeightContainerViewConstType const & weight = stencil.getWeights();

    forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
    {

      // working variables
      real64 localFlux = 0.0;
      real64 localFluxJacobian = 0.0;

      localIndex const er  = seri( iconn, Order::ELEM );
      localIndex const esr = sesri( iconn, Order::ELEM );
      localIndex const ei  = sefi( iconn, Order::ELEM );
      real64 const areaFraction = weight( iconn, Order::ELEM );

      // compute the aquifer influx rate using the pressure influence function and the aquifer props
      real64 dAquiferVolFlux_dPres = 0.0;
      real64 const aquiferVolFlux = aquiferBCWrapper.compute( timeAtBeginningOfStep,
                                                              dt,
                                                              pres[er][esr][ei],
                                                              pres_n[er][esr][ei],
                                                              gravCoef[er][esr][ei],
                                                              areaFraction,
                                                              dAquiferVolFlux_dPres );

      // compute the phase/component aquifer flux
      AquiferBCKernel::compute( aquiferVolFlux,
                                dAquiferVolFlux_dPres,
                                aquiferDens,
                                dens[er][esr][ei][0],
                                dDens_dPres[er][esr][ei][0],
                                dt,
                                localFlux,
                                localFluxJacobian );

      // Add to residual/jacobian
      if( ghostRank[er][esr][ei] < 0 )
      {
        globalIndex const globalRow = dofNumber[er][esr][ei];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( localMatrix.numRows(), localRow );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux );
        localMatrix.addToRow< parallelDeviceAtomic >( localRow,
                                                      &dofNumber[er][esr][ei],
                                                      &localFluxJacobian,
                                                      1 );
      }
    } );
  }

};


} // namespace singlePhaseFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
