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
 * @file DirichletFaceBasedAssemblyKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FVM_DIRICHLETFACEBASEDASSEMBLYKERNEL_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FVM_DIRICHLETFACEBASEDASSEMBLYKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/FVM/FaceBasedAssemblyKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FVM/FaceBasedAssemblyKernelBase.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FVM/MobilityKernel.hpp"


namespace geos
{

namespace singlePhaseFVMKernels
{

/******************************** DirichletFaceBasedAssemblyKernel ********************************/

/**
 * @class DirichFaceBasedAssemblyKernel
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms
 */
template< integer NUM_DOF, typename FLUIDWRAPPER >
class DirichletFaceBasedAssemblyKernel : public FaceBasedAssemblyKernel< NUM_DOF,
                                                                         geos::BoundaryStencilWrapper >
{
public:

  using AbstractBase = singlePhaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;
  using AbstractBase::m_mob;
  using AbstractBase::m_dMob_dPres;
  using AbstractBase::m_dens;
  using AbstractBase::m_dDens_dPres;
  using AbstractBase::m_permeability;
  using AbstractBase::m_dPerm_dPres;

  using AbstractBase::m_localMatrix;
  using AbstractBase::m_localRhs;

  using Base = singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_DOF,
                                                               geos::BoundaryStencilWrapper >;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] faceManager the face manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] dofNumberAccessor
   * @param[in] singlePhaseFlowAccessors
   * @param[in] singlePhaseFluidAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  DirichletFaceBasedAssemblyKernel( geos::globalIndex const rankOffset,
                                    geos::FaceManager const & faceManager,
                                    geos::BoundaryStencilWrapper const & stencilWrapper,
                                    FLUIDWRAPPER const & fluidWrapper,
                                    DofNumberAccessor const & dofNumberAccessor,
                                    SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                    SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                    PermeabilityAccessors const & permeabilityAccessors,
                                    geos::real64 const & dt,
                                    geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
                                    geos::arrayView1d< geos::real64 > const & localRhs )
    : Base( rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_facePres( faceManager.getField< geos::fields::flow::facePressure >() ),
    m_faceGravCoef( faceManager.getField< geos::fields::flow::gravityCoefficient >() ),
    m_fluidWrapper( fluidWrapper )
  { }

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
    StackVariables( geos::localIndex const GEOS_UNUSED_PARAM( size ),
                    geos::localIndex GEOS_UNUSED_PARAM( numElems ) )
    { }

    /// Transmissibility
    geos::real64 transmissibility = 0.0;

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    geos::globalIndex dofColIndices[numDof]{};

    /// Storage for the face local residual
    geos::real64 localFlux[numEqn]{};

    /// Storage for the face local Jacobian
    geos::real64 localFluxJacobian[numEqn][numDof]{};

  };


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( geos::localIndex const iconn,
              StackVariables & stack ) const
  {
    geos::globalIndex const offset =
      m_dofNumber[m_seri( iconn, geos::BoundaryStencil::Order::ELEM )][m_sesri( iconn, geos::BoundaryStencil::Order::ELEM )][m_sei( iconn, geos::BoundaryStencil::Order::ELEM )];

    for( geos::integer jdof = 0; jdof < numDof; ++jdof )
    {
      stack.dofColIndices[jdof] = offset + jdof;
    }
  }

  /**
   * @brief Compute the local Dirichlet face flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = geos::NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( geos::localIndex const iconn,
                    StackVariables & stack,
                    FUNC && compFluxKernelOp = geos::NoOpFunc{} ) const
  {
    using Order = geos::BoundaryStencil::Order;
    geos::localIndex constexpr numElems = geos::BoundaryStencil::maxNumPointsInFlux;

    geos::stackArray1d< geos::real64, numElems > mobility( numElems );
    geos::stackArray1d< geos::real64, numElems > dMobility_dP( numElems );

    geos::localIndex const er = m_seri( iconn, Order::ELEM );
    geos::localIndex const esr = m_sesri( iconn, Order::ELEM );
    geos::localIndex const ei = m_sei( iconn, Order::ELEM );
    geos::localIndex const kf = m_sei( iconn, Order::FACE );

    // Get flow quantities on the elem/face
    geos::real64 faceDens, faceVisc;
    m_fluidWrapper.compute( m_facePres[kf], faceDens, faceVisc );

    mobility[Order::ELEM] = m_mob[er][esr][ei];
    singlePhaseBaseKernels::MobilityKernel::compute( faceDens, faceVisc, mobility[Order::FACE] );

    dMobility_dP[Order::ELEM] = m_dMob_dPres[er][esr][ei];
    dMobility_dP[Order::FACE] = 0.0;

    // Compute average density
    geos::real64 const densMean = 0.5 * ( m_dens[er][esr][ei][0] + faceDens );
    geos::real64 const dDens_dP = 0.5 * m_dDens_dPres[er][esr][ei][0];

    // Evaluate potential difference
    geos::real64 const potDif = ( m_pres[er][esr][ei] - m_facePres[kf] )
                                - densMean * ( m_gravCoef[er][esr][ei] - m_faceGravCoef[kf] );
    geos::real64 const dPotDif_dP = 1.0 - dDens_dP * m_gravCoef[er][esr][ei];

    geos::real64 dTrans_dPerm[3];
    m_stencilWrapper.computeWeights( iconn, m_permeability, stack.transmissibility, dTrans_dPerm );
    geos::real64 const dTrans_dPres = LvArray::tensorOps::AiBi< 3 >( dTrans_dPerm, m_dPerm_dPres[er][esr][ei][0] );

    geos::real64 const f = stack.transmissibility * potDif;
    geos::real64 const dF_dP = stack.transmissibility * dPotDif_dP + dTrans_dPres * potDif;

    // Upwind mobility
    geos::localIndex const k_up = ( potDif >= 0 ) ? Order::ELEM : Order::FACE;
    geos::real64 const mobility_up = mobility[k_up];
    geos::real64 const dMobility_dP_up = dMobility_dP[k_up];

    // call the lambda in the phase loop to allow the reuse of the fluxes and their derivatives
    // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase

    compFluxKernelOp( er, esr, ei, kf, f, dF_dP, mobility_up, dMobility_dP_up );

    // *** end of upwinding

    // Populate local flux vector and derivatives

    stack.localFlux[0] = m_dt * mobility[k_up] * f;
    stack.localFluxJacobian[0][0] = m_dt * ( mobility_up * dF_dP + dMobility_dP_up * f );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = geos::NoOpFunc >
  GEOS_HOST_DEVICE
  void complete( geos::localIndex const iconn,
                 StackVariables & stack,
                 FUNC && assemblyKernelOp = geos::NoOpFunc{} ) const
  {
    using Order = geos::BoundaryStencil::Order;

    geos::localIndex const er = m_seri( iconn, Order::ELEM );
    geos::localIndex const esr = m_sesri( iconn, Order::ELEM );
    geos::localIndex const ei = m_sei( iconn, Order::ELEM );

    if( m_ghostRank[er][esr][ei] < 0 )
    {
      // Add to global residual/jacobian
      geos::globalIndex const dofIndex = m_dofNumber[er][esr][ei];
      geos::localIndex const localRow = LvArray::integerConversion< geos::localIndex >( dofIndex - m_rankOffset );

      RAJA::atomicAdd( geos::parallelDeviceAtomic{}, &m_localRhs[localRow], stack.localFlux[0] );

      m_localMatrix.template addToRowBinarySearchUnsorted< geos::parallelDeviceAtomic >
        ( localRow,
        stack.dofColIndices,
        stack.localFluxJacobian[0],
        numDof );

      assemblyKernelOp( localRow );
    }
  }

protected:

  /// Views on face pressure, temperature, and composition
  geos::arrayView1d< geos::real64 const > const m_facePres;

  /// View on the face gravity coefficient
  geos::arrayView1d< geos::real64 const > const m_faceGravCoef;

  /// Reference to the fluid wrapper
  FLUIDWRAPPER const m_fluidWrapper;

};

/**
 * @class DirichletFaceBasedAssemblyKernelFactory
 */
class DirichletFaceBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] faceManager reference to the face manager
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the boundary stencil wrapper
   * @param[in] fluidBase the single phase fluid constitutive model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( geos::globalIndex const rankOffset,
                   geos::string const & dofKey,
                   geos::string const & solverName,
                   geos::FaceManager const & faceManager,
                   geos::ElementRegionManager const & elemManager,
                   geos::BoundaryStencilWrapper const & stencilWrapper,
                   geos::constitutive::SingleFluidBase & fluidBase,
                   geos::real64 const & dt,
                   geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
                   geos::arrayView1d< geos::real64 > const & localRhs )
  {
    constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
    {
      using FluidType = TYPEOFREF( fluid );
      typename FluidType::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

      geos::integer constexpr NUM_DOF = 1;
      using kernelType = DirichletFaceBasedAssemblyKernel< NUM_DOF, typename FluidType::KernelWrapper >;

      geos::ElementRegionManager::ElementViewAccessor< geos::arrayView1d< geos::globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< geos::globalIndex, 1 >( dofKey );

      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      typename kernelType::SinglePhaseFlowAccessors singlePhaseFlowAccessors( elemManager, solverName );
      typename kernelType::SinglePhaseFluidAccessors singlePhaseFluidAccessors( elemManager, solverName );
      typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

      kernelType kernel( rankOffset,
                         faceManager,
                         stencilWrapper,
                         fluidWrapper,
                         dofNumberAccessor,
                         singlePhaseFlowAccessors,
                         singlePhaseFluidAccessors,
                         permeabilityAccessors,
                         dt,
                         localMatrix,
                         localRhs );

      kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }

};

}

}

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FVM_DIRICHLETFACEBASEDASSEMBLYKERNEL_HPP
