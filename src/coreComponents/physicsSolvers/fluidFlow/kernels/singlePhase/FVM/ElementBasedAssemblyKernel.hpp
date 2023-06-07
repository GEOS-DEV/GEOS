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

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_ELEMENTBASEDASSEMBLYKERNEL_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_ELEMENTBASEDASSEMBLYKERNEL_HPP

#include "physicsSolvers/fluidFlow/fields/FlowSolverBaseFields.hpp"


namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @brief Internal struct to provide no-op defaults used in the inclusion
 *   of lambda functions into kernel component functions.
 * @struct NoOpFunc
 */
struct NoOpFunc
{
  template< typename ... Ts >
  GEOS_HOST_DEVICE
  constexpr void
  operator()( Ts && ... ) const
  { }
};

/**
 * @class ElementBasedAssemblyKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */
template< typename SUBREGION_TYPE, integer NUM_DOF >
class ElementBasedAssemblyKernel
{

public:

  /// Compute time value for the number of degrees of freedom
  static constexpr geos::integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr geos::integer numEqn = NUM_DOF;

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
  ElementBasedAssemblyKernel( geos::globalIndex const rankOffset,
                              geos::string const dofKey,
                              SUBREGION_TYPE const & subRegion,
                              geos::constitutive::SingleFluidBase const & fluid,
                              geos::constitutive::CoupledSolidBase const & solid,
                              geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
                              geos::arrayView1d< geos::real64 > const & localRhs )
    :
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.template getReference< geos::array1d< geos::globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_deltaVolume( subRegion.template getField< geos::fields::flow::deltaVolume >() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_porosityNew( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    m_density_n( fluid.density_n() ),
    m_density( fluid.density() ),
    m_dDensity_dPres( fluid.dDensity_dPressure() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  { }

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    // Pore volume information

    /// Pore volume at time n+1
    geos::real64 poreVolume = 0.0;

    /// Pore volume at the previous converged time step
    geos::real64 poreVolume_n = 0.0;

    /// Derivative of pore volume with respect to pressure
    geos::real64 dPoreVolume_dPres = 0.0;

    // Residual information

    /// Index of the local row corresponding to this element
    geos::localIndex localRow = -1;

    /// Index of the matrix row/column corresponding to the dof in this element
    geos::globalIndex dofIndices[numDof]{};

    /// Storage for the element local residual vector
    geos::real64 localResidual[numEqn]{};

    /// Storage for the element local Jacobian matrix
    geos::real64 localJacobian[numEqn][numDof]{};

  };

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  geos::integer elemGhostRank( geos::localIndex const ei ) const
  { return m_elemGhostRank( ei ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( geos::localIndex const ei,
              StackVariables & stack ) const
  {
    // initialize the pore volume
    stack.poreVolume = ( m_volume[ei] + m_deltaVolume[ei] ) * m_porosityNew[ei][0];
    stack.poreVolume_n = m_volume[ei] * m_porosity_n[ei][0];
    stack.dPoreVolume_dPres = ( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dPres[ei][0];

    // set row index and degrees of freedom indices for this element
    stack.localRow = m_dofNumber[ei] - m_rankOffset;
    for( geos::integer idof = 0; idof < numDof; ++idof )
    {
      stack.dofIndices[idof] = m_dofNumber[ei] + idof;
    }
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeAccumulation( geos::localIndex const ei,
                            StackVariables & stack,
                            FUNC && kernelOp = NoOpFunc{} ) const
  {
    // Residual contribution is mass conservation in the cell
    stack.localResidual[0] = stack.poreVolume * m_density[ei][0] - stack.poreVolume_n * m_density_n[ei][0];

    // Derivative of residual wrt to pressure in the cell
    stack.localJacobian[0][0] = stack.dPoreVolume_dPres * m_density[ei][0] + m_dDensity_dPres[ei][0] * stack.poreVolume;

    // Customize the kernel with this lambda
    kernelOp();
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( geos::localIndex const GEOS_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {
    // add contribution to global residual and jacobian (no need for atomics here)
    m_localMatrix.template addToRow< geos::serialAtomic >( stack.localRow,
                                                           stack.dofIndices,
                                                           stack.localJacobian[0],
                                                           numDof );
    m_localRhs[stack.localRow] += stack.localResidual[0];

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
  launch( geos::localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    geos::forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( geos::localIndex const ei )
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
  geos::globalIndex const m_rankOffset;

  /// View on the dof numbers
  geos::arrayView1d< geos::globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  geos::arrayView1d< geos::integer const > const m_elemGhostRank;

  /// View on the element volumes
  geos::arrayView1d< geos::real64 const > const m_volume;
  geos::arrayView1d< geos::real64 const > const m_deltaVolume;

  /// Views on the porosity
  geos::arrayView2d< geos::real64 const > const m_porosity_n;
  geos::arrayView2d< geos::real64 const > const m_porosityNew;
  geos::arrayView2d< geos::real64 const > const m_dPoro_dPres;

  /// Views on density
  geos::arrayView2d< geos::real64 const > const m_density_n;
  geos::arrayView2d< geos::real64 const > const m_density;
  geos::arrayView2d< geos::real64 const > const m_dDensity_dPres;

  /// View on the local CRS matrix
  geos::CRSMatrixView< geos::real64, geos::globalIndex const > const m_localMatrix;
  /// View on the local RHS
  geos::arrayView1d< geos::real64 > const m_localRhs;

};

/**
 * @class SurfaceElementBasedAssemblyKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation in SurfaceElementSubRegion
 */
class SurfaceElementBasedAssemblyKernel : public ElementBasedAssemblyKernel< geos::SurfaceElementSubRegion, 1 >
{

public:

  using Base = ElementBasedAssemblyKernel< geos::SurfaceElementSubRegion, 1 >;

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
  SurfaceElementBasedAssemblyKernel( geos::globalIndex const rankOffset,
                                     geos::string const dofKey,
                                     geos::SurfaceElementSubRegion const & subRegion,
                                     geos::constitutive::SingleFluidBase const & fluid,
                                     geos::constitutive::CoupledSolidBase const & solid,
                                     geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
                                     geos::arrayView1d< geos::real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs )
#if ALLOW_CREATION_MASS
    ,
    m_creationMass( subRegion.getReference< geos::array1d< geos::real64 > >( geos::SurfaceElementSubRegion::viewKeyStruct::creationMassString() ) )
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
  void computeAccumulation( geos::localIndex const ei,
                            Base::StackVariables & stack ) const
  {
    Base::computeAccumulation( ei, stack, [&]()
    {
#if ALLOW_CREATION_MASS
      if( Base::m_volume[ei] * Base::m_density_n[ei][0] > 1.1 * m_creationMass[ei] )
      {
        stack.localResidual[0] += m_creationMass[ei] * 0.25;
      }
#endif
    } );
  }

protected:

#if ALLOW_CREATION_MASS
  geos::arrayView1d< geos::real64 const > const m_creationMass;
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
  createAndLaunch( geos::globalIndex const rankOffset,
                   geos::string const dofKey,
                   geos::CellElementSubRegion const & subRegion,
                   geos::constitutive::SingleFluidBase const & fluid,
                   geos::constitutive::CoupledSolidBase const & solid,
                   geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
                   geos::arrayView1d< geos::real64 > const & localRhs )
  {
    geos::integer constexpr NUM_DOF = 1;

    ElementBasedAssemblyKernel< geos::CellElementSubRegion, NUM_DOF >
    kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
    ElementBasedAssemblyKernel< geos::CellElementSubRegion, NUM_DOF >::template launch< POLICY >( subRegion.size(), kernel );
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
  createAndLaunch( geos::globalIndex const rankOffset,
                   geos::string const dofKey,
                   geos::SurfaceElementSubRegion const & subRegion,
                   geos::constitutive::SingleFluidBase const & fluid,
                   geos::constitutive::CoupledSolidBase const & solid,
                   geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
                   geos::arrayView1d< geos::real64 > const & localRhs )
  {
    SurfaceElementBasedAssemblyKernel
      kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
    SurfaceElementBasedAssemblyKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

}

}

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_ELEMENTBASEDASSEMBLYKERNEL_HPP
