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
 * @file TeamLaplaceFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_TEAMLAPLACEFEMKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_TEAMLAPLACEFEMKERNELS_HPP_

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "finiteElement/kernelInterface/TeamKernelBase.hpp"
#include "finiteElement/kernelInterface/TeamKernelFunctions.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

//*****************************************************************************
/**
 * @brief Implements kernels for solving Laplace's equation.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### TeamLaplaceFEMKernel Description
 * Implements the KernelBase interface functions required for solving Laplace's
 * equation using on of the finite element kernel application functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the TeamKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class TeamLaplaceFEMKernel :
  public finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                        CONSTITUTIVE_TYPE,
                                        FE_TYPE,
                                        1,
                                        1 >
{
public:
  /// An alias for the base class.
  using Base = finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                              CONSTITUTIVE_TYPE,
                                              FE_TYPE,
                                              1,
                                              1 >;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::TeamKernelBase::TeamKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  TeamLaplaceFEMKernel( NodeManager const & nodeManager,
                        EdgeManager const & edgeManager,
                        FaceManager const & faceManager,
                        localIndex const targetRegionIndex,
                        SUBREGION_TYPE const & elementSubRegion,
                        FE_TYPE const & finiteElementSpace,
                        CONSTITUTIVE_TYPE & inputConstitutiveType,
                        arrayView1d< globalIndex const > const inputDofNumber,
                        globalIndex const rankOffset,
                        CRSMatrixView< real64, globalIndex const > const inputMatrix,
                        arrayView1d< real64 > const inputRhs,
                        string const fieldName ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition() ),
    m_primaryField( nodeManager.template getReference< array1d< real64 > >( fieldName )),
    m_rhs( inputRhs )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
    GEOSX_UNUSED_VAR( inputDofNumber );
    GEOSX_UNUSED_VAR( rankOffset );
    GEOSX_UNUSED_VAR( inputMatrix );
    GEOSX_UNUSED_VAR( inputRhs );

  }

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::TeamKernelBase::StackVariables
   *
   * Adds a stack array for the primary field.
   */
  struct StackVariables : public Base::StackVariables
  {
    /**
     * @brief Constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables( TeamLaplaceFEMKernel const & kernelComponent, LaunchContext & ctx ):
      Base::StackVariables( ctx ),
      kernelComponent( kernelComponent )
    {}

    TeamLaplaceFEMKernel const & kernelComponent;

    static constexpr size_t dim = 3;
    static constexpr size_t num_dofs_mesh_1d = 2; // TODO
    static constexpr size_t num_dofs_1d = 2; // TODO
    static constexpr size_t num_quads_1d = Base::StackVariables::num_quads_1d; // TODO

    // TODO abstract and encapsulate into object
    RAJA_TEAM_SHARED real64 mesh_basis[num_dofs_mesh_1d][num_quads_1d];
    RAJA_TEAM_SHARED real64 mesh_basis_gradient[num_dofs_mesh_1d][num_quads_1d];
    RAJA_TEAM_SHARED real64 basis[num_dofs_1d][num_quads_1d];
    RAJA_TEAM_SHARED real64 basis_gradient[num_dofs_1d][num_quads_1d];
    RAJA_TEAM_SHARED real64 weights[num_quads_1d];
    // TODO take into account batch_size / alias shared buffers / Generalize for non-tensor elements
    RAJA_TEAM_SHARED real64 mesh_nodes[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim];
    RAJA_TEAM_SHARED real64 jacobians[num_quads_1d][num_quads_1d][num_quads_1d][dim][dim];
    RAJA_TEAM_SHARED real64 dofs_in[num_dofs_1d][num_dofs_1d][num_dofs_1d];
    RAJA_TEAM_SHARED real64 q_gradient_values[num_quads_1d][num_quads_1d][num_quads_1d][dim];
    RAJA_TEAM_SHARED real64 Du[num_quads_1d][num_quads_1d][num_quads_1d][dim];
    RAJA_TEAM_SHARED real64 dofs_out[num_dofs_1d][num_dofs_1d][num_dofs_1d];
  };

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void kernelSetup( StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( stack );
    // TODO load/compute the different B and G.
  }

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geosx::finiteElement::TeamKernelBase::setup
   *
   * For the TeamLaplaceFEMKernel implementation, global values from the
   * primaryField, and degree of freedom numbers are placed into element local
   * stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( StackVariables & stack,
              localIndex const k ) const
  {
    Base::setup( stack, k );
    /// Computation of the Jacobians
    readField( stack, stack.kernelComponent.m_X, stack.mesh_nodes );
    interpolateGradientAtQuadraturePoints( stack,
                                           stack.mesh_basis,
                                           stack.mesh_basis_gradient,
                                           stack.mesh_nodes,
                                           stack.jacobians );

    /// Computation of the Gradient of the solution field
    readField( stack, stack.kernelComponent.m_primaryField, stack.dofs_in );
    interpolateGradientAtQuadraturePoints( stack,
                                           stack.basis,
                                           stack.basis_gradient,
                                           stack.dofs_in,
                                           stack.q_gradient_values );
  }

  /**
   * @copydoc geosx::finiteElement::TeamKernelBase::quadraturePointKernel
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( StackVariables & stack,
                              localIndex const quad_x,
                              localIndex const quad_y,
                              localIndex const quad_z ) const
  {
    constexpr int dim = StackVariables::dim;
    /// QFunction for Laplace operator
    // Load q-local gradients
    real64 u[dim];
    for (size_t d = 0; d < dim; d++)
    {
      u[d] = stack.q_gradient_values[quad_x][quad_y][quad_z][d];
    }
    // load q-local jacobian
    real64 J[dim][dim];
    for (size_t i = 0; i < dim; i++)
    {
      for (size_t j = 0; j < dim; j++)
      {
        J[i][j] = stack.jacobians[quad_x][quad_y][quad_z][i][j];
      }
    }
    real64 const detJ = J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2])
                      - J[1][0] * (J[0][1] * J[2][2] - J[2][1] * J[0][2])
                      + J[2][0] * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
    // adj(J_q)
    real64 AdjJ[dim][dim];
    AdjJ[0][0] = (J[1][1] * J[2][2]) - (J[1][2] * J[2][1]);
    AdjJ[0][1] = (J[2][1] * J[0][2]) - (J[0][1] * J[2][2]);
    AdjJ[0][2] = (J[0][1] * J[1][2]) - (J[1][1] * J[0][2]);
    AdjJ[1][0] = (J[2][0] * J[1][2]) - (J[1][0] * J[2][2]);
    AdjJ[1][1] = (J[0][0] * J[2][2]) - (J[0][2] * J[2][0]);
    AdjJ[1][2] = (J[1][0] * J[0][2]) - (J[0][0] * J[1][2]);
    AdjJ[2][0] = (J[1][0] * J[2][1]) - (J[2][0] * J[1][1]);
    AdjJ[2][1] = (J[2][0] * J[0][1]) - (J[0][0] * J[2][1]);
    AdjJ[2][2] = (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]);
    // Compute D_q = w_q * det(J_q) * J_q^-1 * J_q^-T = w_q / det(J_q) * adj(J_q) adj(J_q)^T
    real64 const weight = stack.weights[quad_x] * stack.weights[quad_y] * stack.weights[quad_z];
    real64 D[dim][dim];
    D[0][0] = weight / detJ * (AdjJ[0][0]*AdjJ[0][0] + AdjJ[0][1]*AdjJ[0][1] + AdjJ[0][2]*AdjJ[0][2]);
    D[1][0] = weight / detJ * (AdjJ[0][0]*AdjJ[1][0] + AdjJ[0][1]*AdjJ[1][1] + AdjJ[0][2]*AdjJ[1][2]);
    D[2][0] = weight / detJ * (AdjJ[0][0]*AdjJ[2][0] + AdjJ[0][1]*AdjJ[2][1] + AdjJ[0][2]*AdjJ[2][2]);
    D[0][1] = D[1][0];
    D[1][1] = weight / detJ * (AdjJ[1][0]*AdjJ[1][0] + AdjJ[1][1]*AdjJ[1][1] + AdjJ[1][2]*AdjJ[1][2]);
    D[2][1] = weight / detJ * (AdjJ[1][0]*AdjJ[2][0] + AdjJ[1][1]*AdjJ[2][1] + AdjJ[1][2]*AdjJ[2][2]);
    D[0][2] = D[2][0];
    D[1][2] = D[2][1];
    D[2][2] = weight / detJ * (AdjJ[2][0]*AdjJ[2][0] + AdjJ[2][1]*AdjJ[2][1] + AdjJ[2][2]*AdjJ[2][2]);
    // Compute D*u
    for (size_t d = 0; d < dim; d++)
    {
      stack.Du[quad_x][quad_y][quad_z][d] = D[d][0] * u[0]
                                          + D[d][1] * u[1]
                                          + D[d][2] * u[2];
    }
  }

  /**
   * @copydoc geosx::finiteElement::TeamKernelBase::complete
   *
   * Form element residual from the fully formed element Jacobian dotted with
   * the primary field and map the element local Jacobian/Residual to the
   * global matrix/vector.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( StackVariables & stack ) const
  {
    // Applying gradient of the test functions
    applyGradientTestFunctions( stack, stack.basis, stack.basis_gradient, stack.Du, stack.dofs_out );
    writeField( stack, stack.dofs_out, stack.kernelComponent.m_rhs );
  
    constexpr size_t num_dofs_1d = StackVariables::num_dofs_1d;
    RAJA_TEAM_SHARED real64 maxForce = 0; // TODO take into account batch_size
    // TODO put this into a lambda "iterator" function
    RAJA::expt::loop<thread_x> (stack.ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      RAJA::expt::loop<thread_y> (stack.ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
      {
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          maxForce = fmax( maxForce, fabs( stack.dofs_out[ dof_x ][ dof_y ][ dof_z ] ) ); // TODO make atomic
        }
      } );
    } );
    return maxForce;
  }

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The global primary field array.
  arrayView1d< real64 const > const m_primaryField;

  /// The global residual vector.
  arrayView1d< real64 > const m_rhs;
};

/// The factory used to construct a TeamLaplaceFEMKernel.
using TeamLaplaceFEMKernelFactory = finiteElement::KernelFactory< TeamLaplaceFEMKernel,
                                                                  arrayView1d< globalIndex const > const,
                                                                  globalIndex const,
                                                                  CRSMatrixView< real64, globalIndex const > const,
                                                                  arrayView1d< real64 > const,
                                                                  string const >;

} // namesapce geosx

#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_TEAMLAPLACEFEMKERNELS_HPP_
