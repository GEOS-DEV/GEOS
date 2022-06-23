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
#include "finiteElement/kernelInterface/QuadratureFunctionsHelper.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

template < size_t Order >
class LagrangeBasis;

template <>
class LagrangeBasis<1> : public finiteElement::LagrangeBasis1 { };

//*****************************************************************************
/**
 * @brief Implements kernels for solving Laplace's equation.
 * @copydoc geosx::finiteElement::TeamKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### TeamLaplaceFEMKernel Description
 * Implements the TeamKernelBase interface functions required for solving Laplace's
 * equation using the finite element kernel application functions such as
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

  template < size_t num_dofs_mesh_1d, size_t num_quads_1d, size_t dim, size_t batch_size >
  struct MeshStackVariables
  {
    MeshStackVariables( LaunchContext & ctx )
    {
      using RAJA::RangeSegment;

      // Mesh basis
      GEOSX_STATIC_SHARED real64 s_mesh_basis[num_dofs_mesh_1d][num_quads_1d];
      mesh_basis = &s_mesh_basis;
      loop<thread_z> (ctx, RangeSegment(0, 1), [&] (const int tidz)
      {
        loop<thread_y> (ctx, RangeSegment(0, num_dofs_mesh_1d), [&] (size_t d)
        {
          loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t q)
          {
            GEOSX_UNUSED_VAR( tidz );
            s_mesh_basis[ d ][ q ] =
              LagrangeBasis<num_dofs_mesh_1d-1>::value(
                d, LagrangeBasis<num_quads_1d-1>::parentSupportCoord( q ) );
          } );
        } );
      } );

      // Mesh basis gradient
      GEOSX_STATIC_SHARED real64 s_mesh_basis_gradient[num_dofs_mesh_1d][num_quads_1d];
      mesh_basis_gradient = &s_mesh_basis_gradient;
      loop<thread_z> (ctx, RangeSegment(0, 1), [&] (const int tidz)
      {
        loop<thread_y> (ctx, RangeSegment(0, num_dofs_mesh_1d), [&] (size_t d)
        {
          loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t q)
          {
            GEOSX_UNUSED_VAR( tidz );
            s_mesh_basis_gradient[ d ][ q ] =
              LagrangeBasis<num_dofs_mesh_1d-1>::gradient(
                d, LagrangeBasis<num_quads_1d-1>::parentSupportCoord( q ) );
          } );
        } );
      } );

      size_t const batch_index = GEOSX_THREAD_ID(z);
      // Mesh nodes
      GEOSX_STATIC_SHARED real64 s_mesh_nodes[batch_size][num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim];
      mesh_nodes = &s_mesh_nodes[batch_index];

      // Mesh jacobians
      GEOSX_STATIC_SHARED real64 s_jacobians[batch_size][num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim][dim];
      jacobians = &s_jacobians[batch_index];
    }

    // Mesh basis
    real64 ( * mesh_basis )[num_dofs_mesh_1d][num_quads_1d];
    real64 const ( & getBasis() const )[num_dofs_mesh_1d][num_quads_1d]
    {
      return *mesh_basis;
    }
    real64 ( & getBasis() )[num_dofs_mesh_1d][num_quads_1d]
    {
      return *mesh_basis;
    }

    // Mesh basis gradient
    real64 ( * mesh_basis_gradient )[num_dofs_mesh_1d][num_quads_1d];
    real64 const ( & getBasisGradient() const )[num_dofs_mesh_1d][num_quads_1d]
    {
      return *mesh_basis_gradient;
    }
    real64 ( & getBasisGradient() )[num_dofs_mesh_1d][num_quads_1d]
    {
      return *mesh_basis_gradient;
    }

    // Mesh nodes
    real64 ( * mesh_nodes )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim]; // Could be in registers
    real64 const ( & getNodes() const )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim]
    {
      return *mesh_nodes;
    }
    real64 ( & getNodes() )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim]
    {
      return *mesh_nodes;
    }

    // Mesh jacobians
    real64 ( * jacobians )[num_quads_1d][num_quads_1d][num_quads_1d][dim][dim]; // Can be in registers
    real64 ( & getJacobians() )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim][dim]
    {
      return *jacobians;
    }
    real64 const ( & getJacobians() const )[num_dofs_mesh_1d][num_dofs_mesh_1d][num_dofs_mesh_1d][dim][dim]
    {
      return *jacobians;
    }
  };

  template < size_t num_dofs_1d, size_t num_quads_1d, size_t dim, size_t batch_size >
  struct ElementStackVariables
  {
    ElementStackVariables( LaunchContext & ctx )
    {
      using RAJA::RangeSegment;

      // Element basis
      GEOSX_STATIC_SHARED real64 s_basis[num_dofs_1d][num_quads_1d];
      basis = &s_basis;
      loop<thread_z> (ctx, RangeSegment(0, 1), [&] (const int tidz)
      {
        loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t d)
        {
          loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t q)
          {
            GEOSX_UNUSED_VAR( tidz );
            s_basis[ d ][ q ] =
              LagrangeBasis<num_dofs_1d-1>::value(
                d, LagrangeBasis<num_quads_1d-1>::parentSupportCoord( q ) );
          } );
        } );
      } );

      // Element basis gradient
      GEOSX_STATIC_SHARED real64 s_basis_gradient[num_dofs_1d][num_quads_1d];
      basis_gradient = &s_basis_gradient;
      loop<thread_z> (ctx, RangeSegment(0, 1), [&] (const int tidz)
      {
        loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t d)
        {
          loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t q)
          {
            GEOSX_UNUSED_VAR( tidz );
            s_basis_gradient[ d ][ q ] =
              LagrangeBasis<num_dofs_1d-1>::gradient(
                d, LagrangeBasis<num_quads_1d-1>::parentSupportCoord( q ) );
          } );
        } );
      } );

      size_t const batch_index = GEOSX_THREAD_ID(z);
      // Element input dofs of the primary field
      GEOSX_STATIC_SHARED real64 s_dofs_in[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d];
      dofs_in = &s_dofs_in[batch_index];

      // Element primary field gradients at quadrature points
      GEOSX_STATIC_SHARED real64 s_q_gradient_values[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d][dim];
      q_gradient_values = &s_q_gradient_values[batch_index];

      // Element "geometric factors"
      GEOSX_STATIC_SHARED real64 s_Du[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d][dim];
      Du = &s_Du[batch_index];

      // Element contribution to the residual
      GEOSX_STATIC_SHARED real64 s_dofs_out[batch_size][num_dofs_1d][num_dofs_1d][num_dofs_1d];
      dofs_out = &s_dofs_out[batch_index];
    }

    // Element basis
    real64 ( * basis )[num_dofs_1d][num_quads_1d];
    real64 const ( & getBasis() const )[num_dofs_1d][num_quads_1d]
    {
      return *basis;
    }
    real64 ( & getBasis() )[num_dofs_1d][num_quads_1d]
    {
      return *basis;
    }

    // Element basis gradient
    real64 ( * basis_gradient )[num_dofs_1d][num_quads_1d];
    real64 const ( & getBasisGradient() const )[num_dofs_1d][num_quads_1d]
    {
      return *basis_gradient;
    }
    real64 ( & getBasisGradient() )[num_dofs_1d][num_quads_1d]
    {
      return *basis_gradient;
    }

    // Element input dofs of the primary field
    real64 ( * dofs_in )[num_dofs_1d][num_dofs_1d][num_dofs_1d]; // Could be in registers
    real64 const ( & getDofsIn() const )[num_dofs_1d][num_dofs_1d][num_dofs_1d]
    {
      return *dofs_in;
    }
    real64 ( & getDofsIn() )[num_dofs_1d][num_dofs_1d][num_dofs_1d]
    {
      return *dofs_in;
    }

    // Element primary field gradients at quadrature points
    real64 ( * q_gradient_values )[num_quads_1d][num_quads_1d][num_quads_1d][dim]; // Can be in registers
    real64 const ( & getGradientValues() const )[num_quads_1d][num_quads_1d][num_quads_1d][dim]
    {
      return *q_gradient_values;
    }
    real64 ( & getGradientValues() )[num_quads_1d][num_quads_1d][num_quads_1d][dim]
    {
      return *q_gradient_values;
    }

    // Element "geometric factors"
    real64 ( * Du )[num_quads_1d][num_quads_1d][num_quads_1d][dim]; // Could be in registers
    real64 const ( & getQuadValues() const )[num_quads_1d][num_quads_1d][num_quads_1d][dim]
    {
      return *Du;
    }
    real64 ( & getQuadValues() )[num_quads_1d][num_quads_1d][num_quads_1d][dim]
    {
      return *Du;
    }

    // Element contribution to the residual
    real64 ( * dofs_out )[num_dofs_1d][num_dofs_1d][num_dofs_1d]; // Can be in registers
    real64 const ( & getDofsOut() const )[num_dofs_1d][num_dofs_1d][num_dofs_1d]
    {
      return *dofs_out;
    }
    real64 ( & getDofsOut() )[num_dofs_1d][num_dofs_1d][num_dofs_1d]
    {
      return *dofs_out;
    }
  };

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::TeamKernelBase::StackVariables
   *
   * Adds a stack array for the primary field.
   */
  struct StackVariables : public Base::StackVariables
  {
    static constexpr size_t dim = 3;
    static constexpr size_t num_dofs_mesh_1d = 2; // TODO
    static constexpr size_t num_dofs_1d = 2; // TODO
    using Base::StackVariables::num_quads_1d;
    using Base::StackVariables::batch_size;
  
    /**
     * @brief Constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables( TeamLaplaceFEMKernel const & kernelComponent, LaunchContext & ctx ):
      Base::StackVariables( ctx ),
      kernelComponent( kernelComponent ),
      mesh( ctx ),
      element( ctx )
    {}

    TeamLaplaceFEMKernel const & kernelComponent;

    // TODO alias shared buffers / Generalize for non-tensor elements
    MeshStackVariables< num_dofs_mesh_1d, num_quads_1d, dim, batch_size > mesh;
    ElementStackVariables< num_dofs_1d, num_quads_1d, dim, batch_size > element;
    // TODO abstract and encapsulate into object
    real64 ( * weights )[num_quads_1d];

    /// Shared memory buffers, using buffers allows to avoid using too much shared memory.
    static constexpr size_t buffer_size = num_quads_1d * num_quads_1d * num_quads_1d * dim * dim;
    real64 * shared_mem_buffer_1;
    real64 * shared_mem_buffer_2;
  };

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void kernelSetup( StackVariables & stack ) const
  {
    // Initialize quadrature weights
    GEOSX_STATIC_SHARED real64 s_weights[StackVariables::num_quads_1d];
    stack.weights = &s_weights;
    // TODO generalize
    s_weights[0] = 1.0;
    s_weights[1] = 1.0;

    // "Allocate" shared memory for the buffers.
    constexpr size_t batch_size = StackVariables::batch_size;
    constexpr size_t buffer_size = StackVariables::buffer_size;
    size_t const tidz = GEOSX_THREAD_ID(z);

    GEOSX_STATIC_SHARED real64 shared_buffer_1[batch_size][buffer_size];
    stack.shared_mem_buffer_1 = ( real64 * )&shared_buffer_1[tidz];
    GEOSX_STATIC_SHARED real64 shared_buffer_2[batch_size][buffer_size];
    stack.shared_mem_buffer_2 = ( real64 * )&shared_buffer_2[tidz];
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
              localIndex const element_index ) const
  {
    Base::setup( stack, element_index );
    /// Computation of the Jacobians
    readField( stack, stack.kernelComponent.m_X, stack.mesh.getNodes() );
    interpolateGradientAtQuadraturePoints( stack,
                                           stack.mesh.getBasis(),
                                           stack.mesh.getBasisGradient(),
                                           stack.mesh.getNodes(),
                                           stack.mesh.getJacobians() );

    /// Computation of the Gradient of the solution field
    readField( stack, stack.kernelComponent.m_primaryField, stack.element.getDofsIn() );
    interpolateGradientAtQuadraturePoints( stack,
                                           stack.element.getBasis(),
                                           stack.element.getBasisGradient(),
                                           stack.element.getDofsIn(),
                                           stack.element.getGradientValues() );
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
    /// QFunction for Laplace operator
    constexpr int dim = StackVariables::dim;

    // Load q-local gradients
    real64 grad_u[ dim ];
    qLocalLoad( quad_x, quad_y, quad_z, stack.element.getGradientValues(), grad_u );

    // load q-local jacobian
    real64 J[ dim ][ dim ];
    qLocalLoad( quad_x, quad_y, quad_z, stack.mesh.getJacobians(), J );

    // Compute D_q = w_q * det(J_q) * J_q^-1 * J_q^-T = w_q / det(J_q) * adj(J_q) * adj(J_q)^T
    real64 const weight = (*stack.weights)[ quad_x ] * (*stack.weights)[ quad_y ] * (*stack.weights)[ quad_z ];

    real64 const detJ = determinant( J );

    real64 AdjJ[ dim ][ dim ];
    adjugate( J, AdjJ );

    real64 D[ dim ][ dim ];
    D[0][0] = weight / detJ * (AdjJ[0][0]*AdjJ[0][0] + AdjJ[0][1]*AdjJ[0][1] + AdjJ[0][2]*AdjJ[0][2]);
    D[1][0] = weight / detJ * (AdjJ[0][0]*AdjJ[1][0] + AdjJ[0][1]*AdjJ[1][1] + AdjJ[0][2]*AdjJ[1][2]);
    D[2][0] = weight / detJ * (AdjJ[0][0]*AdjJ[2][0] + AdjJ[0][1]*AdjJ[2][1] + AdjJ[0][2]*AdjJ[2][2]);
    D[0][1] = D[1][0];
    D[1][1] = weight / detJ * (AdjJ[1][0]*AdjJ[1][0] + AdjJ[1][1]*AdjJ[1][1] + AdjJ[1][2]*AdjJ[1][2]);
    D[2][1] = weight / detJ * (AdjJ[1][0]*AdjJ[2][0] + AdjJ[1][1]*AdjJ[2][1] + AdjJ[1][2]*AdjJ[2][2]);
    D[0][2] = D[2][0];
    D[1][2] = D[2][1];
    D[2][2] = weight / detJ * (AdjJ[2][0]*AdjJ[2][0] + AdjJ[2][1]*AdjJ[2][1] + AdjJ[2][2]*AdjJ[2][2]);

    // Compute D_q * grad_u_q
    real64 Du[ dim ];
    for (size_t d = 0; d < dim; d++)
    {
      Du[d] = D[d][0] * grad_u[0]
            + D[d][1] * grad_u[1]
            + D[d][2] * grad_u[2];
    }
    qLocalWrite( quad_x, quad_y, quad_z, Du, stack.element.getQuadValues() );
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
    using RAJA::RangeSegment;

    // Applying gradient of the test functions
    applyGradientTestFunctions( stack,
                                stack.element.getBasis(),
                                stack.element.getBasisGradient(),
                                stack.element.getQuadValues(),
                                stack.element.getDofsOut() );
    writeAddField( stack, stack.element.getDofsOut(), stack.kernelComponent.m_rhs );

    constexpr size_t num_dofs_1d = StackVariables::num_dofs_1d;
    real64 maxForce = 0;
    // TODO put this into a lambda "iterator" function
    auto & dofs_out = stack.element.getDofsOut();
    loop<thread_x> (stack.ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      loop<thread_y> (stack.ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
      {
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          maxForce = fmax( maxForce, fabs( dofs_out[ dof_x ][ dof_y ][ dof_z ] ) );
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
