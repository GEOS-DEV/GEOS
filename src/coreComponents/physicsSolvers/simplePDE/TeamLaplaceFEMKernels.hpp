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

#include "finiteElement/TeamKernelInterface/TeamKernelBase.hpp"
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/TeamKernelFunctions.hpp"
#include "finiteElement/TeamKernelInterface/QuadraturePointKernelFunctions/QuadratureFunctionsHelper.hpp"
#include "finiteElement/TeamKernelInterface/QuadraturePointKernelFunctions/qLocalLoad.hpp"
#include "finiteElement/TeamKernelInterface/QuadraturePointKernelFunctions/qLocalWrite.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/MeshStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/ElementStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/QuadratureWeightsStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/SharedMemStackVariables.hpp"

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
                        CONSTITUTIVE_TYPE & inputConstitutiveType, // end of default args
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
    static constexpr localIndex dim = 3;
    static constexpr localIndex num_dofs_mesh_1d = 2; // TODO
    static constexpr localIndex num_dofs_1d = 2; // TODO
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
      element( ctx ),
      weights( ctx ),
      shared_mem( ctx )
    {}

    TeamLaplaceFEMKernel const & kernelComponent;

    // TODO alias shared buffers / Generalize for non-tensor elements
    MeshStackVariables< num_dofs_mesh_1d, num_quads_1d, dim, batch_size > mesh;
    ElementStackVariables< num_dofs_1d, num_quads_1d, dim, batch_size > element;
    QuadratureWeightsStackVariables< num_quads_1d > weights;

    /// Shared memory buffers, using buffers allows to avoid using too much shared memory.
    static constexpr localIndex buffer_size = num_quads_1d * num_quads_1d * num_quads_1d * dim;
    static constexpr localIndex num_buffers = 2 * dim;
    SharedMemStackVariables< buffer_size, num_buffers, batch_size > shared_mem;
  };

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void kernelSetup( StackVariables & stack ) const
  {
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
                                           stack.mesh.basis.getValuesAtQuadPts(),
                                           stack.mesh.basis.getGradientValuesAtQuadPts(),
                                           stack.mesh.getNodes(),
                                           stack.mesh.getJacobians() );

    /// Computation of the Gradient of the solution field
    readField( stack, stack.kernelComponent.m_primaryField, stack.element.getDofsIn() );
    interpolateGradientAtQuadraturePoints( stack,
                                           stack.element.basis.getValuesAtQuadPts(),
                                           stack.element.basis.getGradientValuesAtQuadPts(),
                                           stack.element.getDofsIn(),
                                           stack.element.getGradientValues() );
  }

  /**
   * @copydoc geosx::finiteElement::TeamKernelBase::quadraturePointKernel
   */
  template < typename QuadraturePointIndex >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( StackVariables & stack,
                              QuadraturePointIndex const & quad_index ) const
  {
    /// QFunction for Laplace operator
    constexpr int dim = StackVariables::dim;

    // Load q-local gradients
    real64 grad_u[ dim ];
    qLocalLoad( quad_index, stack.element.getGradientValues(), grad_u );

    // load q-local jacobian
    real64 J[ dim ][ dim ];
    qLocalLoad( quad_index, stack.mesh.getJacobians(), J );

    // Compute D_q = w_q * det(J_q) * J_q^-1 * J_q^-T = w_q / det(J_q) * adj(J_q) * adj(J_q)^T
    real64 const weight = stack.weights( quad_index );

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
    for (localIndex d = 0; d < dim; d++)
    {
      Du[d] = D[d][0] * grad_u[0]
            + D[d][1] * grad_u[1]
            + D[d][2] * grad_u[2];
    }
    qLocalWrite( quad_index, Du, stack.element.getQuadValues() );
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
                                stack.element.basis.getValuesAtQuadPts(),
                                stack.element.basis.getGradientValuesAtQuadPts(),
                                stack.element.getQuadValues(),
                                stack.element.getDofsOut() );
    writeAddField( stack, stack.element.getDofsOut(), stack.kernelComponent.m_rhs );

    constexpr localIndex num_dofs_1d = StackVariables::num_dofs_1d;
    real64 maxForce = 0;
    // TODO put this into a lambda "iterator" function
    auto & dofs_out = stack.element.getDofsOut();
    loop<thread_x> (stack.ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
    {
      loop<thread_y> (stack.ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
      {
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
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
                                                                  arrayView1d< real64 > const,
                                                                  string const >;

} // namesapce geosx

#endif // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_TEAMLAPLACEFEMKERNELS_HPP_
