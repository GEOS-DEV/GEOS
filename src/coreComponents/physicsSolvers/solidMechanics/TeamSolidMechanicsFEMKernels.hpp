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
 * @file TeamSolidMechanicsFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_TEAMSOLIDMECHANICSFEMKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_TEAMSOLIDMECHANICSFEMKERNELS_HPP_

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "finiteElement/TeamKernelInterface/TeamKernelBase.hpp"
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/TeamKernelFunctions.hpp"
#include "finiteElement/TeamKernelInterface/QuadraturePointKernelFunctions/QuadratureFunctionsHelper.hpp"
#include "finiteElement/TeamKernelInterface/QuadraturePointKernelFunctions/qLocalLoad.hpp"
#include "finiteElement/TeamKernelInterface/QuadraturePointKernelFunctions/qLocalWrite.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/MeshStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/VectorElementStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/QuadratureWeightsStackVariables.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/SharedStackVariables.hpp"
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
 * ### TeamSolidMechanicsFEMKernel Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static equilibrium equations using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the ImplicitKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class TeamSolidMechanicsFEMKernel :
  public finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                        CONSTITUTIVE_TYPE,
                                        FE_TYPE,
                                        3,
                                        3 >
{
public:
  /// An alias for the base class.
  using Base = finiteElement::TeamKernelBase< SUBREGION_TYPE,
                                              CONSTITUTIVE_TYPE,
                                              FE_TYPE,
                                              3,
                                              3 >;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::TeamKernelBase::TeamKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  TeamSolidMechanicsFEMKernel( NodeManager const & nodeManager,
                               EdgeManager const & edgeManager,
                               FaceManager const & faceManager,
                               localIndex const targetRegionIndex,
                               SUBREGION_TYPE const & elementSubRegion,
                               FE_TYPE const & finiteElementSpace,
                               CONSTITUTIVE_TYPE & inputConstitutiveType, // end of default args
                               arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const inputSrc,
                               arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const inputDst ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition() ),
    m_src( inputSrc ),
    m_dst( inputDst )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
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
    StackVariables( TeamSolidMechanicsFEMKernel const & kernelComponent, LaunchContext & ctx ):
      Base::StackVariables( ctx ),
      kernelComponent( kernelComponent ),
      mesh( ctx ),
      element( ctx ),
      weights( ctx ),
      shared_mem( ctx )
    {}

    TeamSolidMechanicsFEMKernel const & kernelComponent;

    // TODO alias shared buffers / Generalize for non-tensor elements
    MeshStackVariables< num_dofs_mesh_1d, num_quads_1d, dim, batch_size > mesh;
    VectorElementStackVariables< num_dofs_1d, num_quads_1d, dim, batch_size > element;
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
   * For the TeamSolidMechanicsFEMKernel implementation, global values from the
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
    readField( stack, stack.kernelComponent.m_src, stack.element.getDofsIn() );
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
    real64 grad_u[ dim ][ dim ];
    qLocalLoad( quad_index, stack.element.getGradientValues(), grad_u );

    // load q-local jacobian
    real64 J[ dim ][ dim ];
    qLocalLoad( quad_index, stack.mesh.getJacobians(), J );

    // Compute D_q = w_q * det(J_q) * J_q^-1 * J_q^-T = w_q / det(J_q) * adj(J_q) * adj(J_q)^T
    real64 const weight = stack.weights( quad_index );

    real64 const detJ = determinant( J );

    real64 AdjJ[ dim ][ dim ];
    adjugate( J, AdjJ );

    // physical gradient of displacement
    real64 grad_phys_u[ dim ][ dim ];
    for (localIndex i = 0; i < dim; i++)
    {
      for (localIndex j = 0; j < dim; j++)
      {
        real64 val = 0.0;
        for (localIndex k = 0; k < dim; k++)
        {
          val += AdjJ[ i ][ k ] * grad_u[ k ][ j ] / detJ;
        }
        grad_phys_u[ i ][ j ] = val;
      }
    }

    // Computation of the strain
    real64 strain[ dim ][ dim ];
    for (localIndex j = 0; j < dim; j++)
    {
      for (localIndex i = 0; i < dim; i++)
      {
        strain[ i ][ j ] = 0.5 * ( grad_phys_u[ i ][ j ] + grad_phys_u[ j ][ i ] );
      }
    }

    // Local material stiffness values
    // real64 D_mat[ 2*dim ][ 2*dim ];
    // qLocalLoad( quad_index, stack.element.getMaterialStiffnes(), D_mat );
    
    // Computation of the stress
    real64 stress[ dim ][ dim ];
    // computeStress< DiscretizationOps >( D_mat, strain, stress ); // Tensor product
    // DiscretizationOps::computeStress( stack.element.getMaterialStiffnes(), strain, stress ); // Tensor product
    computeStress( stack, strain, stress );

    // Apply J^-T
    real64 D[ dim ][ dim ];
    for (localIndex i = 0; i < dim; i++)
    {
      for (localIndex j = 0; j < dim; j++)
      {
        real64 val = 0.0;
        for (localIndex k = 0; k < dim; k++)
        {
          val += AdjJ[ k ][ i ] * stress[ k ][ j ] / detJ;
        }
        D[ i ][ j ] = val;
      }
    }
    qLocalWrite( quad_index, D, stack.element.getQuadValues() );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeStress( StackVariables & stack,
                      real64 const ( &strain )[ 3 ][ 3 ],
                      real64 ( & stress )[ 3 ][ 3 ] ) const
  {
    localIndex const k = stack.element_index;
    real64 symm_strain[6];
    symm_strain[0] = strain[0][0];
    symm_strain[1] = strain[1][1];
    symm_strain[2] = strain[2][2];
    symm_strain[3] = 2 * strain[0][1];
    symm_strain[4] = 2 * strain[0][2];
    symm_strain[5] = 2 * strain[1][2];
    real64 symm_stress[6];
    stack.kernelComponent.m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, 0, symm_strain, symm_stress );
    stress[0][0] = symm_stress[0];
    stress[1][1] = symm_stress[1];
    stress[2][2] = symm_stress[2];
    stress[0][1] = symm_stress[3];
    stress[0][2] = symm_stress[4];
    stress[1][2] = symm_stress[5];
    stress[1][0] = symm_stress[3];
    stress[2][0] = symm_stress[4];
    stress[2][1] = symm_stress[5];
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
    writeAddField( stack, stack.element.getDofsOut(), stack.kernelComponent.m_dst );

    // constexpr localIndex num_dofs_1d = StackVariables::num_dofs_1d;
    // real64 maxForce = 0;
    // // TODO put this into a lambda "iterator" function
    // auto & dofs_out = stack.element.getDofsOut();
    // loop<thread_x> (stack.ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
    // {
    //   loop<thread_y> (stack.ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
    //   {
    //     for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    //     {
    //       maxForce = fmax( maxForce, fabs( dofs_out[ dof_x ][ dof_y ][ dof_z ] ) );
    //     }
    //   } );
    // } );
    // return maxForce;
    return 0.0;
  }

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_src;
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const m_dst;
};

/// The factory used to construct a TeamSolidMechanicsFEMKernel.
using TeamSolidMechanicsFEMKernelFactory = finiteElement::KernelFactory< TeamSolidMechanicsFEMKernel,
                                                                         arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const,
                                                                         arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const >;

} // namesapce geosx

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_TEAMSOLIDMECHANICSFEMKERNELS_HPP_
