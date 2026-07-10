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
 * @file SinglePhasePoromechanicsConformingFracturesALMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALMKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALMKERNELS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsConformingContactKernelsHelper.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"

namespace geos
{

namespace poromechanicsALMKernels
{

/**
 * @brief Kernel to assemble the pressure force contribution to displacement equations
 *        AND the Jacobian derivative dR_u/dP for fully implicit coupling.
 *
 * This kernel is similar to SolidMechanicsConformingPressureContributionKernels::AssemblePressureContribution
 * but also adds the Jacobian terms needed for fully implicit poromechanics coupling.
 * The sign convention matches the existing sequential coupling kernel.
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class AssembleForceResidualDerivativeWrtPressure :
  public solidMechanicsConformingContactKernels::ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                                                                FE_TYPE >
{

public:
  /// Alias for the base class
  using Base = solidMechanicsConformingContactKernels::ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                                                                      FE_TYPE >;

  /// Number of nodes per element
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// The number of displacement dofs per element.
  static constexpr int numUdofs = Base::numUdofs;

  /// The number of bubble dofs per element.
  static constexpr int numBdofs = Base::numBdofs;

  using Base::m_elemsToFaces;
  using Base::m_faceToNodes;
  using Base::m_X;
  using Base::m_finiteElementSpace;
  using Base::m_dofNumber;
  using Base::m_bDofNumber;
  using Base::m_dofRankOffset;
  using Base::m_rotationMatrix;
  using Base::m_matrix;
  using Base::m_rhs;

  /**
   * @brief Constructor
   */
  AssembleForceResidualDerivativeWrtPressure( NodeManager const & nodeManager,
                                               EdgeManager const & edgeManager,
                                               FaceManager const & faceManager,
                                               localIndex const targetRegionIndex,
                                               FaceElementSubRegion & elementSubRegion,
                                               FE_TYPE const & finiteElementSpace,
                                               CONSTITUTIVE_TYPE & inputConstitutiveType,
                                               arrayView1d< globalIndex const > const uDofNumber,
                                               arrayView1d< globalIndex const > const bDofNumber,
                                               globalIndex const rankOffset,
                                               CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                               arrayView1d< real64 > const inputRhs,
                                               real64 const inputDt,
                                               arrayView1d< localIndex const > const & faceElementList,
                                               arrayView1d< globalIndex const > const & pressureDofNumber ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          uDofNumber,
          bDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt,
          faceElementList ),
    m_pressure( elementSubRegion.getField< fields::flow::pressure >().toView() ),
    m_pressureDofNumber( pressureDofNumber )
  {}

  /**
   * @copydoc finiteElement::InterfaceKernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
      dispEqnRowIndices{},
      bEqnRowIndices{},
      pressureColIndex( 0 ),
      unitNormal{},
      localRu{},
      localRb{},
      dRudP{},
      dRbdP{}
    {}

    /// Row indices for displacement equations
    globalIndex dispEqnRowIndices[numUdofs];

    /// Row indices for bubble equations
    globalIndex bEqnRowIndices[numBdofs];

    /// Column index for pressure
    globalIndex pressureColIndex;

    /// Unit normal from rotation matrix
    real64 unitNormal[3];

    /// Local residual for displacement
    real64 localRu[numUdofs];

    /// Local residual for bubble
    real64 localRb[numBdofs];

    /// Jacobian dRu/dP
    real64 dRudP[numUdofs];

    /// Jacobian dRb/dP
    real64 dRbdP[numBdofs];
  };

  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
  }

  /**
   * @brief Setup for each element
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    constexpr int shift = numNodesPerElem * 3;

    int permutation[numNodesPerElem];
    m_finiteElementSpace.getPermutation( permutation );

    localIndex const kf0 = m_elemsToFaces[k][0];
    localIndex const kf1 = m_elemsToFaces[k][1];

    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const kn0 = m_faceToNodes( kf0, a );
      localIndex const kn1 = m_faceToNodes( kf1, a );
      for( int i = 0; i < 3; ++i )
      {
        stack.dispEqnRowIndices[a * 3 + i] = m_dofNumber[kn0] + i - m_dofRankOffset;
        stack.dispEqnRowIndices[shift + a * 3 + i] = m_dofNumber[kn1] + i - m_dofRankOffset;
        stack.X[a][i] = m_X[m_faceToNodes( kf0, permutation[a] )][i];
      }
    }

    for( int i = 0; i < 3; ++i )
    {
      stack.unitNormal[i] = m_rotationMatrix( k, i, 0 );
    }

    for( int i = 0; i < 3; ++i )
    {
      stack.bEqnRowIndices[i] = m_bDofNumber[kf0] + i - m_dofRankOffset;
      stack.bEqnRowIndices[3 + i] = m_bDofNumber[kf1] + i - m_dofRankOffset;
    }

    stack.pressureColIndex = m_pressureDofNumber[k];
  }

  /**
   * @brief Complete - compute residuals and Jacobians
   */
  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    // Compute residual and Jacobian for displacement DOFs:
    // R_u = transpose(Atu) * unitNormal * pressure
    // dR_u/dP = transpose(Atu) * unitNormal
    // Note: This matches the existing AssemblePressureContribution kernel used for sequential coupling.
    LvArray::tensorOps::Ri_eq_AjiBj< numUdofs, 3 >( stack.localRu, stack.localAtu, stack.unitNormal );
    LvArray::tensorOps::Ri_eq_AjiBj< numUdofs, 3 >( stack.dRudP, stack.localAtu, stack.unitNormal );

    // Scale residual by pressure
    LvArray::tensorOps::scale< numUdofs >( stack.localRu, m_pressure[k] );

    // Compute residual and Jacobian for bubble DOFs:
    // R_b = transpose(Atb) * unitNormal * pressure
    // dR_b/dP = transpose(Atb) * unitNormal
    LvArray::tensorOps::Ri_eq_AjiBj< numBdofs, 3 >( stack.localRb, stack.localAtb, stack.unitNormal );
    LvArray::tensorOps::Ri_eq_AjiBj< numBdofs, 3 >( stack.dRbdP, stack.localAtb, stack.unitNormal );

    // Scale residual by pressure
    LvArray::tensorOps::scale< numBdofs >( stack.localRb, m_pressure[k] );

    // Assemble displacement contributions
    for( localIndex i = 0; i < numUdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[i] );

      if( dof < 0 || dof >= m_matrix.numRows() )
        continue;

      // Add residual
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRu[i] );

      // Add Jacobian
      m_matrix.template addToRow< parallelDeviceAtomic >( dof,
                                                          &stack.pressureColIndex,
                                                          &stack.dRudP[i],
                                                          1 );
    }

    // Assemble bubble contributions
    for( localIndex i = 0; i < numBdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.bEqnRowIndices[i] );

      if( dof < 0 || dof >= m_matrix.numRows() )
        continue;

      // Add residual
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRb[i] );

      // Add Jacobian
      m_matrix.template addToRow< parallelDeviceAtomic >( dof,
                                                          &stack.pressureColIndex,
                                                          &stack.dRbdP[i],
                                                          1 );
    }

    return 0.0;
  }

protected:

  /// Pressure array
  arrayView1d< real64 const > const m_pressure;

  /// Pressure DOF numbers
  arrayView1d< globalIndex const > const m_pressureDofNumber;
};

/// Factory for AssembleForceResidualDerivativeWrtPressure kernel
using AssembleForceResidualDerivativeWrtPressureFactory =
  finiteElement::InterfaceKernelFactory< AssembleForceResidualDerivativeWrtPressure,
                                         arrayView1d< globalIndex const > const,
                                         arrayView1d< globalIndex const > const,
                                         globalIndex const,
                                         CRSMatrixView< real64, globalIndex const > const,
                                         arrayView1d< real64 > const,
                                         real64 const,
                                         arrayView1d< localIndex const > const,
                                         arrayView1d< globalIndex const > const >;


/**
 * @brief Kernel to assemble the Jacobian derivative dR_p/du and dR_p/db for fully implicit coupling.
 *
 * This kernel computes the derivative of the flow accumulation residual with respect to
 * nodal displacement and bubble displacement, using the properly computed Atu/Atb matrices.
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class AssembleFluidMassResidualDerivativeWrtDisplacement :
  public solidMechanicsConformingContactKernels::ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                                                                FE_TYPE >
{

public:
  /// Alias for the base class
  using Base = solidMechanicsConformingContactKernels::ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                                                                      FE_TYPE >;

  /// Number of nodes per element
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// The number of displacement dofs per element.
  static constexpr int numUdofs = Base::numUdofs;

  /// The number of bubble dofs per element.
  static constexpr int numBdofs = Base::numBdofs;

  using Base::m_elemsToFaces;
  using Base::m_faceToNodes;
  using Base::m_X;
  using Base::m_finiteElementSpace;
  using Base::m_dofNumber;
  using Base::m_bDofNumber;
  using Base::m_dofRankOffset;
  using Base::m_rotationMatrix;
  using Base::m_matrix;
  using Base::m_rhs;

  /**
   * @brief Constructor
   */
  AssembleFluidMassResidualDerivativeWrtDisplacement( NodeManager const & nodeManager,
                                                       EdgeManager const & edgeManager,
                                                       FaceManager const & faceManager,
                                                       localIndex const targetRegionIndex,
                                                       FaceElementSubRegion & elementSubRegion,
                                                       FE_TYPE const & finiteElementSpace,
                                                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                                                       arrayView1d< globalIndex const > const uDofNumber,
                                                       arrayView1d< globalIndex const > const bDofNumber,
                                                       globalIndex const rankOffset,
                                                       CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                                       arrayView1d< real64 > const inputRhs,
                                                       real64 const inputDt,
                                                       arrayView1d< localIndex const > const & faceElementList,
                                                       arrayView1d< globalIndex const > const & pressureDofNumber,
                                                       arrayView1d< real64 const > const & density,
                                                       arrayView1d< integer const > const & fractureState ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          uDofNumber,
          bDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt,
          faceElementList ),
    m_pressureDofNumber( pressureDofNumber ),
    m_density( density ),
    m_fractureState( fractureState )
  {}

  /**
   * @copydoc finiteElement::InterfaceKernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
      dispColIndices{},
      bColIndices{},
      pressureRowIndex( 0 ),
      unitNormal{},
      dRpdU{},
      dRpdB{}
    {}

    /// Column indices for displacement
    globalIndex dispColIndices[numUdofs];

    /// Column indices for bubble
    globalIndex bColIndices[numBdofs];

    /// Row index for pressure equation
    globalIndex pressureRowIndex;

    /// Unit normal from rotation matrix
    real64 unitNormal[3];

    /// Jacobian dRp/dU
    real64 dRpdU[numUdofs];

    /// Jacobian dRp/dB
    real64 dRpdB[numBdofs];
  };

  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
  }

  /**
   * @brief Setup for each element
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    constexpr int shift = numNodesPerElem * 3;

    int permutation[numNodesPerElem];
    m_finiteElementSpace.getPermutation( permutation );

    localIndex const kf0 = m_elemsToFaces[k][0];
    localIndex const kf1 = m_elemsToFaces[k][1];

    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const kn0 = m_faceToNodes( kf0, a );
      localIndex const kn1 = m_faceToNodes( kf1, a );
      for( int i = 0; i < 3; ++i )
      {
        stack.dispColIndices[a * 3 + i] = m_dofNumber[kn0] + i;
        stack.dispColIndices[shift + a * 3 + i] = m_dofNumber[kn1] + i;
        stack.X[a][i] = m_X[m_faceToNodes( kf0, permutation[a] )][i];
      }
    }

    for( int i = 0; i < 3; ++i )
    {
      stack.unitNormal[i] = m_rotationMatrix( k, i, 0 );
    }

    for( int i = 0; i < 3; ++i )
    {
      stack.bColIndices[i] = m_bDofNumber[kf0] + i;
      stack.bColIndices[3 + i] = m_bDofNumber[kf1] + i;
    }

    stack.pressureRowIndex = m_pressureDofNumber[k] - m_dofRankOffset;
  }

  /**
   * @brief Complete - compute Jacobians dRp/dU and dRp/dB
   */
  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    // Only assemble for open fractures
    if( m_fractureState[k] != fields::contact::FractureState::Open )
    {
      return 0.0;
    }

    // Compute Jacobian for pressure equation w.r.t. displacement:
    // dR_p/du = density * unitNormal^T * Atu = density * area * dAperture/dU
    // This is: dRpdU[j] = sum_i( unitNormal[i] * Atu[i][j] ) * density
    // Use Ri_eq_AjiBj for transpose: R[i] = sum_j(A[j][i] * B[j])
    LvArray::tensorOps::Ri_eq_AjiBj< numUdofs, 3 >( stack.dRpdU, stack.localAtu, stack.unitNormal );
    LvArray::tensorOps::scale< numUdofs >( stack.dRpdU, m_density[k] );

    // Compute Jacobian for pressure equation w.r.t. bubble displacement:
    // dR_p/db = density * unitNormal^T * Atb
    LvArray::tensorOps::Ri_eq_AjiBj< numBdofs, 3 >( stack.dRpdB, stack.localAtb, stack.unitNormal );
    LvArray::tensorOps::scale< numBdofs >( stack.dRpdB, m_density[k] );

    // Assemble into matrix
    localIndex const localRow = LvArray::integerConversion< localIndex >( stack.pressureRowIndex );

    if( localRow >= 0 && localRow < m_matrix.numRows() )
    {
      // Add dRp/dU
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                               stack.dispColIndices,
                                                                               stack.dRpdU,
                                                                               numUdofs );

      // Add dRp/dB
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                               stack.bColIndices,
                                                                               stack.dRpdB,
                                                                               numBdofs );
    }

    return 0.0;
  }

protected:

  /// Pressure DOF numbers
  arrayView1d< globalIndex const > const m_pressureDofNumber;

  /// Fluid density
  arrayView1d< real64 const > const m_density;

  /// Fracture state
  arrayView1d< integer const > const m_fractureState;
};

/// Factory for AssembleFluidMassResidualDerivativeWrtDisplacement kernel
using AssembleFluidMassResidualDerivativeWrtDisplacementFactory =
  finiteElement::InterfaceKernelFactory< AssembleFluidMassResidualDerivativeWrtDisplacement,
                                         arrayView1d< globalIndex const > const,
                                         arrayView1d< globalIndex const > const,
                                         globalIndex const,
                                         CRSMatrixView< real64, globalIndex const > const,
                                         arrayView1d< real64 > const,
                                         real64 const,
                                         arrayView1d< localIndex const > const,
                                         arrayView1d< globalIndex const > const,
                                         arrayView1d< real64 const > const,
                                         arrayView1d< integer const > const >;


/**
 * @brief Kernel to compute d(aperture)/du and d(aperture)/db for each fracture element.
 *
 * This kernel computes the derivative of aperture with respect to displacement and bubble DOFs
 * using the properly computed Atu/Atb matrices. The results are stored for later use in
 * flux derivative computation.
 *
 * d(aperture)/du = (1/area) * unitNormal^T * Atu
 * d(aperture)/db = (1/area) * unitNormal^T * Atb
 *
 * The Atu/Atb matrices encode the sign convention directly: face 0 gets -N*detJ
 * and face 1 gets +N*detJ. This means (1/area) * Atu^T * unitNormal already gives
 * the correct dAperture/dU matching the LM convention (negative for face 0 nodes,
 * positive for face 1 nodes).
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ComputeApertureDerivatives :
  public solidMechanicsConformingContactKernels::ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                                                                FE_TYPE >
{

public:
  /// Alias for the base class
  using Base = solidMechanicsConformingContactKernels::ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                                                                      FE_TYPE >;

  /// Number of nodes per element
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// The number of displacement dofs per element.
  static constexpr int numUdofs = Base::numUdofs;

  /// The number of bubble dofs per element.
  static constexpr int numBdofs = Base::numBdofs;

  using Base::m_elemsToFaces;
  using Base::m_faceToNodes;
  using Base::m_X;
  using Base::m_finiteElementSpace;
  using Base::m_dofNumber;
  using Base::m_bDofNumber;
  using Base::m_dofRankOffset;
  using Base::m_rotationMatrix;
  using Base::m_matrix;
  using Base::m_rhs;

  /**
   * @brief Constructor
   */
  ComputeApertureDerivatives( NodeManager const & nodeManager,
                               EdgeManager const & edgeManager,
                               FaceManager const & faceManager,
                               localIndex const targetRegionIndex,
                               FaceElementSubRegion & elementSubRegion,
                               FE_TYPE const & finiteElementSpace,
                               CONSTITUTIVE_TYPE & inputConstitutiveType,
                               arrayView1d< globalIndex const > const uDofNumber,
                               arrayView1d< globalIndex const > const bDofNumber,
                               globalIndex const rankOffset,
                               CRSMatrixView< real64, globalIndex const > const inputMatrix,
                               arrayView1d< real64 > const inputRhs,
                               real64 const inputDt,
                               arrayView1d< localIndex const > const & faceElementList,
                               arrayView2d< real64 > const & dAperturedU,
                               arrayView2d< real64 > const & dAperturedB ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          uDofNumber,
          bDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt,
          faceElementList ),
    m_area( elementSubRegion.getElementArea().toViewConst() ),
    m_dAperturedU( dAperturedU ),
    m_dAperturedB( dAperturedB )
  {}

  /**
   * @copydoc finiteElement::InterfaceKernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
      unitNormal{},
      localDAperturedU{},
      localDAperturedB{}
    {}

    /// Unit normal from rotation matrix
    real64 unitNormal[3];

    /// Local d(aperture)/dU
    real64 localDAperturedU[numUdofs];

    /// Local d(aperture)/dB
    real64 localDAperturedB[numBdofs];
  };

  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
  }

  /**
   * @brief Setup for each element
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    int permutation[numNodesPerElem];
    m_finiteElementSpace.getPermutation( permutation );

    localIndex const kf0 = m_elemsToFaces[k][0];

    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      for( int i = 0; i < 3; ++i )
      {
        stack.X[a][i] = m_X[m_faceToNodes( kf0, permutation[a] )][i];
      }
    }

    for( int i = 0; i < 3; ++i )
    {
      stack.unitNormal[i] = m_rotationMatrix( k, i, 0 );
    }
  }

  /**
   * @brief Complete - compute and store d(aperture)/dU and d(aperture)/dB
   */
  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    real64 const invArea = 1.0 / m_area[k];

    // Compute d(aperture)/dU = (1/area) * unitNormal^T * Atu
    // This gives: dAperturedU[j] = (1/area) * sum_i(unitNormal[i] * Atu[i][j])
    // Use Ri_eq_AjiBj for transpose: R[i] = sum_j(A[j][i] * B[j])
    LvArray::tensorOps::Ri_eq_AjiBj< numUdofs, 3 >( stack.localDAperturedU, stack.localAtu, stack.unitNormal );
    LvArray::tensorOps::scale< numUdofs >( stack.localDAperturedU, invArea );

    // Compute d(aperture)/dB = (1/area) * unitNormal^T * Atb
    LvArray::tensorOps::Ri_eq_AjiBj< numBdofs, 3 >( stack.localDAperturedB, stack.localAtb, stack.unitNormal );
    LvArray::tensorOps::scale< numBdofs >( stack.localDAperturedB, invArea );

    // Store the results
    for( localIndex i = 0; i < numUdofs; ++i )
    {
      m_dAperturedU( k, i ) = stack.localDAperturedU[i];
    }
    for( localIndex i = 0; i < numBdofs; ++i )
    {
      m_dAperturedB( k, i ) = stack.localDAperturedB[i];
    }

    return 0.0;
  }

protected:

  /// Element area
  arrayView1d< real64 const > const m_area;

  /// Storage for d(aperture)/dU
  arrayView2d< real64 > const m_dAperturedU;

  /// Storage for d(aperture)/dB
  arrayView2d< real64 > const m_dAperturedB;
};

/// Factory for ComputeApertureDerivatives kernel
using ComputeApertureDerivativesFactory =
  finiteElement::InterfaceKernelFactory< ComputeApertureDerivatives,
                                         arrayView1d< globalIndex const > const,
                                         arrayView1d< globalIndex const > const,
                                         globalIndex const,
                                         CRSMatrixView< real64, globalIndex const > const,
                                         arrayView1d< real64 > const,
                                         real64 const,
                                         arrayView1d< localIndex const > const,
                                         arrayView2d< real64 > const,
                                         arrayView2d< real64 > const >;

} // namespace poromechanicsALMKernels


namespace poromechanicsMatrixBubbleKernels
{

/**
 * @brief Kernel to compute the contribution of matrix cell pressure on bubble DOFs
 *        with full Jacobian for fully-implicit poromechanics coupling.
 *
 * This kernel is based on SolidMechanicsPressureFaceBubbleKernels but adds the Jacobian
 * term dR_b/dP_matrix needed for fully implicit coupling.
 *
 * Physics:
 *   R_b = -∫ B_b^T * (biot * p * I) * detJ dΩ
 *   dR_b/dP = -∫ B_b^T * (biot * I) * detJ dΩ
 *
 * where B_b is the strain-displacement matrix for bubble functions.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class MatrixPressureBubbleKernels :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:

  /// Alias for the base class
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  /// Number of nodes per element
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// Compile time value for the number of faces per element.
  static constexpr int numFacesPerElem = FE_TYPE::numFaces;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;
  using Base::m_constitutiveUpdate;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_dt;

  /**
   * @brief Constructor
   */
  MatrixPressureBubbleKernels( NodeManager const & nodeManager,
                               EdgeManager const & edgeManager,
                               FaceManager const & faceManager,
                               localIndex const targetRegionIndex,
                               SUBREGION_TYPE const & elementSubRegion,
                               FE_TYPE const & finiteElementSpace,
                               CONSTITUTIVE_TYPE & inputConstitutiveType,
                               arrayView1d< globalIndex const > const uDofNumber,
                               arrayView1d< globalIndex const > const bDofNumber,
                               globalIndex const rankOffset,
                               CRSMatrixView< real64, globalIndex const > const inputMatrix,
                               arrayView1d< real64 > const inputRhs,
                               real64 const inputDt,
                               string const pressureDofKey,
                               string const fluidModelKey ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          uDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt ),
    m_X( nodeManager.referencePosition()),
    m_bDofNumber( bDofNumber ),
    m_pDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( pressureDofKey ) ),
    m_bubbleElems( elementSubRegion.bubbleElementsList() ),
    m_elemsToFaces( elementSubRegion.faceElementsList() ),
    m_pressure( elementSubRegion.template getField< fields::flow::pressure >().toViewConst() ),
    m_incrBubbleDisp( faceManager.getField< fields::contact::incrementalBubbleDisplacement >().toViewConst() ),
    m_fluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >(
                      elementSubRegion.template getReference< string >( fluidModelKey ) ).density() )
  {}

  //***************************************************************************

  /**
   * @copydoc finiteElement::ImplicitKernelBase::StackVariables
   */
  struct StackVariables
  {
public:

    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;

    /// The number of bubble dofs per element (all faces).
    static constexpr int numBubbleUdofs = numFacesPerElem * 3;

    /// The number of pressure dofs per element.
    static constexpr int numPdofs = 1;

    /**
     * Default constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      bEqnRowIndices{},
      bColIndices{},
      pColIndex( 0 ),
      pRowIndex( 0 ),
      localRb{},
      localdRbdP{},
      localdRpdB{},
      X{ {} },
      pLocal{}
    {}

    /// C-array storage for the element local row degrees of freedom (bubble).
    globalIndex bEqnRowIndices[3];

    /// C-array storage for the (global) bubble column degrees of freedom (for A_pb).
    globalIndex bColIndices[3];

    /// Column index for pressure DOF (for A_bp)
    globalIndex pColIndex;

    /// Row index (local) for pressure DOF (for A_pb)
    globalIndex pRowIndex;

    /// C-array storage for the element local Rb residual vector.
    real64 localRb[numBubbleUdofs];

    /// C-array storage for the element local dRb/dP Jacobian.
    real64 localdRbdP[numBubbleUdofs];

    /// C-array storage for the element local dRp/dB Jacobian (A_pb bulk term).
    real64 localdRpdB[numBubbleUdofs];

    /// local nodal coordinates
    real64 X[ numNodesPerElem ][ 3 ];

    /// local pressure
    real64 pLocal[numPdofs];

  };

  //***************************************************************************

  /**
   * @copydoc ::geos::finiteElement::KernelBase::kernelLaunch
   *
   * @detail it uses the kernelLaunch interface of KernelBase but it only launches the kernel
   * on the set of elements that have bubble dof within the subregion.
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    GEOS_UNUSED_VAR( numElems );

    // Define a RAJA reduction variable to get the maximum residual contribution.
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    forAll< POLICY >( kernelComponent.m_bubbleElems.size(),
                      [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( i, stack );
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( i, q, stack );
      }
      maxResidual.max( kernelComponent.complete( i, stack ) );
    } );

    return maxResidual.get();
  }

  /**
   * @brief Setup for each element
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const kk,
              StackVariables & stack ) const
  {
    localIndex const k = m_bubbleElems[kk];

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
        stack.X[ a ][ i ] = m_X[ localNodeIndex ][ i ];
      }
    }

    localIndex const localFaceIndex = m_elemsToFaces[kk][0];

    for( int i=0; i<3; ++i )
    {
      stack.bEqnRowIndices[i] = m_bDofNumber[localFaceIndex] + i - m_dofRankOffset;
      stack.bColIndices[i] = m_bDofNumber[localFaceIndex] + i;  // global column for A_pb
    }

    stack.pColIndex = m_pDofNumber[k];
    stack.pRowIndex = m_pDofNumber[k] - m_dofRankOffset;
    stack.pLocal[0] = m_pressure( k );
  }

  /**
   * @brief Quadrature point kernel - compute contributions at each quadrature point
   */
  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const kk,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    localIndex const k = m_bubbleElems[kk];

    constexpr int nBubbleUdof = numFacesPerElem*3;

    real64 dBubbleNdX[ numFacesPerElem ][ 3 ];
    // Initialize to zero (needed because calcGradFaceBubbleN may be a placeholder in some FE)
    LvArray::tensorOps::fill< numFacesPerElem, 3 >( dBubbleNdX, 0 );

    real64 const detJ = m_finiteElementSpace.calcGradFaceBubbleN( q, stack.X, dBubbleNdX );

    real64 biotCoefficient;
    m_constitutiveUpdate.getBiotCoefficient( k, biotCoefficient );

    // Assemble the strain-displacement matrix for bubble functions: B_b
    real64 strainBubbleMatrix[6][nBubbleUdof];
    solidMechanicsConformingContactKernelsHelper::
      assembleStrainOperator< 6, nBubbleUdof, numFacesPerElem >( strainBubbleMatrix, dBubbleNdX );

    // Compute biotPressure = -biot * p * I (Voigt notation: [1,1,1,0,0,0])
    real64 biotPressure[6] = {0};
    LvArray::tensorOps::symAddIdentity< 3 >( biotPressure, -biotCoefficient * stack.pLocal[0] );

    // Compute biotIdentity = -biot * I (for Jacobian)
    real64 biotIdentity[6] = {0};
    LvArray::tensorOps::symAddIdentity< 3 >( biotIdentity, -biotCoefficient );

    // Residual contribution: R_b += B_b^T * biotPressure * detJ
    // (note: scaledAdd with -detJ because biotPressure already has negative sign)
    real64 Rb_gauss[nBubbleUdof];
    LvArray::tensorOps::Ri_eq_AjiBj< nBubbleUdof, 6 >( Rb_gauss, strainBubbleMatrix, biotPressure );
    LvArray::tensorOps::scaledAdd< nBubbleUdof >( stack.localRb, Rb_gauss, -detJ );

    // Jacobian contribution: dR_b/dP += B_b^T * biotIdentity * detJ
    real64 dRbdP_gauss[nBubbleUdof];
    LvArray::tensorOps::Ri_eq_AjiBj< nBubbleUdof, 6 >( dRbdP_gauss, strainBubbleMatrix, biotIdentity );
    LvArray::tensorOps::scaledAdd< nBubbleUdof >( stack.localdRbdP, dRbdP_gauss, -detJ );

    // ---- A_pb^Omega : transpose Biot term for the fluid-mass equation ----
    // The bubble mode contributes to the cell volumetric strain, hence to the
    // fluid mass storage.  With dPorosity_dVolStrain = biot (BiotPorosity),
    //   dFluidMassIncrement_dVolStrainIncrement = biot * fluidDensity,
    // so dR_p/db = fluidDensity * dR_b/dp (same geometric factor biot*div(beta),
    // times the fluid density that puts it in mass units). We reuse dRbdP_gauss.
    real64 const rhof = m_fluidDensity( k, 0 );
    LvArray::tensorOps::scaledAdd< nBubbleUdof >( stack.localdRpdB, dRbdP_gauss, -detJ * rhof );
  }

  /**
   * @brief Complete - assemble residual and Jacobian into global system
   */
  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const kk,
                   StackVariables & stack ) const
  {
    localIndex const parentFaceIndex = m_elemsToFaces[kk][1];

    // Extract only the components corresponding to the parent face
    // on which the bubble function was applied.
    real64 localRb[3];
    real64 localdRbdP[3];
    real64 localdRpdB[3];
    for( localIndex i = 0; i < 3; ++i )
    {
      localRb[i] = stack.localRb[parentFaceIndex*3+i];
      localdRbdP[i] = stack.localdRbdP[parentFaceIndex*3+i];
      localdRpdB[i] = stack.localdRpdB[parentFaceIndex*3+i];
    }

    // ---- A_bp : bubble-momentum row (dR_b/dP) ----
    for( localIndex i=0; i < 3; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.bEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() )
        continue;

      // Add residual
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], localRb[i] );

      // Add Jacobian dR_b/dP
      m_matrix.template addToRow< parallelDeviceAtomic >( dof,
                                                          &stack.pColIndex,
                                                          &localdRbdP[i],
                                                          1 );
    }

    // ---- A_pb^Omega : fluid-mass row (dR_p/db) + storage residual ----
    localIndex const pRow = LvArray::integerConversion< localIndex >( stack.pRowIndex );
    if( pRow >= 0 && pRow < m_matrix.numRows() )
    {
      // Residual: R_p += fluidDensity * biot * tr(grad beta) * delta_b (bubble storage),
      // using the incremental bubble displacement of the bubble face.
      localIndex const bFace = m_elemsToFaces[kk][0];
      real64 rp = 0.0;
      for( localIndex i = 0; i < 3; ++i )
      {
        rp += localdRpdB[i] * m_incrBubbleDisp[bFace][i];
      }
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[pRow], rp );

      // Jacobian dR_p/db (three bubble columns of the parent face)
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( pRow,
                                                                              stack.bColIndices,
                                                                              localdRpdB,
                                                                              3 );
    }

    return 0.0;
  }

protected:

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The global degree of freedom number of bubble
  arrayView1d< globalIndex const > const m_bDofNumber;

  /// The global degree of freedom number of pressure
  arrayView1d< globalIndex const > const m_pDofNumber;

  /// The array containing the list of bubble elements.
  arrayView1d< localIndex const > const m_bubbleElems;

  /// The array containing the element to bubble face map.
  arrayView2d< localIndex const > const m_elemsToFaces;

  /// The array containing the pressure of each element.
  arrayView1d< real64 const > const m_pressure;

  /// Incremental bubble displacement (face field) -- used for the fluid-mass storage residual.
  arrayView2d< real64 const > const m_incrBubbleDisp;

  /// Fluid density [elem][q].
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_fluidDensity;

};

/// The factory used to construct MatrixPressureBubbleKernels.
using MatrixPressureBubbleFactory = finiteElement::KernelFactory< MatrixPressureBubbleKernels,
                                                                   arrayView1d< globalIndex const > const,
                                                                   arrayView1d< globalIndex const > const,
                                                                   globalIndex const,
                                                                   CRSMatrixView< real64, globalIndex const > const,
                                                                   arrayView1d< real64 > const,
                                                                   real64 const,
                                                                   string const,
                                                                   string const >;

} // namespace poromechanicsMatrixBubbleKernels

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALMKERNELS_HPP_ */
