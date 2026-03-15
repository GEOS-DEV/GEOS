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
 * @file SolidMechanicsConformingPressureContributionKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSPRESSURECONTRIBUTIONKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSPRESSURECONTRIBUTIONKERNELS_HPP_

#include "SolidMechanicsConformingContactKernelsBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp" // needed to register pressure(_n)

namespace geos
{

namespace solidMechanicsConformingContactKernels
{

/**
 * @brief Implements kernels for Conforming Contact.
 * @copydoc geos::finiteElement::InterfaceKernelBase
 *
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          typename MATRIX_VIEW = DefaultGlobalMatrixView >
class AssemblePressureContribution :
  public ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                       FE_TYPE,
                                       MATRIX_VIEW >
{

public:
  /// Alias for the base class;
  using Base = ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                             FE_TYPE,
                                             MATRIX_VIEW >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
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
   * @copydoc geos::finiteElement::InterfaceKernelBase::InterfaceKernelBase
   */
  AssemblePressureContribution( NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager const & faceManager,
                                localIndex const targetRegionIndex,
                                FaceElementSubRegion & elementSubRegion,
                                FE_TYPE const & finiteElementSpace,
                                CONSTITUTIVE_TYPE & inputConstitutiveType,
                                arrayView1d< globalIndex const > const uDofNumber,
                                arrayView1d< globalIndex const > const bDofNumber,
                                globalIndex const rankOffset,
                                MATRIX_VIEW const inputMatrix,
                                arrayView1d< real64 > const inputRhs,
                                real64 const inputDt,
                                arrayView1d< localIndex const > const & faceElementList ):
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
    m_pressure( elementSubRegion.getField< fields::flow::pressure >().toView() )
  {}

  //***************************************************************************
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
                                       unitNormal{},
                                       localRu{},
                                       localRb{}
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex bEqnRowIndices[numBdofs];

    /// Stack storage for the element local pressure contribution matrix.
    real64 unitNormal[3];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local Rb residual vector.
    real64 localRb[numBdofs];
  };

  //***************************************************************************

  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
  }
  //END_kernelLauncher

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geos::finiteElement::InterfaceKernelBase::setup
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
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const kn0 = m_faceToNodes( kf0, a );
      localIndex const kn1 = m_faceToNodes( kf1, a );
      for( int i=0; i<3; ++i )
      {
        stack.dispEqnRowIndices[a*3+i] = m_dofNumber[kn0]+i-m_dofRankOffset;
        stack.dispEqnRowIndices[shift + a*3+i] = m_dofNumber[kn1]+i-m_dofRankOffset;
        stack.X[ a ][ i ] = m_X[ m_faceToNodes( kf0, permutation[ a ] ) ][ i ];
      }
    }

    for( int i=0; i<3; ++i )
    {
      stack.unitNormal[ i ] = m_rotationMatrix( k, i, 0 );
    }

    for( int i=0; i<3; ++i )
    {
      // need to grab the index.
      stack.bEqnRowIndices[i]   = m_bDofNumber[kf0] + i - m_dofRankOffset;
      stack.bEqnRowIndices[3+i] = m_bDofNumber[kf1] + i - m_dofRankOffset;
    }

  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {

    // transp(Atu)*unitNormal
    LvArray::tensorOps::Ri_eq_AjiBj< numUdofs, 3 >( stack.localRu, stack.localAtu, stack.unitNormal );

    // transp(Atb)*unitNormal
    LvArray::tensorOps::Ri_eq_AjiBj< numBdofs, 3 >( stack.localRb, stack.localAtb, stack.unitNormal );

    // Compute the local residuals
    LvArray::tensorOps::scale< numUdofs >( stack.localRu, m_pressure[k] );

    LvArray::tensorOps::scale< numBdofs >( stack.localRb, m_pressure[k] );

    for( localIndex i=0; i < numUdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRu[i] );
    }

    for( localIndex i=0; i < numBdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.bEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRb[i] );
    }

    return 0.0;
  }

protected:

  arrayView1d< real64 const > const m_pressure;

};

using AssemblePressureContributionFactory =
  finiteElement::InterfaceKernelFactory< AssemblePressureContribution,
                                         arrayView1d< globalIndex const > const,
                                         arrayView1d< globalIndex const > const,
                                         globalIndex const,
                                         DefaultGlobalMatrixView const,
                                         arrayView1d< real64 > const,
                                         real64 const,
                                         arrayView1d< localIndex const > const >;

} // namespace SolidMechanicsConformingContactKernels

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSPRESSURECONTRIBUTIONKERNELS_HPP_ */
