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
 * @file MFLaplaceFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_MFLaplaceFEMKernels_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_MFLaplaceFEMKernels_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
// #include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"


namespace geosx
{

/**
 * @brief Implements kernels for solving the equations of motion using the
 *   explicit Newmark method under the small strain assumption.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### Explicit Small Strain Description
 * Implements the KernelBase interface functions required for explicit time
 * integration of the equations of motion using the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the interface for KernelBase is used, but
 * MFLaplaceFEMKernels only conforms to the interface set by KernelBase, and
 * does not inherit from KernelBase.
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3`.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class MFLaplaceFEMKernels : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                              CONSTITUTIVE_TYPE,
                                                              FE_TYPE,
                                                              1,
                                                              1 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          1,
                                          1 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

//*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   * @param elementListName The name of the entry that holds the list of
   *   elements to be processed during this kernel launch.
   */
  MFLaplaceFEMKernels( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       arrayView1d< real64 const > const inputSrc,
                       arrayView1d< real64 > const inputDst ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition()),
    m_input( inputSrc ),
    m_res( inputDst )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
  }

  //*****************************************************************************
  /**
   * @copydoc geosx::finiteElement::KernelBase::StackVariables
   *
   * ### MFLaplaceFEMKernels Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOSX_HOST_DEVICE
    StackVariables():
      fLocal{ 0.0 },
      varLocal{ 0.0 },
      xLocal()
    {}

    /// C-array stack storage for the element local force
    real64 fLocal[ numNodesPerElem ];

    /// C-array stack storage for element local primary variable values.
    real64 varLocal[ numNodesPerElem ];

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif
  };
  //***************************************************************************


  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    #pragma unroll
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      stack.varLocal[ a ] = m_input( nodeIndex );
      #pragma unroll
      for( int i=0; i<numDofPerTrialSupportPoint; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X( nodeIndex, i );
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### MFLaplaceFEMKernels Description
   * Calculates the shape function derivatives, and the strain tensor. Then
   * calls the constitutive update, and also performs the integration of
   * the stress divergence, rather than using the dedicated component function
   * to allow for some variable reuse.
   */
  GEOSX_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::complete
   *
   * ### MFLaplaceFEMKernels Description
   * Performs the distribution of the nodal force out to the rank local arrays.
   */
  GEOSX_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables const & stack ) const;











  GEOSX_HOST_DEVICE
  void setup( localIndex const k,
              real64 (&xLocal) [ numNodesPerElem ][ 3 ],
              real64 (&varLocal) [ numNodesPerElem ] ) const
{
  #pragma unroll
  for( localIndex a=0; a< numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    varLocal[ a ] = m_input( nodeIndex );
    #pragma unroll
    for( int i=0; i<3; ++i )
    {
      xLocal[ a ][ i ] = m_X( nodeIndex, i );
    }
  }
}

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### MFLaplaceFEMKernels Description
   * Calculates the shape function derivatives, and the strain tensor. Then
   * calls the constitutive update, and also performs the integration of
   * the stress divergence, rather than using the dedicated component function
   * to allow for some variable reuse.
   */
  GEOSX_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const qa,
                              localIndex const qb,
                              localIndex const qc,
                              real64 const (&xLocal)[ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                              real64 const (&varLocal)[ numNodesPerElem ][ numDofPerTrialSupportPoint ],
                              real64 (&fLocal) [ numNodesPerElem ][ numDofPerTrialSupportPoint ] ) const;


  template< int qa, int qb, int qc >
  GEOSX_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              real64 const (&xLocal)[ numNodesPerElem ][ 3 ],
                              real64 const (&varLocal)[ numNodesPerElem ],
                              real64 (&fLocal) [ numNodesPerElem ] ) const
{
  real64 invJ[3][3] = {{0}};
  real64 parentGradVar[3] = {0};

  FE_TYPE::template parentGradient2< qa, qb, qc >( xLocal, varLocal, invJ, parentGradVar);

//  FE_TYPE::template parentGradient< qa, qb, qc >( xLocal, invJ);
  real64 const detJ = LvArray::tensorOps::invert< 3 >( invJ );

//  FE_TYPE::template parentGradient< qa, qb, qc >( varLocal, parentGradVar);
  real64 gradVar[3] = {0};

  #pragma unroll
  for( int j = 0; j < 3; ++j )
  {
    #pragma unroll
    for( int kk = 0; kk < 3; ++kk )
    {
      gradVar[j] = gradVar[j] + parentGradVar[kk] * invJ[kk][j];
    }
    gradVar[j] *= -detJ;
  }
  // real64 strain[6] = {0};
  // strain[0] = gradVar[0][0];
  // strain[1] = gradVar[1][1];
  // strain[2] = gradVar[2][2];
  // strain[3] = gradVar[2][1] + gradVar[1][2];
  // strain[4] = gradVar[2][0] + gradVar[0][2];
  // strain[5] = gradVar[1][0] + gradVar[0][1];


  // real64 stressLocal[ 6 ] = {0};
  // m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, qa+2*qb+4*qc, strain, stressLocal );

  // for( localIndex c = 0; c < 6; ++c )
  // {
  //   stressLocal[ c ] *= -detJ;
  // }

  FE_TYPE::template plusGradNajAij< qa, qb, qc >( invJ, gradVar, fLocal );

}


  /**
   * @copydoc geosx::finiteElement::KernelBase::complete
   *
   * ### MFLaplaceFEMKernels Description
   * Performs the distribution of the nodal force out to the rank local arrays.
   */
  GEOSX_HOST_DEVICE
  real64 complete( localIndex const k,
                   real64 const (&fLocal) [ numNodesPerElem ] ) const
{
  for( localIndex a = 0; a < numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    RAJA::atomicAdd< parallelDeviceAtomic >( &m_res( nodeIndex ), fLocal[ a ] );
  }
  return 0;
}
















  /**
   * @copydoc geosx::finiteElement::KernelBase::kernelLaunch
   *
   * ### MFLaplaceFEMKernels Description
   * Copy of the KernelBase::kernelLaunch function without the exclusion of ghost
   * elements.
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( numElems,
                      [=] GEOSX_DEVICE ( localIndex const k )
    {
      real64 fLocal[ KERNEL_TYPE::numNodesPerElem ] = {0};
      real64 varLocal[ KERNEL_TYPE::numNodesPerElem ];
      real64 xLocal[ KERNEL_TYPE::numNodesPerElem ][ 3 ];

      kernelComponent.setup( k, xLocal, varLocal );
      kernelComponent.template quadraturePointKernel<0, 0, 0>( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel<0, 0, 1>( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel<0, 1, 0>( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel<0, 1, 1>( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel<1, 0, 0>( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel<1, 0, 1>( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel<1, 1, 0>( k, xLocal, varLocal, fLocal );
      kernelComponent.template quadraturePointKernel<1, 1, 1>( k, xLocal, varLocal, fLocal );
      kernelComponent.complete( k, fLocal );

    } );
    return 0;
  }


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal displacement array.
  arrayView1d< real64 const > const m_input;

  /// The array containing the nodal acceleration array, which is used to store
  /// the force.
  arrayView1d< real64 > const m_res;

  /// The list of elements to process for the kernel launch.


};



/// The factory used to construct a MFLaplaceFEMKernels kernel.
using MFLaplaceFEMKernelsFactory = finiteElement::KernelFactory< MFLaplaceFEMKernels,
                                                                 arrayView1d< real64 const > const,
                                                                 arrayView1d< real64 > const >;



} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_MFLaplaceFEMKernels_HPP_
