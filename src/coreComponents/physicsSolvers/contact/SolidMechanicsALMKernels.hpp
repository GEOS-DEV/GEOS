/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsALMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_

#include "finiteElement/kernelInterface/InterfaceKernelBase.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ALMKernelsBase :
  public finiteElement::InterfaceKernelBase< SUBREGION_TYPE,
                                             CONSTITUTIVE_TYPE,
                                             FE_TYPE,
                                             3, 3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::InterfaceKernelBase< SUBREGION_TYPE,
                                                   CONSTITUTIVE_TYPE,
                                                   FE_TYPE, 
                                                   3, 3 >;

  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  ALMKernelsBase( NodeManager const & nodeManager,
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
                  real64 const inputDt, 
                  arrayView1d< localIndex const > const & faceElementsList ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt ),
    m_faceElementsList(faceElementsList)
{}

  struct StackVariables  
  {
    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3 * 2;

  public:
  
    GEOS_HOST_DEVICE
    StackVariables():
      dispEqnRowIndices{ 0 },
      dispColIndices{ 0 },
      localRu{ 0.0 },
      localAutAtu{ { 0.0 } },
      tLocal()
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local AutAtu matrix.
    real64 localAutAtu[numUdofs][numUdofs];

    /// Stack storage for the element local lagrange multiplier vector
    real64 tLocal[3];
  };

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

    forAll< POLICY >( kernelComponent.m_faceElementsList.size(),
                      [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      std::cout << "# QuadPoints: " << numQuadraturePointsPerElem << std::endl;
      std::cout << "# NodesperElem: " << numNodesPerElem << std::endl;

      localIndex k = kernelComponent.m_faceElementsList[i];
      typename KERNEL_TYPE::StackVariables stack;

      std::cout << i << " " << k << std::endl;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
      maxResidual.max( kernelComponent.complete( k, stack ) );
    } );

    return maxResidual.get();
  }
  //END_kernelLauncher

  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR(k, stack);
  }

  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR(k, q, stack);
  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k, stack );

    return 0.0;
  }

protected:

  arrayView1d< localIndex const > const m_faceElementsList;

  arrayView3d< real64 const > const m_rotationMatrix;

};

using ALMFactory = finiteElement::KernelFactory< ALMKernelsBase,
                                                 arrayView1d< globalIndex const > const,
                                                 globalIndex const,
                                                 CRSMatrixView< real64, globalIndex const > const,
                                                 arrayView1d< real64 > const,
                                                 real64 const, 
                                                 arrayView1d< localIndex const > const >;


struct ComputeRotationMatricesKernel
{
  template< typename POLICY >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const > const & faceNormal,
          ArrayOfArraysView< localIndex const > const & elemsToFaces,
          arrayView3d< real64 > const & rotationMatrix )
  {

    GEOS_UNUSED_VAR(faceNormal, elemsToFaces, rotationMatrix);

    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      localIndex const & f0 = elemsToFaces[k][0];
      localIndex const & f1 = elemsToFaces[k][1];

      //stackArray1d< real64, 3 > Nbar( 3 );
      real64 Nbar[3];
      Nbar[0] = faceNormal[f0][0] - faceNormal[f1][0];
      Nbar[1] = faceNormal[f0][1] - faceNormal[f1][1];
      Nbar[2] = faceNormal[f0][2] - faceNormal[f1][2];

      LvArray::tensorOps::normalize< 3 >( Nbar );
      computationalGeometry::RotationMatrix_3D( Nbar, rotationMatrix[k] );

    } );
  }

};

} // namespace SolidMechanicsALMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_ */