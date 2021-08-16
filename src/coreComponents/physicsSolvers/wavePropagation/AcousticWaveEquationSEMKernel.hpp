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
 * @file AcousticWaveEquationSEMKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"


namespace geosx
{

/// Namespace to contain the acoustic wave kernels.
namespace AcousticWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{

  /**
   * @brief Convert a mesh element point coordinate into a coorinate on the reference element
   * @param coords coordinate of the point
   * @param coordsOnRefElem to contain the coordinate computed in the reference element
   * @param elementIndex index of the element containing the coords
   * @param faceNodeIndices array of face of the element
   * @param elemsToNodes map to obtaint global nodes from element index
   * @param X array of mesh nodes coordinates
   * @return true if coords is inside the element num index
   */
  template< typename FE_TYPE >
  GEOSX_HOST_DEVICE
  static bool
  computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                        real64 (& coordsOnRefElem)[3],
                                        localIndex const & elementIndex,
                                        arrayView1d< arrayView1d< localIndex const > const > const faceNodeIndices,
                                        arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X )
  {
    if( computationalGeometry::isPointInsidePolyhedron( X, faceNodeIndices, coords ) )
    {
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      real64 xLocal[numNodesPerElem][3];
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = X( elemsToNodes( elementIndex, a ), i );
        }
      }

      // coordsOnRefElem = invJ*(coords-coordsNode_0)
      localIndex const q = 0;

      real64 invJ[3][3] = {{0}};
      FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

      real64 coordsRef[3] = {0};
      for( localIndex i = 0; i < 3; ++i )
      {
        coordsRef[i] = coords[i] - xLocal[q][i];
      }

      for( localIndex i = 0; i < 3; ++i )
      {
        // Init at (-1,-1,-1) as the origin of the referential elem
        coordsOnRefElem[i] = -1.0;
        for( localIndex j = 0; j < 3; ++j )
        {
          coordsOnRefElem[i] += invJ[i][j]*coordsRef[j];
        }
      }
      return true;
    }
    else
    {
      return false;
    }
  }

  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          localIndex const numNodesPerElem,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView1d< arrayView1d< localIndex > const > const faceNodeIndices,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsLocal,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstants,
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants )
  {

    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      CellBlock::getAllFaceNodes< ElementType::Hexahedron >( k,
                                                             faceNodeIndices,
                                                             elemsToNodes );

      /// loop over all the source that haven't been found yet
      for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
      {
        if( sourceIsLocal[isrc] == 0 )
        {
          real64 const coords[3] = { sourceCoordinates[isrc][0],
                                     sourceCoordinates[isrc][1],
                                     sourceCoordinates[isrc][2] };

          real64 coordsOnRefElem[3]{};
          bool const sourceFound = computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                                    coordsOnRefElem,
                                                                                    k,
                                                                                    faceNodeIndices.toNestedViewConst(),
                                                                                    elemsToNodes,
                                                                                    X );
          if( sourceFound )
          {
            sourceIsLocal[isrc] = 1;
            real64 Ntest[8];
            finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

            for( localIndex a=0; a< numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes[k][a];
              sourceConstants[isrc][a] = Ntest[a];
            }
          }
        }
      } // End loop over all source

      /// loop over all the receiver that haven't been found yet
      for( localIndex ircv = 0; ircv < receiverCoordinates.size( 0 ); ++ircv )
      {
        if( receiverIsLocal[ircv] == 0 )
        {
          real64 const coords[3] = { receiverCoordinates[ircv][0],
                                     receiverCoordinates[ircv][1],
                                     receiverCoordinates[ircv][2] };

          real64 coordsOnRefElem[3]{};
          bool const receiverFound = computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                                      coordsOnRefElem,
                                                                                      k,
                                                                                      faceNodeIndices.toNestedViewConst(),
                                                                                      elemsToNodes,
                                                                                      X );
          if( receiverFound )
          {
            receiverIsLocal[ircv] = 1;

            real64 Ntest[8];
            finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

            for( localIndex a=0; a< numNodesPerElem; ++a )
            {
              receiverNodeIds[ircv][a] = elemsToNodes[k][a];
              receiverConstants[ircv][a] = Ntest[a];
            }
          }
        }
      } // End loop over receiver
    } );
  }
};

template< typename FE_TYPE >
struct MassAndDampingMatrixKernel
{

  MassAndDampingMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  template< typename EXEC_POLICY >
  void
  launch( localIndex const size,
          localIndex const numFacesPerElem,
          localIndex const numNodesPerFace,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView1d< integer const > const facesDomainBoundaryIndicator,
          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
          arrayView2d< real64 const > const faceNormal,
          arrayView1d< real64 const > const velocity,
          arrayView1d< real64 > const mass,
          arrayView1d< real64 > const damping )
  {
    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real64 const invC2 = 1.0 / ( velocity[k] * velocity[k] );
      real64 xLocal[ numNodesPerElem ][ 3 ];
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = X( elemsToNodes( k, a ), i );
        }
      }

      real64 N[ numNodesPerElem ];
      real64 gradN[ numNodesPerElem ][ 3 ];

      for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
      {
        FE_TYPE::calcN( q, N );
        real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

        for( localIndex a = 0; a < numNodesPerElem; ++a )
        {
          // RAJA atomic
          mass[elemsToNodes[k][a]] +=  invC2 * detJ * N[a];
        }
      }

      real64 const alpha = 1.0 / velocity[k];

      for( localIndex kfe = 0; kfe < numFacesPerElem; ++kfe )
      {
        localIndex const iface = elemsToFaces[k][kfe];

        /// Face on the domain boundary and not on free surface
        if( facesDomainBoundaryIndicator[iface] == 1 && freeSurfaceFaceIndicator[iface] != 1 )
        {
          for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );
            real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

            real64 invJ[3][3] = {{ 0 }};
            FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

            for( localIndex a = 0; a < numNodesPerFace; ++a )
            {
              /// compute ds=||detJ*invJ*normalFace_{kfe}||
              real64 tmp[3] = { 0 };
              real64 ds = 0.0;
              for( localIndex i = 0; i < 3; ++i )
              {
                for( localIndex j = 0; j < 3; ++j )
                {
                  tmp[i] += invJ[j][i] * faceNormal[iface][j];
                }
                ds += tmp[i] * tmp[i];
              }
              ds = sqrt( ds );

              localIndex const inode = facesToNodes[iface][a];
              damping[inode] += alpha * detJ * ds * N[a];
            }
          }
        }
      }
    } ); // end loop over element
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};


/**
 * @brief Implements kernels for solving the acoustic wave equations
 *   explicit central FD method and SEM
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### AcousticWaveEquationSEMKernel Description
 * Implements the KernelBase interface functions required for solving
 * the acoustic wave equations using the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1`.
 */


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitAcousticSEM : public finiteElement::KernelBase< SUBREGION_TYPE,
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

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

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
   *   elements to be processed during this kernel launch.
   */
  ExplicitAcousticSEM( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       real64 const dt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_X( nodeManager.referencePosition() ),
    m_p_n( nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >() ),
    m_stiffnessVector( nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector >() ),
    m_dt( dt )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
    GEOSX_UNUSED_VAR( targetRegionIndex );
  }



  //*****************************************************************************
  /**
   * @copydoc geosx::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitAcousticSEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOSX_HOST_DEVICE
    StackVariables():
      xLocal()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
  };
  //***************************************************************************


  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    /// numDofPerTrialSupportPoint = 1
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X[ nodeIndex ][ i ];
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitAcousticSEM Description
   * Calculates stiffness vector
   *
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 gradN[ numNodesPerElem ][ 3 ];

    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, gradN );

    for( localIndex i=0; i<numNodesPerElem; ++i )
    {
      for( localIndex j=0; j<numNodesPerElem; ++j )
      {
        real64 const Rh_ij = detJ * LvArray::tensorOps::AiBi< 3 >( gradN[ i ], gradN[ j ] );

        m_stiffnessVector[m_elemsToNodes[k][i]] += Rh_ij*m_p_n[m_elemsToNodes[k][j]];

      }
    }
  }


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal pressure array.
  arrayView1d< real64 const > const m_p_n;

  /// The array containing the product of the stiffness matrix and the nodal pressure.
  arrayView1d< real64 > const m_stiffnessVector;

  /// The time increment for this time integration step.
  real64 const m_dt;


};


/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitAcousticSEMFactory = finiteElement::KernelFactory< ExplicitAcousticSEM,
                                                                 real64 >;

} // namespace AcousticWaveEquationSEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_
