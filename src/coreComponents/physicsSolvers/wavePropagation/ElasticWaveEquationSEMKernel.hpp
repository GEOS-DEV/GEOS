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
 * @file ElasticWaveEquationSEMKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEMKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"


namespace geosx
{

/// Namespace to contain the elastic wave kernels.
namespace elasticWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{

  /**
   * @brief Convert a mesh element point coordinate into a coordinate on the reference element
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
                                        real64 const (&elemCenter)[3],
                                        arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const elemsToNodes,
                                        arraySlice1d< localIndex const > const elemsToFaces,
                                        ArrayOfArraysView< localIndex const > const & facesToNodes,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                                        real64 (& coordsOnRefElem)[3] )
  {
    bool const isInsidePolyhedron =
      computationalGeometry::isPointInsidePolyhedron( X,
                                                      elemsToFaces,
                                                      facesToNodes,
                                                      elemCenter,
                                                      coords );
    if( isInsidePolyhedron )
    {
      real64 xLocal[FE_TYPE::numNodes][3]{};
      for( localIndex a = 0; a < FE_TYPE::numNodes; ++a )
      {
        LvArray::tensorOps::copy< 3 >( xLocal[a], X[ elemsToNodes[a] ] );
      }

      // coordsOnRefElem = invJ*(coords-coordsNode_0)
      real64 invJ[3][3]{};
      FE_TYPE::invJacobianTransformation( 0, xLocal, invJ );

      for( localIndex i = 0; i < 3; ++i )
      {
        // init at (-1,-1,-1) as the origin of the referential elem
        coordsOnRefElem[i] = -1.0;
        for( localIndex j = 0; j < 3; ++j )
        {
          coordsOnRefElem[i] += invJ[i][j] * (coords[j] - xLocal[0][j]);
        }
      }
      return true;
    }
    return false;

  }

  GEOSX_HOST_DEVICE
  static real64
  evaluateRicker( real64 const & time_n,
                  real64 const & f0,
                  localIndex const & order )
  {
    real64 const o_tpeak = 1.0/f0;
    real64 pulse = 0.0;
    if((time_n <= -0.9*o_tpeak) || (time_n >= 2.9*o_tpeak))
    {
      return pulse;
    }

    constexpr real64 pi = M_PI;
    real64 const lam = (f0*pi)*(f0*pi);

    switch( order )
    {
      case 2:
      {
        pulse = 2.0*lam*(2.0*lam*(time_n-o_tpeak)*(time_n-o_tpeak)-1.0)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
      }
      break;
      case 1:
      {
        pulse = -2.0*lam*(time_n-o_tpeak)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
      }
      break;
      case 0:
      {
        pulse = -(time_n-o_tpeak)*exp( -2*lam*(time_n-o_tpeak)*(time_n-o_tpeak) );
      }
      break;
      default:
        GEOSX_ERROR( "This option is not supported yet, rickerOrder must be 0, 1 or 2" );
    }

    return pulse;
  }
  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] numNodesPerElem number of nodes per element
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] facesToNodes map from faces to nodes
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsLocal flag indicating whether the source is local or not
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceNodeConstants constant part of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverNodeConstants constant part of the receiver term
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          localIndex const numNodesPerElem,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const & facesToNodes,
          arrayView2d< real64 const > const & elemCenter,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsLocal,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstants,
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants,
          arrayView2d< real64 > const sourceValue,
          real64 const dt,
          real64 const timeSourceFrequency,
          localIndex const rickerOrder )
  {

    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      real64 const center[3] = { elemCenter[k][0],
                                 elemCenter[k][1],
                                 elemCenter[k][2] };

      // Step 1: locate the sources, and precompute the source term

      /// loop over all the source that haven't been found yet
      for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
      {
        if( sourceIsLocal[isrc] == 0 )
        {
          real64 const coords[3] = { sourceCoordinates[isrc][0],
                                     sourceCoordinates[isrc][1],
                                     sourceCoordinates[isrc][2] };

          real64 coordsOnRefElem[3]{};
          bool const sourceFound =
            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                             center,
                                                             elemsToNodes[k],
                                                             elemsToFaces[k],
                                                             facesToNodes,
                                                             X,
                                                             coordsOnRefElem );
          if( sourceFound )
          {
            sourceIsLocal[isrc] = 1;
            real64 Ntest[8];
            finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes[k][a];
              sourceConstants[isrc][a] = Ntest[a];
            }

            for( localIndex cycle = 0; cycle < sourceValue.size( 0 ); ++cycle )
            {
              real64 const time = cycle*dt;
              sourceValue[cycle][isrc] = evaluateRicker( time, timeSourceFrequency, rickerOrder );
            }

          }
        }
      } // end loop over all sources


      // Step 2: locate the receivers, and precompute the receiver term

      /// loop over all the receivers that haven't been found yet
      for( localIndex ircv = 0; ircv < receiverCoordinates.size( 0 ); ++ircv )
      {
        if( receiverIsLocal[ircv] == 0 )
        {
          real64 const coords[3] = { receiverCoordinates[ircv][0],
                                     receiverCoordinates[ircv][1],
                                     receiverCoordinates[ircv][2] };

          real64 coordsOnRefElem[3]{};
          bool const receiverFound =
            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                             center,
                                                             elemsToNodes[k],
                                                             elemsToFaces[k],
                                                             facesToNodes,
                                                             X,
                                                             coordsOnRefElem );
          if( receiverFound )
          {
            receiverIsLocal[ircv] = 1;

            real64 Ntest[8];
            finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              receiverNodeIds[ircv][a] = elemsToNodes[k][a];
              receiverConstants[ircv][a] = Ntest[a];
            }
          }
        }
      } // end loop over receivers
    } );

  }
};

template< typename FE_TYPE >
struct MassAndDampingMatrixKernel
{

  MassAndDampingMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the mass and damping matrices
   * @param[in] size the number of cells in the subRegion
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] density cell-wise density
   * @param[out] mass diagonal of the mass matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
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
          arrayView1d< real64 const > const density,
          arrayView1d< real64 > const velocityVp,
          arrayView1d< real64 > const velocityVs,
          arrayView1d< real64 > const damping_x,
          arrayView1d< real64 > const damping_y,
          arrayView1d< real64 > const damping_z,
          arrayView1d< real64 > const mass )
  {
    forAll< EXEC_POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
  

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
          real64 const localIncrement = density[k] * detJ * N[a];
          RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes[k][a]], localIncrement );
        }
      }

      for( localIndex kfe = 0; kfe < numFacesPerElem; ++kfe )
      {
        localIndex const iface = elemsToFaces[k][kfe];

        if( facesDomainBoundaryIndicator[iface] == 1 && freeSurfaceFaceIndicator[iface] != 1 )
        {
          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );
            real64 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

            real64 invJ[3][3]={{0}};
            FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

            for( localIndex a=0; a < numNodesPerFace; ++a )
            {
              /// compute ds=||detJ*invJ*normalFace_{kfe}||
              real64 tmp[3]={0};
              real64 ds = 0.0;
              for( localIndex b=0; b<3; ++b )
              {
                for( localIndex j = 0; j < 3; ++j )
                {
                  tmp[b] += invJ[j][b]*faceNormal[iface][j];
                }
                ds +=tmp[b]*tmp[b];
              }
              ds = std::sqrt( ds );

              localIndex numNodeGl = facesToNodes[iface][a];

              // Damping in x=xpos and x=xneg direction
              if( faceNormal[iface][0] < 0.0 || faceNormal[iface][0] > 0.0)
              {
                real64 const alpha_x = density[k] * velocityVp[k];
                real64 const alpha_y = density[k] * velocityVs[k];
                real64 const alpha_z = density[k] * velocityVs[k];
                damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
                damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
                damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];
              }
              // Damping in y=ypos and y=yneg direction
              if( faceNormal[iface][1] < 0.0 || faceNormal[iface][1] > 0.0)
              {
                real64 const alpha_x = density[k] * velocityVs[k];
                real64 const alpha_y = density[k] * velocityVp[k];
                real64 const alpha_z = density[k] * velocityVs[k];
                damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
                damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
                damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];
              }
              // Damping in z=zpos and z=zneg direction
              if( faceNormal[iface][2] < 0.0 || faceNormal[iface][2] > 0.0)
              {
                real64 const alpha_x = density[k] * velocityVs[k];
                real64 const alpha_y = density[k] * velocityVs[k];
                real64 const alpha_z = density[k] * velocityVp[k];
                damping_x[numNodeGl] += alpha_x*detJ*ds*N[a];
                damping_y[numNodeGl] += alpha_y*detJ*ds*N[a];
                damping_z[numNodeGl] += alpha_z*detJ*ds*N[a];
              }
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
class ExplicitElasticSEM : public finiteElement::KernelBase< SUBREGION_TYPE,
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
   *   elements to be processed during this kernel launch.
   */
  ExplicitElasticSEM( NodeManager & nodeManager,
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
    m_ux_n( nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_n >() ),
    m_uy_n( nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_n >() ),
    m_uz_n( nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_n >() ),
    m_stiffnessVector_x( nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector_x >() ),
    m_stiffnessVector_y( nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector_y >() ),
    m_stiffnessVector_z( nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector_z >() ),
    m_density( elementSubRegion.template getExtrinsicData< extrinsicMeshData::MediumDensity >() ),
    m_velocityVp( elementSubRegion.template getExtrinsicData< extrinsicMeshData::MediumVelocityVp >() ),
    m_velocityVs( elementSubRegion.template getExtrinsicData< extrinsicMeshData::MediumVelocityVs >() ),
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
    real64 xLocal[ numNodesPerElem ][ 3 ]{};   
    real64 mu=0;
    real64 lambda=0;
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

    real64 sigmaxx[numNodesPerElem] = {{0.0}};
    real64 sigmaxy[numNodesPerElem] = {{0.0}};
    real64 sigmaxz[numNodesPerElem] = {{0.0}};
    real64 sigmayy[numNodesPerElem] = {{0.0}};
    real64 sigmayx[numNodesPerElem] = {{0.0}};
    real64 sigmayz[numNodesPerElem] = {{0.0}};
    real64 sigmazz[numNodesPerElem] = {{0.0}};
    real64 sigmazx[numNodesPerElem] = {{0.0}};
    real64 sigmazy[numNodesPerElem] = {{0.0}};

    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, gradN );

    for( localIndex i=0; i<numNodesPerElem; ++i )
    {
      sigmaxx[i] = detJ * ((stack.lambda + 2*stack.mu) * gradN[i][0]*m_ux_n[m_elemsToNodes[k][i]] + stack.lambda * ( gradN[i][1]*m_uy_n[m_elemsToNodes[k][i]]+ gradN[i][2]*m_uz_n[m_elemsToNodes[k][i]]));
      sigmayy[i] = detJ * ((stack.lambda + 2*stack.mu) * gradN[i][0]*m_uy_n[m_elemsToNodes[k][i]] + stack.lambda * ( gradN[i][1]*m_ux_n[m_elemsToNodes[k][i]]+ gradN[i][2]*m_uz_n[m_elemsToNodes[k][i]]));
      sigmazz[i] = detJ * ((stack.lambda + 2*stack.mu) * gradN[i][0]*m_uz_n[m_elemsToNodes[k][i]] + stack.lambda * ( gradN[i][1]*m_ux_n[m_elemsToNodes[k][i]]+ gradN[i][2]*m_uy_n[m_elemsToNodes[k][i]]));
      sigmaxy[i] = detJ * (stack.mu * ( gradN[i][1]*m_ux_n[m_elemsToNodes[k][i]]+ gradN[i][0]*m_uy_n[m_elemsToNodes[k][i]]));
      sigmaxz[i] = detJ * (stack.mu * ( gradN[i][2]*m_ux_n[m_elemsToNodes[k][i]]+ gradN[i][0]*m_uz_n[m_elemsToNodes[k][i]]));
      sigmayz[i] = detJ * (stack.mu * ( gradN[i][2]*m_uy_n[m_elemsToNodes[k][i]]+ gradN[i][1]*m_uz_n[m_elemsToNodes[k][i]]));

      sigmayx[i] = sigmaxy[i];
      sigmazy[i] = sigmayz[i];
      sigmazx[i] = sigmaxz[i];
    }


    for( localIndex i=0; i<numNodesPerElem; ++i )
    {
      for( localIndex j=0; j<numNodesPerElem; ++j )
      {
        // real64 const Rxx_ij = detJ* ((stack.lambda+2.0*stack.mu)*gradN[j][0]*gradN[i][0] + stack.mu * gradN[j][1]*gradN[i][1] + stack.mu * gradN[j][2]*gradN[i][2]);
        // real64 const Ryy_ij = detJ* ((stack.lambda+2.0*stack.mu)*gradN[j][1]*gradN[i][1] + stack.mu * gradN[j][0]*gradN[i][0] + stack.mu * gradN[j][2]*gradN[i][2]);
        // real64 const Rzz_ij = detJ*((stack.lambda+2.0*stack.mu)*gradN[j][2]*gradN[i][2] + stack.mu * gradN[j][1]*gradN[i][1] + stack.mu * gradN[j][0]*gradN[i][0]);
        // real64 const Rxy_ij =  detJ*(stack.mu * gradN[j][1]*gradN[i][0] + stack.lambda * gradN[j][0]*gradN[i][1]);
        // real64 const Ryx_ij =  detJ*(stack.mu * gradN[j][0]*gradN[i][1] + stack.lambda * gradN[j][1]*gradN[i][0]);
        // real64 const Rxz_ij =  detJ*(stack.mu * gradN[j][2]*gradN[i][0] + stack.lambda * gradN[j][0]*gradN[i][2]);
        // real64 const Rzx_ij =  detJ*(stack.mu * gradN[j][0]*gradN[i][2] + stack.lambda * gradN[j][2]*gradN[i][0]);
        // real64 const Ryz_ij =  detJ*(stack.mu * gradN[j][2]*gradN[i][1] + stack.lambda * gradN[j][1]*gradN[i][2]);
        // real64 const Rzy_ij =  detJ*(stack.mu * gradN[j][1]*gradN[i][2] + stack.lambda * gradN[j][2]*gradN[i][1]);

        real64 const Rxx_ij = gradN[i][0];
        real64 const Rxz_ij = gradN[i][2];
        real64 const Rxy_ij = gradN[i][1];
        real64 const Ryy_ij = gradN[i][1];
        real64 const Ryx_ij = gradN[i][0];
        real64 const Ryz_ij = gradN[i][2];
        real64 const Rzx_ij = gradN[i][0];
        real64 const Rzy_ij = gradN[i][1];
        real64 const Rzz_ij = gradN[i][2];

        real64 const localIncrement_x = Rxx_ij * sigmaxx[j] + Rxy_ij*sigmaxy[j] + Rxz_ij*sigmaxz[j];
        real64 const localIncrement_y = Ryx_ij * sigmayx[j] + Ryy_ij*sigmayy[j] + Ryz_ij*sigmayz[j];
        real64 const localIncrement_z = Rzx_ij * sigmazx[j] + Rzy_ij*sigmazy[j] + Rzz_ij*sigmazz[j];

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_x[m_elemsToNodes[k][i]], localIncrement_x );
        RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_y[m_elemsToNodes[k][i]], localIncrement_y );
        RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_z[m_elemsToNodes[k][i]], localIncrement_z );
      }

    }

  }


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal displacement array in x direction.
  arrayView1d< real64 > const m_ux_n;

  /// The array containing the nodal displacement array in y direction.
  arrayView1d< real64 > const m_uy_n;

  /// The array containing the nodal displacement array in z direction.
  arrayView1d< real64 > const m_uz_n;

  /// The array containing the product of the stiffness matrix and the nodal pressure.
  arrayView1d< real64 > const m_stiffnessVector_x;

  /// The array containing the product of the stiffness matrix and the nodal pressure.
  arrayView1d< real64 > const m_stiffnessVector_y;

  /// The array containing the product of the stiffness matrix and the nodal pressure.
  arrayView1d< real64 > const m_stiffnessVector_z;

  /// The array containing the density of the medium
  arrayView1d< real64 const > const m_density;

  /// The array containing the P-wavespeed
  arrayView1d< real64 const > const m_velocityVp;

  /// The array containing the S-wavespeed
  arrayView1d< real64 const > const m_velocityVs;

  /// The time increment for this time integration step.
  real64 const m_dt;


};


/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitElasticSEMFactory = finiteElement::KernelFactory< ExplicitElasticSEM,
                                                                real64 >;

} // namespace ElasticWaveEquationSEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEMKERNEL_HPP_
