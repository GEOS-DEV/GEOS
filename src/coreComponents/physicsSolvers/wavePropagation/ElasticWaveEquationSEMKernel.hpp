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
   * @brief Check if the source point is inside an element or not
   * @param numFacesPerElem number of face on an element
   * @param elemCenter array containing the center of the elements
   * @param faceNormal array containing the normal of all faces
   * @param faceCenter array containing the center of all faces
   * @param elemsToFaces map to get the global faces from element index and local face index
   * @param coords coordinate of the point
   * @return true if coords is inside the element
   */

  GEOSX_HOST_DEVICE
  static bool
  locateSourceElement( real64 const numFacesPerElem,
                       real64 const (&elemCenter)[3],
                       arrayView2d< real64 const > const faceNormal,
                       arrayView2d< real64 const > const faceCenter,
                       arraySlice1d< localIndex const > const elemsToFaces,
                       real64 const (&coords)[3] )
  {
    //Loop over the element faces
    real64 tmpVector[3]{};
    for( localIndex kfe = 0; kfe < numFacesPerElem; ++kfe )
    {

      localIndex const iface = elemsToFaces[kfe];
      real64 faceCenterOnFace[3] = {faceCenter[iface][0],
                                    faceCenter[iface][1],
                                    faceCenter[iface][2]};
      real64 faceNormalOnFace[3] = {faceNormal[iface][0],
                                    faceNormal[iface][1],
                                    faceNormal[iface][2]};

      //Test to make sure if the normal is outwardly directed
      LvArray::tensorOps::copy< 3 >( tmpVector, faceCenterOnFace );
      LvArray::tensorOps::subtract< 3 >( tmpVector, elemCenter );
      if( LvArray::tensorOps::AiBi< 3 >( tmpVector, faceNormalOnFace ) < 0.0 )
      {
        LvArray::tensorOps::scale< 3 >( faceNormalOnFace, -1 );
      }

      // compute the vector face center to query point
      LvArray::tensorOps::subtract< 3 >( faceCenterOnFace, coords );
      localIndex const s = computationalGeometry::sign( LvArray::tensorOps::AiBi< 3 >( faceNormalOnFace, faceCenterOnFace ));

      // all dot products should be non-negative (we enforce outward normals)
      if( s < 0 )
      {
        return false;
      }

    }
    return true;
  }

  /**
   * @brief Convert a mesh element point coordinate into a coordinate on the reference element
   * @tparam FE_TYPE finite element type
   * @param[in] coords coordinate of the point
   * @param[in] elemsToNodes map to obtaint global nodes from element index
   * @param[in] X array of mesh nodes coordinates
   * @param[out] coordsOnRefElem to contain the coordinate computed in the reference element
   */
  template< typename FE_TYPE >
  GEOSX_HOST_DEVICE
  static void
  computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                        arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const elemsToNodes,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                                        real64 (& coordsOnRefElem)[3] )
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
   * @param[in] numFacesPerElem number of face on an element
   * @param[in] X coordinates of the nodes
   * @param[in] elemGhostRank array containing the ghost rank
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] faceNormal array containing the normal of all faces
   * @param[in] faceCenter array containing the center of all faces
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[in] dt time-step
   * @param[in] timeSourceFrequency Peak frequency of the source
   * @param[in] rickerOrder Order of the Ricker wavelet
   * @param[out] sourceIsLocal flag indicating whether the source is local or not
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceConstantsx constant part of the source terms in x-direction
   * @param[out] sourceConstantsy constant part of the source terms in y-direction
   * @param[out] sourceConstantsz constant part of the source terms in z-direction
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverNodeConstants constant part of the receiver term
   * @param[out] sourceValue array containing the value of the time dependent source (Ricker for e.g)
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          localIndex const numFacesPerElem,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView1d< integer const > const elemGhostRank,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          arrayView2d< real64 const > const & elemCenter,
          arrayView2d< real64 const > const faceNormal,
          arrayView2d< real64 const > const faceCenter,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsLocal,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstantsx,
          arrayView2d< real64 > const sourceConstantsy,
          arrayView2d< real64 > const sourceConstantsz,
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

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

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

          real64 xLocal[numNodesPerElem][3];

          for( localIndex a=0; a< numNodesPerElem; ++a )
          {
            for( localIndex i=0; i<3; ++i )
            {
              xLocal[a][i] = X( elemsToNodes( k, a ), i );
            }
          }


          bool const sourceFound =
            locateSourceElement( numFacesPerElem,
                                 center,
                                 faceNormal,
                                 faceCenter,
                                 elemsToFaces[k],
                                 coords );

          if( sourceFound && elemGhostRank[k] < 0 )
          {
            real64 coordsOnRefElem[3]{};


            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                             elemsToNodes[k],
                                                             X,
                                                             coordsOnRefElem );
            sourceIsLocal[isrc] = 1;

            //Compute source coefficients: this generate a P-wave and an "unwanted" S-wave. It is classical in the case of the elastic wave
            // equation at order 2, the S-wave can be attenuated by refining the mesh or get to high order
            //However, we will propably use elastic wave at 1st order for the FWI case.
            for( localIndex c=0; c<2; ++c )
            {
              for( localIndex b=0; b<2; ++b )
              {
                for( localIndex a=0; a<2; ++a )
                {
                  real64 const Grad[3] = { finiteElement::LagrangeBasis1::gradient( a, coordsOnRefElem[0] )*
                                           finiteElement::LagrangeBasis1::value( b, coordsOnRefElem[1] )*
                                           finiteElement::LagrangeBasis1::value( c, coordsOnRefElem[2] ),
                                           finiteElement::LagrangeBasis1::value( a, coordsOnRefElem[0] )*
                                           finiteElement::LagrangeBasis1::gradient( b, coordsOnRefElem[1] )*
                                           finiteElement::LagrangeBasis1::value( c, coordsOnRefElem[2] ),
                                           finiteElement::LagrangeBasis1::value( a, coordsOnRefElem[0] )*
                                           finiteElement::LagrangeBasis1::value( b, coordsOnRefElem[1] )*
                                           finiteElement::LagrangeBasis1::gradient( c, coordsOnRefElem[2] )};

                  localIndex const nodeIndex = finiteElement::LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );

                  real64 invJ[3][3]={{0}};
                  FE_TYPE::invJacobianTransformation( nodeIndex, xLocal, invJ );
                  sourceNodeIds[isrc][nodeIndex] = elemsToNodes[k][nodeIndex];
                  sourceConstantsx[isrc][nodeIndex] = Grad[0] * invJ[0][0] + Grad[1] * invJ[0][1] + Grad[2] * invJ[0][2];
                  sourceConstantsy[isrc][nodeIndex] = Grad[0] * invJ[1][0] + Grad[1] * invJ[1][1] + Grad[2] * invJ[1][2];
                  sourceConstantsz[isrc][nodeIndex] = Grad[0] * invJ[2][0] + Grad[1] * invJ[2][1] + Grad[2] * invJ[2][2];
                }
              }
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
            locateSourceElement( numFacesPerElem,
                                 center,
                                 faceNormal,
                                 faceCenter,
                                 elemsToFaces[k],
                                 coords );

          if( receiverFound && elemGhostRank[k] < 0 )
          {
            computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                             elemsToNodes[k],
                                                             X,
                                                             coordsOnRefElem );

            receiverIsLocal[ircv] = 1;

            real64 Ntest[numNodesPerElem];
            //finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

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
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] numFacesPerElem number of faces per element
   * @param[in] numNodesPerFace number of nodes per face
   * @param[in] X coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] facesToNodes map from face to nodes
   * @param[in] facesDomainBoundaryIndicator flag equal to 1 if the face is on the boundary, and to 0 otherwise
   * @param[in] freeSurfaceFaceIndicator flag equal to 1 if the face is on the free surface, and to 0 otherwise
   * @param[in] faceNormal normal vectors at the faces
   * @param[in] density cell-wise density
   * @param[in] velocityVp cell-wise P-wavespeed
   * @param[in] velocityVp cell-wise S-wavespeed
   * @param[out] dampingx diagonal of the damping matrix (x-part)
   * @param[out] dampingy diagonal of the damping matrix (y-part)
   * @param[out] dampingz diagonal of the damping matrix (z-part)
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
          arrayView1d< real64 > const dampingx,
          arrayView1d< real64 > const dampingy,
          arrayView1d< real64 > const dampingz,
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

              real64 const alphax = density[k] * (velocityVp[k]*(faceNormal[iface][0]*faceNormal[iface][0]) + velocityVs[k]*(faceNormal[iface][1]*faceNormal[iface][1]) +
                                                  velocityVs[k]*(faceNormal[iface][2]*faceNormal[iface][2]) );

              real64 const alphay = density[k] * (velocityVs[k]*(faceNormal[iface][0]*faceNormal[iface][0]) + velocityVp[k]*(faceNormal[iface][1]*faceNormal[iface][1]) +
                                                  velocityVs[k]*(faceNormal[iface][2]*faceNormal[iface][2]) );

              real64 const alphaz = density[k] * (velocityVs[k]*(faceNormal[iface][0]*faceNormal[iface][0]) + velocityVs[k]*(faceNormal[iface][1]*faceNormal[iface][1]) +
                                                  velocityVp[k]*(faceNormal[iface][2]*faceNormal[iface][2]) );

              RAJA::atomicAdd< ATOMIC_POLICY >( &dampingx[numNodeGl], alphax*detJ*ds*N[a] );
              RAJA::atomicAdd< ATOMIC_POLICY >( &dampingy[numNodeGl], alphay*detJ*ds*N[a] );
              RAJA::atomicAdd< ATOMIC_POLICY >( &dampingz[numNodeGl], alphaz*detJ*ds*N[a] );

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
 * @brief Implements kernels for solving the elastic wave equations
 *   explicit central FD method and SEM
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### ElasticWaveEquationSEMKernel Description
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
    m_stiffnessVectorx( nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVectorx >() ),
    m_stiffnessVectory( nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVectory >() ),
    m_stiffnessVectorz( nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVectorz >() ),
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
   * ### ExplicitElasticSEM Description
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
    stack.mu = m_density[k] * m_velocityVs[k] * m_velocityVs[k];
    stack.lambda = m_density[k] *m_velocityVp[k] * m_velocityVp[k] - 2.0*stack.mu;
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitElasticSEM Description
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
        real64 const Rxx_ij = detJ* ((stack.lambda+2.0*stack.mu)*gradN[j][0]*gradN[i][0] + stack.mu * gradN[j][1]*gradN[i][1] + stack.mu * gradN[j][2]*gradN[i][2]);
        real64 const Ryy_ij = detJ* ((stack.lambda+2.0*stack.mu)*gradN[j][1]*gradN[i][1] + stack.mu * gradN[j][0]*gradN[i][0] + stack.mu * gradN[j][2]*gradN[i][2]);
        real64 const Rzz_ij = detJ*((stack.lambda+2.0*stack.mu)*gradN[j][2]*gradN[i][2] + stack.mu * gradN[j][1]*gradN[i][1] + stack.mu * gradN[j][0]*gradN[i][0]);
        real64 const Rxy_ij =  detJ*(stack.mu * gradN[j][1]*gradN[i][0] + stack.lambda * gradN[j][0]*gradN[i][1]);
        real64 const Ryx_ij =  detJ*(stack.mu * gradN[j][0]*gradN[i][1] + stack.lambda * gradN[j][1]*gradN[i][0]);
        real64 const Rxz_ij =  detJ*(stack.mu * gradN[j][2]*gradN[i][0] + stack.lambda * gradN[j][0]*gradN[i][2]);
        real64 const Rzx_ij =  detJ*(stack.mu * gradN[j][0]*gradN[i][2] + stack.lambda * gradN[j][2]*gradN[i][0]);
        real64 const Ryz_ij =  detJ*(stack.mu * gradN[j][2]*gradN[i][1] + stack.lambda * gradN[j][1]*gradN[i][2]);
        real64 const Rzy_ij =  detJ*(stack.mu * gradN[j][1]*gradN[i][2] + stack.lambda * gradN[j][2]*gradN[i][1]);

        real64 const localIncrementx = (Rxx_ij * m_ux_n[m_elemsToNodes[k][j]] + Rxy_ij*m_uy_n[m_elemsToNodes[k][j]] + Rxz_ij*m_uz_n[m_elemsToNodes[k][j]]);
        real64 const localIncrementy = (Ryx_ij * m_ux_n[m_elemsToNodes[k][j]] + Ryy_ij*m_uy_n[m_elemsToNodes[k][j]] + Ryz_ij*m_uz_n[m_elemsToNodes[k][j]]);
        real64 const localIncrementz = (Rzx_ij * m_ux_n[m_elemsToNodes[k][j]] + Rzy_ij*m_uy_n[m_elemsToNodes[k][j]] + Rzz_ij*m_uz_n[m_elemsToNodes[k][j]]);

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVectorx[m_elemsToNodes[k][i]], localIncrementx );
        RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVectory[m_elemsToNodes[k][i]], localIncrementy );
        RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVectorz[m_elemsToNodes[k][i]], localIncrementz );
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
  arrayView1d< real64 > const m_stiffnessVectorx;

  /// The array containing the product of the stiffness matrix and the nodal pressure.
  arrayView1d< real64 > const m_stiffnessVectory;

  /// The array containing the product of the stiffness matrix and the nodal pressure.
  arrayView1d< real64 > const m_stiffnessVectorz;

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
