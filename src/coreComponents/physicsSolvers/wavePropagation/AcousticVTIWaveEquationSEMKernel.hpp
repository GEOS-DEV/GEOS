/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOS Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file AcousticVTIWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "WaveSolverBaseFields.hpp"
#include "WaveSolverUtils.hpp"


namespace geos
{

/// Namespace to contain the acoustic wave kernels.
namespace acousticVTIWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{
  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] numFacesPerElem number of faces per element
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemGhostRank the ghost ranks
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] facesToNodes map from faces to nodes
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceNodeConstants constant part of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverNodeConstants constant part of the receiver term
   * @param[out] sourceValue the value of the source
   * @param[in] dt the time step size
   * @param[in] timeSourceFrequency the time frequency of the source
   * @param[in] timeSourceDelay the time delay of the source
   * @param[in] rickerOrder the order of the ricker
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          localIndex const numFacesPerElem,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView1d< integer const > const elemGhostRank,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          arrayView2d< real64 const > const & elemCenter,
          arrayView2d< real64 const > const faceNormal,
          arrayView2d< real64 const > const faceCenter,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsAccessible,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstants,
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants,
          arrayView2d< real32 > const sourceValue,
          real64 const dt,
          real32 const timeSourceFrequency,
          real32 const timeSourceDelay,
          localIndex const rickerOrder )
  {
    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      real64 const center[3] = { elemCenter[k][0],
                                 elemCenter[k][1],
                                 elemCenter[k][2] };

      // Step 1: locate the sources, and precompute the source term

      /// loop over all the source that haven't been found yet
      for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
      {
        if( sourceIsAccessible[isrc] == 0 )
        {
          real64 const coords[3] = { sourceCoordinates[isrc][0],
                                     sourceCoordinates[isrc][1],
                                     sourceCoordinates[isrc][2] };

          bool const sourceFound =
            WaveSolverUtils::locateSourceElement( numFacesPerElem,
                                                  center,
                                                  faceNormal,
                                                  faceCenter,
                                                  elemsToFaces[k],
                                                  coords );
          if( sourceFound )
          {
            real64 coordsOnRefElem[3]{};


            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              elemsToNodes[k],
                                                                              nodeCoords,
                                                                              coordsOnRefElem );

            sourceIsAccessible[isrc] = 1;
            real64 Ntest[numNodesPerElem];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes[k][a];
              sourceConstants[isrc][a] = Ntest[a];
            }

            for( localIndex cycle = 0; cycle < sourceValue.size( 0 ); ++cycle )
            {
              sourceValue[cycle][isrc] = WaveSolverUtils::evaluateRicker( cycle * dt, timeSourceFrequency, timeSourceDelay, rickerOrder );
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
            WaveSolverUtils::locateSourceElement( numFacesPerElem,
                                                  center,
                                                  faceNormal,
                                                  faceCenter,
                                                  elemsToFaces[k],
                                                  coords );

          if( receiverFound && elemGhostRank[k] < 0 )
          {
            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              elemsToNodes[k],
                                                                              nodeCoords,
                                                                              coordsOnRefElem );

            receiverIsLocal[ircv] = 1;

            real64 Ntest[numNodesPerElem];
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
struct MassMatrixKernel
{

  MassMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the mass matrices
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] velocity cell-wise velocity
   * @param[in] density cell-wise density
   * @param[out] mass diagonal of the mass matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView1d< real32 const > const velocity,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 > const mass )

  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
    {

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real32 const invC2 = 1.0 / pow( velocity[e], 2 ) / density[e];
      real64 xLocal[ numNodesPerElem ][ 3 ];
      for( localIndex a = 0; a < numNodesPerElem; ++a )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = nodeCoords( elemsToNodes( e, a ), i );
        }
      }

      for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
      {
        real32 const localIncrement = invC2 * m_finiteElement.computeMassTerm( q, xLocal );
        RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes( e, q )], localIncrement );
      }
    } ); // end loop over element
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

template< typename FE_TYPE >
struct DampingMatrixKernel
{

  DampingMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the damping matrices
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToFaces map from faces to elements
   * @param[in] facesToNodes map from face to nodes
   * @param[in] facesDomainBoundaryIndicator flag equal to 1 if the face is on the boundary, and to 0 otherwise
   * @param[in] freeSurfaceFaceIndicator flag equal to 1 if the face is on the free surface, and to 0 otherwise
   * @param[in] lateralSurfaceFaceIndicator flag equal to 1 if the face is on a lateral surface, and to 0 otherwise
   * @param[in] bottomSurfaceFaceIndicator flag equal to 1 if the face is on the bottom surface, and to 0 otherwise
   * @param[in] velocity cell-wise velocity
   * @param[in] density cell-wise density
   * @param[in] vti_epsilon cell-wise Thomsen epsilon parameter
   * @param[in] vti_delta cell-wise Thomsen delta parameter
   * @param[in] vti_sigma cell-wise sigma parameter in Fletcher's equation
   * @param[out] damping_p diagonal of the damping matrix for quantities in p in p-equation
   * @param[out] damping_q diagonal of the damping matrix for quantities in q in q-equation
   * @param[out] damping_pq diagonal of the damping matrix for quantities in q in p-equation
   * @param[out] damping_qp diagonal of the damping matrix for quantities in p in q-equation
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView1d< integer const > const facesDomainBoundaryIndicator,
          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
          arrayView1d< localIndex const > const lateralSurfaceFaceIndicator,
          arrayView1d< localIndex const > const bottomSurfaceFaceIndicator,
          arrayView1d< real32 const > const velocity,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 const > const vti_epsilon,
          arrayView1d< real32 const > const vti_delta,
          arrayView1d< real32 const > const vti_sigma,
          arrayView1d< real32 > const damping_p,
          arrayView1d< real32 > const damping_q,
          arrayView1d< real32 > const damping_pq,
          arrayView1d< real32 > const damping_qp )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
    {
      for( localIndex i = 0; i < elemsToFaces.size( 1 ); ++i )
      {
        localIndex const f = elemsToFaces( e, i );
        // face on the domain boundary and not on free surface
        if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
        {
          constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
          real64 xLocal[ numNodesPerFace ][ 3 ];
          for( localIndex a = 0; a < numNodesPerFace; ++a )
          {
            for( localIndex d = 0; d < 3; ++d )
            {
              xLocal[a][d] = nodeCoords( facesToNodes( f, a ), d );
            }
          }
          real32 vti_f = 1 - (vti_epsilon[e] - vti_delta[e]) / vti_sigma[e];
          if( lateralSurfaceFaceIndicator[f] == 1 )
          {
            // ABC coefficients
            real32 alpha = 1.0 / (velocity[e] * density[e] * sqrt( 1+2*vti_epsilon[e] ));
            // VTI coefficients
            real32 vti_p_xy = 0;
            real32 vti_q_xy = 0;
            real32 vti_qp_xy= 0;

            vti_p_xy  = (1+2*vti_epsilon[e]);
            vti_q_xy  = -(vti_f - 1);
            vti_qp_xy = (vti_f+2*vti_delta[e]);
            for( localIndex q = 0; q < numNodesPerFace; ++q )
            {
              real32 const aux = m_finiteElement.computeDampingTerm( q, xLocal );
              real32 const localIncrement_p = alpha * vti_p_xy  * aux;
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping_p[facesToNodes( f, q )], localIncrement_p );

              real32 const localIncrement_q = alpha*vti_q_xy * aux;
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping_q[facesToNodes( f, q )], localIncrement_q );

              real32 const localIncrement_qp = alpha*vti_qp_xy * aux;
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qp[facesToNodes( f, q )], localIncrement_qp );
            }
          }
          if( bottomSurfaceFaceIndicator[f] == 1 )
          {
            // ABC coefficients updated to fit horizontal velocity
            real32 alpha = 1.0 / (velocity[e] * density[e]);
            // VTI coefficients
            real32 vti_p_z  = 0;
            real32 vti_pq_z = 0;
            real32 vti_q_z  = 0;
            vti_p_z  = -(vti_f - 1);
            vti_pq_z = vti_f;
            vti_q_z  = 1;
            for( localIndex q = 0; q < numNodesPerFace; ++q )
            {
              real32 const aux = m_finiteElement.computeDampingTerm( q, xLocal );
              real32 const localIncrement_p = alpha * vti_p_z * aux;
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping_p[facesToNodes( f, q )], localIncrement_p );

              real32 const localIncrement_pq = alpha*vti_pq_z * aux;
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pq[facesToNodes( f, q )], localIncrement_pq );

              real32 const localIncrement_q = alpha * vti_q_z * aux;
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping_q[facesToNodes( f, q )], localIncrement_q );
            }
          }
        }
      }
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};



/**
 * @brief Implements kernels for solving the acoustic wave equations in a VTI medium
 *   explicit central FD method and SEM
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### AcousticVTIWaveEquationSEMKernel Description
 * Implements the KernelBase interface functions required for solving
 * the acoustic wave equations using the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1`.
 */


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitAcousticVTISEM : public finiteElement::KernelBase< SUBREGION_TYPE,
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
   * @copydoc geos::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   *   elements to be processed during this kernel launch.
   */
  ExplicitAcousticVTISEM( NodeManager & nodeManager,
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
    m_nodeCoords( nodeManager.getField< fields::referencePosition32 >() ),
    m_p_n( nodeManager.getField< fields::wavesolverfields::Pressure_p_n >() ),
    m_q_n( nodeManager.getField< fields::wavesolverfields::Pressure_q_n >() ),
    m_stiffnessVector_p( nodeManager.getField< fields::wavesolverfields::StiffnessVector_p >() ),
    m_stiffnessVector_q( nodeManager.getField< fields::wavesolverfields::StiffnessVector_q >() ),
    m_density( elementSubRegion.template getField< fields::wavesolverfields::AcousticDensity >() ),
    m_vti_epsilon( elementSubRegion.template getField< fields::wavesolverfields::AcousticEpsilon >() ),
    m_vti_delta( elementSubRegion.template getField< fields::wavesolverfields::AcousticDelta >() ),
    m_vti_sigma( elementSubRegion.template getField< fields::wavesolverfields::AcousticSigma >() ),
    m_dt( dt )
  {
    GEOS_UNUSED_VAR( edgeManager );
    GEOS_UNUSED_VAR( faceManager );
    GEOS_UNUSED_VAR( targetRegionIndex );
  }

  //*****************************************************************************
  /**
   * @copydoc geos::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitAcousticVTISEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
  };
  //***************************************************************************


  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    /// numDofPerTrialSupportPoint = 1
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_nodeCoords[ nodeIndex ][ i ];
      }
    }
  }

  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitAcousticVTISEM Description
   * Calculates stiffness vector
   *
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Pseudo Stiffness xy
    m_finiteElementSpace.template computeStiffnessxyTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 vti_f = 1 - (m_vti_epsilon[k] - m_vti_delta[k]) / m_vti_sigma[k];
      real32 const localIncrement_p = val / m_density[k] *(-1-2*m_vti_epsilon[k])*m_p_n[m_elemsToNodes[k][j]];
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_p[m_elemsToNodes[k][i]], localIncrement_p );
      real32 const localIncrement_q = val / m_density[k] * ((-2*m_vti_delta[k]-vti_f)*m_p_n[m_elemsToNodes[k][j]] +(vti_f-1)*m_q_n[m_elemsToNodes[k][j]]);
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_q[m_elemsToNodes[k][i]], localIncrement_q );
    } );

    // Pseudo-Stiffness z

    m_finiteElementSpace.template computeStiffnesszTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 vti_f = 1 - (m_vti_epsilon[k] - m_vti_delta[k]) / m_vti_sigma[k];
      real32 const localIncrement_p = val / m_density[k] * ((vti_f-1)*m_p_n[m_elemsToNodes[k][j]] - vti_f*m_q_n[m_elemsToNodes[k][j]]);
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_p[m_elemsToNodes[k][i]], localIncrement_p );

      real32 const localIncrement_q = -val  / m_density[k] * m_q_n[m_elemsToNodes[k][j]];
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_q[m_elemsToNodes[k][i]], localIncrement_q );
    } );
  }

protected:
  /// The array containing the nodal position array.
  arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const m_nodeCoords;

  /// The array containing the nodal pressure array.
  arrayView1d< real32 const > const m_p_n;

  /// The array containing the nodal auxiliary variable array.
  arrayView1d< real32 const > const m_q_n;

  /// The array containing the product of the stiffness matrix and the nodal pressure for the equation in p.
  arrayView1d< real32 > const m_stiffnessVector_p;

  /// The array containing the product of the stiffness matrix and the nodal pressure for the equation in q.
  arrayView1d< real32 > const m_stiffnessVector_q;

  /// The array containing the density parameter.
  arrayView1d< real32 const > const m_density;

  /// The array containing the epsilon Thomsen parameter.
  arrayView1d< real32 const > const m_vti_epsilon;

  /// The array containing the delta Thomsen parameter.
  arrayView1d< real32 const > const m_vti_delta;

  /// The array containing the sigma parameter.
  arrayView1d< real32 const > const m_vti_sigma;

  /// The time increment for this time integration step.
  real64 const m_dt;


};



/// The factory used to construct a ExplicitAcousticVTIWaveEquation kernel.
using ExplicitAcousticVTISEMFactory = finiteElement::KernelFactory< ExplicitAcousticVTISEM,
                                                                    real64 >;

} // namespace acousticVTIWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIWAVEEQUATIONSEMKERNEL_HPP_
