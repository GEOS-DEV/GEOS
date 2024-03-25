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
 * @file AcousticWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "WaveSolverUtils.hpp"
#if !defined( GEOS_USE_HIP )
#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#endif
#include "AcousticFields.hpp"

namespace geos
{

using namespace fields;

/// Namespace to contain the acoustic wave kernels.
namespace acousticWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{

  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
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
              sourceNodeIds[isrc][a] = elemsToNodes( k, a );
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
              receiverNodeIds[ircv][a] = elemsToNodes( k, a );
              receiverConstants[ircv][a] = Ntest[a];
            }
          }
        }
      } // end loop over receivers

    } );

  }
};

/**
 * @brief Implements kernels for solving the acoustic wave equations
 *   explicit central FD method and SEM
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### AcousticWaveEquationSEMKernel Description
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
    m_nodeCoords( nodeManager.getField< fields::referencePosition32 >() ),
    m_p_n( nodeManager.getField< acousticfields::Pressure_n >() ),
    m_stiffnessVector( nodeManager.getField< acousticfields::StiffnessVector >() ),
    m_density( elementSubRegion.template getField< acousticfields::AcousticDensity >() ),
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
   * ### ExplicitAcousticSEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal(),
      stiffnessVectorLocal()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ 8 ][ 3 ];
    real32 stiffnessVectorLocal[ numNodesPerElem ]{};
    real32 invDensity;
  };
  //***************************************************************************


  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    stack.invDensity = 1./m_density[k];
    for( localIndex a=0; a< 8; a++ )
    {
      localIndex const nodeIndex =  m_elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_nodeCoords[ nodeIndex ][ i ];
      }
    }
  }

  /**
   * @copydoc geos::finiteElement::KernelBase::complete
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    for( int i=0; i<numNodesPerElem; i++ )
    {
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector[m_elemsToNodes( k, i )], stack.stiffnessVectorLocal[i] );
    }
    return 0;
  }


  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitAcousticSEM Description
   * Calculates stiffness vector
   *
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    m_finiteElementSpace.template computeStiffnessTerm( q, stack.xLocal, [&] ( const int i, const int j, const real64 val )
    {
      real32 const localIncrement = stack.invDensity*val*m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal[ i ] += localIncrement;
    } );
  }

protected:
  /// The array containing the nodal position array.
  arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const m_nodeCoords;

  /// The array containing the nodal pressure array.
  arrayView1d< real32 const > const m_p_n;

  /// The array containing the product of the stiffness matrix and the nodal pressure.
  arrayView1d< real32 > const m_stiffnessVector;

  /// The array containing the cell-wise density
  arrayView1d< real32 const > const m_density;

  /// The time increment for this time integration step.
  real64 const m_dt;


};



/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitAcousticSEMFactory = finiteElement::KernelFactory< ExplicitAcousticSEM,
                                                                 real64 >;


} // namespace acousticWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_
