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

template< typename FE_TYPE >
struct ComputeTimeStep
{

  ComputeTimeStep( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Compute timestep using power iteration method
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  real64
  launch( localIndex const sizeElem,
          localIndex const sizeNode,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView1d< real32 > const mass )
  {

    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

    real64 const epsilon = 0.00001;
    localIndex const nIterMax = 10000;
    localIndex numberIter = 0;
    localIndex counter = 0;
    real64 lambdaNew = 0.0;

    array1d< real32 > const p( sizeNode );
    array1d< real32 > const pAux( sizeNode );

    arrayView1d< real32 > const pView = p;
    arrayView1d< real32 > const pAuxView = pAux;

    //Randomize p values
    srand( time( NULL ));
    for( localIndex a = 0; a < sizeNode; ++a )
    {
      pView[a] = (real64)rand()/(real64) RAND_MAX;
    }

    //Step 1: Normalize randomized pressure
    real64 normP= 0.0;
    WaveSolverUtils::dotProduct( sizeNode, pView, pView, normP );

    forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      pView[a]/= sqrt( normP );
    } );

    //Step 2: Initial iteration of (M^{-1}K)p
    pAuxView.zero();
    forAll< EXEC_POLICY >( sizeElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      real64 xLocal[ 8 ][ 3 ];
      for( localIndex a = 0; a < 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = nodeCoords[nodeIndex][i];
        }
      }
      for( localIndex q = 0; q < numNodesPerElem; ++q )
      {
        m_finiteElement.computeStiffnessTerm( q, xLocal, [&] ( const int i, const int j, const real64 val )
        {
          real32 const localIncrement = val*pView[elemsToNodes[k][j]];
          RAJA::atomicAdd< parallelDeviceAtomic >( &pAuxView[elemsToNodes[k][i]], localIncrement );
        } );
      }

    } );

    forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      pAuxView[a]/= mass[a];
    } );
    real64 lambdaOld = lambdaNew;

    //Compute lambdaNew using two dotProducts
    real64 dotProductPPaux = 0.0;
    normP = 0.0;
    WaveSolverUtils::dotProduct( sizeNode, pView, pAuxView, dotProductPPaux );
    WaveSolverUtils::dotProduct( sizeNode, pView, pView, normP );

    lambdaNew = dotProductPPaux/normP;

    real64 normPaux = 0.0;
    WaveSolverUtils::dotProduct( sizeNode, pAuxView, pAuxView, normPaux );
    forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      pView[a] = pAuxView[a]/( normPaux );
    } );

    //Step 3: Do previous algorithm until we found the max eigenvalues
    do
    {
      pAuxView.zero();
      forAll< EXEC_POLICY >( sizeElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
      {

        real64 xLocal[ 8 ][ 3 ];
        for( localIndex a = 0; a < 8; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
          for( localIndex i = 0; i < 3; ++i )
          {
            xLocal[a][i] = nodeCoords[nodeIndex][i];
          }
        }
        for( localIndex q = 0; q < numNodesPerElem; ++q )
        {
          m_finiteElement.computeStiffnessTerm( q, xLocal, [&] ( const int i, const int j, const real64 val )
          {
            real32 const localIncrement = val*pView[elemsToNodes[k][j]];
            RAJA::atomicAdd< parallelDeviceAtomic >( &pAuxView[elemsToNodes[k][i]], localIncrement );
          } );
        }

      } );

      forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        pAuxView[a]/= mass[a];
      } );

      lambdaOld = lambdaNew;

      dotProductPPaux = 0.0;
      normP=0.0;
      WaveSolverUtils::dotProduct( sizeNode, pView, pAuxView, dotProductPPaux );
      WaveSolverUtils::dotProduct( sizeNode, pView, pView, normP );

      lambdaNew = dotProductPPaux/normP;

      normPaux = 0.0;
      WaveSolverUtils::dotProduct( sizeNode, pAuxView, pAuxView, normPaux );

      forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        pView[a] = pAuxView[a]/( normPaux );
      } );

      if( abs( lambdaNew-lambdaOld )/abs( lambdaNew )<= epsilon )
      {
        counter++;
      }
      else
      {
        counter=0;
      }

      numberIter++;


    }
    while (counter < 10 && numberIter < nIterMax);

    GEOS_THROW_IF( numberIter> nIterMax, "Power Iteration algorithm does not converge", std::runtime_error );

    real64 dt = 1.99/sqrt( abs( lambdaNew ));

    return dt;

  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
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
