/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsLagrangianSSLEKernels.hpp
 */

#pragma once

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"
#include "finiteElement/FiniteElementShapeFunctionKernel.hpp"
#include "finiteElement/Kinematics.h"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace SolidMechanicsLagrangianSSLEKernels
{


/**
 * @struct Structure to wrap templated function that implements the explicit time integration kernels.
 */
struct ExplicitKernel
{
#if defined(GEOSX_USE_CUDA)
  #define CALCFEMSHAPE
#endif
  // If UPDATE_STRESS is undef, then stress is not updated at all.
//  #define UPDATE_STRESS 1 // uses total displacement to and adds material stress state to integral for nodalforces.
  #define UPDATE_STRESS 2 // uses velocity*dt and updates material stress state.

  /**
   * @brief Launch of the element processing kernel for explicit time integration.
   * @tparam NUM_NODES_PER_ELEM The number of nodes/dof per element.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam CONSTITUTIVE_TYPE the type of the constitutive relation that is being used.
   * @param constitutiveRelation A pointer to the constitutive relation that is being used.
   * @param elementList The list of elements to be processed
   * @param elemsToNodes The map from the elements to the nodes that form that element.
   * @param dNdX The derivatives of the shape functions wrt the reference configuration.
   * @param detJ The determinant of the transformation matrix (Jacobian) to the parent element.
   * @param u The nodal array of total displacements.
   * @param vel The nodal array of velocity.
   * @param acc The nodal array of force/acceleration.
   * @param stress The stress at each element quadrature point.
   * @param dt The timestep
   * @return The achieved timestep.
   */
  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const constitutiveRelation,
          SortedArrayView< localIndex const > const & elementList,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView4d< real64 const > const & dNdX,
          arrayView2d< real64 const > const & detJ,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u,
          arrayView2d< real64 const, nodes::VELOCITY_USD > const & vel,
          arrayView2d< real64, nodes::ACCELERATION_USD > const & acc,
          real64 const dt )
  {
    GEOSX_MARK_FUNCTION;

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelUpdates();

#if defined(CALCFEMSHAPE)
    GEOSX_UNUSED_VAR( dNdX );
    GEOSX_UNUSED_VAR( detJ );
#else
    GEOSX_UNUSED_VAR( X );
#endif

#if UPDATE_STRESS == 2
    GEOSX_UNUSED_VAR( u );
#else
    GEOSX_UNUSED_VAR( vel );
#endif

    using KERNEL_POLICY = parallelDevicePolicy< 32 >;
    RAJA::forall< KERNEL_POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, elementList.size() ),
                                   [=] GEOSX_DEVICE ( localIndex const index )
    {
      localIndex const k = elementList[ index ];

      real64 fLocal[ NUM_NODES_PER_ELEM ][ 3 ] = {{0}};
      real64 varLocal[ NUM_NODES_PER_ELEM ][ 3 ];

#if defined(CALCFEMSHAPE)
      real64 xLocal[ 8 ][ 3 ];
#endif

      for( localIndex a=0; a< NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
#if defined(CALCFEMSHAPE)
        LvArray::tensorOps::copy< 3 >( xLocal[ a ], X[ nodeIndex ] );
#endif

#if UPDATE_STRESS==2
        LvArray::tensorOps::scaledCopy< 3 >( varLocal[ a ], vel[ nodeIndex ], dt );
#else
        LvArray::tensorOps::copy< 3 >( varLocal[ a ], u[ nodeIndex ] );
#endif
      }

      //Compute Quadrature
      for( localIndex q = 0; q < NUM_QUADRATURE_POINTS; ++q )
      {

#if defined(CALCFEMSHAPE)
        real64 dNdX[ 8 ][ 3 ];
        real64 const detJ = FiniteElementShapeKernel::shapeFunctionDerivatives( q, xLocal, dNdX );
        #define DNDX dNdX
        #define DETJ detJ
#else //defined(CALCFEMSHAPE)
        #define DNDX dNdX[k][q]
        #define DETJ detJ( k, q )
#endif //defined(CALCFEMSHAPE)

        real64 strain[6] = {0};
        for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
        {
          strain[ 0 ] = strain[ 0 ] + DNDX[ a ][ 0 ] * varLocal[ a ][ 0 ];
          strain[ 1 ] = strain[ 1 ] + DNDX[ a ][ 1 ] * varLocal[ a ][ 1 ];
          strain[ 2 ] = strain[ 2 ] + DNDX[ a ][ 2 ] * varLocal[ a ][ 2 ];
          strain[ 3 ] = strain[ 3 ] + DNDX[ a ][ 2 ] * varLocal[ a ][ 1 ] + DNDX[ a ][ 1 ] * varLocal[ a ][ 2 ];
          strain[ 4 ] = strain[ 4 ] + DNDX[ a ][ 2 ] * varLocal[ a ][ 0 ] + DNDX[ a ][ 0 ] * varLocal[ a ][ 2 ];
          strain[ 5 ] = strain[ 5 ] + DNDX[ a ][ 1 ] * varLocal[ a ][ 0 ] + DNDX[ a ][ 0 ] * varLocal[ a ][ 1 ];
        }

#if UPDATE_STRESS == 2
        constitutive.SmallStrain( k, q, strain );
#else
        real64 stressLocal[ 6 ] = {0};
        constitutive.SmallStrainNoState( k, strain, stressLocal );
#endif

#if UPDATE_STRESS == 2
        LvArray::tensorOps::scaledCopy< 6 >( strain, constitutive.m_stress[ k ][ q ], -DETJ );
#elif UPDATE_STRESS == 1
        LvArray::tensorOps::copy< 6 >( strain, stressLocal );
        LvArray::tensorOps::add< 6 >( strain, constitutive.m_stress[ k ][ q ] );
        LvArray::tensorOps::scale< 6 >( strain, -DETJ );
#else
        LvArray::tensorOps::scaledCopy< 6 >( strain, stressLocal, -DETJ );
#endif

        for( localIndex a=0; a< NUM_NODES_PER_ELEM; ++a )
        {
          LvArray::tensorOps::plusSymAijBj< 3 >( fLocal[ a ], strain, DNDX[ a ] );
        }
      }    //quadrature loop

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        for( int b = 0; b < 3; ++b )
        {
          RAJA::atomicAdd< parallelDeviceAtomic >( &acc( nodeIndex, b ), fLocal[ a ][ b ] );
        }
      }
    } );
    return dt;
  }

#undef CALCFEMSHAPE
#undef DNDX
#undef DETJ
#undef UPDATE_STRESS

};

/**
 * @struct Structure to wrap templated function that implements the implicit time integration kernel.
 */
struct ImplicitKernel
{
// #if defined(GEOSX_USE_CUDA)
  #define CALCFEMSHAPE
// #endif

  template< int NUM_NODES_PER_ELEM, int NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const constitutiveRelation,
          localIndex const numElems,
          real64 const GEOSX_UNUSED_PARAM( dt ),
          arrayView4d< real64 const > const & _dNdX,
          arrayView2d< real64 const > const & _detJ,
          FiniteElementBase const * const GEOSX_UNUSED_PARAM( fe ),
          arrayView1d< integer const > const & GEOSX_UNUSED_PARAM( elemGhostRank ),
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView1d< globalIndex const > const & globalDofNumber,
          globalIndex const dofRankOffset,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & _X,
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & GEOSX_UNUSED_PARAM( disp ),
          arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & uhat,
          arrayView1d< R1Tensor const > const & GEOSX_UNUSED_PARAM( vtilde ),
          arrayView1d< R1Tensor const > const & GEOSX_UNUSED_PARAM( uhattilde ),
          arrayView2d< real64 const > const & density,
          arrayView1d< real64 const > const & fluidPressure,
          arrayView1d< real64 const > const & deltaFluidPressure,
          real64 const biotCoefficient,
          TimeIntegrationOption const GEOSX_UNUSED_PARAM( tiOption ),
          real64 const GEOSX_UNUSED_PARAM( stiffnessDamping ),
          real64 const GEOSX_UNUSED_PARAM( massDamping ),
          real64 const GEOSX_UNUSED_PARAM( newmarkBeta ),
          real64 const GEOSX_UNUSED_PARAM( newmarkGamma ),
          R1Tensor const & gravityVector,
          CRSMatrixView< real64, globalIndex const > const & matrix,
          arrayView1d< real64 > const & rhs )
  {
    GEOSX_MARK_FUNCTION;
    constexpr int NDIM = 3;

    // if the following is not static, then gcc8.1 gives a "error: use of 'this' in a constant expression"
    static constexpr int ndof = NDIM * NUM_NODES_PER_ELEM;

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelUpdates();

    arrayView3d< real64 const, solid::STRESS_USD > const & stress = constitutiveRelation->getStress();

  #if defined(CALCFEMSHAPE)
    GEOSX_UNUSED_VAR( _dNdX );
    GEOSX_UNUSED_VAR( _detJ );
    auto const & X = _X;
  #else
    GEOSX_UNUSED_VAR( _X );
    auto const & dNdX = _dNdX;
    auto const & detJ = _detJ;
  #endif

    RAJA::ReduceMax< parallelDeviceReduce, real64 > maxForce( 0 );
    RAJA::forall< parallelDevicePolicy< 32 > >( RAJA::TypedRangeSegment< localIndex >( 0, numElems ),
                                                [=] GEOSX_DEVICE ( localIndex const k )
    {
      globalIndex elementLocalDofIndex[ ndof ];
      real64 R[ ndof ] = { 0 };
      real64 dRdU[ ndof ][ ndof ] = {{ 0 }};

      R1Tensor uhat_local[NUM_NODES_PER_ELEM];

      real64 c[6][6];
      constitutive.getElasticStiffness( k, c );

    #if defined(CALCFEMSHAPE)
      real64 xLocal[ NUM_NODES_PER_ELEM ][ 3 ];
    #endif

      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const localNodeIndex = elemsToNodes( k, a );
        for( int i=0; i<NDIM; ++i )
        {
          elementLocalDofIndex[ a * NDIM + i] = globalDofNumber[localNodeIndex]+i;
        }
      }

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        uhat_local[ a ] = uhat[ nodeIndex ];
      #if defined(CALCFEMSHAPE)
        for( int i = 0; i < NDIM; ++i )
        { xLocal[ a ][ i ] = X[ nodeIndex ][ i ]; }
      #endif
      }

      for( int q = 0; q < NUM_QUADRATURE_POINTS; ++q )
      {

      #if defined(CALCFEMSHAPE)
        real64 dNdX[ 8 ][ 3 ];
        real64 const detJ_kq = FiniteElementShapeKernel::shapeFunctionDerivatives( q, xLocal, dNdX );
        #define DNDX dNdX
        #define DETJ detJ_kq
      #else //defined(CALCFEMSHAPE)
        #define DNDX dNdX[k][q]
        #define DETJ detJ( k, q )
      #endif //defined(CALCFEMSHAPE)

        real64 stress0[ 6 ] = LVARRAY_TENSOROPS_INIT_LOCAL_6( DETJ * stress[ k ][ q ] );
        if( !fluidPressure.empty() )
        {
          LvArray::tensorOps::symAddIdentity< 3 >( stress0, -DETJ * biotCoefficient * (fluidPressure[ k ] + deltaFluidPressure[ k ]) );
        }

        for( integer a = 0; a < NUM_NODES_PER_ELEM; ++a )
        {
          real64 temp[ 3 ];
          LvArray::tensorOps::symAijBj< 3 >( temp, stress0, DNDX[ a ] );

          maxForce.max( LvArray::tensorOps::maxAbsoluteEntry< 3 >( temp ) );

          R[ a * NDIM + 0 ] -= temp[ 0 ];
          R[ a * NDIM + 1 ] -= temp[ 1 ];
          R[ a * NDIM + 2 ] -= temp[ 2 ];

          for( int b = 0; b < NUM_NODES_PER_ELEM; ++b )
          {
            dRdU[ a * NDIM + 0 ][ b * NDIM + 0 ] -= ( c[ 0 ][ 0 ] * DNDX[ a ][ 0 ] * DNDX[ b ][ 0 ] +
                                                      c[ 5 ][ 5 ] * DNDX[ a ][ 1 ] * DNDX[ b ][ 1 ] +
                                                      c[ 4 ][ 4 ] * DNDX[ a ][ 2 ] * DNDX[ b ][ 2 ] ) * DETJ;
            dRdU[ a * NDIM + 0 ][ b * NDIM + 1 ] -= ( c[ 5 ][ 5 ] * DNDX[ a ][ 1 ] * DNDX[ b ][ 0 ] +
                                                      c[ 0 ][ 1 ] * DNDX[ a ][ 0 ] * DNDX[ b ][ 1 ] ) * DETJ;
            dRdU[ a * NDIM + 0 ][ b * NDIM + 2 ] -= ( c[ 4 ][ 4 ] * DNDX[ a ][ 2 ] * DNDX[ b ][ 0 ] +
                                                      c[ 0 ][ 2 ] * DNDX[ a ][ 0 ] * DNDX[ b ][ 2 ] ) * DETJ;

            dRdU[ a * NDIM + 1 ][ b * NDIM + 0 ] -= ( c[ 0 ][ 1 ] * DNDX[ a ][ 1 ] * DNDX[ b ][ 0 ] +
                                                      c[ 5 ][ 5 ] * DNDX[ a ][ 0 ] * DNDX[ b ][ 1 ] ) * DETJ;
            dRdU[ a * NDIM + 1 ][ b * NDIM + 1 ] -= ( c[ 5 ][ 5 ] * DNDX[ a ][ 0 ] * DNDX[ b ][ 0 ] +
                                                      c[ 1 ][ 1 ] * DNDX[ a ][ 1 ] * DNDX[ b ][ 1 ] +
                                                      c[ 3 ][ 3 ] * DNDX[ a ][ 2 ] * DNDX[ b ][ 2 ] ) * DETJ;
            dRdU[ a * NDIM + 1 ][ b * NDIM + 2 ] -= ( c[ 3 ][ 3 ] * DNDX[ a ][ 2 ] * DNDX[ b ][ 1 ] +
                                                      c[ 1 ][ 2 ] * DNDX[ a ][ 1 ] * DNDX[ b ][ 2 ] ) * DETJ;

            dRdU[ a * NDIM + 2 ][ b * NDIM + 0 ] -= ( c[ 0 ][ 2 ] * DNDX[ a ][ 2 ] * DNDX[ b ][ 0 ] +
                                                      c[ 4 ][ 4 ] * DNDX[ a ][ 0 ] * DNDX[ b ][ 2 ] ) * DETJ;
            dRdU[ a * NDIM + 2 ][ b * NDIM + 1 ] -= ( c[ 1 ][ 2 ] * DNDX[ a ][ 2 ] * DNDX[ b ][ 1 ] +
                                                      c[ 3 ][ 3 ] * DNDX[ a ][ 1 ] * DNDX[ b ][ 2 ] ) * DETJ;
            dRdU[ a * NDIM + 2 ][ b * NDIM + 2 ] -= ( c[ 4 ][ 4 ] * DNDX[ a ][ 0 ] * DNDX[ b ][ 0 ] +
                                                      c[ 3 ][ 3 ] * DNDX[ a ][ 1 ] * DNDX[ b ][ 1 ] +
                                                      c[ 2 ][ 2 ] * DNDX[ a ][ 2 ] * DNDX[ b ][ 2 ] ) * DETJ;
          }
        }

        R1Tensor gravityForce = gravityVector;
        gravityForce *= DETJ * density( k, q );
        R[ q * NDIM + 0 ] += gravityForce[ 0 ];
        R[ q * NDIM + 1 ] += gravityForce[ 1 ];
        R[ q * NDIM + 2 ] += gravityForce[ 2 ];
      }

      // TODO It is simpler to do this...try it.
      //  dRdU.Multiply(dof_np1,R);
      for( int a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        for( int b = 0; b < NUM_NODES_PER_ELEM; ++b )
        {
          for( int i = 0; i < NDIM; ++i )
          {
            for( int j = 0; j < NDIM; ++j )
            {
              R[ a * NDIM + i ] += dRdU[ a * NDIM + i ][ b * NDIM + j ] * uhat_local[ b ][ j ];
            }
          }
        }
      }

      for( int localNode = 0; localNode < NUM_NODES_PER_ELEM; ++localNode )
      {
        for( int dim = 0; dim < NDIM; ++dim )
        {
          localIndex const dof = LvArray::integerConversion< localIndex >( elementLocalDofIndex[ NDIM * localNode + dim ] - dofRankOffset );
          if( dof < 0 || dof >= matrix.numRows() ) continue;
          matrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                       elementLocalDofIndex,
                                                                       dRdU[ NDIM * localNode + dim ],
                                                                       NUM_NODES_PER_ELEM * NDIM );

          maxForce.max( fabs( R[ NDIM * localNode + dim ] ) );
          RAJA::atomicAdd< parallelDeviceAtomic >( &rhs[ dof ], R[ NDIM * localNode + dim ] );
        }
      }
    } );

    return maxForce.get();
  }

  #undef CALCFEMSHAPE
  #undef DNDX
  #undef DETJ
};

} // namespace SolidMechanicsLagrangianSSLEKernels

} // namespace geosx
