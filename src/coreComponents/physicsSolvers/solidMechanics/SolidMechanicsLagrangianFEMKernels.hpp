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
 * @file SolidMechanicsLagrangianFEMKernels.hpp
 */

#pragma once

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"
#include "finiteElement/elementFormulations/TrilinearHexahedronShapeFunctionKernel.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "TimeIntegrationOption.hpp"

namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

inline void velocityUpdate( arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
                            arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
                            real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  localIndex const N = acceleration.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOSX_DEVICE ( localIndex const i )
  {
    LvArray::tensorOps::scaledAdd< 3 >( velocity[ i ], acceleration[ i ], dt );
    LvArray::tensorOps::fill< 3 >( acceleration[ i ], 0 );
  } );
}

inline void velocityUpdate( arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
                            arrayView1d< real64 const > const & mass,
                            arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
                            real64 const dt,
                            SortedArrayView< localIndex const > const & indices )
{
  GEOSX_MARK_FUNCTION;

  forAll< parallelDevicePolicy<> >( indices.size(), [=] GEOSX_DEVICE ( localIndex const i )
  {
    localIndex const a = indices[ i ];
    LvArray::tensorOps::scale< 3 >( acceleration[ a ], 1.0 / mass[ a ] );
    LvArray::tensorOps::scaledAdd< 3 >( velocity[ a ], acceleration[ a ], dt );
  } );
}

inline void displacementUpdate( arrayView2d< real64 const, nodes::VELOCITY_USD > const & velocity,
                                arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & uhat,
                                arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u,
                                real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  localIndex const N = velocity.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOSX_DEVICE ( localIndex const i )
  {
    LvArray::tensorOps::scaledCopy< 3 >( uhat[ i ], velocity[ i ], dt );
    LvArray::tensorOps::add< 3 >( u[ i ], uhat[ i ] );
  } );
}


/**
 * @brief Function to select which templated kernel function to call.
 * @tparam KERNELWRAPPER A struct or class that contains the following method
 *  "Launch<NUM_NODES_PER_ELEM, NUM_QUADRATURE_POINTS, CONSTITUTIVE_TYPE>( CONSTITUTIVE_TYPE *, PARAMS... )"
 * @tparam PARAMS Variadic parameter pack to pass arguments to Launch function.
 * @param NUM_NODES_PER_ELEM The number of nodes in an element.
 * @param NUM_QUADRATURE_POINTS The number of quadrature points in an element.
 * @param params Variadic parameter list to hold all parameters that are forwarded to the kernel function.
 * @return Depends on the kernel.
 */
template< typename KERNELWRAPPER, typename ... PARAMS >
inline real64
ElementKernelLaunchSelector( localIndex NUM_NODES_PER_ELEM,
                             localIndex NUM_QUADRATURE_POINTS,
                             constitutive::ConstitutiveBase * const constitutiveRelation,
                             PARAMS && ... params )
{
  real64 rval = 0;

  using namespace constitutive;

  ConstitutivePassThru< SolidBase >::Execute( constitutiveRelation,
                                              [&]( auto * const constitutive )
  {
    rval = finiteElementLaunchDispatch< KERNELWRAPPER >( NUM_NODES_PER_ELEM, NUM_QUADRATURE_POINTS, &constitutive, std::forward< PARAMS >( params )... );
  } );
  return rval;
}

/**
 * @struct Structure to wrap templated function that implements the explicit time integration kernel.
 */
struct ExplicitKernel
{

#if defined(GEOSX_USE_CUDA)
  #define CALCFEMSHAPE
#endif


  template< int N, int USD >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static
  void Integrate( arraySlice1d< real64 const, USD > const & fieldVar,
  #if defined(CALCFEMSHAPE)
                  real64 const (&dNdX)[ N ][ 3 ],
  #else
                  arraySlice2d< real64 const > const & dNdX,
  #endif
                  real64 const detJ,
                  real64 const detF,
                  real64 const ( &fInv )[ 3 ][ 3 ],
                  real64 ( & result )[ N ][ 3 ] )
  {
    GEOSX_ASSERT_EQ( fieldVar.size(), 6 );

    real64 const integrationFactor = -detJ * detF;

    real64 P[ 3 ][ 3 ];
    LvArray::tensorOps::symAikBjk< 3 >( P, fieldVar, fInv );
    LvArray::tensorOps::scale< 3, 3 >( P, integrationFactor );

    for( int a = 0; a < N; ++a )    // loop through all shape functions in element
    {
      LvArray::tensorOps::plusAijBj< 3, 3 >( result[ a ], P, dNdX[ a ] );
    }
  }

  /**
   * @brief Launch of the element processing kernel for explicit time integration.
   * @tparam NUM_NODES_PER_ELEM The number of nodes/dof per element.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam CONSTITUTIVE_TYPE the type of the constitutive relation that is being used.
   * @param A pointer to the constitutive relation that is being used.
   * @param elementList The list of elements to be processed
   * @param elemsToNodes The map from the elements to the nodes that form that element.
   * @param dNdX The derivatives of the shape functions wrt the reference configuration.
   * @param detJ The determinant of the transformation matrix (Jacobian) to the parent element.
   * @param u The nodal array of total displacements.
   * @param vel The nodal array of velocity.
   * @param acc The nodal array of force/acceleration.
   * @param stress The stress at each element quadrature point
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


#if defined(CALCFEMSHAPE)
    GEOSX_UNUSED_VAR( dNdX );
    GEOSX_UNUSED_VAR( detJ );
#else
    GEOSX_UNUSED_VAR( X );
#endif



    typename CONSTITUTIVE_TYPE::KernelWrapper constitutive = constitutiveRelation->createKernelUpdates();

    using KERNEL_POLICY = parallelDevicePolicy< 32 >;
    RAJA::forall< KERNEL_POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, elementList.size() ),
                                   [=] GEOSX_DEVICE ( localIndex const index )
    {
      localIndex const k = elementList[ index ];

      real64 v_local[ NUM_NODES_PER_ELEM ][ 3 ];
      real64 u_local[ NUM_NODES_PER_ELEM ][ 3 ];
      real64 f_local[ NUM_NODES_PER_ELEM ][ 3 ] = { { 0 } };
#if defined(CALCFEMSHAPE)
      real64 X_local[ NUM_NODES_PER_ELEM ][ 3 ];
#endif

      for( localIndex a=0; a< NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        LvArray::tensorOps::copy< 3 >( v_local[ a ], vel[ nodeIndex ] );
        LvArray::tensorOps::copy< 3 >( u_local[ a ], u[ nodeIndex ] );
#if defined(CALCFEMSHAPE)
        LvArray::tensorOps::copy< 3 >( X_local[ a ], X[ nodeIndex ] );
#endif
      }

      //Compute Quadrature
      for( localIndex q = 0; q<NUM_QUADRATURE_POINTS; ++q )
      {
#if defined(CALCFEMSHAPE)
        real64 dNdX[ 8 ][ 3 ];
        real64 const detJ = TrilinearHexahedronShapeFunctionKernel::shapeFunctionDerivatives( q, X_local, dNdX );
        #define DNDX dNdX
        #define DETJ detJ
#else
        #define DNDX dNdX[k][q]
        #define DETJ detJ( k, q )
#endif
        real64 dUhatdX[ 3 ][ 3 ], dUdX[ 3 ][ 3 ];
        CalculateGradients< NUM_NODES_PER_ELEM >( dUhatdX, dUdX, v_local, u_local, DNDX );
        LvArray::tensorOps::scale< 3, 3 >( dUhatdX, dt );

        real64 F[ 3 ][ 3 ], Ldt[ 3 ][ 3 ], fInv[ 3 ][ 3 ];

        // calculate du/dX
        LvArray::tensorOps::scaledCopy< 3, 3 >( F, dUhatdX, 0.5 );
        LvArray::tensorOps::add< 3, 3 >( F, dUdX );
        LvArray::tensorOps::addIdentity< 3 >( F, 1.0 );
        LvArray::tensorOps::invert< 3 >( fInv, F );

        // chain rule: calculate dv/du = dv/dX * dX/du
        LvArray::tensorOps::AikBkj< 3, 3, 3 >( Ldt, dUhatdX, fInv );

        // calculate gradient (end of step)
        LvArray::tensorOps::copy< 3, 3 >( F, dUhatdX );
        LvArray::tensorOps::add< 3, 3 >( F, dUdX );
        LvArray::tensorOps::addIdentity< 3 >( F, 1.0 );
        real64 const detF = LvArray::tensorOps::invert< 3 >( fInv, F );

        real64 Rot[ 3 ][ 3 ];
        real64 Dadt[ 6 ];
        HughesWinget( Rot, Dadt, Ldt );

        constitutive.HypoElastic( k, q, Dadt, Rot );

        Integrate< NUM_NODES_PER_ELEM >( constitutive.m_stress[k][q].toSliceConst(),
                                         DNDX,
                                         DETJ,
                                         detF,
                                         fInv,
                                         f_local );
      }    //quadrature loop

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        RAJA::atomicAdd< parallelDeviceAtomic >( &acc( nodeIndex, 0 ), f_local[ a ][ 0 ] );
        RAJA::atomicAdd< parallelDeviceAtomic >( &acc( nodeIndex, 1 ), f_local[ a ][ 1 ] );
        RAJA::atomicAdd< parallelDeviceAtomic >( &acc( nodeIndex, 2 ), f_local[ a ][ 2 ] );
      }

    } );

    return dt;
  }

#undef CALCFEMSHAPE
#undef DNDX
#undef DETJ


  static inline real64
  CalculateSingleNodalForce( localIndex const k,
                             localIndex const targetNode,
                             localIndex const numQuadraturePoints,
                             arrayView4d< real64 const > const & dNdX,
                             arrayView2d< real64 const > const & detJ,
                             arrayView3d< real64 const, solid::STRESS_USD > const & stress,
                             R1Tensor & force )
  {
    GEOSX_MARK_FUNCTION;
    localIndex const & a = targetNode;

    //Compute Quadrature
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      force[ 0 ] -= ( stress( k, q, 0 ) * dNdX( k, q, a, 0 ) +
                      stress( k, q, 5 ) * dNdX( k, q, a, 1 ) +
                      stress( k, q, 4 ) * dNdX( k, q, a, 2 ) ) * detJ( k, q );
      force[ 1 ] -= ( stress( k, q, 5 ) * dNdX( k, q, a, 0 ) +
                      stress( k, q, 1 ) * dNdX( k, q, a, 1 ) +
                      stress( k, q, 3 ) * dNdX( k, q, a, 2 ) ) * detJ( k, q );
      force[ 2 ] -= ( stress( k, q, 4 ) * dNdX( k, q, a, 0 ) +
                      stress( k, q, 3 ) * dNdX( k, q, a, 1 ) +
                      stress( k, q, 2 ) * dNdX( k, q, a, 2 ) ) * detJ( k, q );

    }//quadrature loop

    return 0;
  }

};

struct ImplicitKernel
{

  /**
   * @brief Launch of the element processing kernel for implicit time integration.
   * @tparam NUM_NODES_PER_ELEM The number of nodes/dof per element.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam CONSTITUTIVE_TYPE the type of the constitutive relation that is being used.
   * @param constitutiveRelation A pointer to the constitutive relation that is being used.
   * @param numElems The number of elements the kernel will process.
   * @param dt The timestep.
   * @param dNdX The derivatives of the shape functions wrt the reference configuration.
   * @param detJ The determinant of the transformation matrix (Jacobian) to the parent element.
   * @param fe A pointer to the finite element class used in this kernel.
   * @param elemGhostRank An array containing the values of the owning ranks for ghost elements.
   * @param elemsToNodes The map from the elements to the nodes that form that element.
   * @param globalDofNumber The map from localIndex to the globalDOF number.
   * @param disp The array of total displacements.
   * @param uhat The array of incremental displacements (displacement for this step).
   * @param vtilde The array for the velocity predictor.
   * @param uhattilde The array for the incremental displacement predictor.
   * @param density The array containing the density
   * @param fluidPressure Array containing element fluid pressure at the beginning of the step.
   * @param deltaFluidPressure Array containing the change in element fluid pressure over this step.
   * @param biotCoefficient The biotCoefficient used to calculate effective stress.
   * @param tiOption The time integration option used for the integration.
   * @param stiffnessDamping The stiffness damping coefficient for the Newmark method assuming Rayleigh damping.
   * @param massDamping The mass damping coefficient for the Newmark method assuming Rayleigh damping.
   * @param newmarkBeta The value of \beta in the Newmark update.
   * @param newmarkGamma The value of \gamma in the Newmark update.
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix sparse matrix containing the derivatives of the residual wrt displacement
   * @param rhs parallel vector containing the global residual
   * @return The maximum nodal force contribution from all elements.
   */
  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const GEOSX_UNUSED_PARAM( constitutiveRelation ),
          localIndex const GEOSX_UNUSED_PARAM( numElems ),
          real64 const GEOSX_UNUSED_PARAM( dt ),
          arrayView4d< real64 const > const & GEOSX_UNUSED_PARAM( dNdX ),
          arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( detJ ),
          FiniteElementBase const * const GEOSX_UNUSED_PARAM( fe ),
          arrayView1d< integer const > const & GEOSX_UNUSED_PARAM( elemGhostRank ),
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & GEOSX_UNUSED_PARAM( elemsToNodes ),
          arrayView1d< globalIndex const > const & GEOSX_UNUSED_PARAM( globalDofNumber ),
          globalIndex const GEOSX_UNUSED_PARAM( dofRankOffset ),
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & GEOSX_UNUSED_PARAM( X ),
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & GEOSX_UNUSED_PARAM( disp ),
          arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & GEOSX_UNUSED_PARAM( uhat ),
          arrayView1d< R1Tensor const > const & GEOSX_UNUSED_PARAM( vtilde ),
          arrayView1d< R1Tensor const > const & GEOSX_UNUSED_PARAM( uhattilde ),
          arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( density ),
          arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( fluidPressure ),
          arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( deltaFluidPressure ),
          real64 const GEOSX_UNUSED_PARAM( biotCoefficient ),
          TimeIntegrationOption const GEOSX_UNUSED_PARAM( tiOption ),
          real64 const GEOSX_UNUSED_PARAM( stiffnessDamping ),
          real64 const GEOSX_UNUSED_PARAM( massDamping ),
          real64 const GEOSX_UNUSED_PARAM( newmarkBeta ),
          real64 const GEOSX_UNUSED_PARAM( newmarkGamma ),
          R1Tensor const & GEOSX_UNUSED_PARAM( gravityVector ),
          CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( matrix ),
          arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( rhs ) )
  {
    GEOSX_ERROR( "SolidMechanicsLagrangianFEM::CRSImplicitKernel::Launch() not implemented" );
    return 0;
  }
};


} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
