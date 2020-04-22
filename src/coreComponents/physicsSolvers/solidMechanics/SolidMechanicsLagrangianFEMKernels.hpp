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
#include "constitutive/ConstitutiveBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/solid/solidSelector.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "finiteElement/Kinematics.h"
#include "common/TimingMacros.hpp"

namespace geosx
{

/**
 * @enum timeIntegrationOption
 *
 * The options for time integration
 */
enum class timeIntegrationOption : int
{
  QuasiStatic,    //!< QuasiStatic
  ImplicitDynamic,//!< ImplicitDynamic
  ExplicitDynamic //!< ExplicitDynamic
};

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
    for( int j = 0; j < 3; ++j )
    {
      velocity( i, j ) += dt * acceleration( i, j );
      acceleration( i, j ) = 0;
    }
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
    for( int j = 0; j < 3; ++j )
    {
      acceleration( a, j ) /= mass[ a ];
      velocity( a, j ) += dt * acceleration( a, j );
    }
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
    for( int j = 0; j < 3; ++j )
    {
      uhat( i, j ) = velocity( i, j ) * dt;
      u( i, j ) += uhat( i, j );
    }
  } );
}

template< int N >
inline void Integrate( const R2SymTensor & fieldvar,
                       arraySlice1d< R1Tensor const > const & dNdX,
                       real64 const & detJ,
                       real64 const & detF,
                       const R2Tensor & fInv,
                       R1Tensor * GEOSX_RESTRICT const result )
{
  real64 const integrationFactor = detJ * detF;

  R2Tensor P;
  P.AijBkj( fieldvar, fInv );
  P *= integrationFactor;

  for( int a=0; a<N; ++a )    // loop through all shape functions in element
  {
    result[a].minusAijBj( P, dNdX[a] );
  }
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

  constitutive::constitutiveUpdatePassThru( constitutiveRelation, [&]( auto & constitutive )
  {
    using CONSTITUTIVE_TYPE = TYPEOFREF( constitutive );
    if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==8 )
    {
      rval = KERNELWRAPPER::template Launch< 8, 8, CONSTITUTIVE_TYPE >( &constitutive, std::forward< PARAMS >( params )... );
    }
    else if( NUM_NODES_PER_ELEM==4 && NUM_QUADRATURE_POINTS==1 )
    {
      rval = KERNELWRAPPER::template Launch< 4, 1, CONSTITUTIVE_TYPE >( &constitutive, std::forward< PARAMS >( params )... );
    }
  } );
  return rval;
}

/**
 * @struct Structure to wrap templated function that implements the explicit time integration kernel.
 */
struct ExplicitKernel
{
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
          arrayView3d< R1Tensor const > const & dNdX,
          arrayView2d< real64 const > const & detJ,
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u,
          arrayView2d< real64 const, nodes::VELOCITY_USD > const & vel,
          arrayView2d< real64, nodes::ACCELERATION_USD > const & acc,
          real64 const dt )
  {

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelWrapper();

    forAll< serialPolicy >( elementList.size(), [=] ( localIndex const i )
    {
      localIndex const k = elementList[ i ];
      R1Tensor v_local[NUM_NODES_PER_ELEM];
      R1Tensor u_local[NUM_NODES_PER_ELEM];
      R1Tensor f_local[NUM_NODES_PER_ELEM];

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        u_local[ a ] = u[ nodeIndex ];
        v_local[ a ] = vel[ nodeIndex ];
      }

      //Compute Quadrature
      for( localIndex q = 0; q<NUM_QUADRATURE_POINTS; ++q )
      {
        R2Tensor dUhatdX, dUdX;
        CalculateGradients< NUM_NODES_PER_ELEM >( dUhatdX, dUdX, v_local, u_local, dNdX[k][q] );
        dUhatdX *= dt;

        R2Tensor F, Ldt, fInv;

        // calculate du/dX
        F = dUhatdX;
        F *= 0.5;
        F += dUdX;
        F.PlusIdentity( 1.0 );
        fInv.Inverse( F );

        // chain rule: calculate dv/du = dv/dX * dX/du
        Ldt.AijBjk( dUhatdX, fInv );

        // calculate gradient (end of step)
        F = dUhatdX;
        F += dUdX;
        F.PlusIdentity( 1.0 );
        real64 detF = F.Det();
        fInv.Inverse( F );


        R2Tensor Rot;
        R2SymTensor Dadt;
        HughesWinget( Rot, Dadt, Ldt );

        constitutive.HypoElastic( k, q, Dadt.Data(), Rot );

        Integrate< NUM_NODES_PER_ELEM >( constitutive.m_stress[k][q], dNdX[k][q], detJ[k][q], detF, fInv, f_local );
      }//quadrature loop

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        RAJA::atomicAdd< serialAtomic >( &acc( nodeIndex, 0 ), f_local[ a ][ 0 ] );
        RAJA::atomicAdd< serialAtomic >( &acc( nodeIndex, 1 ), f_local[ a ][ 1 ] );
        RAJA::atomicAdd< serialAtomic >( &acc( nodeIndex, 2 ), f_local[ a ][ 2 ] );
      }

    } );

    return dt;
  }


  static inline real64
  CalculateSingleNodalForce( localIndex const k,
                             localIndex const targetNode,
                             localIndex const numQuadraturePoints,
                             arrayView3d< R1Tensor const > const & dNdX,
                             arrayView2d< real64 const > const & detJ,
                             arrayView3d< real64 const, solid::STRESS_USD > const & stress,
                             R1Tensor & force )
  {
    GEOSX_MARK_FUNCTION;
    localIndex const & a = targetNode;

    //Compute Quadrature
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      force[ 0 ] -= ( stress( k, q, 0 ) * dNdX( k, q, a )[ 0 ] +
                      stress( k, q, 5 ) * dNdX( k, q, a )[ 1 ] +
                      stress( k, q, 4 ) * dNdX( k, q, a )[ 2 ] ) * detJ( k, q );
      force[ 1 ] -= ( stress( k, q, 5 ) * dNdX( k, q, a )[ 0 ] +
                      stress( k, q, 1 ) * dNdX( k, q, a )[ 1 ] +
                      stress( k, q, 3 ) * dNdX( k, q, a )[ 2 ] ) * detJ( k, q );
      force[ 2 ] -= ( stress( k, q, 4 ) * dNdX( k, q, a )[ 0 ] +
                      stress( k, q, 3 ) * dNdX( k, q, a )[ 1 ] +
                      stress( k, q, 2 ) * dNdX( k, q, a )[ 2 ] ) * detJ( k, q );

    }//quadrature loop

    return 0;
  }

};

/**
 * @struct Structure to wrap templated function that implements the implicit time integration kernel.
 */
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
          arrayView3d< R1Tensor const > const & GEOSX_UNUSED_PARAM( dNdX ),
          arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( detJ ),
          FiniteElementBase const * const GEOSX_UNUSED_PARAM( fe ),
          arrayView1d< integer const > const & GEOSX_UNUSED_PARAM( elemGhostRank ),
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & GEOSX_UNUSED_PARAM( elemsToNodes ),
          arrayView1d< globalIndex const > const & GEOSX_UNUSED_PARAM( globalDofNumber ),
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & GEOSX_UNUSED_PARAM( disp ),
          arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & GEOSX_UNUSED_PARAM( uhat ),
          arrayView1d< R1Tensor const > const & GEOSX_UNUSED_PARAM( vtilde ),
          arrayView1d< R1Tensor const > const & GEOSX_UNUSED_PARAM( uhattilde ),
          arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( density ),
          arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( fluidPressure ),
          arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( deltaFluidPressure ),
          real64 const GEOSX_UNUSED_PARAM( biotCoefficient ),
          timeIntegrationOption const GEOSX_UNUSED_PARAM( tiOption ),
          real64 const GEOSX_UNUSED_PARAM( stiffnessDamping ),
          real64 const GEOSX_UNUSED_PARAM( massDamping ),
          real64 const GEOSX_UNUSED_PARAM( newmarkBeta ),
          real64 const GEOSX_UNUSED_PARAM( newmarkGamma ),
          R1Tensor const & GEOSX_UNUSED_PARAM( gravityVector ),
          DofManager const * const GEOSX_UNUSED_PARAM( dofManager ),
          ParallelMatrix * const GEOSX_UNUSED_PARAM( matrix ),
          ParallelVector * const GEOSX_UNUSED_PARAM( rhs ) )
  {
    GEOSX_ERROR( "SolidMechanicsLagrangianFEM::ImplicitElementKernelWrapper::Launch() not implemented" );
    return 0;
  }

};

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
