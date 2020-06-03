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

struct StressCalculationKernel
{
  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const constitutiveRelation,
          localIndex const numElems,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView4d< real64 const > const & dNdX,
          arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( detJ ),
          arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & uhat )
  {
    GEOSX_MARK_FUNCTION;

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelWrapper();

    arrayView3d< real64, solid::STRESS_USD > const & stress = constitutiveRelation->getStress();


//    using KERNEL_POLICY = parallelDevicePolicy< 256 >;
    using KERNEL_POLICY = parallelHostPolicy;
    RAJA::forall< KERNEL_POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, numElems ),
                                   [&] ( localIndex const k )
    {
      real64 uhat_local[ NUM_NODES_PER_ELEM ][ 3 ];
      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        LvArray::tensorOps::copy< 3 >( uhat_local[ a ], uhat[ elemsToNodes( k, a ) ] );
      }

      real64 c[ 6 ][ 6 ];
      constitutive.GetStiffness( k, c );

      //Compute Quadrature
      for( localIndex q = 0; q < NUM_QUADRATURE_POINTS; ++q )
      {
        for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
        {
          real64 const v0_x_dNdXa0 = uhat_local[ a ][ 0 ] * dNdX[ k ][ q ][ a ][ 0 ];
          real64 const v1_x_dNdXa1 = uhat_local[ a ][ 1 ] * dNdX[ k ][ q ][ a ][ 1 ];
          real64 const v2_x_dNdXa2 = uhat_local[ a ][ 2 ] * dNdX[ k ][ q ][ a ][ 2 ];

          stress( k, q, 0 ) += ( v0_x_dNdXa0 * c[ 0 ][ 0 ] + v1_x_dNdXa1 * c[ 0 ][ 1 ] + v2_x_dNdXa2*c[ 0 ][ 2 ] );
          stress( k, q, 1 ) += ( v0_x_dNdXa0 * c[ 1 ][ 0 ] + v1_x_dNdXa1 * c[ 1 ][ 1 ] + v2_x_dNdXa2*c[ 1 ][ 2 ] );
          stress( k, q, 2 ) += ( v0_x_dNdXa0 * c[ 2 ][ 0 ] + v1_x_dNdXa1 * c[ 2 ][ 1 ] + v2_x_dNdXa2*c[ 2 ][ 2 ] );
          stress( k, q, 3 ) += ( uhat_local[ a ][ 2 ] * dNdX[ k ][ q ][ a ][ 1 ] + uhat_local[ a ][ 1 ] * dNdX[ k ][ q ][ a ][ 2 ] ) * c[ 3 ][ 3 ];
          stress( k, q, 4 ) += ( uhat_local[ a ][ 2 ] * dNdX[ k ][ q ][ a ][ 0 ] + uhat_local[ a ][ 0 ] * dNdX[ k ][ q ][ a ][ 2 ] ) * c[ 4 ][ 4 ];
          stress( k, q, 5 ) += ( uhat_local[ a ][ 1 ] * dNdX[ k ][ q ][ a ][ 0 ] + uhat_local[ a ][ 0 ] * dNdX[ k ][ q ][ a ][ 1 ] ) * c[ 5 ][ 5 ];
        }
      }//quadrature loop

    } );

    return 0.0;
  }
};

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
          LvArray::SortedArrayView< localIndex const, localIndex > const & elementList,
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

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelWrapper();

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
      real64 xLocal[ NUM_NODES_PER_ELEM ][ 3 ];
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
  template< int NUM_NODES_PER_ELEM, int NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const constitutiveRelation,
          localIndex const numElems,
          real64 const dt,
          arrayView4d< real64 const > const & dNdX,
          arrayView2d< real64 const > const & detJ,
          FiniteElementBase const * const fe,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView1d< globalIndex const > const & globalDofNumber,
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & disp,
          arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & uhat,
          arrayView1d< R1Tensor const > const & vtilde,
          arrayView1d< R1Tensor const > const & uhattilde,
          arrayView2d< real64 const > const & density,
          arrayView1d< real64 const > const & fluidPressure,
          arrayView1d< real64 const > const & deltaFluidPressure,
          real64 const biotCoefficient,
          TimeIntegrationOption const tiOption,
          real64 const stiffnessDamping,
          real64 const massDamping,
          real64 const newmarkBeta,
          real64 const newmarkGamma,
          R1Tensor const & gravityVector,
          DofManager const * const GEOSX_UNUSED_PARAM( dofManager ),
          ParallelMatrix * const matrix,
          ParallelVector * const rhs )
  {
    GEOSX_MARK_FUNCTION;
    constexpr int dim = 3;

    // if the following is not static, then gcc8.1 gives a "error: use of 'this' in a constant expression"
    static constexpr int ndof = dim * NUM_NODES_PER_ELEM;
    RAJA::ReduceMax< serialReduce, double > maxForce( 0 );

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelWrapper();

    arrayView3d< real64 const, solid::STRESS_USD > const & stress = constitutiveRelation->getStress();

    RAJA::forall< serialPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, numElems ),
                                  [=] ( localIndex const k )
    {
      stackArray1d< globalIndex, ndof >       elementLocalDofIndex( ndof );
      stackArray1d< real64, ndof >            R( ndof );
      stackArray2d< real64, ndof *ndof >  dRdU( ndof, ndof );
      stackArray1d< real64, ndof >       element_dof_np1( ndof );

      stackArray1d< real64, ndof > R_InertiaMassDamping( ndof );
      stackArray2d< real64, ndof *ndof > dRdU_InertiaMassDamping( ndof, ndof );
      stackArray1d< real64, ndof > R_StiffnessDamping( ndof );
      stackArray2d< real64, ndof *ndof > dRdU_StiffnessDamping( ndof, ndof );

      R1Tensor u_local[NUM_NODES_PER_ELEM];
      R1Tensor uhat_local[NUM_NODES_PER_ELEM];
      R1Tensor vtilde_local[NUM_NODES_PER_ELEM];
      R1Tensor uhattilde_local[NUM_NODES_PER_ELEM];

      dRdU = 0.0;
      R = 0.0;

      dRdU_InertiaMassDamping = 0.0;
      R_InertiaMassDamping = 0.0;
      dRdU_StiffnessDamping = 0.0;
      R_StiffnessDamping = 0.0;

      real64 c[6][6];
      constitutive.GetStiffness( k, c );

      if( elemGhostRank[k] < 0 )
      {
        for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
        {

          localIndex localNodeIndex = elemsToNodes[k][a];

          for( int i=0; i<dim; ++i )
          {
            elementLocalDofIndex[static_cast< int >(a)*dim+i] = globalDofNumber[localNodeIndex]+i;

            // TODO must add last solution estimate for this to be valid
            element_dof_np1( static_cast< int >(a)*dim+i ) = disp[localNodeIndex][i];
          }
        }

        if( tiOption == TimeIntegrationOption::ImplicitDynamic )
        {
          GEOSX_ERROR( "Option not supported" );
          for( localIndex i = 0; i < NUM_NODES_PER_ELEM; ++i )
          {
            localIndex const nodeID = elemsToNodes( k, i );
            u_local[ i ] = disp[ nodeID ];
            uhat_local[ i ] = uhat[ nodeID ];
            vtilde_local[ i ] = vtilde[ i ];
            uhattilde_local[ i ] = uhattilde[ i ];
          }
        }
        else
        {
          for( localIndex i = 0; i < NUM_NODES_PER_ELEM; ++i )
          {
            localIndex const nodeID = elemsToNodes( k, i );
            u_local[ i ] = disp[ nodeID ];
            uhat_local[ i ] = uhat[ nodeID ];
          }
        }



        R1Tensor dNdXa;
        R1Tensor dNdXb;

        for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
        {
          const realT detJq = detJ[k][q];
          std::vector< double > const & N = fe->values( q );

          for( integer a=0; a<NUM_NODES_PER_ELEM; ++a )
          {
            //      realT const * const dNdXa = dNdX(q,a).Data();
            dNdXa = dNdX[k][q][a];

            for( integer b=0; b<NUM_NODES_PER_ELEM; ++b )
            {
              //        realT const * const dNdXb = dNdX(q,b).Data();
              dNdXb = dNdX[k][q][b];

              dRdU( a*dim+0, b*dim+0 ) -= ( c[0][0]*dNdXa[0]*dNdXb[0] + c[5][5]*dNdXa[1]*dNdXb[1] + c[4][4]*dNdXa[2]*dNdXb[2] ) * detJq;
              dRdU( a*dim+0, b*dim+1 ) -= ( c[5][5]*dNdXa[1]*dNdXb[0] + c[0][1]*dNdXa[0]*dNdXb[1] ) * detJq;
              dRdU( a*dim+0, b*dim+2 ) -= ( c[4][4]*dNdXa[2]*dNdXb[0] + c[0][2]*dNdXa[0]*dNdXb[2] ) * detJq;

              dRdU( a*dim+1, b*dim+0 ) -= ( c[0][1]*dNdXa[1]*dNdXb[0] + c[5][5]*dNdXa[0]*dNdXb[1] ) * detJq;
              dRdU( a*dim+1, b*dim+1 ) -= ( c[5][5]*dNdXa[0]*dNdXb[0] + c[1][1]*dNdXa[1]*dNdXb[1] + c[3][3]*dNdXa[2]*dNdXb[2] ) * detJq;
              dRdU( a*dim+1, b*dim+2 ) -= ( c[3][3]*dNdXa[2]*dNdXb[1] + c[1][2]*dNdXa[1]*dNdXb[2] ) * detJq;

              dRdU( a*dim+2, b*dim+0 ) -= ( c[0][2]*dNdXa[2]*dNdXb[0] + c[4][4]*dNdXa[0]*dNdXb[2] ) * detJq;
              dRdU( a*dim+2, b*dim+1 ) -= ( c[1][2]*dNdXa[2]*dNdXb[1] + c[3][3]*dNdXa[1]*dNdXb[2] ) * detJq;
              dRdU( a*dim+2, b*dim+2 ) -= ( c[4][4]*dNdXa[0]*dNdXb[0] + c[3][3]*dNdXa[1]*dNdXb[1] + c[2][2]*dNdXa[2]*dNdXb[2] ) * detJq;

              if( tiOption == TimeIntegrationOption::ImplicitDynamic )
              {

                real64 integrationFactor = density( k, q ) * N[a] * N[b] * detJq;
                real64 temp1 = ( massDamping * newmarkGamma/( newmarkBeta * dt ) + 1.0 / ( newmarkBeta * dt * dt ) )* integrationFactor;

                for( int i=0; i<dim; ++i )
                {
                  realT const acc = 1.0 / ( newmarkBeta * dt * dt ) * ( uhat_local[b][i] - uhattilde_local[b][i] );
                  realT const vel = vtilde_local[b][i] + newmarkGamma/( newmarkBeta * dt ) *( uhat_local[b][i] - uhattilde_local[b][i] );

                  dRdU_InertiaMassDamping( a*dim+i, b*dim+i ) -= temp1;
                  R_InertiaMassDamping( a*dim+i ) -= ( massDamping * vel + acc ) * integrationFactor;
                }
              }
            }
          }
        }

        for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
        {
          const realT detJq = detJ( k, q );
          real64 stress0[ 6 ] = LVARRAY_TENSOROPS_INIT_LOCAL_6( detJq * stress[ k ][ q ] );
          if( !fluidPressure.empty() )
          {
            LvArray::tensorOps::addIdentityToSymmetric< 3 >( stress0, -detJq * biotCoefficient * (fluidPressure[ k ] + deltaFluidPressure[ k ]) );
          }

          for( integer a=0; a<NUM_NODES_PER_ELEM; ++a )
          {
            real64 temp[ 3 ];
            LvArray::tensorOps::symAijBj< 3 >( temp, stress0, dNdX[ k ][ q ][ a ] );

            maxForce.max( LvArray::tensorOps::maxAbsoluteEntry< 3 >( temp ) );

            R( a*dim+0 ) -= temp[0];
            R( a*dim+1 ) -= temp[1];
            R( a*dim+2 ) -= temp[2];
          }

          R1Tensor gravityForce = gravityVector;
          gravityForce *= detJq * density( k, q );
          R( q*dim+0 ) += gravityForce[0];
          R( q*dim+1 ) += gravityForce[1];
          R( q*dim+2 ) += gravityForce[2];
        }


        // TODO It is simpler to do this...try it.
        //  dRdU.Multiply(dof_np1,R);
        for( integer a=0; a<NUM_NODES_PER_ELEM; ++a )
        {
          realT nodeForce = 0;
          for( integer b=0; b<NUM_NODES_PER_ELEM; ++b )
          {
            for( int i=0; i<dim; ++i )
            {
              for( int j=0; j<dim; ++j )
              {
                R( a*dim+i ) += dRdU( a*dim+i, b*dim+j ) * uhat_local[b][j];
              }
            }

            if( tiOption == TimeIntegrationOption::ImplicitDynamic )
            {
              for( int i=0; i<dim; ++i )
              {
                for( int j=0; j<dim; ++j )
                {
                  R_StiffnessDamping( a*dim+i ) += stiffnessDamping *
                                                   dRdU( a*dim+i,
                                                         b*dim+j ) *
                                                   ( vtilde_local[b][j] + newmarkGamma/(newmarkBeta * dt)*(uhat_local[b][j]-uhattilde_local[b][j]) );
                }
              }
            }

          }

          nodeForce = std::max( std::max( R( a*dim+0 ), R( a*dim+1 ) ), R( a*dim+2 ) );
          maxForce.max( fabs( nodeForce ) );
        }


        if( tiOption == TimeIntegrationOption::ImplicitDynamic )
        {
          GEOSX_ERROR( "NOT IMPLEMENTED" );
//          dRdU_StiffnessDamping = dRdU;
//          dRdU_StiffnessDamping.Scale( stiffnessDamping * newmarkGamma / ( newmarkBeta * dt ) );
//
//          dRdU += dRdU_InertiaMassDamping;
//          dRdU += dRdU_StiffnessDamping;
//          R    += R_InertiaMassDamping;
//          R    += R_StiffnessDamping;
        }

        // TODO remove local epetra objects, remove use of unwrapped()
        matrix->add( elementLocalDofIndex.data(), elementLocalDofIndex.data(), dRdU.data(), ndof, ndof );
        rhs->add( elementLocalDofIndex.data(), R.data(), ndof );
      }
    } );

    return maxForce.get();
  }
};

} // namespace SolidMechanicsLagrangianSSLEKernels

} // namespace geosx
