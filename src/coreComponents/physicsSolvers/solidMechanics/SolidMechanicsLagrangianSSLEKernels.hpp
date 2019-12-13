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
#include "SolidMechanicsLagrangianFEMKernels.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

#include "RAJA/RAJA.hpp"

#include "finiteElement/Kinematics.h"
#include "common/TimingMacros.hpp"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

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
          arrayView2d< localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM > const & elemsToNodes,
          arrayView3d< R1Tensor const> const & dNdX,
          arrayView2d<real64 const> const & GEOSX_UNUSED_ARG( detJ ),
          arrayView1d< R1Tensor const > const & u )
  {
    GEOSX_MARK_FUNCTION;

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelWrapper();

    arrayView2d< R2SymTensor > const & stress = constitutiveRelation->getStress();

    using KERNEL_POLICY = parallelDevicePolicy< 256 >;
    RAJA::forall< KERNEL_POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, numElems ),
                                   GEOSX_DEVICE_LAMBDA ( localIndex const k )
    {
      real64 u_local[ NUM_NODES_PER_ELEM ][ 3 ];

      for ( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        for ( int b = 0; b < 3; ++b )
        {
          u_local[ a ][ b ] = u[ elemsToNodes[ k ][ a ] ][ b ];
        }
      }

      real64 c[ 6 ][ 6 ];
      constitutive.GetStiffness( k, c );

      //Compute Quadrature
      for ( localIndex q = 0; q < NUM_QUADRATURE_POINTS; ++q )
      {
        real64 * const restrict p_stress = stress[ k ][ q ].Data();
        for ( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
        {
          real64 const v0_x_dNdXa0 = u_local[ a ][ 0 ] * dNdX[ k ][ q ][ a ][ 0 ];
          real64 const v1_x_dNdXa1 = u_local[ a ][ 1 ] * dNdX[ k ][ q ][ a ][ 1 ];
          real64 const v2_x_dNdXa2 = u_local[ a ][ 2 ] * dNdX[ k ][ q ][ a ][ 2 ];

          p_stress[ 0 ] += ( v0_x_dNdXa0 * c[ 0 ][ 0 ] + v1_x_dNdXa1 * c[ 0 ][ 1 ] + v2_x_dNdXa2*c[ 0 ][ 2 ] ) ;
          p_stress[ 2 ] += ( v0_x_dNdXa0 * c[ 1 ][ 0 ] + v1_x_dNdXa1 * c[ 1 ][ 1 ] + v2_x_dNdXa2*c[ 1 ][ 2 ] ) ;
          p_stress[ 5 ] += ( v0_x_dNdXa0 * c[ 2 ][ 0 ] + v1_x_dNdXa1 * c[ 2 ][ 1 ] + v2_x_dNdXa2*c[ 2 ][ 2 ] ) ;
          p_stress[ 4 ] += ( u_local[ a ][ 2 ] * dNdX[ k ][ q ][ a ][ 1 ] + u_local[ a ][ 1 ] * dNdX[ k ][ q ][ a ][ 2 ] ) * c[ 3 ][ 3 ] ;
          p_stress[ 3 ] += ( u_local[ a ][ 2 ] * dNdX[ k ][ q ][ a ][ 0 ] + u_local[ a ][ 0 ] * dNdX[ k ][ q ][ a ][ 2 ] ) * c[ 4 ][ 4 ] ;
          p_stress[ 1 ] += ( u_local[ a ][ 1 ] * dNdX[ k ][ q ][ a ][ 0 ] + u_local[ a ][ 0 ] * dNdX[ k ][ q ][ a ][ 1 ] ) * c[ 5 ][ 5 ] ;
        }
      }//quadrature loop

    });

    return 0.0;
  }
};

/**
 * @struct Structure to wrap templated function that implements the explicit time integration kernels.
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
   * @param meanStress The mean stress at each element quadrature point
   * @param devStress The deviator stress at each element quadrature point.
   * @param dt The timestep
   * @return The achieved timestep.
   */
  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const constitutiveRelation,
          LvArray::SortedArrayView<localIndex const, localIndex> const & elementList,
          arrayView2d<localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM> const & elemsToNodes,
          arrayView3d< R1Tensor const> const & dNdX,
          arrayView2d<real64 const> const & detJ,
          arrayView1d<R1Tensor const> const & GEOSX_UNUSED_ARG( u ),
          arrayView1d<R1Tensor const> const & vel,
          arrayView1d<R1Tensor> const & acc,
          arrayView2d<R2SymTensor> const & stress,
          real64 const dt )
  {
    GEOSX_MARK_FUNCTION;

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelWrapper();

    using KERNEL_POLICY = parallelDevicePolicy< 256 >;
    RAJA::forall< KERNEL_POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, elementList.size() ),
                                   GEOSX_DEVICE_LAMBDA ( localIndex const i )
    {
      localIndex const k = elementList[ i ];

      real64 v_local[ NUM_NODES_PER_ELEM ][ 3 ];
      real64 f_local[ NUM_NODES_PER_ELEM ][ 3 ] = {};

      for ( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        for ( int b = 0; b < 3; ++b )
        {
          v_local[ a ][ b ] = vel[ elemsToNodes[ k ][ a ] ][ b ];
        }
      }

      real64 c[ 6 ][ 6 ];
      constitutive.GetStiffness( k, c );

      //Compute Quadrature
      for ( localIndex q = 0; q < NUM_QUADRATURE_POINTS; ++q )
      {
        real64 p_stress[ 6 ] = { 0 };
        for ( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
        {
          real64 const v0_x_dNdXa0 = v_local[ a ][ 0 ] * dNdX[ k ][ q ][ a ][ 0 ];
          real64 const v1_x_dNdXa1 = v_local[ a ][ 1 ] * dNdX[ k ][ q ][ a ][ 1 ];
          real64 const v2_x_dNdXa2 = v_local[ a ][ 2 ] * dNdX[ k ][ q ][ a ][ 2 ];

          p_stress[ 0 ] += ( v0_x_dNdXa0 * c[ 0 ][ 0 ] + v1_x_dNdXa1 * c[ 0 ][ 1 ] + v2_x_dNdXa2*c[ 0 ][ 2 ] ) * dt;
          p_stress[ 1 ] += ( v0_x_dNdXa0 * c[ 1 ][ 0 ] + v1_x_dNdXa1 * c[ 1 ][ 1 ] + v2_x_dNdXa2*c[ 1 ][ 2 ] ) * dt;
          p_stress[ 2 ] += ( v0_x_dNdXa0 * c[ 2 ][ 0 ] + v1_x_dNdXa1 * c[ 2 ][ 1 ] + v2_x_dNdXa2*c[ 2 ][ 2 ] ) * dt;
          p_stress[ 3 ] += ( v_local[ a ][ 2 ] * dNdX[ k ][ q ][ a ][ 1 ] + v_local[ a ][ 1 ] * dNdX[ k ][ q ][ a ][ 2 ] ) * c[ 3 ][ 3 ] * dt;
          p_stress[ 4 ] += ( v_local[ a ][ 2 ] * dNdX[ k ][ q ][ a ][ 0 ] + v_local[ a ][ 0 ] * dNdX[ k ][ q ][ a ][ 2 ] ) * c[ 4 ][ 4 ] * dt;
          p_stress[ 5 ] += ( v_local[ a ][ 1 ] * dNdX[ k ][ q ][ a ][ 0 ] + v_local[ a ][ 0 ] * dNdX[ k ][ q ][ a ][ 1 ] ) * c[ 5 ][ 5 ] * dt;
        }

        real64 * const restrict p_Stress = stress[ k ][ q ].Data();
        p_Stress[ 0 ] += p_stress[ 0 ];
        p_Stress[ 2 ] += p_stress[ 1 ];
        p_Stress[ 5 ] += p_stress[ 2 ];
        p_Stress[ 4 ] += p_stress[ 3 ];
        p_Stress[ 3 ] += p_stress[ 4 ];
        p_Stress[ 1 ] += p_stress[ 5 ];

        for ( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
        {
          f_local[ a ][ 0 ] -= ( p_Stress[ 1 ] * dNdX[ k ][ q ][ a ][ 1 ]
                               + p_Stress[ 3 ] * dNdX[ k ][ q ][ a ][ 2 ]
                               + p_Stress[ 0 ] * dNdX[ k ][ q ][ a ][ 0 ] ) * detJ[ k ][ q ];
          f_local[ a ][ 1 ] -= ( p_Stress[ 1 ] * dNdX[ k ][ q ][ a ][ 0 ]
                               + p_Stress[ 4 ] * dNdX[ k ][ q ][ a ][ 2 ]
                               + p_Stress[ 2 ] * dNdX[ k ][ q ][ a ][ 1 ] ) * detJ[ k ][ q ];
          f_local[ a ][ 2 ] -= ( p_Stress[ 3 ] * dNdX[ k ][ q ][ a ][ 0 ]
                               + p_Stress[ 4 ] * dNdX[ k ][ q ][ a ][ 1 ]
                               + p_Stress[ 5 ] * dNdX[ k ][ q ][ a ][ 2 ] ) * detJ[ k ][ q ];
        }
      }//quadrature loop

      for ( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        for ( int b = 0; b < 3; ++b )
        {
          RAJA::atomicAdd<RAJA::auto_atomic>( &acc[ elemsToNodes[ k ][ a ] ][ b ], f_local[ a ][ b ] );
        }
      }
    });

    return dt;
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
  template< int NUM_NODES_PER_ELEM, int NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const constitutiveRelation,
          localIndex const numElems,
          real64 const dt,
          arrayView3d<R1Tensor const> const & dNdX,
          arrayView2d<real64 const > const& detJ,
          FiniteElementBase const * const fe,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView2d< localIndex const, CellBlock::NODE_MAP_UNIT_STRIDE_DIM > const & elemsToNodes,
          arrayView1d< globalIndex const > const & globalDofNumber,
          arrayView1d< R1Tensor const > const & disp,
          arrayView1d< R1Tensor const > const & uhat,
          arrayView1d< R1Tensor const > const & vtilde,
          arrayView1d< R1Tensor const > const & uhattilde,
          arrayView1d< real64 const > const & density,
          arrayView1d< real64 const > const & fluidPressure,
          arrayView1d< real64 const > const & deltaFluidPressure,
          arrayView1d< real64 const > const & biotCoefficient,
          timeIntegrationOption const tiOption,
          real64 const stiffnessDamping,
          real64 const massDamping,
          real64 const newmarkBeta,
          real64 const newmarkGamma,
          DofManager const * const GEOSX_UNUSED_ARG( dofManager ),
          ParallelMatrix * const matrix,
          ParallelVector * const rhs )
  {
    GEOSX_MARK_FUNCTION;
    constexpr int dim = 3;
    RAJA::ReduceMax< serialReduce, double > maxForce( 0 );

    typename CONSTITUTIVE_TYPE::KernelWrapper const & constitutive = constitutiveRelation->createKernelWrapper();

    arrayView2d<R2SymTensor const> const & stress = constitutiveRelation->getStress();

    RAJA::forall< serialPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, numElems ),
                                  GEOSX_LAMBDA ( localIndex const k )
    {
      Epetra_LongLongSerialDenseVector elementLocalDofIndex( dim * NUM_NODES_PER_ELEM );
      Epetra_SerialDenseVector         R                   ( dim * NUM_NODES_PER_ELEM );
      Epetra_SerialDenseMatrix         dRdU                ( dim * NUM_NODES_PER_ELEM,
                                                             dim * NUM_NODES_PER_ELEM );
      Epetra_SerialDenseVector         element_dof_np1     ( dim * NUM_NODES_PER_ELEM );

      Epetra_SerialDenseVector R_InertiaMassDamping(R);
      Epetra_SerialDenseMatrix dRdU_InertiaMassDamping(dRdU);
      Epetra_SerialDenseVector R_StiffnessDamping(R);
      Epetra_SerialDenseMatrix dRdU_StiffnessDamping(dRdU);

      R1Tensor u_local[NUM_NODES_PER_ELEM];
      R1Tensor uhat_local[NUM_NODES_PER_ELEM];
      R1Tensor vtilde_local[NUM_NODES_PER_ELEM];
      R1Tensor uhattilde_local[NUM_NODES_PER_ELEM];

      dRdU.Scale(0);
      R.Scale(0);

      dRdU_InertiaMassDamping.Scale(0);
      R_InertiaMassDamping.Scale(0);
      dRdU_StiffnessDamping.Scale(0);
      R_StiffnessDamping.Scale(0);

      real64 c[6][6];
      constitutive.GetStiffness( k, c );

      if(elemGhostRank[k] < 0)
      {
        for( localIndex a=0 ; a<NUM_NODES_PER_ELEM ; ++a)
        {

          localIndex localNodeIndex = elemsToNodes[k][a];

          for( int i=0 ; i<dim ; ++i )
          {
            elementLocalDofIndex[static_cast<int>(a)*dim+i] = globalDofNumber[localNodeIndex]+i;

            // TODO must add last solution estimate for this to be valid
            element_dof_np1(static_cast<int>(a)*dim+i) = disp[localNodeIndex][i];
          }
        }

        if( tiOption == timeIntegrationOption::ImplicitDynamic )
        {
          GEOSX_ERROR("Option not supported");
          CopyGlobalToLocal< NUM_NODES_PER_ELEM, R1Tensor>( elemsToNodes[k],
                                      disp, uhat, vtilde, uhattilde,
                                      u_local, uhat_local, vtilde_local, uhattilde_local );
        }
        else
        {
          CopyGlobalToLocal<NUM_NODES_PER_ELEM,R1Tensor>( elemsToNodes[k], disp, uhat, u_local, uhat_local );
        }




        R1Tensor dNdXa;
        R1Tensor dNdXb;


        for( integer q=0 ; q<NUM_QUADRATURE_POINTS ; ++q )
        {
          const realT detJq = detJ[k][q];
          std::vector<double> const & N = fe->values(q);

          for( integer a=0 ; a<NUM_NODES_PER_ELEM ; ++a )
          {
      //      realT const * const dNdXa = dNdX(q,a).Data();
            dNdXa = dNdX[k][q][a];

            for( integer b=0 ; b<NUM_NODES_PER_ELEM ; ++b )
            {
      //        realT const * const dNdXb = dNdX(q,b).Data();
              dNdXb = dNdX[k][q][b];

              dRdU(a*dim+0,b*dim+0) -= ( c[0][0]*dNdXa[0]*dNdXb[0] + c[5][5]*dNdXa[1]*dNdXb[1] + c[4][4]*dNdXa[2]*dNdXb[2] ) * detJq;
              dRdU(a*dim+0,b*dim+1) -= ( c[5][5]*dNdXa[1]*dNdXb[0] + c[0][1]*dNdXa[0]*dNdXb[1] ) * detJq;
              dRdU(a*dim+0,b*dim+2) -= ( c[4][4]*dNdXa[2]*dNdXb[0] + c[0][2]*dNdXa[0]*dNdXb[2] ) * detJq;

              dRdU(a*dim+1,b*dim+0) -= ( c[0][1]*dNdXa[1]*dNdXb[0] + c[5][5]*dNdXa[0]*dNdXb[1] ) * detJq;
              dRdU(a*dim+1,b*dim+1) -= ( c[5][5]*dNdXa[0]*dNdXb[0] + c[1][1]*dNdXa[1]*dNdXb[1] + c[3][3]*dNdXa[2]*dNdXb[2] ) * detJq;
              dRdU(a*dim+1,b*dim+2) -= ( c[3][3]*dNdXa[2]*dNdXb[1] + c[1][2]*dNdXa[1]*dNdXb[2] ) * detJq;

              dRdU(a*dim+2,b*dim+0) -= ( c[0][2]*dNdXa[2]*dNdXb[0] + c[4][4]*dNdXa[0]*dNdXb[2] ) * detJq;
              dRdU(a*dim+2,b*dim+1) -= ( c[1][2]*dNdXa[2]*dNdXb[1] + c[3][3]*dNdXa[1]*dNdXb[2] ) * detJq;
              dRdU(a*dim+2,b*dim+2) -= ( c[4][4]*dNdXa[0]*dNdXb[0] + c[3][3]*dNdXa[1]*dNdXb[1] + c[2][2]*dNdXa[2]*dNdXb[2] ) * detJq;

              if( tiOption == timeIntegrationOption::ImplicitDynamic )
              {

                real64 integrationFactor = density[k] * N[a] * N[b] * detJq;
                real64 temp1 = ( massDamping * newmarkGamma/( newmarkBeta * dt ) + 1.0 / ( newmarkBeta * dt * dt ) )* integrationFactor;

                for( int i=0 ; i<dim ; ++i )
                {
                  realT const acc = 1.0 / ( newmarkBeta * dt * dt ) * ( uhat_local[b][i] - uhattilde_local[b][i] );
                  realT const vel = vtilde_local[b][i] + newmarkGamma/( newmarkBeta * dt ) *( uhat_local[b][i] - uhattilde_local[b][i] );

                  dRdU_InertiaMassDamping(a*dim+i,b*dim+i) -= temp1;
                  R_InertiaMassDamping(a*dim+i) -= ( massDamping * vel + acc ) * integrationFactor;
                }
              }
            }
          }
        }



          R1Tensor temp;
          for( integer q=0 ; q<NUM_QUADRATURE_POINTS ; ++q )
          {
            R2SymTensor referenceStress = stress(k,q);
            if( !fluidPressure.empty() )
            {
              referenceStress.PlusIdentity( - biotCoefficient[0] * (fluidPressure[k] + deltaFluidPressure[k]));
            }

            const realT detJq = detJ[k][q];
            R2SymTensor stress0 = referenceStress;
            stress0 *= detJq;
            for( integer a=0 ; a<NUM_NODES_PER_ELEM ; ++a )
            {
              dNdXa = dNdX[k][q][a];

              temp.AijBj(stress0,dNdXa);
              realT maxF = temp.MaxVal();
              maxForce.max( maxF );

              R(a*dim+0) -= temp[0];
              R(a*dim+1) -= temp[1];
              R(a*dim+2) -= temp[2];
            }
          }


      // TODO It is simpler to do this...try it.
      //  dRdU.Multiply(dof_np1,R);
        for( integer a=0 ; a<NUM_NODES_PER_ELEM ; ++a )
        {
          realT nodeForce = 0;
          for( integer b=0 ; b<NUM_NODES_PER_ELEM ; ++b )
          {
            for( int i=0 ; i<dim ; ++i )
            {
              for( int j=0 ; j<dim ; ++j )
              {
                R(a*dim+i) += dRdU(a*dim+i,b*dim+j) * uhat_local[b][j];
              }
            }

            if( tiOption == timeIntegrationOption::ImplicitDynamic )
            {
              for( int i=0 ; i<dim ; ++i )
              {
                for( int j=0 ; j<dim ; ++j )
                {
                  R_StiffnessDamping(a*dim+i) += stiffnessDamping * dRdU(a*dim+i,b*dim+j) * ( vtilde_local[b][j] + newmarkGamma/(newmarkBeta * dt)*(uhat_local[b][j]-uhattilde_local[b][j]) );
                }
              }
            }

          }

          nodeForce = std::max( std::max( R(a*dim+0), R(a*dim+1) ),  R(a*dim+2) );
          maxForce.max( fabs( nodeForce ) );
        }


        if( tiOption == timeIntegrationOption::ImplicitDynamic )
        {
          dRdU_StiffnessDamping = dRdU;
          dRdU_StiffnessDamping.Scale( stiffnessDamping * newmarkGamma / ( newmarkBeta * dt ) );

          dRdU += dRdU_InertiaMassDamping;
          dRdU += dRdU_StiffnessDamping;
          R    += R_InertiaMassDamping;
          R    += R_StiffnessDamping;
        }

        // TODO remove local epetra objects, remove use of unwrappedPointer()
        matrix->unwrappedPointer()->SumIntoGlobalValues( elementLocalDofIndex, dRdU);
        rhs->unwrappedPointer()->SumIntoGlobalValues( elementLocalDofIndex, R);
      }
    });

    return maxForce.get();
  }
};

} // namespace SolidMechanicsLagrangianSSLEKernels

} // namespace geosx
