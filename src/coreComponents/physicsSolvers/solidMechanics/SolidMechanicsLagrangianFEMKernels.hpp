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
#include "finiteElement/FiniteElementShapeFunctionKernel.hpp"
#include "finiteElement/kernelInterface/RegionLoopSparsity.hpp"
#include "finiteElement/Kinematics.h"
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
    using CONSTITUTIVE_TYPE = TYPEOFPTR( constitutive );
    if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==8 )
    {
      rval = KERNELWRAPPER::template Launch< 8, 8, CONSTITUTIVE_TYPE >( constitutive,
                                                                        std::forward< PARAMS >( params )... );
    }
    else if( NUM_NODES_PER_ELEM==4 && NUM_QUADRATURE_POINTS==1 )
    {
      rval = KERNELWRAPPER::template Launch< 4, 1, CONSTITUTIVE_TYPE >( constitutive,
                                                                        std::forward< PARAMS >( params )... );
    }
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
                  real64 const (&dNdX)[8][3],
  #else
                  arraySlice1d< R1Tensor const > const & dNdX,
  #endif
                  real64 const detJ,
                  real64 const detF,
                  R2Tensor const & fInv,
                  R1Tensor (& result)[N] )
  {
    GEOSX_ASSERT_EQ( fieldVar.size(), 6 );

    real64 const integrationFactor = detJ * detF;

    real64 P[ 3 ][ 3 ];
    P[ 0 ][ 0 ] = ( fieldVar[ 0 ] * fInv( 0, 0 ) + fieldVar[ 5 ] * fInv( 0, 1 ) + fieldVar[ 4 ] * fInv( 0, 2 ) ) * integrationFactor;
    P[ 0 ][ 1 ] = ( fieldVar[ 0 ] * fInv( 1, 0 ) + fieldVar[ 5 ] * fInv( 1, 1 ) + fieldVar[ 4 ] * fInv( 1, 2 ) ) * integrationFactor;
    P[ 0 ][ 2 ] = ( fieldVar[ 0 ] * fInv( 2, 0 ) + fieldVar[ 5 ] * fInv( 2, 1 ) + fieldVar[ 4 ] * fInv( 2, 2 ) ) * integrationFactor;

    P[ 1 ][ 0 ] = ( fieldVar[ 5 ] * fInv( 0, 0 ) + fieldVar[ 1 ] * fInv( 0, 1 ) + fieldVar[ 3 ] * fInv( 0, 2 ) ) * integrationFactor;
    P[ 1 ][ 1 ] = ( fieldVar[ 5 ] * fInv( 1, 0 ) + fieldVar[ 1 ] * fInv( 1, 1 ) + fieldVar[ 3 ] * fInv( 1, 2 ) ) * integrationFactor;
    P[ 1 ][ 2 ] = ( fieldVar[ 5 ] * fInv( 2, 0 ) + fieldVar[ 1 ] * fInv( 2, 1 ) + fieldVar[ 3 ] * fInv( 2, 2 ) ) * integrationFactor;

    P[ 2 ][ 0 ] = ( fieldVar[ 4 ] * fInv( 0, 0 ) + fieldVar[ 3 ] * fInv( 0, 1 ) + fieldVar[ 2 ] * fInv( 0, 2 ) ) * integrationFactor;
    P[ 2 ][ 1 ] = ( fieldVar[ 4 ] * fInv( 1, 0 ) + fieldVar[ 3 ] * fInv( 1, 1 ) + fieldVar[ 2 ] * fInv( 1, 2 ) ) * integrationFactor;
    P[ 2 ][ 2 ] = ( fieldVar[ 4 ] * fInv( 2, 0 ) + fieldVar[ 3 ] * fInv( 2, 1 ) + fieldVar[ 2 ] * fInv( 2, 2 ) ) * integrationFactor;

    for( int a=0; a<N; ++a )    // loop through all shape functions in element
    {
      result[a][0] -= P[ 0 ][ 0 ] * dNdX[ a ][ 0 ] + P[ 0 ][ 1 ] * dNdX[ a ][ 1 ] + P[ 0 ][ 2 ] * dNdX[ a ][ 2 ];
      result[a][1] -= P[ 1 ][ 0 ] * dNdX[ a ][ 0 ] + P[ 1 ][ 1 ] * dNdX[ a ][ 1 ] + P[ 1 ][ 2 ] * dNdX[ a ][ 2 ];
      result[a][2] -= P[ 2 ][ 0 ] * dNdX[ a ][ 0 ] + P[ 2 ][ 1 ] * dNdX[ a ][ 1 ] + P[ 2 ][ 2 ] * dNdX[ a ][ 2 ];
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
          LvArray::SortedArrayView< localIndex const, localIndex > const & elementList,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView3d< R1Tensor const > const & dNdX,
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



    typename CONSTITUTIVE_TYPE::KernelWrapper constitutive = constitutiveRelation->createKernelWrapper();

    using KERNEL_POLICY = parallelDevicePolicy< 32 >;
    RAJA::forall< KERNEL_POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, elementList.size() ),
                                   [=] GEOSX_DEVICE ( localIndex const index )
    {
      localIndex const k = elementList[ index ];

      R1Tensor v_local[NUM_NODES_PER_ELEM];
      R1Tensor u_local[NUM_NODES_PER_ELEM];
      R1Tensor f_local[NUM_NODES_PER_ELEM];
#if defined(CALCFEMSHAPE)
      real64 X_local[8][3];
#endif
      for( localIndex a=0; a< NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        for( int i=0; i<3; ++i )
        {
#if defined(CALCFEMSHAPE)
          X_local[ a ][ i ] = X[ nodeIndex ][ i ];
#endif
          u_local[ a ][ i ] = u[ nodeIndex ][ i ];
          v_local[ a ][ i ] = vel[ nodeIndex ][ i ];
        }
      }

      //Compute Quadrature
      for( localIndex q = 0; q<NUM_QUADRATURE_POINTS; ++q )
      {
#if defined(CALCFEMSHAPE)
        real64 dNdX[ 8 ][ 3 ];
        real64 const detJ = FiniteElementShapeKernel::shapeFunctionDerivatives( q, X_local, dNdX );
#define DNDX dNdX
#define DETJ detJ
#else
#define DNDX dNdX[k][q]
#define DETJ detJ( k, q )
#endif
        R2Tensor dUhatdX, dUdX;
        CalculateGradients< NUM_NODES_PER_ELEM >( dUhatdX, dUdX, v_local, u_local, DNDX );
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
        real64 const detF = F.Det();
        fInv.Inverse( F );


        R2Tensor Rot;
        R2SymTensor Dadt;
        HughesWinget( Rot, Dadt, Ldt );

        constitutive.HypoElastic( k, q, Dadt.Data(), Rot );

        Integrate< NUM_NODES_PER_ELEM >( constitutive.m_stress[k][q].toSliceConst(),
                                         DNDX,
                                         DETJ,
                                         detF,
                                         fInv,
                                         f_local );
// );
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
          TimeIntegrationOption const GEOSX_UNUSED_PARAM( tiOption ),
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

class QuasiStatic : public finiteElement::RegionLoopSparsity
{
public:
  using Base = finiteElement::RegionLoopSparsity;
  static constexpr int numTestDofPerSP = 3;
  static constexpr int numTrialDofPerSP = 3;

  struct Parameters : public Base::Parameters
  {
    Parameters( real64 const inputGravityVector[3] ):
      m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] }
    {}

    real64 const m_gravityVector[3];
  };


  template< int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
            int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >
  struct StackVariables : Base::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM*numTestDofPerSP,
                                                NUM_TRIAL_SUPPORT_POINTS_PER_ELEM*numTrialDofPerSP >
  {
public:
    using StackVariablesBase = Base::StackVariables< NUM_TEST_SUPPORT_POINTS_PER_ELEM*numTestDofPerSP,
                                                     NUM_TRIAL_SUPPORT_POINTS_PER_ELEM*numTrialDofPerSP >;
    using StackVariablesBase::numRows;
    using StackVariablesBase::numCols;
    static constexpr int numNodes = NUM_TEST_SUPPORT_POINTS_PER_ELEM;


    GEOSX_HOST_DEVICE
    StackVariables():
      StackVariablesBase(),
      u_local(),
      uhat_local(),
      constitutiveStiffness{ {0.0} }
    {}

    R1Tensor u_local[numNodes];
    R1Tensor uhat_local[numNodes];
    real64 constitutiveStiffness[6][6];
  };

  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  using SparsityKernels = Base::Kernels< SUBREGION_TYPE,
                                         CONSTITUTIVE_TYPE,
                                         NUM_NODES_PER_ELEM,
                                         NUM_NODES_PER_ELEM,
                                         numTestDofPerSP,
                                         numTrialDofPerSP >;

  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  class Kernels : public Base::Kernels< SUBREGION_TYPE,
                                        CONSTITUTIVE_TYPE,
                                        NUM_NODES_PER_ELEM,
                                        NUM_NODES_PER_ELEM,
                                        numTestDofPerSP,
                                        numTrialDofPerSP >
  {
public:
    using KernelBase = Base::Kernels< SUBREGION_TYPE,
                                      CONSTITUTIVE_TYPE,
                                      NUM_NODES_PER_ELEM,
                                      NUM_NODES_PER_ELEM,
                                      numTestDofPerSP,
                                      numTrialDofPerSP >;

    static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;

    using KernelBase::m_dofNumber;
    using KernelBase::m_matrix;
    using KernelBase::m_rhs;
    using KernelBase::elemsToNodes;
    using KernelBase::constitutiveUpdate;
    using KernelBase::m_finiteElementSpace;

    using StackVars = StackVariables< numNodesPerElem,
                                      numNodesPerElem >;



    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
             ParallelMatrix & inputMatrix,
             ParallelVector & inputRhs,
             NodeManager const & nodeManager,
             SUBREGION_TYPE const & elementSubRegion,
             FiniteElementBase const * const finiteElementSpace,
             CONSTITUTIVE_TYPE * const inputConstitutiveType,
             Parameters const & parameters ):
      KernelBase( inputDofNumber,
                  inputMatrix,
                  inputRhs,
                  nodeManager,
                  elementSubRegion,
                  finiteElementSpace,
                  inputConstitutiveType,
                  parameters ),
      m_disp( nodeManager.totalDisplacement()),
      m_uhat( nodeManager.incrementalDisplacement()),
      dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
      detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) )//,
    {}

    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;
    arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

    arrayView3d< R1Tensor const > const dNdX;
    arrayView2d< real64 const > const detJ;


    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void preKernel( localIndex const k,
                    STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const localNodeIndex = elemsToNodes( k, a );

        stack.u_local[ a ] = m_disp[ localNodeIndex ];
        stack.uhat_local[ a ] = m_uhat[ localNodeIndex ];

        for( int i=0; i<3; ++i )
        {
          stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
          stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        }
      }

    }

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void updateKernel( localIndex const k,
                       localIndex const q,
                       STACK_VARIABLE_TYPE & stack ) const
    {
      real64 strainInc[6] = {0};
      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        strainInc[0] = strainInc[0] + dNdX( k, q, a )[0] * stack.uhat_local[a][0];
        strainInc[1] = strainInc[1] + dNdX( k, q, a )[1] * stack.uhat_local[a][1];
        strainInc[2] = strainInc[2] + dNdX( k, q, a )[2] * stack.uhat_local[a][2];
        strainInc[3] = strainInc[3] + dNdX( k, q, a )[2] * stack.uhat_local[a][1] +
                       dNdX( k, q, a )[1] * stack.uhat_local[a][2];
        strainInc[4] = strainInc[4] + dNdX( k, q, a )[2] * stack.uhat_local[a][0] +
                       dNdX( k, q, a )[0] * stack.uhat_local[a][2];
        strainInc[5] = strainInc[5] + dNdX( k, q, a )[1] * stack.uhat_local[a][0] +
                       dNdX( k, q, a )[0] * stack.uhat_local[a][1];
      }

      constitutiveUpdate.SmallStrain( k, q, strainInc );

      GEOSX_UNUSED_VAR( q )
      constitutiveUpdate.GetStiffness( k, stack.constitutiveStiffness );
    }


    template< typename PARAMETERS_TYPE,
              typename STACK_VARIABLE_TYPE/*,
              typename DYNAMICS_LAMBDA = nvstd::function< void( localIndex, localIndex) >*/ >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void stiffnessKernel( localIndex const k,
                          localIndex const q,
                          PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                          STACK_VARIABLE_TYPE & stack/*,
                          DYNAMICS_LAMBDA && dynamicsTerms = [] GEOSX_DEVICE ( localIndex, localIndex){}*/ ) const
    {
      for( localIndex a=0; a<NUM_NODES_PER_ELEM; ++a )
      {
        for( localIndex b=0; b<NUM_NODES_PER_ELEM; ++b )
        {
          real64 const (&c)[6][6] = stack.constitutiveStiffness;
          stack.localJacobian[ a*3+0 ][ b*3+0 ] -= ( c[0][0]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                     c[5][5]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                     c[4][4]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+0 ][ b*3+1 ] -= ( c[5][5]*dNdX( k, q, a )[1]*dNdX( k, q, b )[0] +
                                                     c[0][1]*dNdX( k, q, a )[0]*dNdX( k, q, b )[1] ) * detJ( k, q );

          stack.localJacobian[ a*3+0 ][ b*3+2 ] -= ( c[4][4]*dNdX( k, q, a )[2]*dNdX( k, q, b )[0] +
                                                     c[0][2]*dNdX( k, q, a )[0]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+1 ][ b*3+1 ] -= ( c[5][5]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                     c[1][1]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                     c[3][3]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+1 ][ b*3+0 ] -= ( c[0][1]*dNdX( k, q, a )[1]*dNdX( k, q, b )[0] +
                                                     c[5][5]*dNdX( k, q, a )[0]*dNdX( k, q, b )[1] ) * detJ( k, q );

          stack.localJacobian[ a*3+1 ][ b*3+2 ] -= ( c[3][3]*dNdX( k, q, a )[2]*dNdX( k, q, b )[1] +
                                                     c[1][2]*dNdX( k, q, a )[1]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+2 ][ b*3+0 ] -= ( c[0][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[0] +
                                                     c[4][4]*dNdX( k, q, a )[0]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+2 ][ b*3+1 ] -= ( c[1][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[1] +
                                                     c[3][3]*dNdX( k, q, a )[1]*dNdX( k, q, b )[2] ) * detJ( k, q );

          stack.localJacobian[ a*3+2 ][ b*3+2 ] -= ( c[4][4]*dNdX( k, q, a )[0]*dNdX( k, q, b )[0] +
                                                     c[3][3]*dNdX( k, q, a )[1]*dNdX( k, q, b )[1] +
                                                     c[2][2]*dNdX( k, q, a )[2]*dNdX( k, q, b )[2] ) * detJ( k, q );

//          dynamicsTerms( a, b );
        }
      }
    }

    template< typename PARAMETERS_TYPE,
              typename STACK_VARIABLE_TYPE/*,
              typename DYNAMICS_LAMBDA = STD_FUNCTION< void( real64 * ) >*/ >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void integrationKernel( localIndex const k,
                            localIndex const q,
                            PARAMETERS_TYPE const & parameters,
                            STACK_VARIABLE_TYPE & stack/*,
                            DYNAMICS_LAMBDA && stressModifier = [] GEOSX_DEVICE ( real64 * ) {}*/ ) const
    {
      real64 stress[6] = { constitutiveUpdate.m_stress( k, q, 0 ),
                           constitutiveUpdate.m_stress( k, q, 1 ),
                           constitutiveUpdate.m_stress( k, q, 2 ),
                           constitutiveUpdate.m_stress( k, q, 3 ),
                           constitutiveUpdate.m_stress( k, q, 4 ),
                           constitutiveUpdate.m_stress( k, q, 5 ) };

//      stressModifier( stress );

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        stack.localResidual[ a * 3 + 0 ] -= ( stress[ 0 ] * dNdX( k, q, a )[ 0 ] +
                                              stress[ 5 ] * dNdX( k, q, a )[ 1 ] +
                                              stress[ 4 ] * dNdX( k, q, a )[ 2 ] -
                                              parameters.m_gravityVector[0] ) * detJ( k, q );
        stack.localResidual[ a * 3 + 1 ] -= ( stress[ 5 ] * dNdX( k, q, a )[ 0 ] +
                                              stress[ 1 ] * dNdX( k, q, a )[ 1 ] +
                                              stress[ 3 ] * dNdX( k, q, a )[ 2 ] -
                                              parameters.m_gravityVector[1] ) * detJ( k, q );
        stack.localResidual[ a * 3 + 2 ] -= ( stress[ 4 ] * dNdX( k, q, a )[ 0 ] +
                                              stress[ 3 ] * dNdX( k, q, a )[ 1 ] +
                                              stress[ 2 ] * dNdX( k, q, a )[ 2 ] -
                                              parameters.m_gravityVector[2] ) * detJ( k, q );
      }
    }

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    //GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 postKernel( localIndex const GEOSX_UNUSED_PARAM( k ),
                       PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                       STACK_VARIABLE_TYPE & stack ) const
    {
      real64 meanForce = 0;
      for( localIndex a=0; a<stack.numRows; ++a )
      {
//        RAJA::atomicMax< RAJA::auto_atomic >( &meanForce, stack.localResidual[a] );
        meanForce = std::max( meanForce, stack.localResidual[a] );
//                meanForce += fabs( stack.localResidual[a] );
      }
//            meanForce /= stack.ndof;

      m_matrix.add( stack.localRowDofIndex,
                    stack.localColDofIndex,
                    &(stack.localJacobian[0][0]),
                    stack.numRows,
                    stack.numCols );

      m_rhs.add( stack.localRowDofIndex,
                 stack.localResidual,
                 stack.numRows );

      return meanForce;
    }

  };
};

//
//
//class ImplicitNewmark
//{
//public:
//  using Base = QuasiStatic;
//
//  struct Parameters : public Base::Parameters
//  {
//    Parameters( real64 const inputGravityVector[3],
//                real64 const inputNewmarkGamma,
//                real64 const inputNewmarkBeta,
//                real64 const inputMassDamping,
//                real64 const inputStiffnessDamping,
//                real64 const inputDt ):
//    Base::Parameters( inputGravityVector ),
//      newmarkGamma(inputNewmarkGamma),
//      newmarkBeta(inputNewmarkBeta),
//      massDamping(inputMassDamping),
//      stiffnessDamping(inputStiffnessDamping),
//      dt(inputDt)
//    {}
//
//    real64 const newmarkGamma;
//    real64 const newmarkBeta;
//    real64 const massDamping;
//    real64 const stiffnessDamping;
//    real64 const dt;
//  };
//
//  template< int NUM_NODES_PER_ELEM, int NUM_DOF_PER_NODE >
//  struct StackVariables : Base::StackVariables< NUM_NODES_PER_ELEM, NUM_DOF_PER_NODE >
//  {
//public:
//    using StackVariablesBase = Base::StackVariables< NUM_NODES_PER_ELEM, NUM_DOF_PER_NODE >;
//
//    using StackVariablesBase::numNodesPerElem;
//    using StackVariablesBase::numDofPerNode;
//    using StackVariablesBase::ndof;
//
////      GEOSX_HOST_DEVICE
//    StackVariables():
//      StackVariablesBase(),
//      dRdU_InertiaMassDamping{ {0.0} },
//      vtilde_local{ { 0.0, 0.0, 0.0} },
//      uhattilde_local{ { 0.0, 0.0, 0.0} }
//    {}
//
//    real64 dRdU_InertiaMassDamping[ ndof ][ ndof ];
//
//    R1Tensor vtilde_local[NUM_NODES_PER_ELEM];
//    R1Tensor uhattilde_local[NUM_NODES_PER_ELEM];
//  };
//
//  template< typename SUBREGION_TYPE,
//            typename CONSTITUTIVE_TYPE >
//  class Kernels : public Base::Kernels< SUBREGION_TYPE, CONSTITUTIVE_TYPE >
//  {
//  public:
//
//    using KernelBase = Base::Kernels<SUBREGION_TYPE,CONSTITUTIVE_TYPE>;
//    using KernelBase::m_dofNumber;
//    using KernelBase::m_matrix;
//    using KernelBase::m_rhs;
//    using KernelBase::elemsToNodes;
//    using KernelBase::constitutiveUpdate;
//    using KernelBase::m_disp;
//    using KernelBase::m_uhat;
//    using KernelBase::dNdX;
//    using KernelBase::detJ;
//    using KernelBase::m_finiteElementSpace;
//
//
//
//    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
//             ParallelMatrix & inputMatrix,
//             ParallelVector & inputRhs,
//             NodeManager const & nodeManager,
//             SUBREGION_TYPE const & elementSubRegion,
//             FiniteElementBase const * const finiteElementSpace,
//             CONSTITUTIVE_TYPE & constitutiveModel,
//             Base::Parameters const & parameters ):
//      KernelBase( inputDofNumber,
//                 inputMatrix,
//                 inputRhs,
//                 nodeManager,
//                 elementSubRegion,
//                 finiteElementSpace,
//                 constitutiveModel,
//                 parameters ),
//      m_vtilde(nodeManager.totalDisplacement()),
//      m_uhattilde(nodeManager.totalDisplacement()),
//      m_density(constitutiveModel.getDensity())
//    {}
//
//    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_vtilde;
//    arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhattilde;
//    arrayView2d< real64 const > const m_density;
//
//    template< typename STACK_VARIABLE_TYPE >
//  //    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    void preKernel( localIndex const k,
//                    STACK_VARIABLE_TYPE & stack ) const
//    {
//      for( localIndex a=0; a<STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
//      {
//        localIndex const localNodeIndex = elemsToNodes( k, a );
//
//        stack.u_local[ a ] = m_disp[ localNodeIndex ];
//        stack.uhat_local[ a ] = m_uhat[ localNodeIndex ];
//        stack.vtilde_local[ a ] = m_vtilde[ localNodeIndex ];
//        stack.uhattilde_local[ a ] = m_uhattilde[ localNodeIndex ];
//
//        for( int i=0; i<3; ++i )
//        {
//          stack.elementLocalDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
//        }
//      }
//
//    }
//
//    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
//  //    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    void stiffnessKernel( localIndex const k,
//                          localIndex const q,
//                          PARAMETERS_TYPE const & parameters,
//                          STACK_VARIABLE_TYPE & stack ) const
//    {
//
//      std::vector< double > const & N = m_finiteElementSpace->values( q );
//
////      real64 N[STACK_VARIABLE_TYPE::numNodesPerElem];
//      real64 const & massDamping = parameters.massDamping;
//      real64 const & newmarkGamma = parameters.newmarkGamma;
//      real64 const & newmarkBeta = parameters.newmarkBeta;
//      real64 const & dt = parameters.dt;
//
//      KernelBase::stiffnessKernel( k, q, parameters, stack, [&]( localIndex const a, localIndex const b )
//      {
//        real64 integrationFactor = m_density( k, q ) * N[a] * N[b] * detJ(k,q);
//        real64 temp1 = ( massDamping * newmarkGamma/( newmarkBeta * dt ) + 1.0 / ( newmarkBeta * dt * dt ) )*
// integrationFactor;
//
//        constexpr int nsdof = STACK_VARIABLE_TYPE::numDofPerNode;
//        for( int i=0; i<nsdof; ++i )
//        {
//          realT const acc = 1.0 / ( newmarkBeta * dt * dt ) * ( stack.uhat_local[b][i] - stack.uhattilde_local[b][i]
// );
//          realT const vel = stack.vtilde_local[b][i] + newmarkGamma/( newmarkBeta * dt ) *( stack.uhat_local[b][i] -
// stack.uhattilde_local[b][i] );
//
//          stack.dRdU_InertiaMassDamping[ a*nsdof+i][ b*nsdof+i ] -= temp1;
//          stack.localResidual[ a*nsdof+i ] -= ( massDamping * vel + acc ) * integrationFactor;
//        }
//      } );
//    }
//
//    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
//  //    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    real64 postKernel( PARAMETERS_TYPE const & parameters,
//                       STACK_VARIABLE_TYPE & stack ) const
//    {
//      constexpr int nsdof = STACK_VARIABLE_TYPE::numDofPerNode;
//      real64 const & stiffnessDamping = parameters.stiffnessDamping;
//      real64 const & newmarkGamma     = parameters.newmarkGamma;
//      real64 const & newmarkBeta      = parameters.newmarkBeta;
//      real64 const & dt               = parameters.dt;
//
//      for( localIndex a=0; a<STACK_VARIABLE_TYPE::numNodesPerElem; ++a )
//      {
//        for( localIndex b=0; b<STACK_VARIABLE_TYPE::numNodesPerElem; ++b )
//        {
//          for( int i=0; i<nsdof; ++i )
//          {
//            for( int j=0; j<nsdof; ++j )
//            {
//              stack.localResidual[ a*nsdof+i ] += stiffnessDamping * stack.localJacobian[ a*nsdof+i][ b*nsdof+j ] *
//                                               ( stack.vtilde_local[b][j] + newmarkGamma/(newmarkBeta *
// dt)*(stack.uhat_local[b][j]-stack.uhattilde_local[b][j]) );
//
//              stack.localJacobian[a*nsdof+i][b*nsdof+j] += stack.localJacobian[a][b] * (1.0 + stiffnessDamping *
// newmarkGamma / ( newmarkBeta * dt ) ) +
//                                           stack.dRdU_InertiaMassDamping[ a ][ b ] ;
//            }
//          }
//        }
//      }
//
//      for( localIndex a=0 ; a<STACK_VARIABLE_TYPE::ndof ; ++a )
//      {
//        for( localIndex b=0 ; b<STACK_VARIABLE_TYPE::ndof ; ++b )
//        {
//          stack.localJacobian[a][b] += stack.localJacobian[a][b] * (1.0 + stiffnessDamping * newmarkGamma / (
// newmarkBeta * dt ) ) +
//                                       stack.dRdU_InertiaMassDamping[ a ][ b ] ;
//        }
//      }
//
//      return KernelBase::postKernel( parameters, stack );
//    }
//  };
//
//};




//class SmallStrainFracturePenaltyContact
//{
//public:
//  using Base = physicsLoopInterface::RegionLoop;
//  static constexpr int maxNumNodesPerFace = 4;
//  static constexpr int numTestDofPerSP = 3;
//  static constexpr int numTrialDofPerSP = 3;
//
//
//  template< int ,
//            int  >
//  struct StackVariables : Base::StackVariables< maxNumNodesPerFace*numTestDofPerSP*2,
//                                                maxNumNodesPerFace*numTrialDofPerSP*2 >
//  {
//public:
//    using StackVariablesBase = Base::StackVariables< maxNumNodesPerFace*numTestDofPerSP*2,
//                                                     maxNumNodesPerFace*numTrialDofPerSP*2 >;
//    using StackVariablesBase::numRows;
//    using StackVariablesBase::numCols;
//
////      GEOSX_HOST_DEVICE
//    StackVariables():
//      StackVariablesBase()
//    {}
//  };
//
//  template< typename SUBREGION_TYPE,
//            typename CONSTITUTIVE_TYPE,
//            int NUM_NODES_PER_ELEM,
//            int >
//  using SparsityKernels = Base::Kernels< SUBREGION_TYPE,
//                                         CONSTITUTIVE_TYPE,
//                                         NUM_NODES_PER_ELEM,
//                                         NUM_NODES_PER_ELEM,
//                                         numTestDofPerSP,
//                                         numTrialDofPerSP >
//  ;
//
//  template< typename SUBREGION_TYPE,
//            typename CONSTITUTIVE_TYPE,
//            int NUM_NODES_PER_ELEM,
//            int >
//  class Kernels : public Base::Kernels< SUBREGION_TYPE,
//                                        CONSTITUTIVE_TYPE,
//                                        NUM_NODES_PER_ELEM,
//                                        NUM_NODES_PER_ELEM,
//                                        numTestDofPerSP,
//                                        numTrialDofPerSP >
//  {
//public:
//    using KernelBase = Base::Kernels< SUBREGION_TYPE,
//                                      CONSTITUTIVE_TYPE,
//                                      NUM_NODES_PER_ELEM,
//                                      NUM_NODES_PER_ELEM,
//                                      numTestDofPerSP,
//                                      numTrialDofPerSP >;
//
//    static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;
//
//    using KernelBase::m_dofNumber;
//    using KernelBase::m_matrix;
//    using KernelBase::m_rhs;
//    using KernelBase::elemsToNodes;
//    using KernelBase::constitutiveUpdate;
//    using KernelBase::m_finiteElementSpace;
//
//    using StackVars = StackVariables< numNodesPerElem,
//                                      numNodesPerElem >;
//
//
//
//    Kernels( arrayView1d< globalIndex const > const & inputDofNumber,
//             ParallelMatrix & inputMatrix,
//             ParallelVector & inputRhs,
//             NodeManager const & nodeManager,
//             SUBREGION_TYPE const & elementSubRegion,
//             FiniteElementBase const * const finiteElementSpace,
//             CONSTITUTIVE_TYPE & inputConstitutiveType,
//             Parameters const & GEOSX_UNUSED_PARAM( parameters ) )://,
//      KernelBase( inputDofNumber,
//                  inputMatrix,
//                  inputRhs,
//                  nodeManager,
//                  elementSubRegion,
//                  finiteElementSpace,
//                  inputConstitutiveType,
//                  inputConstitutiveType.createKernelWrapper() ),
//      m_faceNormal( faceManager->faceNormal() ),
//      m_facesToNodes( faceManager->nodeList() ),
//      area( subRegion.getElementArea() ),
//      elemsToFaces( subRegion.faceList() ),
//
//      dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
//      detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) )//,
//    {}
//
//    arrayView1d< R1Tensor const > const m_faceNormal;
//    ArrayOfArraysView< localIndex const > const m_facesToNodes;
//    arrayView1d< real64 > const area;
//    arrayView2d< localIndex const > const elemsToFaces;
//
//
//
//    template< typename STACK_VARIABLE_TYPE >
//    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    void preKernel( localIndex const k,
//                    STACK_VARIABLE_TYPE & stack ) const
//    {
//      R1Tensor Nbar = m_faceNormal[elemsToFaces[kfe][0]];
//      Nbar -= m_faceNormal[elemsToFaces[kfe][1]];
//      Nbar.Normalize();
//
//      localIndex const kf0 = elemsToFaces[kfe][0];
//      localIndex const kf1 = elemsToFaces[kfe][1];
//      localIndex const numNodesPerFace=m_facesToNodes.sizeOfArray( kf0 );
//      real64 const Ja = area[kfe] / numNodesPerFace;
//
//      stackArray1d< globalIndex, maxDofPerElem > rowDOF( numNodesPerFace*3*2 );
//      stackArray1d< real64, maxDofPerElem > nodeRHS( numNodesPerFace*3*2 );
//      stackArray2d< real64, maxDofPerElem *maxDofPerElem > dRdP( numNodesPerFace*3*2, numNodesPerFace*3*2 );
//
//      for( localIndex a=0; a<numNodesPerFace; ++a )
//      {
//        localIndex const node0 = facesToNodes[kf0][a];
//        localIndex const node1 = facesToNodes[kf1][ a==0 ? a : numNodesPerFace-a ];
//
//        for( int i=0; i<3; ++i )
//        {
//          rowDOF[3*a+i]                     = nodeDofNumber[node0]+i;
//          rowDOF[3*(numNodesPerFace + a)+i] = nodeDofNumber[node1]+i;
//        }
//
//        R1Tensor gap = u[node1];
//        gap -= u[node0];
//        real64 const gapNormal = Dot( gap, Nbar );
//
//
//        if( gapNormal < 0 )
//        {
//          R1Tensor penaltyForce = Nbar;
//          penaltyForce *= -contactStiffness * gapNormal * Ja;
//          for( int i=0; i<3; ++i )
//          {
//            fc[node0] -= penaltyForce;
//            fc[node1] += penaltyForce;
//            nodeRHS[3*a+i]                     -= penaltyForce[i];
//            nodeRHS[3*(numNodesPerFace + a)+i] += penaltyForce[i];
//
//            dRdP( 3*a+i, 3*a+i )                                         -= contactStiffness * Ja * Nbar[i] * Nbar[i];
//            dRdP( 3*a+i, 3*(numNodesPerFace + a)+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
//            dRdP( 3*(numNodesPerFace + a)+i, 3*a+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
//            dRdP( 3*(numNodesPerFace + a)+i, 3*(numNodesPerFace + a)+i ) -= contactStiffness * Ja * Nbar[i] * Nbar[i];
//          }
//        }
//      }
//
//      m_residual.add( rowDOF, nodeRHS );
//      m_matrix.add( rowDOF, rowDOF, dRdP );
//
//
//    }
//
//    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
//    GEOSX_HOST_DEVICE
//    GEOSX_FORCE_INLINE
//    real64 postKernel( PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
//                       STACK_VARIABLE_TYPE & stack ) const
//    {
//      real64 meanForce = 0;
//      for( localIndex a=0; a<stack.numRows; ++a )
//      {
//        meanForce = std::max( meanForce, stack.localResidual[a] );
//        //        meanForce += fabs( stack.localResidual[a] );
//      }
//      //      meanForce /= stack.ndof;
//
//
//      m_matrix.add( stack.localRowDofIndex,
//                    stack.localColDofIndex,
//                    &(stack.localJacobian[0][0]),
//                    stack.numRows,
//                    stack.numCols );
//
//      m_rhs.add( stack.localRowDofIndex,
//                 stack.localResidual,
//                 stack.numRows );
//
//      return meanForce;
//    }
//  };
//};


class ExplicitSmallStrain
{
//#if defined(GEOSX_USE_CUDA)
  #define CALCFEMSHAPE
//#endif
  // If UPDATE_STRESS is undef, then stress is not updated at all.
//  #define UPDATE_STRESS 1 // uses total displacement to and adds material stress state to integral for nodalforces.
#define UPDATE_STRESS 2 // uses velocity*dt and updates material stress state.

public:
  using Base = finiteElement::RegionLoop;
  static constexpr int numTestDofPerSP = 3;
  static constexpr int numTrialDofPerSP = 3;

  struct Parameters : public Base::Parameters
  {

    Parameters( real64 const dt,
                real64 const inputGravityVector[3] ):
      m_dt( dt ),
      m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] }
    {}

    real64 const m_dt;
    real64 const m_gravityVector[3];
  };


  template< int NUM_TEST_SUPPORT_POINTS_PER_ELEM,
            int NUM_TRIAL_SUPPORT_POINTS_PER_ELEM >
  struct StackVariables
  {
public:

    static constexpr int numNodes = NUM_TEST_SUPPORT_POINTS_PER_ELEM;


    GEOSX_HOST_DEVICE
    StackVariables():
      fLocal{ { 0.0} },
      varLocal{ {0.0} }
#if defined(CALCFEMSHAPE)
      ,
      xLocal(),
      dNdX(),
      detJ()
#endif
    {}

    real64 fLocal[ numNodes ][ 3 ];
    real64 varLocal[ numNodes ][ 3 ];
#if defined(CALCFEMSHAPE)
    real64 xLocal[ 8 ][ 3 ];
    real64 dNdX[ 8 ][ 3 ];
    real64 detJ;
#endif
  };

  template< typename SUBREGION_TYPE,
            typename CONSTITUTIVE_TYPE,
            int NUM_NODES_PER_ELEM,
            int >
  class Kernels
  {
public:

    static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;


    using StackVars = StackVariables< numNodesPerElem,
                                      numNodesPerElem >;



    Kernels( arrayView1d< globalIndex const > const & ,
             ParallelMatrix & ,
             ParallelVector & ,
             NodeManager & nodeManager,
             SUBREGION_TYPE const & elementSubRegion,
             FiniteElementBase const * const finiteElementSpace,
             CONSTITUTIVE_TYPE * const inputConstitutiveType,
             Parameters const & parameters )://,
      elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
      elemGhostRank( elementSubRegion.ghostRank() ),
      constitutiveUpdate( inputConstitutiveType->createKernelWrapper() ),
      m_finiteElementSpace( finiteElementSpace ),
#if !defined(CALCFEMSHAPE)
      dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
      detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) ),
#endif
      X(nodeManager.referencePosition()),
      u(nodeManager.totalDisplacement()),
      vel(nodeManager.velocity()),
      acc(nodeManager.acceleration()),
      m_dt(parameters.m_dt)
    {}

    typename SUBREGION_TYPE::NodeMapType::base_type::ViewTypeConst const elemsToNodes;
    arrayView1d< integer const > const elemGhostRank;
    typename CONSTITUTIVE_TYPE::KernelWrapper const constitutiveUpdate;
    FiniteElementBase const * m_finiteElementSpace;
#if !defined(CALCFEMSHAPE)
    arrayView3d< R1Tensor const > const dNdX;
    arrayView2d< real64 const > const detJ;
#endif
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X;
    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const u;
    arrayView2d< real64 const, nodes::VELOCITY_USD > const vel;
    arrayView2d< real64, nodes::ACCELERATION_USD > const acc;
    real64 const m_dt;


    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void preKernel( localIndex const k,
                    STACK_VARIABLE_TYPE & stack ) const
    {
      for( localIndex a=0; a< NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        for( int i=0; i<3; ++i )
        {
#if defined(CALCFEMSHAPE)
          stack.xLocal[ a ][ i ] = X[ nodeIndex ][ i ];
#endif

#if UPDATE_STRESS==2
          stack.varLocal[ a ][ i ] = vel[ nodeIndex ][ i ] * m_dt;
#else
          stack.varLocal[ a ][ i ] = u[ nodeIndex ][ i ];
#endif
        }
      }
    }

    template< typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void updateKernel( localIndex const k,
                       localIndex const q,
                       STACK_VARIABLE_TYPE & stack ) const
    {

#if defined(CALCFEMSHAPE)
        real64 dNdX[ 8 ][ 3 ];
        real64 const detJ = FiniteElementShapeKernel::shapeFunctionDerivatives( q, stack.xLocal, dNdX );
  #define DNDX dNdX
  #define DETJ detJ
#else //defined(CALCFEMSHAPE)
  #define DNDX dNdX[k][q]
  #define DETJ detJ( k, q )
#endif //defined(CALCFEMSHAPE)

        real64 stressLocal[ 6 ] = {0};
        real64 strain[6] = {0};
        for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
        {
          strain[0] = strain[0] + DNDX[ a ][0] * stack.varLocal[ a ][0];
          strain[1] = strain[1] + DNDX[ a ][1] * stack.varLocal[ a ][1];
          strain[2] = strain[2] + DNDX[ a ][2] * stack.varLocal[ a ][2];
          strain[3] = strain[3] + DNDX[ a ][2] * stack.varLocal[ a ][1] + DNDX[ a ][1] * stack.varLocal[ a ][2];
          strain[4] = strain[4] + DNDX[ a ][2] * stack.varLocal[ a ][0] + DNDX[ a ][0] * stack.varLocal[ a ][2];
          strain[5] = strain[5] + DNDX[ a ][1] * stack.varLocal[ a ][0] + DNDX[ a ][0] * stack.varLocal[ a ][1];
        }

#if UPDATE_STRESS == 2
        constitutiveUpdate.SmallStrain( k, q, strain );
#else
        constitutiveUpdate.SmallStrainNoState( k, strain, stressLocal );
#endif

        for( localIndex c = 0; c < 6; ++c )
        {
#if UPDATE_STRESS == 2
          stressLocal[ c ] =  constitutiveUpdate.m_stress( k, q, c ) * (-DETJ);
#elif UPDATE_STRESS == 1
          stressLocal[ c ] = ( stressLocal[ c ] + constitutiveUpdate.m_stress( k, q, c ) ) *(-DETJ);
#else
          stressLocal[ c ] *= -DETJ;
#endif
        }


        for( localIndex a=0; a< NUM_NODES_PER_ELEM; ++a )
        {
          stack.fLocal[ a ][ 0 ] = stack.fLocal[ a ][ 0 ] + ( stressLocal[ 0 ] * DNDX[ a ][ 0 ] + stressLocal[ 5 ] * DNDX[ a ][ 1 ] + stressLocal[ 4 ] * DNDX[ a ][ 2 ] );
          stack.fLocal[ a ][ 1 ] = stack.fLocal[ a ][ 1 ] + ( stressLocal[ 5 ] * DNDX[ a ][ 0 ] + stressLocal[ 1 ] * DNDX[ a ][ 1 ] + stressLocal[ 3 ] * DNDX[ a ][ 2 ] );
          stack.fLocal[ a ][ 2 ] = stack.fLocal[ a ][ 2 ] + ( stressLocal[ 4 ] * DNDX[ a ][ 0 ] + stressLocal[ 3 ] * DNDX[ a ][ 1 ] + stressLocal[ 2 ] * DNDX[ a ][ 2 ] );
        }
    }

    template< typename PARAMETERS_TYPE,
              typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void stiffnessKernel( localIndex const ,
                          localIndex const ,
                          PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                          STACK_VARIABLE_TYPE & ) const
    {
    }

    template< typename PARAMETERS_TYPE,
              typename STACK_VARIABLE_TYPE>
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void integrationKernel( localIndex const GEOSX_UNUSED_PARAM(k),
                            localIndex const GEOSX_UNUSED_PARAM(q),
                            PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM(parameters),
                            STACK_VARIABLE_TYPE & GEOSX_UNUSED_PARAM(stack) ) const
    {
    }

    template< typename PARAMETERS_TYPE, typename STACK_VARIABLE_TYPE >
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    real64 postKernel( localIndex const k,
                       PARAMETERS_TYPE const & GEOSX_UNUSED_PARAM( parameters ),
                       STACK_VARIABLE_TYPE const & stack ) const
    {
      real64 meanForce = 0;

      for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, a );
        for( int b = 0; b < 3; ++b )
        {
          RAJA::atomicAdd< parallelDeviceAtomic >( &acc( nodeIndex, b ), stack.fLocal[ a ][ b ] );
        }
      }

      return meanForce;
    }

    template< typename POLICY,
              int NUM_QUADRATURE_POINTS,
              typename STACK_VARIABLES,
              typename PARAMETERS_TYPE,
              typename KERNEL_CLASS >
    static
    real64 Launch( localIndex const numElems,
                   PARAMETERS_TYPE const & parameters,
                   KERNEL_CLASS const & kernelClass )
    {
      return finiteElement::
             RegionLoop::
             Kernels< SUBREGION_TYPE,
                      CONSTITUTIVE_TYPE,
                      NUM_NODES_PER_ELEM,
                      NUM_NODES_PER_ELEM,
                      numTestDofPerSP,
                      numTrialDofPerSP>::template Launch< POLICY,
                                                          NUM_QUADRATURE_POINTS,
                                                          STACK_VARIABLES,
                                                          PARAMETERS_TYPE,
                                                          KERNEL_CLASS>( numElems, parameters, kernelClass );
    }

//    template< typename POLICY,
//              int NUM_QUADRATURE_POINTS,
//              typename STACK_VARIABLES,
//              typename PARAMETERS_TYPE,
//              typename KERNEL_CLASS >
//    static
//    real64 Launch( localIndex const numElems,
//                   PARAMETERS_TYPE const & parameters,
//                   KERNEL_CLASS const & kernelClass )
//    {
//      GEOSX_MARK_FUNCTION;
//      RAJA::ReduceMax< typename ReducePolicy<POLICY>::type, real64 > maxResidual( 0 );
//
//      forAll< POLICY >( numElems,
//                        [=] GEOSX_DEVICE ( localIndex const k )
//      {
//        STACK_VARIABLES stack;
//
//        kernelClass.preKernel( k, stack );
//        for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
//        {
//          kernelClass.updateKernel( k, q, stack );
//
//          kernelClass.stiffnessKernel( k, q, parameters, stack );
//
//          kernelClass.integrationKernel( k, q, parameters, stack );
//        }
//        kernelClass.postKernel( k, parameters, stack );
//      } );
//      return maxResidual.get();
//    }
//
  };
#undef CALCFEMSHAPE
#undef DNDX
#undef DETJ
#undef UPDATE_STRESS

};

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
