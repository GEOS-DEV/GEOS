/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#include "constitutive/Solid/solidSelector.hpp"
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

#if defined(__CUDACC__)
  using KERNEL_POLICY = RAJA::cuda_exec< 256 >;
#elif defined(GEOSX_USE_OPENMP)
  using KERNEL_POLICY = RAJA::omp_parallel_for_exec;
#else
  using KERNEL_POLICY = RAJA::loop_exec;
#endif

inline void velocityUpdate( arrayView1d<R1Tensor> const & acceleration,
                            arrayView1d<R1Tensor> const & velocity,
                            real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  RAJA::forall< KERNEL_POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, acceleration.size() ),
                                 GEOSX_DEVICE_LAMBDA ( localIndex const i )
  {
    for (int j = 0; j < 3; ++j)
    {
      velocity[ i ][ j ] += dt * acceleration[ i ][ j ];
      acceleration[ i ][ j ] = 0;
    } 
  });
}

inline void velocityUpdate( arrayView1d<R1Tensor> const & acceleration,
                            arrayView1d<real64 const> const & mass, 
                            arrayView1d<R1Tensor> const & velocity,
                            real64 const dt,
                            LvArray::SortedArrayView<localIndex const, localIndex> const & indices )
{
  GEOSX_MARK_FUNCTION;

  RAJA::forall< KERNEL_POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, indices.size() ),
                                 GEOSX_DEVICE_LAMBDA ( localIndex const i )
  {
    localIndex const a = indices[ i ];
    for (int j = 0; j < 3; ++j)
    {
      acceleration[ a ][ j ] /= mass[ a ];
      velocity[ a ][ j ] += dt * acceleration[ a ][ j ];
    }
  });
}

inline void displacementUpdate( arrayView1d<R1Tensor const> const & velocity,
                                arrayView1d<R1Tensor> const & uhat,
                                arrayView1d<R1Tensor> const & u,
                                real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  RAJA::forall< KERNEL_POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, velocity.size() ),
                                 GEOSX_DEVICE_LAMBDA ( localIndex const i )
  {
    for (int j = 0; j < 3; ++j)
    {
      uhat[ i ][ j ] = velocity[ i ][ j ] * dt;
      u[ i ][ j ] += uhat[ i ][ j ];
    }
  });
}

template< int N >
inline void Integrate( const R2SymTensor& fieldvar,
                       arraySlice1d<R1Tensor const> const & dNdX,
                       real64 const& detJ,
                       real64 const& detF,
                       const R2Tensor& fInv,
                       R1Tensor * restrict const result )
{
  real64 const integrationFactor = detJ * detF;

  R2Tensor P;
  P.AijBkj( fieldvar, fInv );
  P *= integrationFactor;

  for( int a=0 ; a<N ; ++a )  // loop through all shape functions in element
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
template< typename KERNELWRAPPER, typename ... PARAMS>
inline real64
ElementKernelLaunchSelector( localIndex NUM_NODES_PER_ELEM,
                             localIndex NUM_QUADRATURE_POINTS,
                             constitutive::ConstitutiveBase * const constitutiveRelation,
                             PARAMS&& ... params)
{
  real64 rval = 0;

  constitutive::constitutiveUpdatePassThru( constitutiveRelation, [&]( auto & constitutive )
  {
    using CONSTITUTIVE_TYPE = TYPEOFREF( constitutive );
    if( NUM_NODES_PER_ELEM==8 && NUM_QUADRATURE_POINTS==8 )
    {
      rval = KERNELWRAPPER::template Launch<8,8, CONSTITUTIVE_TYPE>( &constitutive, std::forward<PARAMS>(params)... );
    }
    else if( NUM_NODES_PER_ELEM==4 && NUM_QUADRATURE_POINTS==1 )
    {
      rval = KERNELWRAPPER::template Launch<4,1, CONSTITUTIVE_TYPE>( &constitutive, std::forward<PARAMS>(params)... );
    }
  });
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
   * @param meanStress The mean stress at each element quadrature point
   * @param devStress The deviator stress at each element quadrature point.
   * @param dt The timestep
   * @return The achieved timestep.
   */
  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const constitutiveRelation,
          LvArray::SortedArrayView<localIndex const, localIndex> const & elementList,
          arrayView2d<localIndex const> const & elemsToNodes,
          arrayView3d< R1Tensor const> const & dNdX,
          arrayView2d<real64 const> const & detJ,
          arrayView1d<R1Tensor const> const & u,
          arrayView1d<R1Tensor const> const & vel,
          arrayView1d<R1Tensor> const & acc,
          arrayView2d<real64> const & meanStress,
          arrayView2d<R2SymTensor> const & devStress,
          real64 const dt )
  {
    forall_in_set<serialPolicy>( elementList.values(),
                              elementList.size(),
                              GEOSX_LAMBDA ( localIndex k) mutable
    {
      R1Tensor v_local[NUM_NODES_PER_ELEM];
      R1Tensor u_local[NUM_NODES_PER_ELEM];
      R1Tensor f_local[NUM_NODES_PER_ELEM];

      CopyGlobalToLocal<NUM_NODES_PER_ELEM,R1Tensor>( elemsToNodes[k],
                                                      u, vel,
                                                      u_local, v_local );

      //Compute Quadrature
      for( localIndex q = 0 ; q<NUM_QUADRATURE_POINTS ; ++q)
      {
        R2Tensor dUhatdX, dUdX;
        CalculateGradients<NUM_NODES_PER_ELEM>( dUhatdX, dUdX, v_local, u_local, dNdX[k][q]);
        dUhatdX *= dt;

        R2Tensor F,Ldt, fInv;

        // calculate du/dX
        F = dUhatdX;
        F *= 0.5;
        F += dUdX;
        F.PlusIdentity(1.0);
        fInv.Inverse(F);

        // chain rule: calculate dv/du = dv/dX * dX/du
        Ldt.AijBjk(dUhatdX, fInv);

        // calculate gradient (end of step)
        F = dUhatdX;
        F += dUdX;
        F.PlusIdentity(1.0);
        real64 detF = F.Det();
        fInv.Inverse(F);


        R2Tensor Rot;
        R2SymTensor Dadt;
        HughesWinget(Rot, Dadt, Ldt);

        constitutiveRelation->StateUpdatePoint( k, q, Dadt, Rot, 0);

        R2SymTensor TotalStress = devStress[k][q];
        TotalStress.PlusIdentity( meanStress[k][q] );

        Integrate<NUM_NODES_PER_ELEM>( TotalStress, dNdX[k][q], detJ[k][q], detF, fInv, f_local );
      }//quadrature loop


      AddLocalToGlobal<NUM_NODES_PER_ELEM>( elemsToNodes[k], f_local, acc );
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
  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS, typename CONSTITUTIVE_TYPE >
  static inline real64
  Launch( CONSTITUTIVE_TYPE * const constitutiveRelation,
          localIndex const numElems,
          real64 const dt,
          arrayView3d<R1Tensor const> const & dNdX,
          arrayView2d<real64 const > const& detJ,
          FiniteElementBase const * const fe,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView2d< localIndex const > const & elemsToNodes,
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
          DofManager const * const dofManager,
          ParallelMatrix * const matrix,
          ParallelVector * const rhs )
  {
    GEOS_ERROR("SolidMechanicsLagrangianFEM::ImplicitElementKernelWrapper::Launch() not implemented");
    return 0;
  }

};

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx
