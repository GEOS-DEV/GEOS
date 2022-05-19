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
 * @file SinglePhaseBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "functions/TableFunction.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"

namespace geosx
{

namespace singlePhaseBaseKernels
{

/******************************** MobilityKernel ********************************/

struct MobilityKernel
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  compute( real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & visc,
           real64 const & dVisc_dPres,
           real64 & mob,
           real64 & dMob_dPres )
  {
    mob = dens / visc;
    dMob_dPres = dDens_dPres / visc - mob / visc * dVisc_dPres;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  compute( real64 const & dens,
           real64 const & visc,
           real64 & mob )
  {
    mob = dens / visc;
  }

  template< typename POLICY >
  static void launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & dDens_dPres,
                      arrayView2d< real64 const > const & visc,
                      arrayView2d< real64 const > const & dVisc_dPres,
                      arrayView1d< real64 > const & mob,
                      arrayView1d< real64 > const & dMob_dPres )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      compute( dens[a][0],
               dDens_dPres[a][0],
               visc[a][0],
               dVisc_dPres[a][0],
               mob[a],
               dMob_dPres[a] );
    } );
  }

  template< typename POLICY >
  static void launch( SortedArrayView< localIndex const > targetSet,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & dDens_dPres,
                      arrayView2d< real64 const > const & visc,
                      arrayView2d< real64 const > const & dVisc_dPres,
                      arrayView1d< real64 > const & mob,
                      arrayView1d< real64 > const & dMob_dPres )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      compute( dens[a][0],
               dDens_dPres[a][0],
               visc[a][0],
               dVisc_dPres[a][0],
               mob[a],
               dMob_dPres[a] );
    } );
  }

  template< typename POLICY >
  static void launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & visc,
                      arrayView1d< real64 > const & mob )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      compute( dens[a][0],
               visc[a][0],
               mob[a] );
    } );
  }

  template< typename POLICY >
  static void launch( SortedArrayView< localIndex const > targetSet,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & visc,
                      arrayView1d< real64 > const & mob )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      compute( dens[a][0],
               visc[a][0],
               mob[a] );
    } );
  }
};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @brief Internal struct to provide no-op defaults used in the inclusion
 *   of lambda functions into kernel component functions.
 * @struct NoOpFunc
 */
struct NoOpFunc
{
  template< typename ... Ts >
  GEOSX_HOST_DEVICE
  constexpr void
  operator()( Ts && ... ) const {}
};

/**
 * @class ElementBasedAssemblyKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */
template< typename SUBREGION_TYPE >
class ElementBasedAssemblyKernel
{

public:

  /**
   * @brief Constructor
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              string const dofKey,
                              SUBREGION_TYPE const & subRegion,
                              constitutive::SingleFluidBase const & fluid,
                              constitutive::CoupledSolidBase const & solid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.template getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_deltaVolume( subRegion.template getExtrinsicData< extrinsicMeshData::flow::deltaVolume >() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_porosityNew( solid.getPorosity() ),
    m_dPoro_dPres( solid.getDporosity_dPressure() ),
    m_density_n( fluid.density_n() ),
    m_density( fluid.density() ),
    m_dDensity_dPres( fluid.dDensity_dPressure() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    // Pore volume information

    /// Pore volume at time n+1
    real64 poreVolume = 0.0;

    /// Pore volume at the previous converged time step
    real64 poreVolume_n = 0.0;

    /// Derivative of pore volume with respect to pressure
    real64 dPoreVolume_dPres = 0.0;

    // Residual information

    /// Index of the matrix row/column corresponding to the dof in this element
    globalIndex dofNumber;

    /// Storage for the element local residual vector
    real64 localResidual;

    /// Storage for the element local Jacobian matrix
    real64 localJacobian;

  };

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOSX_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    // initialize the pore volume
    stack.poreVolume = ( m_volume[ei] + m_deltaVolume[ei] ) * m_porosityNew[ei][0];
    stack.poreVolume_n = m_volume[ei] * m_porosity_n[ei][0];
    stack.dPoreVolume_dPres = ( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dPres[ei][0];

    // set degree of freedom index for this element
    stack.dofNumber = m_dofNumber[ei];
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOSX_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack,
                            FUNC && kernelOp = NoOpFunc{} ) const
  {
    // Residual contribution is mass conservation in the cell
    stack.localResidual = stack.poreVolume * m_density[ei][0] - stack.poreVolume_n * m_density_n[ei][0];

    // Derivative of residual wrt to pressure in the cell
    stack.localJacobian = stack.dPoreVolume_dPres * m_density[ei][0] + m_dDensity_dPres[ei][0] * stack.poreVolume;

    // Customize the kernel with this lambda
    kernelOp();
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void complete( localIndex const GEOSX_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {
    localIndex const localElemDof = stack.dofNumber - m_rankOffset;

    // add contribution to global residual and jacobian (no need for atomics here)
    m_localMatrix.addToRow< serialAtomic >( localElemDof, &stack.dofNumber, &stack.localJacobian, 1 );
    m_localRhs[localElemDof] += stack.localResidual;
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.elemGhostRank( ei ) >= 0 )
      {
        return;
      }

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.computeAccumulation( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;
  arrayView1d< real64 const > const m_deltaVolume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosity_n;
  arrayView2d< real64 const > const m_porosityNew;
  arrayView2d< real64 const > const m_dPoro_dPres;

  /// Views on density
  arrayView2d< real64 const > const m_density_n;
  arrayView2d< real64 const > const m_density;
  arrayView2d< real64 const > const m_dDensity_dPres;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

};

/**
 * @class SurfaceElementBasedAssemblyKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation in SurfaceElementSubRegion
 */
class SurfaceElementBasedAssemblyKernel : public ElementBasedAssemblyKernel< SurfaceElementSubRegion >
{

public:

  using Base = ElementBasedAssemblyKernel< SurfaceElementSubRegion >;

  /**
   * @brief Constructor
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  SurfaceElementBasedAssemblyKernel( globalIndex const rankOffset,
                                     string const dofKey,
                                     SurfaceElementSubRegion const & subRegion,
                                     constitutive::SingleFluidBase const & fluid,
                                     constitutive::CoupledSolidBase const & solid,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs )
#if ALLOW_CREATION_MASS
    , m_creationMass( subRegion.getReference< array1d< real64 > >( SurfaceElementSubRegion::viewKeyStruct::creationMassString() ) )
#endif
  {
#if !defined(ALLOW_CREATION_MASS)
    static_assert( true, "must have ALLOW_CREATION_MASS defined" );
#endif
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOSX_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            Base::StackVariables & stack ) const
  {
    Base::computeAccumulation( ei, stack, [&] ()
    {
#if ALLOW_CREATION_MASS
      if( Base::m_volume[ei] * Base::m_density_n[ei][0] > 1.1 * m_creationMass[ei] )
      {
        stack.localResidual += m_creationMass[ei] * 0.25;
      }
#endif
    } );
  }

protected:

#if ALLOW_CREATION_MASS
  arrayView1d< real64 const > const m_creationMass;
#endif

};

/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   CellElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    ElementBasedAssemblyKernel< CellElementSubRegion >
    kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
    ElementBasedAssemblyKernel< CellElementSubRegion >::template launch< POLICY >( subRegion.size(), kernel );
  }

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   SurfaceElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    SurfaceElementBasedAssemblyKernel
      kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
    SurfaceElementBasedAssemblyKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{
  template< typename POLICY >
  static void launch( arrayView1d< real64 const > const & localResidual,
                      globalIndex const rankOffset,
                      arrayView1d< globalIndex const > const & presDofNumber,
                      arrayView1d< integer const > const & ghostRank,
                      arrayView1d< real64 const > const & volume,
                      arrayView2d< real64 const > const & dens_n,
                      arrayView2d< real64 const > const & poro_n,
                      real64 * localResidualNorm )
  {
    RAJA::ReduceSum< ReducePolicy< POLICY >, real64 > localSum( 0.0 );
    RAJA::ReduceSum< ReducePolicy< POLICY >, real64 > normSum( 0.0 );
    RAJA::ReduceSum< ReducePolicy< POLICY >, localIndex > count( 0 );

    forAll< POLICY >( presDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      if( ghostRank[a] < 0 )
      {
        localIndex const lid = presDofNumber[a] - rankOffset;
        real64 const val = localResidual[lid];
        localSum += val * val;
        normSum += poro_n[a][0] * dens_n[a][0] * volume[a];
        count += 1;
      }
    } );

    localResidualNorm[0] += localSum.get();
    localResidualNorm[1] += normSum.get();
    localResidualNorm[2] += count.get();
  }
};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY >
  static localIndex launch( arrayView1d< real64 const > const & localSolution,
                            globalIndex const rankOffset,
                            arrayView1d< globalIndex const > const & presDofNumber,
                            arrayView1d< integer const > const & ghostRank,
                            arrayView1d< real64 const > const & pres,
                            real64 const scalingFactor )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, localIndex > minVal( 1 );

    forAll< POLICY >( presDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 && presDofNumber[ei] >= 0 )
      {
        localIndex const lid = presDofNumber[ei] - rankOffset;
        real64 const newPres = pres[ei] + scalingFactor * localSolution[lid];

        if( newPres < 0.0 )
        {
          minVal.min( 0 );
        }
      }

    } );
    return minVal.get();
  }

};

/******************************** HydrostaticPressureKernel ********************************/

struct HydrostaticPressureKernel
{

  template< typename FLUID_WRAPPER >
  static bool
  computeHydrostaticPressure( integer const maxNumEquilIterations,
                              real64 const & equilTolerance,
                              real64 const (&gravVector)[ 3 ],
                              FLUID_WRAPPER fluidWrapper,
                              TableFunction::KernelWrapper tempTableWrapper, 
                              real64 const & refElevation,
                              real64 const & refPres,
                              real64 const & refDens,
                              real64 const & newElevation,
                              real64 & newPres,
                              real64 & newDens )
  {
    // Step 1: guess the pressure with the refDensity

    real64 const gravCoef = gravVector[2] * ( refElevation - newElevation );
    real64 const temp = tempTableWrapper.compute( &newElevation ); 
    real64 pres0 = refPres - refDens * gravCoef;
    real64 pres1 = 0.0;

    // Step 2: compute the mass density at this elevation using the guess, and update pressure

    real64 dens = 0.0;
    real64 visc = 0.0;
    fluidWrapper.compute( pres0, temp, dens, visc );
    pres1 = refPres - 0.5 * ( refDens + dens ) * gravCoef;

    // Step 3: fixed-point iteration until convergence

    bool equilHasConverged = false;
    for( localIndex eqIter = 0; eqIter < maxNumEquilIterations; ++eqIter )
    {

      // check convergence
      equilHasConverged = ( LvArray::math::abs( pres0 - pres1 ) < equilTolerance );
      pres0 = pres1;

      // if converged, move on
      if( equilHasConverged )
      {
        break;
      }

      // compute the density at this elevation using the previous pressure, and compute the new pressure
      fluidWrapper.compute( pres0, temp, dens, visc );
      pres1 = refPres - 0.5 * ( refDens + dens ) * gravCoef;
    }

    // Step 4: save the hydrostatic pressure and the corresponding density

    newPres = pres1;
    newDens = dens;

    return equilHasConverged;
  }


  template< typename FLUID_WRAPPER >
  static bool
  launch( localIndex const size,
          integer const maxNumEquilIterations,
          real64 const equilTolerance,
          real64 const (&gravVector)[ 3 ],
          real64 const & minElevation,
          real64 const & elevationIncrement,
          real64 const & datumElevation,
          real64 const & datumPres,
          FLUID_WRAPPER fluidWrapper,
          TableFunction::KernelWrapper tempTableWrapper,
          arrayView1d< arrayView1d< real64 > const > elevationValues, 
          arrayView1d< real64 > pressureValues )
  {
    bool hasConverged = true;

    // Step 1: compute the mass density at the datum elevation

    real64 const datumTemp = tempTableWrapper.compute( &datumElevation ); 

    real64 datumDens = 0.0;
    real64 datumVisc = 0.0;
    
    fluidWrapper.compute( datumPres,
                          datumTemp, 
                          datumDens,
                          datumVisc );

    // Step 2: find the closest elevation to datumElevation

    forAll< parallelHostPolicy >( size, [=] ( localIndex const i )
    {
      real64 const elevation = minElevation + i * elevationIncrement;
      elevationValues[0][i] = elevation;
    } );
    integer const iRef = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                 elevationValues[0].size(),
                                                                 datumElevation );


    // Step 3: compute the mass density and pressure at the reference elevation

    array1d< real64 > dens( pressureValues.size() );

    bool const hasConvergedRef =
      computeHydrostaticPressure( maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  tempTableWrapper, 
                                  datumElevation,
                                  datumPres,
                                  datumDens,
                                  elevationValues[0][iRef],
                                  pressureValues[iRef],
                                  dens[iRef] );
    if( !hasConvergedRef )
    {
      return false;
    }


    // Step 4: for each elevation above the reference elevation, compute the pressure

    localIndex const numEntriesAboveRef = size - iRef - 1;
    forAll< serialPolicy >( numEntriesAboveRef, [=, &hasConverged] ( localIndex const i )
    {
      bool const hasConvergedAboveRef =
        computeHydrostaticPressure( maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    tempTableWrapper, 
                                    elevationValues[0][iRef+i],
                                    pressureValues[iRef+i],
                                    dens[iRef+i],
                                    elevationValues[0][iRef+i+1],
                                    pressureValues[iRef+i+1],
                                    dens[iRef+i+1] );
      if( !hasConvergedAboveRef )
      {
        hasConverged = false;
      }


    } );

    // Step 5: for each elevation below the reference elevation, compute the pressure

    localIndex const numEntriesBelowRef = iRef;
    forAll< serialPolicy >( numEntriesBelowRef, [=, &hasConverged] ( localIndex const i )
    {
      bool const hasConvergedBelowRef =
        computeHydrostaticPressure( maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    tempTableWrapper, 
                                    elevationValues[0][iRef-i],
                                    pressureValues[iRef-i],
                                    dens[iRef-i],
                                    elevationValues[0][iRef-i-1],
                                    pressureValues[iRef-i-1],
                                    dens[iRef-i-1] );
      if( !hasConvergedBelowRef )
      {
        hasConverged = false;
      }
    } );

    return hasConverged;
  }
};


} // namespace singlePhaseBaseKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP
