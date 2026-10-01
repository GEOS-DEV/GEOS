/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseWellConstraintKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLCONSTRAINTKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLCONSTRAINTKERNELS_HPP

#include "codingUtilities/Utilities.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"


#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellBHPConstraints.hpp"
#include "physicsSolvers/fluidFlow/wells/WellVolumeRateConstraint.hpp"


namespace geos
{

namespace singlePhaseWellConstraintKernels
{

/******************************** ControlEquationHelper ********************************/


template< integer IS_THERMAL >
struct ConstraintHelper
{
  GEOS_HOST_DEVICE
  static void setTemperatureDerivative( real64 * const dControlEqn,
                                        real64 const value )
  {
    if constexpr ( IS_THERMAL )
    {
      dControlEqn[ singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >::dT ] = value;
    }
  }

  template< BHPConstraintTypeId I >
  static void assembleConstraintEquation( real64 const & time_n,
                                          WellControls & wellControls,
                                          BHPConstraint< I > & constraint,
                                          WellElementSubRegion const & subRegion,
                                          string const & wellDofKey,
                                          localIndex const & rankOffset,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    // subRegion data
    localIndex const iwelemRef = subRegion.getTopWellElementIndex();
    arrayView1d< globalIndex const > const wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    constitutive::SingleFluidBase & fluidSeparator =  wellControls.getSingleFluidSeparator();
    arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dDensity = fluidSeparator.dDensity();

    arrayView1d< real64 const > const wellElemGravCoef = subRegion.getField< fields::well::gravityCoefficient >();

    // setup row/column indices for constraint equation
    using ROFFSET_WJ = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;
    using COFFSET_WJ = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
    using Deriv = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

    // constraint data
    real64 const targetBHP = constraint.getConstraintValue( time_n );
    real64 const refGravCoef = constraint.getReferenceGravityCoef();

    // current constraint value
    real64 const currentBHP =
      wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );

    // The separator is updated on the host. Copy its small set of scalar
    // derivatives into the device closure and keep the matrix and RHS on device.
    real64 const dDensity_dP = dDensity[iwelemRef][0][Deriv::dP];
    real64 dDensity_dT = 0.0;
    if constexpr ( IS_THERMAL )
    {
      dDensity_dT = dDensity[iwelemRef][0][Deriv::dT];
    }

    forAll< parallelDevicePolicy<> >( 1, [=] GEOS_HOST_DEVICE ( localIndex const )
    {
      globalIndex const dofNumber = wellElemDofNumber[iwelemRef];
      localIndex const eqnRowIndex = LvArray::integerConversion< localIndex >( dofNumber + ROFFSET_WJ::CONTROL - rankOffset );
      globalIndex dofColIndices[COFFSET_WJ::nDer]{};
      for( integer i = 0; i < COFFSET_WJ::nDer; ++i )
      {
        dofColIndices[ i ] = dofNumber + i;
      }

      real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
      real64 dControlEqn[2+IS_THERMAL]{};
      dControlEqn[COFFSET_WJ::dP] = 1.0 + dDensity_dP * diffGravCoef;
      setTemperatureDerivative( dControlEqn, dDensity_dT * diffGravCoef );

      RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], currentBHP - targetBHP );
      localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                        dofColIndices,
                                                                        dControlEqn,
                                                                        COFFSET_WJ::nDer );
    } );
  }
  template< template< typename U > class T, typename U=VolumeRateConstraint >
  static void assembleConstraintEquation( real64 const & time_n,
                                          WellControls & wellControls,
                                          T< VolumeRateConstraint > & constraint,
                                          WellElementSubRegion const & subRegion,
                                          string const & wellDofKey,
                                          localIndex const & rankOffset,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    // subRegion data

    localIndex const iwelemRef = subRegion.getTopWellElementIndex();
    arrayView1d< globalIndex const > const wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // setup row/column indices for constraint equation
    using ROFFSET_WJ = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;
    using COFFSET_WJ = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
    using Deriv = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

    // fluid data
    constitutive::SingleFluidBase & fluidSeparator =  wellControls.getSingleFluidSeparator();
    arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const density = fluidSeparator.density();
    arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const dDensity = fluidSeparator.dDensity();

    // constraint data
    real64 const targetVolRate = constraint.getConstraintValue( time_n );

    // current constraint value
    real64 const currentVolRate =
      wellControls.getReference< real64 >( WellControls::viewKeyStruct::currentVolRateString() );

    integer const useSurfaceConditions = wellControls.useSurfaceConditions();

    // The separator is updated on the host. Copy its small set of scalar
    // properties into the device closure and keep the matrix and RHS on device.
    real64 const densityRef = density[iwelemRef][0];
    real64 const dDensity_dP = dDensity[iwelemRef][0][Deriv::dP];
    real64 dDensity_dT = 0.0;
    if constexpr ( IS_THERMAL )
    {
      dDensity_dT = dDensity[iwelemRef][0][Deriv::dT];
    }

    forAll< parallelDevicePolicy<> >( 1, [=] GEOS_HOST_DEVICE ( localIndex const )
    {
      globalIndex const dofNumber = wellElemDofNumber[iwelemRef];
      localIndex const eqnRowIndex = LvArray::integerConversion< localIndex >( dofNumber + ROFFSET_WJ::CONTROL - rankOffset );
      globalIndex dofColIndices[COFFSET_WJ::nDer]{};
      for( integer i = 0; i < COFFSET_WJ::nDer; ++i )
      {
        dofColIndices[ i ] = dofNumber + i;
      }

      real64 const densInv = 1.0 / densityRef;
      real64 dControlEqn[2+IS_THERMAL]{};
      dControlEqn[COFFSET_WJ::dP] = -( useSurfaceConditions == 0 ) * dDensity_dP * currentVolRate * densInv;
      dControlEqn[COFFSET_WJ::dQ] = densInv;
      setTemperatureDerivative( dControlEqn,
                                -( useSurfaceConditions == 0 ) * dDensity_dT * currentVolRate * densInv );

      RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], currentVolRate - targetVolRate );
      localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                        dofColIndices,
                                                                        dControlEqn,
                                                                        COFFSET_WJ::nDer );
    } );
  }
};

} // end namespace wellConstraintKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINTKERNELS_HPP
