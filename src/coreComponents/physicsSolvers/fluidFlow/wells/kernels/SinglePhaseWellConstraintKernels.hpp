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
  static void assembleConstraintEquation( real64 const & time_n,
                                          WellControls & wellControls,
                                          BHPConstraint & constraint,
                                          WellElementSubRegion const & subRegion,
                                          string const & wellDofKey,
                                          localIndex const & rankOffset,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
  {
    // subRegion data
    localIndex const iwelemRef = subRegion.getTopWellElementIndex();
    arrayView1d< globalIndex const > const & wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & pres = subRegion.getField< fields::well::pressure >();

    constitutive::SingleFluidBase & fluidSeparator =  wellControls.getSingleFluidSeparator();
    arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & density = fluidSeparator.density();
    arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dDensity = fluidSeparator.dDensity();

    arrayView1d< real64 const > const wellElemGravCoef = subRegion.getField< fields::well::gravityCoefficient >();

    // setup row/column indices for constraint equation
    using ROFFSET_WJ = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;
    using COFFSET_WJ = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
    using Deriv = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

    localIndex const eqnRowIndex = wellElemDofNumber[iwelemRef] + ROFFSET_WJ::CONTROL - rankOffset;
    globalIndex dofColIndices[COFFSET_WJ::nDer]{};
    for( integer i = 0; i < COFFSET_WJ::nDer; ++i )
    {
      dofColIndices[ i ] = wellElemDofNumber[iwelemRef] + i;
    }
    // constraint data
    real64 const & targetBHP = constraint.getConstraintValue( time_n );
    real64 const & refGravCoef = constraint.getReferenceGravityCoef();

    // current constraint value
    real64 const & currentBHP =
      wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );

    // residual
    real64 controlEqn = currentBHP - targetBHP;

    // setup Jacobian terms
    real64 dControlEqn[2+IS_THERMAL]{};

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [pres,
                                density,
                                dDensity,
                                wellElemGravCoef,
                                &dControlEqn,
                                &iwelemRef,
                                &refGravCoef] ( localIndex const )
    {
      real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
      dControlEqn[COFFSET_WJ::dP] =   1.0 + dDensity[iwelemRef][0][Deriv::dP] *diffGravCoef;
      if constexpr ( IS_THERMAL )
      {
        dControlEqn[COFFSET_WJ::dT] =  dDensity[iwelemRef][0][Deriv::dT] * diffGravCoef;
      }
    } );

    // add solver matrices
    localRhs[eqnRowIndex] += controlEqn;
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                              dofColIndices,
                                                              dControlEqn,
                                                              COFFSET_WJ::nDer );
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
    arrayView1d< globalIndex const > const & wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );


    // setup row/column indices for constraint equation
    using ROFFSET_WJ = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;
    using COFFSET_WJ = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
    using Deriv = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

    localIndex const eqnRowIndex = wellElemDofNumber[iwelemRef] + ROFFSET_WJ::CONTROL - rankOffset;
    globalIndex dofColIndices[COFFSET_WJ::nDer]{};
    for( integer i = 0; i < COFFSET_WJ::nDer; ++i )
    {
      dofColIndices[ i ] = wellElemDofNumber[iwelemRef] + i;
    }

    // fluid data
    constitutive::SingleFluidBase & fluidSeparator =  wellControls.getSingleFluidSeparator();
    arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & density = fluidSeparator.density();
    arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dDensity = fluidSeparator.dDensity();

    // constraint data
    real64 const & targetVolRate = constraint.getConstraintValue( time_n );

    // current constraint value
    real64 & currentVolRate =
      wellControls.getReference< real64 >( WellControls::viewKeyStruct::currentVolRateString() );

    integer const useSurfaceConditions = wellControls.useSurfaceConditions();

    // residual
    real64 controlEqn =  currentVolRate - targetVolRate;

    // setup Jacobian terms
    real64 dControlEqn[2+IS_THERMAL]{};

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [currentVolRate,
                                density,
                                dDensity,
                                &dControlEqn,
                                &useSurfaceConditions,
                                &iwelemRef] ( localIndex const )
    {
      // compute the inverse of the total density and derivatives
      real64 const densInv = 1.0 / density[iwelemRef][0];

      dControlEqn[COFFSET_WJ::dP] = -( useSurfaceConditions ==  0 ) * dDensity[iwelemRef][0][Deriv::dP] * currentVolRate * densInv;
      dControlEqn[COFFSET_WJ::dQ] = densInv;
      if constexpr ( IS_THERMAL )
      {
        dControlEqn[COFFSET_WJ::dT] = -( useSurfaceConditions ==  0 ) * dDensity[iwelemRef][0][Deriv::dT] * currentVolRate * densInv;
      }

    } );

    // add solver matrices
    localRhs[eqnRowIndex] += controlEqn;
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                              dofColIndices,
                                                              dControlEqn,
                                                              COFFSET_WJ::nDer );
  }
};

} // end namespace wellConstraintKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONSTRAINTKERNELS_HPP
