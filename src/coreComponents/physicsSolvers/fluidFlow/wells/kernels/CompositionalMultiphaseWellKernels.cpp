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
 * @file CompositionalMultiphaseWellKernels.cpp
 */

#include "CompositionalMultiphaseWellKernels.hpp"

//#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
// TODO: move keys to WellControls
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"

namespace geos
{

namespace compositionalMultiphaseWellKernels
{

using namespace constitutive;


/******************************** PressureRelationKernel ********************************/

template< integer NC, integer IS_THERMAL >
GEOS_HOST_DEVICE
void
PressureRelationKernel::
  compute( real64 const & gravCoef,
           real64 const & gravCoefNext,
           real64 const & pres,
           real64 const & presNext,
           real64 const & totalMassDens,
           real64 const & totalMassDensNext,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDensNext,
           real64 & localPresRel,
           real64 ( & localPresRelJacobian )[2*(NC+1 + IS_THERMAL)] )
{
  // local working variables and arrays
  real64 dAvgMassDens_dCompCurrent[NC]{};
  real64 dAvgMassDens_dCompNext[NC]{};

  // compute the average density at the interface between well elements
  real64 const avgMassDens = 0.5 * ( totalMassDensNext + totalMassDens );
  real64 const dAvgMassDens_dPresNext    = 0.5 * dTotalMassDensNext[Deriv::dP];
  real64 const dAvgMassDens_dPresCurrent = 0.5 * dTotalMassDens[Deriv::dP];
  for( integer ic = 0; ic < NC; ++ic )
  {
    dAvgMassDens_dCompNext[ic]    = 0.5 * dTotalMassDensNext[Deriv::dC+ic];
    dAvgMassDens_dCompCurrent[ic] = 0.5 * dTotalMassDens[Deriv::dC+ic];
  }

  // compute depth diff times acceleration
  real64 const gravD = gravCoefNext - gravCoef;

  // TODO: add friction and acceleration terms

  localPresRel = ( presNext - pres - avgMassDens * gravD );

  // localPresRelJacbain contains dP, dC and potentially dT derivatives for neighboring well elements
  // TAG::NEXT is 1, CURRENT is 0 , not sure why indexes are setup as below
  localPresRelJacobian[TAG::NEXT *(NC+1+IS_THERMAL)]    = ( 1 - dAvgMassDens_dPresNext * gravD );
  localPresRelJacobian[TAG::CURRENT *(NC+1+IS_THERMAL)] = ( -1 - dAvgMassDens_dPresCurrent * gravD );

  for( integer ic = 0; ic < NC; ++ic )
  {
    localPresRelJacobian[TAG::NEXT *(NC+1+IS_THERMAL) + ic+1]    = -dAvgMassDens_dCompNext[ic] * gravD;
    localPresRelJacobian[TAG::CURRENT *(NC+1+IS_THERMAL) + ic+1] = -dAvgMassDens_dCompCurrent[ic] * gravD;
  }
  if constexpr ( IS_THERMAL )
  {
    localPresRelJacobian[TAG::NEXT *(NC+1+IS_THERMAL)+NC+1]    =  -0.5 * dTotalMassDensNext[Deriv::dT]* gravD;
    localPresRelJacobian[TAG::CURRENT *(NC+1+IS_THERMAL)+NC+1] = -0.5 * dTotalMassDens[Deriv::dT]* gravD;
  }
}

template< integer NC, integer IS_THERMAL >
void
PressureRelationKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< integer const > const elemStatus,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & wellElemTotalMassDens,
          arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens,
          bool & controlHasSwitched,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  using COFFSET_WJ = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;

  RAJA::ReduceMax< parallelDeviceReduce, localIndex > switchControl( 0 );

  // loop over the well elements to compute the pressure relations between well elements
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    if( elemStatus[iwelem] ==  WellElementSubRegion::WellElemStatus::CLOSED )
    {
      return;
    }

    localIndex const iwelemNext = nextWellElemIndex[iwelem];

    if( iwelemNext >= 0 )   // if iwelemNext >= 0, form momentum equation
    {

      real64 localPresRel = 0;
      real64 localPresRelJacobian[2*(NC+1+IS_THERMAL)]{};

      compute< NC, IS_THERMAL >(
        wellElemGravCoef[iwelem],
        wellElemGravCoef[iwelemNext],
        wellElemPressure[iwelem],
        wellElemPressure[iwelemNext],
        wellElemTotalMassDens[iwelem],
        wellElemTotalMassDens[iwelemNext],
        dWellElemTotalMassDens[iwelem],
        dWellElemTotalMassDens[iwelemNext],
        localPresRel,
        localPresRelJacobian );


      // local working variables and arrays
      globalIndex dofColIndices[2*(NC+1+IS_THERMAL)];

      globalIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::CONTROL - rankOffset;
      dofColIndices[TAG::NEXT *(NC+1+IS_THERMAL)]    = wellElemDofNumber[iwelemNext] + COFFSET_WJ::dP;
      dofColIndices[TAG::CURRENT *(NC+1+IS_THERMAL)] = wellElemDofNumber[iwelem] + COFFSET_WJ::dP;

      for( integer ic = 0; ic < NC; ++ic )
      {
        dofColIndices[TAG::NEXT *(NC+1+IS_THERMAL) + ic+1]    = wellElemDofNumber[iwelemNext] + COFFSET_WJ::dC + ic;
        dofColIndices[TAG::CURRENT *(NC+1+IS_THERMAL) + ic+1] = wellElemDofNumber[iwelem] + COFFSET_WJ::dC + ic;
      }
      if constexpr ( IS_THERMAL )
      {
        dofColIndices[TAG::NEXT *(NC+1+IS_THERMAL)+NC+1]    = wellElemDofNumber[iwelemNext] + COFFSET_WJ::dT;
        dofColIndices[TAG::CURRENT *(NC+1+IS_THERMAL)+NC+1] = wellElemDofNumber[iwelem] + COFFSET_WJ::dT;
      }
      if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                          dofColIndices,
                                                                          localPresRelJacobian,
                                                                          2 * (NC+1+IS_THERMAL) );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], localPresRel );
      }
    }
  } );
  controlHasSwitched = ( switchControl.get() == 1 );
}

#define INST_PressureRelationKernel( NC, IS_THERMAL ) \
  template \
  void PressureRelationKernel:: \
    launch< NC, IS_THERMAL >( localIndex const size, \
                              globalIndex const rankOffset, \
                              arrayView1d< integer const > const elemStatus, \
                              arrayView1d< globalIndex const > const & wellElemDofNumber, \
                              arrayView1d< real64 const > const & wellElemGravCoef, \
                              arrayView1d< localIndex const > const & nextWellElemIndex, \
                              arrayView1d< real64 const > const & wellElemPressure, \
                              arrayView1d< real64 const > const & wellElemTotalMassDens, \
                              arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens, \
                              bool & controlHasSwitched, \
                              CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                              arrayView1d< real64 > const & localRhs )

INST_PressureRelationKernel( 1, 0 );
INST_PressureRelationKernel( 1, 1 );
INST_PressureRelationKernel( 2, 0 );
INST_PressureRelationKernel( 2, 1 );
INST_PressureRelationKernel( 3, 0 );
INST_PressureRelationKernel( 3, 1 );
INST_PressureRelationKernel( 4, 0 );
INST_PressureRelationKernel( 4, 1 );
INST_PressureRelationKernel( 5, 0 );
INST_PressureRelationKernel( 5, 1 );

void
PresTempCompFracInitializationKernel::
  launch( localIndex const perforationSize,
          localIndex const subRegionSize,
          integer const numComps,
          integer const numPhases,
          WellControls const & wellControls,
          real64 const & currentTime,
          ElementViewConst< arrayView1d< real64 const > > const & resPres,
          ElementViewConst< arrayView1d< real64 const > > const & resTemp,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const & resCompDens,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & resPhaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseMassDens,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< integer const > const & perfStatus,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 > const & wellElemPres,
          arrayView1d< real64 > const & wellElemTemp,
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompFrac )
{
  integer constexpr MAX_NUM_COMP = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

  real64 const targetBHP = wellControls.getTargetBHP( currentTime );
  real64 const refWellElemGravCoef = wellControls.getReferenceGravityCoef();
  real64 const initialPresCoef = wellControls.getInitialPressureCoefficient();
  WellControls::Control const currentControl = wellControls.getControl();
  bool const isProducer = wellControls.isProducer();



  // Step 1: we loop over all the perforations on this rank to compute the following quantities:
  //   - Sum of total mass densities over the perforated reservoir elements
  //   - Sum of the temperatures over the perforated reservoir elements
  //   - Sum of the component fractions over the perforated reservoir elements
  // In passing, we save the min gravCoef difference between the reference depth and the perforation depth
  // Note that we use gravCoef instead of depth for the (unlikely) case in which the gravityVector is not aligned with z

  RAJA::ReduceSum< parallelDeviceReduce, integer > numOpenPerfs( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumTotalMassDens( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumTemp( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumCompFrac[MAX_NUM_COMP]{};
  RAJA::ReduceMin< parallelDeviceReduce, real64 > localMinGravCoefDiff( 1e9 );

  forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
  {
    if( perfStatus[iperf] )
    {
      numOpenPerfs += 1;
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];

      // save the min gravCoef difference between the reference depth and the perforation depth (times g)
      localMinGravCoefDiff.min( LvArray::math::abs( refWellElemGravCoef - perfGravCoef[iperf] ) );

      // increment the temperature
      sumTemp += resTemp[er][esr][ei];

      // increment the total mass density
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        sumTotalMassDens += resPhaseVolFrac[er][esr][ei][ip] * resPhaseMassDens[er][esr][ei][0][ip];
      }

      // increment the component fractions
      real64 perfTotalDens = 0.0;
      for( integer ic = 0; ic < numComps; ++ic )
      {
        perfTotalDens += resCompDens[er][esr][ei][ic];
      }
      for( integer ic = 0; ic < numComps; ++ic )
      {
        sumCompFrac[ic] += resCompDens[er][esr][ei][ic] / perfTotalDens;
      }
    }
  } );
  real64 const minGravCoefDiff = MpiWrapper::min( localMinGravCoefDiff.get() );

  integer totalOpenPerfs = MpiWrapper::sum( numOpenPerfs.get() );

  // Step 2: we assign average quantities over the well (i.e., over all the ranks)
  // For composition and temperature, we make a distinction between injection and production

  // for total mass density, we always use the values of the perforated reservoir elements, even for injectors
  real64 const avgTotalMassDens = MpiWrapper::sum( sumTotalMassDens.get() ) / totalOpenPerfs;

  stackArray1d< real64, MAX_NUM_COMP > avgCompFrac( numComps );
  real64 avgTemp = 0;

  // for a producer, we use the temperature and component fractions from the reservoir
  if( isProducer )
  {
    // use average temperature from reservoir
    avgTemp = MpiWrapper::sum( sumTemp.get() ) / totalOpenPerfs;

    // use average comp frac from reservoir
    for( integer ic = 0; ic < numComps; ++ic )
    {
      avgCompFrac[ic] = MpiWrapper::sum( sumCompFrac[ic].get() ) / totalOpenPerfs;
    }
  }
  // for an injector, we use the injection stream values
  else
  {
    // use temperature from injection stream
    avgTemp = wellControls.getInjectionTemperature();

    // use comp frac from injection stream
    for( integer ic = 0; ic < numComps; ++ic )
    {
      avgCompFrac[ic] = wellControls.getInjectionStream()[ic];
    }
  }



  // Step 3: we compute the approximate pressure at the reference depth
  // We make a distinction between pressure-controlled wells and rate-controlled wells

  real64 refPres = 0.0;

  // if the well is controlled by pressure, initialize the reference pressure at the target pressure
  if( currentControl == WellControls::Control::BHP )
  {
    refPres = targetBHP;
  }
  // if the well is controlled by rate, initialize the reference pressure using the pressure at the closest perforation
  else
  {
    RAJA::ReduceMin< parallelDeviceReduce, real64 > localRefPres( 1e9 );
    real64 const alpha = ( isProducer ) ? 1 - initialPresCoef : 1 + initialPresCoef;

    forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
    {
      if( perfStatus[iperf] )
      {
        // get the reservoir (sub)region and element indices
        localIndex const er = resElementRegion[iperf];
        localIndex const esr = resElementSubRegion[iperf];
        localIndex const ei = resElementIndex[iperf];

        // get the perforation pressure and save the estimated reference pressure
        real64 const gravCoefDiff = LvArray::math::abs( refWellElemGravCoef - perfGravCoef[iperf] );
        if( isZero( gravCoefDiff - minGravCoefDiff ) )
        {
          localRefPres.min( alpha * resPres[er][esr][ei] + avgTotalMassDens * ( refWellElemGravCoef - perfGravCoef[iperf] ) );
        }
      }
    } );
    refPres = MpiWrapper::min( localRefPres.get() );
  }



  // Step 4: we are ready to assign the primary variables on the well elements:
  //  - pressure: hydrostatic pressure using our crude approximation of the total mass density
  //  - temperature: uniform, using the average temperature computed above
  //  - component fraction: uniform, using the average component fraction computed above

  RAJA::ReduceMax< parallelDeviceReduce, integer > foundNegativeTemp( 0 );
  RAJA::ReduceMax< parallelDeviceReduce, integer > foundNegativePres( 0 );
  RAJA::ReduceMax< parallelDeviceReduce, integer > foundInconsistentCompFrac( 0 );

  auto const avgCompFracView = avgCompFrac.toViewConst();

  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    wellElemPres[iwelem] = refPres + avgTotalMassDens * ( wellElemGravCoef[iwelem] - refWellElemGravCoef );
    wellElemTemp[iwelem] = avgTemp;

    real64 sumCompFracForCheck = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      wellElemCompFrac[iwelem][ic] = avgCompFracView[ic];
      sumCompFracForCheck += wellElemCompFrac[iwelem][ic];
    }

    if( wellElemPres[iwelem] <= 0 )
    {
      foundNegativePres.max( 1 );
    }
    if( wellElemTemp[iwelem] <= 0 )
    {
      foundNegativeTemp.max( 1 );
    }
    if( !isZero( sumCompFracForCheck - 1.0, constitutive::MultiFluidConstants::minForSpeciesPresence ) )
    {
      foundInconsistentCompFrac.max( 1 );
    }

  } );


  GEOS_THROW_IF( foundNegativePres.get() == 1,
                 wellControls.getDataContext() << "Invalid well initialization, negative pressure was found.",
                 InputError, wellControls.getDataContext() );
  GEOS_THROW_IF( foundNegativeTemp.get() == 1,
                 wellControls.getDataContext() << "Invalid well initialization, negative temperature was found.",
                 InputError, wellControls.getDataContext() );
  GEOS_THROW_IF( foundInconsistentCompFrac.get() == 1,
                 wellControls.getDataContext() << "Invalid well initialization, inconsistent component fractions were found.",
                 InputError, wellControls.getDataContext() );


}

/******************************** CompDensInitializationKernel ********************************/

void
CompDensInitializationKernel::
  launch( localIndex const subRegionSize,
          integer const numComponents,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & wellElemTotalDens,
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens )
{
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    for( integer ic = 0; ic < numComponents; ++ic )
    {
      wellElemCompDens[iwelem][ic] = wellElemCompFrac[iwelem][ic] * wellElemTotalDens[iwelem][0];
    }
  } );
}

/******************************** RateInitializationKernel ********************************/

void
RateInitializationKernel::
  launch( localIndex const subRegionSize,
          WellControls const & wellControls,
          real64 const & time,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens,
          arrayView1d< real64 > const & connRate )
{
  if( wellControls.isProducer() )
  {
    // Use use defined control type to set initial connection rates
    WellConstraintBase const *   constraint = wellControls.getCurrentConstraint();
    real64 const constraintVal = constraint->getConstraintValue( time );
    ConstraintTypeId const controlType = constraint->getControl();
    if( controlType == ConstraintTypeId::PHASEVOLRATE )
    {
      integer const targetPhaseIndex = wellControls.getConstraintPhaseIndex();

      forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {
        connRate[iwelem] = constraintVal * phaseDens[iwelem][0][targetPhaseIndex];
      } );
    }
    else if( controlType == ConstraintTypeId::TOTALVOLRATE )
    {
      forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {
        connRate[iwelem] = LvArray::math::max( 0.1 * constraintVal * totalDens[iwelem][0], -1e3 );
      } );
    }
    else if( controlType == ConstraintTypeId::MASSRATE )
    {
      forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {
        connRate[iwelem] = constraintVal;
      } );
    }
    else if( controlType == ConstraintTypeId::BHP )
    {
      // this assumes phase control present
      integer const targetPhaseIndex = wellControls.getConstraintPhaseIndex();
      std::vector< WellConstraintBase * >  const constraints = wellControls.getProdRateConstraints();
      // Use first rate constraint to set initial connection rates
      real64 const rateVal = constraints[0]->getConstraintValue( time );
      forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {

        connRate[iwelem] = LvArray::math::max( 0.1 * rateVal * phaseDens[iwelem][0][targetPhaseIndex], -1e3 );
      } );
    }
  }
  else
  {
    // Use use defined control type to set initial connection rates
    WellConstraintBase const *   constraint = wellControls.getCurrentConstraint();
    real64 const constraintVal = constraint->getConstraintValue( time );
    ConstraintTypeId const controlType = constraint->getControl();
    if( controlType == ConstraintTypeId::PHASEVOLRATE )
    {
      integer const targetPhaseIndex =   wellControls.getConstraintPhaseIndex();

      forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {
        connRate[iwelem] = LvArray::math::max( 0.1 * constraintVal * phaseDens[iwelem][0][targetPhaseIndex], 1e3 );
      } );
    }
    else if( controlType == ConstraintTypeId::TOTALVOLRATE )
    {
      forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {
        connRate[iwelem] =  constraintVal * totalDens[iwelem][0];
      } );
    }
    else if( controlType == ConstraintTypeId::MASSRATE )
    {
      forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {
        connRate[iwelem] = constraintVal;
      } );
    }
    else if( controlType == ConstraintTypeId::BHP )
    {
      std::vector< WellConstraintBase * >  const constraints = wellControls.getInjRateConstraints();
      // Use first rate constraint to set initial connection rates
      real64 const rateVal = constraints[0]->getConstraintValue( time );
      forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {
        connRate[iwelem] = LvArray::math::min( 0.1 * rateVal * totalDens[iwelem][0], 1e3 );
      } );
    }
  }
}

} // end namespace compositionalMultiphaseWellKernels

} // end namespace geos
