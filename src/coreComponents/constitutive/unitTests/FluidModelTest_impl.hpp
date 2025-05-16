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
 * @file FluidModelTest_impl.hpp
 */

#ifndef GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_IMPL_HPP_
#define GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_IMPL_HPP_

#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

#include "codingUtilities/UnitTestUtilities.hpp"

namespace geos
{
namespace testing
{

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::FluidModelTest():
  m_parent( "parent", m_node )
{
  createFunctionManager();
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
typename FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::FluidModel *
FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::getFluid( string const & name )
{
  return m_parent.getGroupPointer< FLUID_TYPE >( name );
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
template< typename LAMBDA >
typename FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::FluidModel *
FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::createFluid( string const & name, LAMBDA && function )
{
  m_parent.resize( 1 );
  auto & fluid = m_parent.registerGroup< FLUID_TYPE >( name );
  function( fluid );
  fluid.postInputInitializationRecursive();
  m_parent.initialize();
  m_parent.initializePostInitialConditions();
  return &fluid;
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::
testValuesAgainstPreviousImplementation( FluidModel * fluid,
                                         TestPoint const & testPoint,
                                         TestResult const & expectedValues,
                                         real64 const relTol,
                                         real64 const absTol )
{
  integer constexpr size = 1;
  m_parent.resize( size );
  fluid->allocateConstitutiveData( m_parent, 1 );

  string_array const & phaseNames = fluid->phaseNames();

  real64 const pressure = std::get< 0 >( testPoint );
  real64 const temperature = std::get< 1 >( testPoint );
  auto const composition = std::get< 2 >( testPoint );

  array2d< real64, compflow::LAYOUT_PHASE > compositionArray( size, numComp );
  for( integer ic = 0; ic < NUM_COMP; ++ic )
  {
    compositionArray( 0, ic ) = composition[ic];
  }

  FluidWrapper fluidWrapper = fluid->createKernelWrapper();

  auto const compositionView = compositionArray.toViewConst();

  forAll< parallelHostPolicy >( size, [fluidWrapper,
                                       pressure,
                                       temperature,
                                       compositionView]
                                GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
    {
      fluidWrapper.update( k, q, pressure, temperature, compositionView[k] );
    }
  } );

  // Bring everything back to host
  forAll< serialPolicy >( 1, [fluidWrapper] ( localIndex const )
  {
    /* Do nothing */
  } );

  auto compareValues = [&]( string const & name, real64 const calculated, real64 const expected )
  {
    testing::checkRelativeError( calculated,
                                 expected,
                                 relTol,
                                 absTol,
                                 GEOS_FMT(
                                   "\n{} failed.\n"
                                   "Pressure: {}, Temperature: {} Composition: {}.\n"
                                   "Calculated: {}.\n"
                                   "Expacted: {}\n"
                                   "Difference: {}",
                                   name,
                                   pressure, temperature, toString( compositionView[0].toSliceConst() ),
                                   calculated, expected,
                                   expected-calculated ));
  };

  bool const isThermal = fluid->isThermal();

  for( integer phaseIndex = 0; phaseIndex < numPhase; phaseIndex++ )
  {
    real64 const phaseFraction = fluid->phaseFraction()( 0, 0, phaseIndex );
    real64 const expectedPhaseFraction = std::get< 0 >( expectedValues )[phaseIndex];
    compareValues( GEOS_FMT( "Phase fraction ({})", phaseNames[phaseIndex] ), phaseFraction, expectedPhaseFraction );

    real64 const phaseDensity = fluid->phaseDensity()( 0, 0, phaseIndex );
    real64 const expectedPhaseDensity = std::get< 1 >( expectedValues )[phaseIndex];
    compareValues( GEOS_FMT( "Phase density ({})", phaseNames[phaseIndex] ), phaseDensity, expectedPhaseDensity );

    real64 const phaseMassDensity = fluid->phaseMassDensity()( 0, 0, phaseIndex );
    real64 const expectedPhaseMassDensity = std::get< 2 >( expectedValues )[phaseIndex];
    compareValues( GEOS_FMT( "Phase mass density ({})", phaseNames[phaseIndex] ), phaseMassDensity, expectedPhaseMassDensity );

    real64 const phaseViscosity = fluid->phaseViscosity()( 0, 0, phaseIndex );
    real64 const expectedPhaseViscosity = std::get< 3 >( expectedValues )[phaseIndex];
    compareValues( GEOS_FMT( "Phase viscosity ({})", phaseNames[phaseIndex] ), phaseViscosity, expectedPhaseViscosity );

    if( isThermal )
    {
      real64 const phaseEnthalpy = fluid->phaseEnthalpy()( 0, 0, phaseIndex );
      real64 const expectedPhaseEnthalpy = std::get< 4 >( expectedValues )[phaseIndex];
      compareValues( GEOS_FMT( "Phase enthalpy ({})", phaseNames[phaseIndex] ), phaseEnthalpy, expectedPhaseEnthalpy );

      real64 const phaseInternalEnergy = fluid->phaseInternalEnergy()( 0, 0, phaseIndex );
      real64 const expectedPhaseInternalEnergy = std::get< 5 >( expectedValues )[phaseIndex];
      compareValues( GEOS_FMT( "Phase internal energy ({})", phaseNames[phaseIndex] ), phaseInternalEnergy, expectedPhaseInternalEnergy );
    }
  }
  real64 const totalDensity = fluid->totalDensity()( 0, 0 );
  real64 const expectedTotalDensity = std::get< 6 >( expectedValues );
  compareValues( "Total density", totalDensity, expectedTotalDensity );
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
template< integer NDIM, typename ... INDICES, integer USD1, integer USD2, integer USD3, typename >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::testDerivatives( string const propName,
                                                                         string const testPoint,
                                                                         ArrayView< real64 const, NDIM, USD1 > const & valueArray,
                                                                         ArrayView< real64 const, NDIM+1, USD2 > const & derivArray,
                                                                         ArraySlice< real64 const, 1, USD3 > const & displacements,
                                                                         real64 const valueScale,
                                                                         string_array const & dofNames,
                                                                         real64 const relTol,
                                                                         real64 const absTol,
                                                                         INDICES const ... indices )
{
  integer const numberOfDof = dofNames.size();
  real64 const invScale = 1.0 / valueScale;

  for( integer idof = 0; idof < numberOfDof; idof++ )
  {
    real64 const analyticalDerivative = derivArray( 0, 0, indices ..., idof ) * invScale;
    real64 numericalDerivative = 0.0;

    real64 centreValue = valueArray( 0, 0, indices ... ) * invScale;
    real64 rightValue = valueArray( 2*idof+1, 0, indices ... ) * invScale;
    real64 leftValue = valueArray( 2*idof+2, 0, indices ... ) * invScale;

    real64 const dVr = displacements[2*idof+1];
    real64 const dVl = displacements[2*idof+2];

    // Calculate the left, central and right derivative.
    // The selected numerical derivative will be the one nearest the analytical derivative
    real64 minError = LvArray::NumericLimits< real64 >::max;
    for( real64 const derivative : { (rightValue - centreValue)/dVr,
                                     (rightValue - leftValue)/(dVr-dVl),
                                     (leftValue - centreValue)/dVl} )
    {
      real64 const error = LvArray::math::abs( analyticalDerivative - derivative );
      if( error < minError )
      {
        minError = error;
        numericalDerivative = derivative;
      }
    }

    testing::checkRelativeError( analyticalDerivative,
                                 numericalDerivative,
                                 relTol,
                                 absTol,
                                 GEOS_FMT(
                                   "\nd{}/d{} failed.\n"
                                   "Test Point: {}.\n"
                                   "Values: [{}, {}, {}].\n"
                                   "Analytical derivative: {}.\n"
                                   "Numerical derivative: {}",
                                   propName, dofNames[idof],
                                   testPoint,
                                   valueScale*leftValue, valueScale*centreValue, valueScale*rightValue,
                                   valueScale*analyticalDerivative, valueScale*numericalDerivative ));
  }
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::testNumericalDerivatives( FluidModel * fluid,
                                                                                  TestPoint const & data,
                                                                                  real64 const perturbationLevel,
                                                                                  real64 const relTol,
                                                                                  real64 const absTol )
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  // Number of execution points actual value plus left and right pertubations for each variable
  integer constexpr size = 2*numDof + 1;
  m_parent.resize( size );
  fluid->allocateConstitutiveData( m_parent, 1 );

  bool const isThermal = fluid->isThermal();

  array1d< real64 > pressureArray( size );
  array1d< real64 > temperatureArray( size );
  array2d< real64, compflow::LAYOUT_PHASE > compositionArray( size, numComp );
  array1d< real64 > deltaArray( size );

  // Name the degrees of freedom
  string_array dofNames( numDof );
  string_array const & phaseNames = fluid->phaseNames();
  string_array const & componentNames = fluid->componentNames();

  auto const & [pressure, temperature, composition] = data;
  for( integer i = 0; i < size; ++i )
  {
    pressureArray[i] = pressure;
    temperatureArray[i] = temperature;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      compositionArray( i, ic ) = composition[ic];
    }
  }

  string const testValues = GEOS_FMT( "Pressure: {}, Temperature: {}, Composition: {}",
                                      pressureArray[0],
                                      temperatureArray[0],
                                      toString( compositionArray[0].toSliceConst() ) );

  real64 const dP = perturbationLevel * (pressure + perturbationLevel);
  deltaArray[2*Deriv::dP+1] = dP;
  deltaArray[2*Deriv::dP+2] = -dP;
  pressureArray[2*Deriv::dP+1] += deltaArray[2*Deriv::dP+1];
  pressureArray[2*Deriv::dP+2] += deltaArray[2*Deriv::dP+2];
  dofNames[Deriv::dP] = "pressure";

  real64 const dT = perturbationLevel * (temperature + perturbationLevel);
  deltaArray[2*Deriv::dT+1] = dT;
  deltaArray[2*Deriv::dT+2] = -dT;
  temperatureArray[2*Deriv::dT+1] += deltaArray[2*Deriv::dT+1];
  temperatureArray[2*Deriv::dT+2] += deltaArray[2*Deriv::dT+2];
  dofNames[Deriv::dT] = "temperature";

  for( integer ic = 0; ic < numComp; ++ic )
  {
    integer const idof = Deriv::dC+ic;
    real64 const dz = LvArray::math::max( 1.0e-7, perturbationLevel * ( composition[ic] + perturbationLevel ) );
    deltaArray[2*idof+1] = dz;
    deltaArray[2*idof+2] = -dz;
    compositionArray( 2*idof+1, ic ) += deltaArray[2*idof+1];
    compositionArray( 2*idof+2, ic ) += deltaArray[2*idof+2];
    compositionArray( 2*idof+1, ic ) = LvArray::math::min( compositionArray( 2*idof+1, ic ), 1.0 );
    compositionArray( 2*idof+2, ic ) = LvArray::math::max( compositionArray( 2*idof+2, ic ), 0.0 );
    dofNames[idof] = GEOS_FMT( "z({})", componentNames[ic] );
  }

  auto fluidBase = fluid->deliverClone( fluid->getName()+"Copy", &fluid->getParent());
  FluidModel * fluidCopy = dynamicCast< FluidModel * >( fluidBase.get() );

  auto fluidWrapper = fluidCopy->createKernelWrapper();

  auto const pressureView = pressureArray.toViewConst();
  auto const temperatureView = temperatureArray.toViewConst();
  auto const compositionView = compositionArray.toViewConst();

  forAll< parallelHostPolicy >( size, [fluidWrapper,
                                       pressureView,
                                       temperatureView,
                                       compositionView]
                                GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
    {
      fluidWrapper.update( k, q, pressureView[k], temperatureView[k], compositionView[k] );
    }
  } );

  // Bring everything back to host
  forAll< serialPolicy >( 1, [fluidWrapper] ( localIndex const )
  {
    /* Do nothing */
  } );

  // Enthalpy and internal energy values can be quite large and prone to round-off errors when
  // performing finite difference derivatives. We will therefore scale the values using this
  // reference enthalpy value to bring them to a more manageable scale.
  real64 constexpr referenceEnthalpy = 5.0584e5;

  for( integer phaseIndex = 0; phaseIndex < numPhase; phaseIndex++ )
  {
    testDerivatives( GEOS_FMT( "Phase fraction ({})", phaseNames[phaseIndex] ),
                     testValues,
                     fluidCopy->phaseFraction().toViewConst(),
                     fluidCopy->dPhaseFraction().toViewConst(),
                     deltaArray.toSliceConst(),
                     1.0,
                     dofNames,
                     relTol,
                     absTol,
                     phaseIndex );
    testDerivatives( GEOS_FMT( "Phase density ({})", phaseNames[phaseIndex] ),
                     testValues,
                     fluidCopy->phaseDensity().toViewConst(),
                     fluidCopy->dPhaseDensity().toViewConst(),
                     deltaArray.toSliceConst(),
                     1.0,
                     dofNames,
                     relTol,
                     absTol,
                     phaseIndex );
    testDerivatives( GEOS_FMT( "Phase mass density ({})", phaseNames[phaseIndex] ),
                     testValues,
                     fluidCopy->phaseMassDensity().toViewConst(),
                     fluidCopy->dPhaseMassDensity().toViewConst(),
                     deltaArray.toSliceConst(),
                     1.0,
                     dofNames,
                     relTol,
                     absTol,
                     phaseIndex );
    testDerivatives( GEOS_FMT( "Phase viscosity ({})", phaseNames[phaseIndex] ),
                     testValues,
                     fluidCopy->phaseViscosity().toViewConst(),
                     fluidCopy->dPhaseViscosity().toViewConst(),
                     deltaArray.toSliceConst(),
                     1.0,
                     dofNames,
                     relTol,
                     absTol,
                     phaseIndex );
    if( isThermal )
    {
      testDerivatives( GEOS_FMT( "Phase enthalpy ({})", phaseNames[phaseIndex] ),
                       testValues,
                       fluidCopy->phaseEnthalpy().toViewConst(),
                       fluidCopy->dPhaseEnthalpy().toViewConst(),
                       deltaArray.toSliceConst(),
                       referenceEnthalpy,
                       dofNames,
                       relTol,
                       absTol,
                       phaseIndex );
      testDerivatives( GEOS_FMT( "Phase internal energy ({})", phaseNames[phaseIndex] ),
                       testValues,
                       fluidCopy->phaseInternalEnergy().toViewConst(),
                       fluidCopy->dPhaseInternalEnergy().toViewConst(),
                       deltaArray.toSliceConst(),
                       referenceEnthalpy,
                       dofNames,
                       relTol,
                       absTol,
                       phaseIndex );
    }
    for( integer compIndex = 0; compIndex < numComp; compIndex++ )
    {
      testDerivatives( GEOS_FMT( "Phase composition ({} {})", phaseNames[phaseIndex], componentNames[compIndex] ),
                       testValues,
                       fluidCopy->phaseCompFraction().toViewConst(),
                       fluidCopy->dPhaseCompFraction().toViewConst(),
                       deltaArray.toSliceConst(),
                       1.0,
                       dofNames,
                       relTol,
                       absTol,
                       phaseIndex,
                       compIndex );
    }
  }
  testDerivatives( "Total density",
                   testValues,
                   fluidCopy->totalDensity().toViewConst(),
                   fluidCopy->dTotalDensity().toViewConst(),
                   deltaArray.toSliceConst(),
                   1.0,
                   dofNames,
                   relTolerance,
                   absTolerance );
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::writeTableToFile( string const & fileName, char const * content )
{
  std::ofstream os( fileName );
  ASSERT_TRUE( os.is_open() );
  os << content;
  os.close();
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::removeFile( string const & fileName )
{
  int const ret = std::remove( fileName.c_str() );
  ASSERT_TRUE( ret == 0 );
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
TableFunction const *
FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::createTable( string const & tableName,
                                                                real64_array const & coordinates,
                                                                real64_array const & values )
{
  array1d< real64_array > coordinateWrapper;
  coordinateWrapper.emplace_back( coordinates );
  return createTable( tableName, coordinateWrapper, values );
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
TableFunction const *
FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::createTable( string const & tableName,
                                                                array1d< real64_array > const & coordinates,
                                                                real64_array const & values )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction * table = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
  table->setTableCoordinates( coordinates );
  table->setTableValues( values );
  table->reInitializeFunction();

  table->setInterpolationMethod( TableFunction::InterpolationType::Linear );

  return table;
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
template< typename ARRAY, typename LIST >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::fill( ARRAY & array, LIST const & data )
{
  for( auto const & v : data )
  {
    array.emplace_back( v );
  }
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
template< typename ARRAY, typename LIST >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::populate( ARRAY & array, LIST const & data )
{
  integer const n = array.size();
  auto it = std::begin( data );
  for( integer i = 0; i <= n; ++i, ++it )
  {
    array[i] = *it;
  }
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
template< integer USD >
string FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::toString( arraySlice1d< real64 const, USD > const & array )
{
  std::ostringstream os;
  os << array;
  return os.str();
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::populateLinearScale( array1d< real64 > & array,
                                                                             real64 const x0, real64 const x1, integer const n )
{
  real64 const dx = (x1 - x0)/n;
  for( integer i = 0; i <= n; ++i )
  {
    array.emplace_back( x0 + i*dx );
  }
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::populateLogScale( array1d< real64 > & array,
                                                                          real64 const x0, real64 const x1, integer const n )
{
  real64 const dx = LvArray::math::exp( LvArray::math::log( x1/x0 ) / n );
  real64 x = x0;
  for( integer i = 0; i <= n; ++i, x *= dx )
  {
    array.emplace_back( x );
  }
}

template< typename FLUID_TYPE, integer NUM_COMP, integer NUM_PHASE >
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::createFunctionManager()
{
  // TODO: FunctionManager is a singleton
  string const name = FunctionManager::catalogName();
  auto functionManager = std::make_unique< FunctionManager >( name, &m_parent );
  m_parent.registerGroup( name, std::move( functionManager ) );
}

} // namespace testing

} //namespace geos

#endif // GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_IMPL_HPP_
