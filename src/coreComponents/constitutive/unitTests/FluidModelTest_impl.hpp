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

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

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
  return m_parent.getGroupPointer< typename FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::FluidModel >( name );
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
void FluidModelTest< FLUID_TYPE, NUM_COMP, NUM_PHASE >::testNumericalDerivatives( FluidModel * fluid,
                                                                                  TestPoint const & data )
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  integer const size = 2*numDof + 1;
  m_parent.resize( size );
  fluid->allocateConstitutiveData( m_parent, 1 );
  std::cout << "Data Size " << fluid->phaseMassDensity().size() << std::endl;

  array1d< real64 > pressureArray( size );
  array1d< real64 > temperatureArray( size );
  array2d< real64 > compositionArray( size, numComp );
  array1d< real64 > deltaArray( size );

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

  real64 constexpr dP = 1.0e2;
  deltaArray[2*Deriv::dP+1] = dP;
  deltaArray[2*Deriv::dP+2] = -dP;
  pressureArray[2*Deriv::dP+1] += deltaArray[2*Deriv::dP+1];
  pressureArray[2*Deriv::dP+2] += deltaArray[2*Deriv::dP+2];

  real64 constexpr dT = 1.0e-3;
  deltaArray[2*Deriv::dT+1] = dT;
  deltaArray[2*Deriv::dT+2] = -dT;
  temperatureArray[2*Deriv::dT+1] += deltaArray[2*Deriv::dT+1];
  temperatureArray[2*Deriv::dT+2] += deltaArray[2*Deriv::dT+2];

  real64 constexpr dz = 1.0e-7;
  for( integer ic = 0; ic < numComp; ++ic )
  {
    integer const idof = Deriv::dC+ic;
    deltaArray[2*idof+1] = dz;
    deltaArray[2*idof+2] = -dz;
    compositionArray( 2*idof+1, ic ) += deltaArray[2*idof+1];
    compositionArray( 2*idof+2, ic ) += deltaArray[2*idof+2];
    compositionArray( 2*idof+1, ic ) = LvArray::math::min( compositionArray( 2*idof+1, ic ), 1.0 );
    compositionArray( 2*idof+2, ic ) = LvArray::math::max( compositionArray( 2*idof+2, ic ), 0.0 );
  }

  FluidWrapper fluidWrapper = fluid->createKernelWrapper();

  forAll< serialPolicy >( size, [fluidWrapper,
                                 pres=pressureArray.toViewConst(),
                                 temp=temperatureArray.toViewConst(),
                                 comp=compositionArray.toViewConst()]
                          GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
    {
      fluidWrapper.update( k, q, pres[k], temp[k], comp[k] );
    }
  } );

  auto const & phaseFraction = fluid->phaseFraction();
  std::cout << "FRACTION: " << phaseFraction( 0, 0, 0 ) << " " << phaseFraction( 0, 0, 1 ) << std::endl;
  auto const & phaseCompFraction = fluid->phaseCompFraction();
  std::cout << "XMF: " << phaseCompFraction[0][0][0] << std::endl;
  std::cout << "YMF: " << phaseCompFraction[0][0][1] << std::endl;
  auto const & phaseDensity = fluid->phaseMassDensity();
  auto const & dPhaseDensity = fluid->dPhaseMassDensity();
  for( integer ip = 0; ip < numPhase; ip++ )
  {
    real64 const centreValue = phaseDensity( 0, 0, ip );
    for( integer idof = 0; idof < numDof; idof++ )
    {
      real64 const analyticDerivative = dPhaseDensity( 0, 0, ip, idof );
      real64 const dv = 1.0 / (analyticDerivative + 1.0e-16);
      real64 const rightValue = phaseDensity( 2*idof+1, 0, ip );
      real64 const leftValue = phaseDensity( 2*idof+2, 0, ip );
      real64 const dVr = deltaArray[2*idof+1];
      real64 const dVl = deltaArray[2*idof+2];
      real64 const rightDerivative = (rightValue - centreValue)/dVr;
      real64 const leftDerivative = (leftValue - centreValue)/dVl;
      real64 const centreDerivative = (rightValue - leftValue)/(dVr-dVl);
      std::cout << "A: " << analyticDerivative << "\n| "
                << std::scientific << std::setprecision( 6 ) << std::setw( 14 ) << leftValue  << " | "
                << std::scientific << std::setprecision( 6 ) << std::setw( 14 ) << centreValue  << " | "
                << std::scientific << std::setprecision( 6 ) << std::setw( 14 ) << rightValue  << "\n| "
                << std::scientific << std::setprecision( 6 ) << std::setw( 14 ) << leftDerivative  << " | "
                << std::scientific << std::setprecision( 6 ) << std::setw( 14 ) << centreDerivative  << " | "
                << std::scientific << std::setprecision( 6 ) << std::setw( 14 ) << rightDerivative  << "\n| "
                << std::scientific << std::setprecision( 6 ) << std::setw( 14 ) << LvArray::math::abs( dv*(leftDerivative-analyticDerivative))  << " | "
                << std::scientific << std::setprecision( 6 ) << std::setw( 14 ) << LvArray::math::abs( dv*(centreDerivative-analyticDerivative))  << " | "
                << std::scientific << std::setprecision( 6 ) << std::setw( 14 ) << LvArray::math::abs( dv*(rightDerivative-analyticDerivative)) << "\n"
                << "---------------------------------------------------------"
                << std::endl;
    }
    std::cout << "===========================================================\n";
  }
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
  string const name = FunctionManager::catalogName();
  auto functionManager = std::make_unique< FunctionManager >( name, &m_parent );
  m_parent.registerGroup( name, std::move( functionManager ) );
}

} // namespace testing

} //namespace geos

#endif // GEOS_CORECOMPONENTS_CONSTITUTIVE_UNITTESTS_FLUIDMODELTEST_IMPL_HPP_
