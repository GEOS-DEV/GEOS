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
 * @file HeatCapacityCoefficients.cpp
 */

#include "HeatCapacityCoefficients.hpp"
#include "ComponentProperties.hpp"
#include "EquationOfState.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "dataRepository/InputFlags.hpp"
#include "common/PhysicsConstants.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

HeatCapacityCoefficients::HeatCapacityCoefficients( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

std::unique_ptr< ModelParameters >
HeatCapacityCoefficients::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< HeatCapacityCoefficients >() != nullptr )
  {
    return parameters;
  }
  parameters = EquationOfState::create( std::move( parameters ) );
  return std::make_unique< HeatCapacityCoefficients >( std::move( parameters ) );
}

void HeatCapacityCoefficients::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::enthalpyReferenceTemperatureString(), &m_referenceTemperature ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_referenceTemperature ).
    setDescription( "The reference temperature for enthalpy calculation" );

  fluid->registerWrapper( viewKeyStruct::referenceEnthalpyString(), &m_referenceEnthalpy ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "The enthalpy of each component in each of the phases at the reference temperature." );

  fluid->registerWrapper( viewKeyStruct::componentHeatCapacityCoefficientsString(), &m_coefficients ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "The polynomial coefficients for the specific heat capacity of each component." );

  // Register extra wrappers to enable auto-cloning
  fluid->registerWrapper( "heatCapacityPhaseTypes", &m_phaseTypes )
    .setSizedFromParent( 0 )
    .setRestartFlags( dataRepository::RestartFlags::NO_WRITE );
}

void HeatCapacityCoefficients::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                            ComponentProperties const & componentProperties )
{
  integer const numPhases = fluid->numFluidPhases();
  integer const numComps = fluid->numFluidComponents();

  // If the reference enthalpies are given, then there must be as many as there are components
  if( m_referenceEnthalpy.empty() )
  {
    m_referenceEnthalpy.resize( numPhases, numComps );
    m_referenceEnthalpy.zero();
  }
  else
  {
    integer const dim0 = m_referenceEnthalpy.size( 0 );
    integer const dim1 = m_referenceEnthalpy.size( 1 );

    GEOS_THROW_IF_NE_MSG( dim0, numPhases,
                          GEOS_FMT( "{}: '{}' the first dimension must be equal to the number of phases {}",
                                    fluid->getFullName(),
                                    viewKeyStruct::referenceEnthalpyString(),
                                    numPhases ),
                          InputError );
    GEOS_THROW_IF_NE_MSG( dim1, numComps,
                          GEOS_FMT( "{}: '{}' the second dimension must be equal to the number of components {}",
                                    fluid->getFullName(),
                                    viewKeyStruct::referenceEnthalpyString(),
                                    numComps ),
                          InputError );
  }

  // Reference enthalpies must be ordered
  m_phaseTypes.resize( numPhases );
  string_array const & phaseNames = fluid->phaseNames();
  string_array const & componentNames = fluid->componentNames();

  integer ipL = -1, ipV = -1, ipA = -1;
  for( integer ip = 0; ip < numPhases; ip++ )
  {
    PhaseType const phaseType = getPhaseTypeFromName( phaseNames[ip] );
    m_phaseTypes[ip] = static_cast< integer >(phaseType);
    switch( phaseType )
    {
      case PhaseType::LIQUID:
        ipL = ip;
        break;
      case PhaseType::VAPOUR:
        ipV = ip;
        break;
      case PhaseType::AQUEOUS:
        ipA = ip;
        break;
    }
  }
  // Enthalpy must be Oil < Water < Gas
  stdVector< std::pair< integer, integer > > pairs;
  if( 0 <= ipL && 0 <= ipV )
  {
    pairs.emplace_back( ipL, ipV );
  }
  if( 0 <= ipL && 0 <= ipA )
  {
    pairs.emplace_back( ipL, ipA );
  }
  if( 0 <= ipA && 0 <= ipV )
  {
    pairs.emplace_back( ipA, ipV );
  }
  for( auto const & [ip1, ip2] : pairs )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const phase1Enthalpy = m_referenceEnthalpy( ip1, ic );
      real64 const phase2Enthalpy = m_referenceEnthalpy( ip2, ic );
      GEOS_THROW_IF_GT_MSG( phase1Enthalpy, phase2Enthalpy,
                            GEOS_FMT( "{}: '{}' for component {}, the {} reference enthalpy {} must not "
                                      "be greater than the {} reference enthalpy {}.",
                                      fluid->getFullName(),
                                      viewKeyStruct::referenceEnthalpyString(),
                                      componentNames[ic], phaseNames[ip1], phase1Enthalpy, phaseNames[ip2], phase2Enthalpy ),
                            InputError );
    }
  }

  {
    integer const dim0 = m_coefficients.size( 0 );
    integer const dim1 = m_coefficients.size( 1 );

    // First dimension must be equal to number of components
    GEOS_THROW_IF_NE_MSG( dim0, numComps,
                          GEOS_FMT( "{}: '{}' the first dimension must be equal to the number of components {}",
                                    fluid->getFullName(),
                                    viewKeyStruct::componentHeatCapacityCoefficientsString(),
                                    numComps ),
                          InputError );
    // Second dimension must be equal to 5
    GEOS_THROW_IF_NE_MSG( dim1, 5,
                          GEOS_FMT( "{}: '{}' the second dimension must be equal 5",
                                    fluid->getFullName(),
                                    viewKeyStruct::componentHeatCapacityCoefficientsString() ),
                          InputError );
  }

  // Reference temperature must not be negative
  GEOS_THROW_IF_LE_MSG( m_referenceTemperature, 0.0,
                        GEOS_FMT( "{}: '{}' the reference temperature {} must not be zero or negative.",
                                  fluid->getFullName(),
                                  viewKeyStruct::enthalpyReferenceTemperatureString(),
                                  m_referenceTemperature ),
                        InputError );

  // Determine temperature range
  constexpr real64 zeroC = constants::zeroDegreesCelsiusInKelvin;
  real64 minTemperature = LvArray::math::min( zeroC, m_referenceTemperature );
  real64 maxTemperature = LvArray::math::max( zeroC, m_referenceTemperature );
  auto const criticalTemperature = componentProperties.getComponentCriticalTemperature();
  for( integer ic = 0; ic < numComps; ++ic )
  {
    minTemperature = LvArray::math::min( minTemperature, criticalTemperature[ic] );
    maxTemperature = LvArray::math::max( maxTemperature, criticalTemperature[ic] );
  }

  // Extend interval by 10% in each direction
  real64 const dt = LvArray::math::max( maxTemperature - minTemperature, 100.0 );
  minTemperature -= 0.1 * dt;
  maxTemperature += 0.1 * dt;

  // Transform to reference temperature space
  minTemperature -= m_referenceTemperature;
  maxTemperature -= m_referenceTemperature;

  real64 negT = 0.0;
  real64 negHT = 0.0;

  for( integer ic = 0; ic < numComps; ++ic )
  {
    bool isPositive = isPolynomialPositive( m_coefficients[ic], minTemperature, maxTemperature, negT, negHT );
    GEOS_THROW_IF( !isPositive,
                   GEOS_FMT( "{}: '{}' coefficients for component {} ({}) give a zero or negative heat "
                             "capacity of {} at temperature difference {} corresponding to a temperature {}.",
                             fluid->getFullName(),
                             viewKeyStruct::componentHeatCapacityCoefficientsString(),
                             componentNames[ic], m_coefficients[ic], negHT, negT, negT + m_referenceTemperature ),
                   InputError );
  }
}

/*
 * @brief Evaluates whether a fourth-order polynomial remains non-negative across a specified interval.
 *
 * @details
 * This method verifies that the polynomial h(T) does not fall below a safety threshold
 * (epsilon) for all T in the range [T0, T1].
 * The algorithm combines an integration-like stepping mechanism with an adaptively computed
 * local Lipschitz constant to determine safe step sizes.
 * At any current temperature, the maximum possible rate of decrease (the local Lipschitz constant)
 * is calculated by finding the absolute maximum of the polynomial's derivative on the remaining
 * interval [currentT, T1]. Because the derivative is a cubic function, its maximums can only
 * exist at the boundaries of the interval or at the roots of the second derivative (which is a
 * quadratic equation).
 * Once this maximum absolute derivative is found, a safe forward step is calculated as
 * deltaT = (h(currentT) - epsilon) / localLipschitz. This step mathematically guarantees
 * that the function cannot drop below epsilon within the step.
 * If the function evaluates to less than epsilon at any point, the iteration halts and
 * updates the reference parameters with the failing coordinates.
 *
 * @param a Array of 5 coefficients representing the polynomial, ordered from T^0 to T^4.
 * @param T0 Start of the evaluation range for temperature T.
 * @param T1 End of the evaluation range for temperature T.
 * @param[out] T Reference output parameter that will store the failing temperature if a drop occurs.
 * @param[out] hT Reference output parameter that will store the failing polynomial value.
 * @return true if the polynomial remains strictly >= epsilon, false otherwise.
 */
bool HeatCapacityCoefficients::isPolynomialPositive( arraySlice1d< real64 const > const a,
                                                     real64 const T0, real64 const T1,
                                                     real64 & T, real64 & hT )
{
  constexpr real64 epsilon = MultiFluidConstants::epsilon;
  constexpr real64 minStep = 1.0e-13;

  real64 currentT = T0;
  real64 currentH;

  // Precalculate the coefficients of the second derivative h''(T)
  real64 const cA = 12.0 * a[4];
  real64 const cB = 6.0 * a[3];
  real64 const cC = 2.0 * a[2];

  // Precalculate the roots of the second derivative
  real64 roots[2]{};
  integer numRoots = 0;

  if( LvArray::math::abs( cA ) < epsilon )
  {
    if( !(LvArray::math::abs( cB ) < epsilon))
    {
      roots[numRoots++] = -cC / cB;
    }
  }
  else
  {
    real64 const discriminant = cB * cB - 4.0 * cA * cC;
    if( discriminant >= 0.0 )
    {
      real64 const sqrtDisc = LvArray::math::sqrt( discriminant );
      roots[numRoots++] = (-cB - sqrtDisc) / (2.0 * cA);
      roots[numRoots++] = (-cB + sqrtDisc) / (2.0 * cA);
    }
  }

  // Lambda function to compute the absolute value of the first derivative
  auto const getAbsDeriv = [&]( real64 const t ) -> real64 {
    real64 const deriv = a[1] + t * (2.0 * a[2] + t * (3.0 * a[3] + t * 4.0 * a[4]));
    return LvArray::math::abs( deriv );
  };

  while( true )
  {
    // Evaluate the polynomial using Horner's method for stability and performance
    currentH = a[0] + currentT * (a[1] + currentT * (a[2] + currentT * (a[3] + currentT * a[4])));

    // If the function falls below the designated epsilon threshold, fail immediately
    if( currentH < epsilon )
    {
      T = currentT;
      hT = currentH;
      return false;
    }

    // If we successfully reached or surpassed the end of the interval, return success
    if( currentT >= T1 )
    {
      T = currentT;
      hT = currentH;
      return true;
    }

    // Compute the local Lipschitz constant (max absolute derivative) over [currentT, T1]
    real64 const valCur = getAbsDeriv( currentT );
    real64 const valEnd = getAbsDeriv( T1 );

    real64 maxAbsDeriv = LvArray::math::max( valCur, valEnd );

    for( integer i = 0; i < numRoots; ++i )
    {
      if( roots[i] > currentT && roots[i] < T1 )
      {
        maxAbsDeriv = LvArray::math::max( maxAbsDeriv, getAbsDeriv( roots[i] ));
      }
    }

    // If the derivative is zero across the entire remainder, the function is flat and safe
    if( LvArray::math::abs( maxAbsDeriv ) < epsilon )
    {
      T = T1;
      hT = currentH;
      return true;
    }

    // Calculate the maximally safe local step size
    real64 const calculatedStep = (currentH - epsilon) / maxAbsDeriv;

    // Ensure forward progress even when the derivative is extremely large near epsilon
    real64 const actualStep = LvArray::math::max( calculatedStep, minStep );

    currentT += actualStep;

    // Ensure we do not overshoot the boundary on the final iteration
    currentT = LvArray::math::min( currentT, T1 );
  }
}

} // end namespace compositional
} // end namespace constitutive
} // end namespace geos
