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
 * @file FlashParameters.cpp
 */

#include "FlashParameters.hpp"
#include "ComponentProperties.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

FlashParameters::FlashParameters( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) ),
  m_continuousParameters( NUM_FLOAT ),
  m_discreteParameters( NUM_INTEGER )
{
  m_continuousParameters[STABILITY_THRESHOLD] = -1.0e-8;
  m_continuousParameters[STABILITY_TOLERANCE] = 1.0e-8;
  m_continuousParameters[FLASH_TOLERANCE] = 1.0e-8;
  m_continuousParameters[NEGATIVE_FLASH_TOLERANCE] = -1.0e-1;

  m_discreteParameters[STABILITY_MAX_ITERATIONS] = 300;
  m_discreteParameters[FLASH_MAX_ITERATIONS] = 300;
}

std::unique_ptr< ModelParameters > FlashParameters::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< FlashParameters >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< FlashParameters >( std::move( parameters ) );
}

void FlashParameters::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::stabilityThresholdString(), &m_continuousParameters[STABILITY_THRESHOLD] ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_continuousParameters[STABILITY_THRESHOLD] ).
    setDescription( "The value of the tangent plane distance below which a mixture is determined to be unstable" );

  fluid->registerWrapper( viewKeyStruct::stabilityToleranceString(), &m_continuousParameters[STABILITY_TOLERANCE] ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_continuousParameters[STABILITY_TOLERANCE] ).
    setDescription( "The tolerance to use to determine convergence to a stationary point in the stability test" );

  fluid->registerWrapper( viewKeyStruct::stabilityMaxIterationsString(), &m_discreteParameters[STABILITY_MAX_ITERATIONS] ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_discreteParameters[STABILITY_MAX_ITERATIONS] ).
    setDescription( "The maximum number of successive substitution iterations used during the stability test. "
                    "A zero or negative value skips the stability test" );

  fluid->registerWrapper( viewKeyStruct::flashMaxIterationsString(), &m_discreteParameters[FLASH_MAX_ITERATIONS] ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_discreteParameters[FLASH_MAX_ITERATIONS] ).
    setDescription( "The maximum number of iterations used during the flash" );

  fluid->registerWrapper( viewKeyStruct::flashToleranceString(), &m_continuousParameters[FLASH_TOLERANCE] ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_continuousParameters[FLASH_TOLERANCE] ).
    setDescription( "The tolerance to use to determine convergece of the flash" );
}

void FlashParameters::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                   ComponentProperties const & componentProperties )
{
  GEOS_UNUSED_VAR( fluid );
  GEOS_UNUSED_VAR( componentProperties );

  auto const checkLowerBound = [&]( auto const & value, auto const & bound, string const & attribute )
  {
    GEOS_THROW_IF_LT_MSG( value, bound,
                          GEOS_FMT( "invalid number of value in attribute '{}'. Should be greater than {}",
                                    bound, attribute ),
                          InputError, fluid->getDataContext() );
  };

  auto const checkUpperBound = [&]( auto const & value, auto const & bound, string const & attribute )
  {
    GEOS_THROW_IF_GT_MSG( value, bound,
                          GEOS_FMT( "invalid number of value in attribute '{}'. Should be less than {}",
                                    bound, attribute ),
                          InputError, fluid->getDataContext());
  };

  real64 constexpr epsilon = MultiFluidConstants::epsilon;

  // Stability threshold should be negative in the range [-10^{-2}, -10^{-14}]
  real64 const stabilityThreshold = m_continuousParameters[STABILITY_THRESHOLD];
  checkUpperBound( stabilityThreshold, -1.0e-14, viewKeyStruct::stabilityThresholdString() );
  checkLowerBound( stabilityThreshold, -1.0e-2, viewKeyStruct::stabilityThresholdString() );

  // Stability tolerance should be positive
  real64 const stabilityTolerance = m_continuousParameters[STABILITY_TOLERANCE];
  checkLowerBound( stabilityTolerance, epsilon, viewKeyStruct::stabilityThresholdString() );

  // Flash tolerance should be positive
  real64 const flashTolerance = m_continuousParameters[FLASH_TOLERANCE];
  checkLowerBound( flashTolerance, epsilon, viewKeyStruct::flashToleranceString() );

  // Max number of flash iterations should be at least 1
  integer const flashMaxIterations = m_discreteParameters[FLASH_MAX_ITERATIONS];
  checkLowerBound( flashMaxIterations, 1, viewKeyStruct::flashMaxIterationsString() );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos
