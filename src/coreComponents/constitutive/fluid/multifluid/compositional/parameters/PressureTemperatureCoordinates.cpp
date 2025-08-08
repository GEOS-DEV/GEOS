/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PressureTemperatureCoordinates.hpp
 */

#include "PressureTemperatureCoordinates.hpp"

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

PressureTemperatureCoordinates::PressureTemperatureCoordinates( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

std::unique_ptr< ModelParameters > PressureTemperatureCoordinates::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< PressureTemperatureCoordinates >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< PressureTemperatureCoordinates >( std::move( parameters ) );
}

void PressureTemperatureCoordinates::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::pressureCoordinatesString(), &m_pressureCoordinates ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "List of pressure values for interpolation of function values." );

  fluid->registerWrapper( viewKeyStruct::temperatureCoordinatesString(), &m_temperatureCoordinates ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "List of temperature values for interpolation of function values." );
}

void PressureTemperatureCoordinates::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                                  ComponentProperties const & componentProperties )
{
  GEOS_UNUSED_VAR( componentProperties );

  if( !m_pressureCoordinates.empty())
  {
    // If coordinates are provided then there must be at least 2
    GEOS_THROW_IF_LT_MSG( m_pressureCoordinates.size(), 2,
                          GEOS_FMT( "{}: invalid number of pressure coordinates provided in {}. "
                                    "At least 2 values must be provided", fluid->getFullName(), viewKeyStruct::pressureCoordinatesString() ),
                          InputError );

    // Values must be strictly increasing
    GEOS_THROW_IF( !isStrictlyIncreasing( m_pressureCoordinates.toSliceConst()),
                   GEOS_FMT( "{}: invalid values of pressure coordinates provided in {}. "
                             "Values must be strictly increasing.", fluid->getFullName(), viewKeyStruct::pressureCoordinatesString() ),
                   InputError );
  }

  if( !m_temperatureCoordinates.empty())
  {
    // If coordinates are provided then there must be at least 2
    GEOS_THROW_IF_LT_MSG( m_temperatureCoordinates.size(), 2,
                          GEOS_FMT( "{}: invalid number of temperature coordinates provided in {}. "
                                    "At least 2 values must be provided", fluid->getFullName(), viewKeyStruct::temperatureCoordinatesString() ),
                          InputError );

    // Values must be strictly increasing
    GEOS_THROW_IF( !isStrictlyIncreasing( m_temperatureCoordinates.toSliceConst()),
                   GEOS_FMT( "{}: invalid values of temperature coordinates provided in {}. "
                             "Values must be strictly increasing.", fluid->getFullName(), viewKeyStruct::temperatureCoordinatesString() ),
                   InputError );
  }
}

bool PressureTemperatureCoordinates::isStrictlyIncreasing( arraySlice1d< real64 const > const & array )
{
  localIndex const size = array.size();
  GEOS_ASSERT( 1 < size );
  real64 constexpr epsilon = MultiFluidConstants::epsilon;
  for( localIndex i = 1; i < size; ++i )
  {
    if( array[i] - array[i-1] < epsilon )
    {
      return false;
    }
  }
  return true;
}

std::pair< integer, integer >
PressureTemperatureCoordinates::getVariableIndices( FunctionBase const * function )
{
  auto const & inputVarNames = function->getWrapper< string_array >( dataRepository::keys::inputVarNames ).reference();
  if( inputVarNames.size() == 2 && inputVarNames[0] == "temperature" )
  {
    return {1, 0};
  }
  return {0, 1};
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos
