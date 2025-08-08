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
 * @file PhillipsBrineViscosity.cpp
 */

#include "PhillipsBrineViscosity.hpp"

#include "constitutive/fluid/multifluid/compositional/parameters/PressureTemperatureCoordinates.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/BrineSalinity.hpp"

#include "constitutive/fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp"

#include "common/Units.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

PhillipsBrineViscosity::PhillipsBrineViscosity( string const & name,
                                                ComponentProperties const & componentProperties,
                                                integer const phaseIndex,
                                                ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties )
{
  GEOS_UNUSED_VAR( phaseIndex );
  m_brineViscosityTable = makeViscosityTable( name,
                                              componentProperties,
                                              modelParameters );
}

PhillipsBrineViscosityUpdate::PhillipsBrineViscosityUpdate( TableFunction const & brineViscosityTable ):
  m_brineViscosityTable( brineViscosityTable.createKernelWrapper() )
{}

PhillipsBrineViscosity::KernelWrapper
PhillipsBrineViscosity::createKernelWrapper() const
{
  return KernelWrapper( *m_brineViscosityTable );
}

std::unique_ptr< ModelParameters >
PhillipsBrineViscosity::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  std::unique_ptr< ModelParameters > params = std::move( parameters );
  params = PressureTemperatureCoordinates::create( std::move( params ) );
  params = BrineSalinity::create( std::move( params ) );
  return params;
}

TableFunction const *
PhillipsBrineViscosity::makeViscosityTable( string const & name,
                                            ComponentProperties const & componentProperties,
                                            ModelParameters const & modelParameters )
{
  GEOS_UNUSED_VAR( name );
  GEOS_UNUSED_VAR( componentProperties );

  FunctionManager & functionManager = FunctionManager::getInstance();

  string const tableName = name + "_brine_viscosity_table";

  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    TableFunction * const viscosityTable = functionManager.getGroupPointer< TableFunction >( tableName );
    viscosityTable->initializeFunction();
    viscosityTable->setDimUnits( { units::Temperature } );
    viscosityTable->setValueUnits( units::Viscosity );
    return viscosityTable;
  }

  string const functionName = "compositional_phillips_viscosity";

  TableFunction const * pureWaterViscosityTable =
    PVTProps::PureWaterProperties::makeSaturationViscosityTable( functionName, functionManager );

  BrineSalinity const * brineSalinity = modelParameters.get< BrineSalinity >();
  real64 const salinity = brineSalinity->m_salinity;

  PressureTemperatureCoordinates const * coordinates = modelParameters.get< PressureTemperatureCoordinates >();
  auto const temperatureCoordinates = coordinates->m_temperatureCoordinates.toSliceConst();
  localIndex const nTemperatures = temperatureCoordinates.size();
  GEOS_THROW_IF_LE_MSG( nTemperatures, 0,
                        GEOS_FMT( "{}: Failed to determine temperature points for Phillips brine viscosity interpolation. "
                                  "Provide values for {}.",
                                  name,
                                  PressureTemperatureCoordinates::viewKeyStruct::temperatureCoordinatesString() ),
                        InputError );

  // These coefficients come from Phillips et al. (1981), equation (1), pages 5-6
  constexpr real64 a = 0.0816;
  constexpr real64 b = 0.0122;
  constexpr real64 c = 0.000128;
  constexpr real64 d = 0.000629;
  constexpr real64 k = -0.7;

  // Compute the model coefficients
  real64 const m_coef0 = 1.0 + salinity * (a + salinity * (b + salinity * c));
  real64 const m_coef1 =  d * (1.0 - exp( k * salinity ));

  array1d< real64 > viscosity( nTemperatures );
  for( localIndex j = 0; j < nTemperatures; ++j )
  {
    real64 const temperature = units::convertKToC( temperatureCoordinates[j] );
    real64 const pureWaterViscosity = pureWaterViscosityTable->evaluate( &temperature );
    real64 const viscMultiplier = m_coef0 + m_coef1 * temperature;
    viscosity[j] = pureWaterViscosity * viscMultiplier;
  }

  // Construct the actual table object
  array1d< real64_array > tableCoords( 1 );
  tableCoords[0].resize( nTemperatures );
  for( localIndex j = 0; j < nTemperatures; ++j )
  {
    tableCoords[0][j] = temperatureCoordinates[j];
  }
  TableFunction * const viscosityTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
  viscosityTable->setTableCoordinates( tableCoords, { units::Temperature } );
  viscosityTable->setTableValues( viscosity, units::Viscosity );
  viscosityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
  return viscosityTable;
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
