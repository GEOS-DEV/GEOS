/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FenghourCO2Viscosity.cpp
 */

#include "constitutive/fluid/PVTFunctions/NewFenghourCO2ViscosityFunction.hpp"

#include "constitutive/fluid/PVTFunctions/NewSpanWagnerCO2DensityFunction.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "managers/GeosxState.hpp"

namespace geosx
{

namespace PVTProps
{

FenghourCO2Viscosity::FenghourCO2Viscosity( array1d< string > const & inputPara,
                                            array1d< string > const & componentNames,
                                            array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( inputPara[1],
                   componentNames,
                   componentMolarWeight )
{
  makeTable( inputPara );
}

void FenghourCO2Viscosity::makeTable( array1d< string > const & inputPara )
{
  array1d< array1d< real64 > > coordinates;
  coordinates.resize( 2 );

  real64 TStart = -1.0;
  real64 TEnd = -1.0;
  real64 dT = -1.0;
  real64 PStart = -1.0;
  real64 PEnd = -1.0;
  real64 dP = -1.0;

  GEOSX_ERROR_IF( inputPara.size() < 8, "Invalid FenghourCO2Viscosity input!" );

  try
  {
    PStart = stod( inputPara[2] );
    PEnd = stod( inputPara[3] );
    dP = stod( inputPara[4] );

    TStart = stod( inputPara[5] );
    TEnd = stod( inputPara[6] );
    dT = stod( inputPara[7] );
  }
  catch( const std::invalid_argument & e )
  {
    GEOSX_ERROR( "Invalid FenghourCO2Viscosity argument:" + string( e.what()) );
  }

  real64 P = PStart;
  while( P <= PEnd )
  {
    coordinates[0].emplace_back( P );
    P += dP;
  }

  real64 T = TStart;
  while( T <= TEnd )
  {
    coordinates[1].emplace_back( T );
    T += dT;
  }

  localIndex nP = coordinates[0].size();
  localIndex nT = coordinates[1].size();
  array1d< real64 > valuesDensity( nP * nT );
  array1d< real64 > valuesViscosity( nP * nT );
  SpanWagnerCO2Density::calculateCO2Density( coordinates, valuesDensity );
  calculateCO2Viscosity( coordinates, valuesDensity, valuesViscosity );

  FunctionManager & functionManager = getGlobalState().getFunctionManager();
  m_CO2ViscosityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", "CO2ViscosityTable" ) );
  m_CO2ViscosityTable->setTableCoordinates( coordinates );
  m_CO2ViscosityTable->setTableValues( valuesViscosity );
  m_CO2ViscosityTable->reInitializeFunction();
  m_CO2ViscosityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
}

void FenghourCO2Viscosity::fenghourCO2ViscosityFunction( real64 const & temperatureCent,
                                                         real64 const & density,
                                                         real64 & viscosity )
{
  constexpr real64 espar = 251.196;
  constexpr real64 esparInv = 1.0 / espar;
  // coefficients from Table (1) of Fenghour and Wakeham (1998)
  constexpr real64 aa[5] = { 0.235156, -0.491266, 5.211155e-2, 5.347906e-2, -1.537102e-2 };
  // coefficients from Table (3) of Fenghour and Wakeham (1998)
  constexpr real64 d11 =  0.4071119e-2;
  constexpr real64 d21 =  0.7198037e-4;
  constexpr real64 d64 =  0.2411697e-16;
  constexpr real64 d81 =  0.2971072e-22;
  constexpr real64 d82 = -0.1627888e-22;
  // we neglect critical viscosity
  constexpr real64 vcrit = 0.0;

  // temperature in Kelvin
  real64 const temperatureKelvin = temperatureCent + 273.15;
  // equation (5) of Fenghour and Wakeham (1998)
  real64 const Tred = temperatureKelvin * esparInv;
  real64 const x = log( Tred );
  // equation (4) of Fenghour and Wakeham (1998)
  real64 const lnGfun = aa[0] + x * (aa[1] + x * (aa[2] + x *(aa[3] + x * aa[4])));
  real64 const GfunInv = exp( -lnGfun );
  // equation (3) of Fenghour and Wakeham (1998)
  real64 const vlimit = 1.00697 * sqrt( temperatureKelvin ) * GfunInv;

  // equation (8) of Fenghour and Wakeham (1998)
  real64 const d2 = density * density;
  real64 const vxcess = density * (d11 + density * (d21 + d2*d2*(d64 / (Tred*Tred*Tred) + d2*(d81 + d82/Tred))));

  // equation (1) of Fenghour and Wakeham (1998)
  viscosity = 1e-6 * (vlimit + vxcess + vcrit);
}

void FenghourCO2Viscosity::calculateCO2Viscosity( array1d< array1d< real64 > > const & coordinates,
                                                  array1d< real64 > const & valuesDensity,
                                                  array1d< real64 > const & valuesViscosity )
{
  localIndex const numPressures = coordinates[0].size();
  localIndex const numTemperatures = coordinates[1].size();

  for( localIndex i = 0; i < numPressures; ++i )
  {
    for( localIndex j = 0; j < numTemperatures; ++j )
    {
      fenghourCO2ViscosityFunction( coordinates[1][j],
                                    valuesDensity[j*numPressures+i],
                                    valuesViscosity[j*numPressures+i] );
    }
  }
}

FenghourCO2Viscosity::KernelWrapper FenghourCO2Viscosity::createKernelWrapper()
{
  return KernelWrapper( m_componentNames,
                        m_componentMolarWeight,
                        m_CO2ViscosityTable );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, FenghourCO2Viscosity, array1d< string > const &, array1d< string > const &, array1d< real64 > const & )

} // end namespace PVTProps

} // end namespace geosx
