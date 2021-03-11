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
 * @file PVTFunctionHelpers.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONHELPERS_HPP
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONHELPERS_HPP

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class PTTableCoordinates
{
public:
  PTTableCoordinates()
  { coords.resize( 2 ); }

  localIndex nPressures() const { return coords[0].size(); }
  localIndex nTemperatures() const { return coords[1].size(); }

  void appendPressure( const real64 & pres ) { coords[0].emplace_back( pres ); }
  void appendTemperature( const real64 & temp ) { coords[1].emplace_back( temp ); }

  real64 const & getPressure( localIndex i ) const { return coords[0][i]; }
  real64 const & getTemperature( localIndex i ) const { return coords[1][i]; }

  array1d< array1d< real64 > > const & getCoords() const { return coords; }

private:
  array1d< array1d< real64 > > coords;
};

struct PVTFunctionHelpers
{

  /**
   * @brief Look for the expectedNames in the names provided by the user
   * @param[in] inputNames phase or component names provided by the user
   * @param[in] expectedNames expected names that can be accepted by GEOSX
   * @return the index of the phase or component that was found
   */
  template< typename InputRange, typename ExpectedRange >
  static localIndex
  findName( InputRange const & inputNames,
            ExpectedRange const & expectedNames )
  {
    using std::begin;
    using std::end;
    auto it = std::find_first_of( begin( inputNames ), end( inputNames ), begin( expectedNames ), end( expectedNames ) );
    localIndex const id = std::distance( begin( inputNames ), it );
    return id;
  }

  static real64
  initializePropertyTableWithSalinity( string_array const & inputParameters,
                                       PTTableCoordinates & tableCoords )
  {
    real64 TStart = -1.0;
    real64 TEnd = -1.0;
    real64 dT = -1.0;
    real64 PStart = -1.0;
    real64 PEnd = -1.0;
    real64 dP = -1.0;
    real64 salinity = -1.0;

    GEOSX_ERROR_IF( inputParameters.size() != 9, "Invalid property input!" );

    try
    {
      PStart = stod( inputParameters[2] );
      PEnd = stod( inputParameters[3] );
      dP = stod( inputParameters[4] );

      TStart = stod( inputParameters[5] );
      TEnd = stod( inputParameters[6] );
      dT = stod( inputParameters[7] );

      salinity = stod( inputParameters[8] );
    }
    catch( const std::invalid_argument & e )
    {
      GEOSX_ERROR( "Invalid BrineCO2Density argument:" + string( e.what()) );
    }

    for( real64 P = PStart; P <= PEnd; P += dP )
    {
      tableCoords.appendPressure( P );
    }
    for( real64 T = TStart; T <= TEnd; T += dT )
    {
      tableCoords.appendTemperature( T );
    }
    return salinity;
  }

  static void
  initializePropertyTable( string_array const & inputParameters,
                           PTTableCoordinates & tableCoords )
  {
    real64 TStart = -1.0;
    real64 TEnd = -1.0;
    real64 dT = -1.0;
    real64 PStart = -1.0;
    real64 PEnd = -1.0;
    real64 dP = -1.0;

    GEOSX_ERROR_IF( inputParameters.size() != 8, "Invalid property input!" );

    try
    {
      PStart = stod( inputParameters[2] );
      PEnd = stod( inputParameters[3] );
      dP = stod( inputParameters[4] );

      TStart = stod( inputParameters[5] );
      TEnd = stod( inputParameters[6] );
      dT = stod( inputParameters[7] );
    }
    catch( const std::invalid_argument & e )
    {
      GEOSX_ERROR( "Invalid property argument:" + string( e.what()) );
    }

    for( real64 P = PStart; P <= PEnd; P += dP )
    {
      tableCoords.appendPressure( P );
    }
    for( real64 T = TStart; T <= TEnd; T += dT )
    {
      tableCoords.appendTemperature( T );
    }
  }

};

} // namespace PVTProps

} // namespace constitutive

} // namespace geosx


#endif
