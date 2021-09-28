/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PVTFunctionHelpers.hpp
 */

#include "common/DataTypes.hpp"

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONHELPERS_HPP
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONHELPERS_HPP

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

/**
 * @class BlackOilTables
 *
 * A class to read and store the PVT tables from file
 */
class BlackOilTables
{
public:

  /**
   * @brief Read a table from file and check its row length
   * @param[in] fileName the name of the file
   * @param[in] minRowLength the expected minimum row length (3 for water, 4 for gas and oil)
   * @param[out] data the data from the table
   */
  static void
  readTable( string const & fileName,
             integer minRowLength,
             array1d< array1d< real64 > > & data );

  /**
   * @brief Read all the raw data from file, and copy the raw data into the phase tables
   * @param[in] ipOil the index of the oil phase
   * @param[in] ipWater the index of the water phase
   * @param[in] ipGas the index of the gas phase
   * @param[in] phases the array of phase indices
   * @param[in] tableFileNames the names of the table files
   */
  void
  buildAllTables( localIndex const ipOil,
                  localIndex const ipWater,
                  localIndex const ipGas,
                  array1d< integer > const & phases,
                  path_array const & tableFileNames );

  /**
   * @brief Getter for the oil table
   * @return the oil table
   */
  array1d< array1d< real64 > > const &
  getOilTable() const { return m_oilTable; }

  /**
   * @brief Getter for the gas table
   * @return the gas table
   */
  array1d< array1d< real64 > > const &
  getGasTable() const { return m_gasTable; }

  /**
   * @brief Getter for the water table
   * @return the water table
   */
  array1d< real64 > const &
  getWaterTable() const { return m_waterTable; }

private:

  /**
   * @brief Call function readTable for each phase
   * @param[in] ipWater the index of the water phase
   * @param[in] phases the array of phase indices
   * @param[in] tableFileNames the names of the table files
   * @return the raw data for all phases
   */
  array1d< array1d< array1d< real64 > > >
  readAllTables( localIndex const ipWater,
                 array1d< integer > const & phases,
                 path_array const & tableFileNames ) const;

  /// The oil phase PVT table read from file
  array1d< array1d< real64 > > m_oilTable;

  /// The gas phase PVT table read from file
  array1d< array1d< real64 > > m_gasTable;

  /// The water phase PVT table read from file
  array1d< real64 > m_waterTable;

};


class PTTableCoordinates
{
public:
  PTTableCoordinates()
  { coords.resize( 2 ); }

  localIndex nPressures() const { return coords[coordType::PRES].size(); }
  localIndex nTemperatures() const { return coords[coordType::TEMP].size(); }

  void appendPressure( const real64 & pres ) { coords[coordType::PRES].emplace_back( pres ); }
  void appendTemperature( const real64 & temp ) { coords[coordType::TEMP].emplace_back( temp ); }

  real64 const & getPressure( localIndex i ) const { return coords[coordType::PRES][i]; }
  real64 const & getTemperature( localIndex i ) const { return coords[coordType::TEMP][i]; }

  array1d< array1d< real64 > > const & getCoords() const { return coords; }

private:

  struct coordType
  {
    static constexpr integer PRES = 0;
    static constexpr integer TEMP = 1;
  };

  array1d< array1d< real64 > > coords;
};

namespace PVTFunctionHelpers
{

/**
 * @brief Look for the expectedNames in the names provided by the user
 * @tparam InputRange type of the input range
 * @tparam ExpectedRange type of the expected range
 * @param[in] input phase or component names provided by the user
 * @param[in] expected expected names that can be accepted by GEOSX
 * @return the index of the phase or component that was found
 */
template< typename InputRange, typename ExpectedRange >
inline integer
findName( InputRange const & input,
          ExpectedRange const & expected,
          string const & attribute )
{
  using std::begin;
  using std::end;
  auto const it = std::find_first_of( begin( input ), end( input ), begin( expected ), end( expected ) );
  GEOSX_THROW_IF( it == end( input ),
                  GEOSX_FMT( "Name '{}' not found in `{}`.\nExpected one of: {}.\nInput provided: {}.",
                             *begin( expected ), attribute,
                             stringutilities::join( begin( expected ), end( expected ), ", " ),
                             stringutilities::join( begin( input ), end( input ) ) ),
                  InputError );
  return static_cast< integer >( std::distance( begin( input ), it ) );
}

/**
 * @brief Populate the coordinate table with pressure and temperature
 * @param[in] inputParameters the strings reads in the file provided by the user
 * @param[inout] tableCoords the (p,T) coordinates of the table
 */
inline void
initializePropertyTable( string_array const & inputParameters,
                         PTTableCoordinates & tableCoords )

{
  GEOSX_THROW_IF( inputParameters.size() < 8,
                  "Invalid property input!",
                  InputError );

  try
  {
    real64 const PStart = stod( inputParameters[2] );
    real64 const PEnd = stod( inputParameters[3] );
    real64 const dP = stod( inputParameters[4] );

    GEOSX_THROW_IF( PStart >= PEnd, "PStart must be strictly smaller than PEnd",
                    InputError );

    constexpr real64 T_K = 273.15;

    real64 const TStart = stod( inputParameters[5] ) - T_K;
    real64 const TEnd = stod( inputParameters[6] )- T_K;
    real64 const dT = stod( inputParameters[7] );

    GEOSX_THROW_IF( TStart < 10, "Temperature must be in Kelvin and must be larger than 283.15 K",
                    InputError );
    GEOSX_THROW_IF( TEnd > 350, "Temperature must be in Kelvin and must be smaller than 623.15 K",
                    InputError );
    GEOSX_THROW_IF( TStart >= TEnd, "TStart must be strictly smaller than TEnd",
                    InputError );

    for( real64 P = PStart; P <= PEnd; P += dP )
    {
      tableCoords.appendPressure( P );
    }
    for( real64 T = TStart; T <= TEnd; T += dT )
    {
      tableCoords.appendTemperature( T );
    }
  }
  catch( const std::invalid_argument & e )
  {
    GEOSX_THROW( "Invalid property argument:" + string( e.what() ), InputError );
  }
}

} // namespace PVTFunctionHelpers

} // namespace PVTProps

} // namespace constitutive

} // namespace geosx


#endif
