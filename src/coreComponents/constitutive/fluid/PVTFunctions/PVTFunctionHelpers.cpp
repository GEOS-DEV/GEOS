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
 * @file PVTFunctionHelpers.cpp
 */

#include "codingUtilities/StringUtilities.hpp"
#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "LvArray/src/sortedArrayManipulation.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

void
BlackOilTables::readTable( string const & fileName,
                           array1d< array1d< real64 > > & data,
                           integer minRowLength )
{
  std::ifstream is( fileName );
  GEOSX_ASSERT_GE_MSG( is.is_open(), true, "Could not open file: " + fileName );

  constexpr std::streamsize bufSize = 256;
  char buf[bufSize];

  // Read line-by-line until eof
  while( is.getline( buf, bufSize ) )
  {
    string str( buf );

    // Remove whitespace and end-of-line characters, if any
    stringutilities::trim( str );

    // Remove # and -- (Eclipse-style) comments
    stringutilities::removeStringAndFollowingContentFromLine( "#", str );
    stringutilities::removeStringAndFollowingContentFromLine( "--", str );

    // Skip empty or comment-only strings
    if( str.empty() )
    {
      continue;
    }

    // Add and read a new line entry
    data.emplace_back( 0 );
    stringutilities::fromStringTo( str, data.back() );

    // Remove line entry of no data read
    if( data.back().empty() )
    {
      data.pop_back();
    }
  }

  is.close();

  for( localIndex i = 0; i < data.size(); ++i )
  {
    GEOSX_ERROR_IF( data[i].size() < minRowLength,
                    "Too few entries in a row of table " + fileName
                    + ", minimum " + std::to_string( minRowLength ) + " required" );
  }
}

array1d< array1d< array1d< real64 > > >
BlackOilTables::readAllTables( localIndex const ipWater,
                               array1d< integer > const & phaseTypes,
                               path_array const & tableFileNames ) const
{
  localIndex const numPhases = phaseTypes.size();
  array1d< array1d< array1d< real64 > > > result( numPhases );

  for( localIndex i = 0; i < numPhases; ++i )
  {
    localIndex const rowMinLength = ( phaseTypes[i] == ipWater ) ? 4 : 3;
    readTable( tableFileNames[i], result[i], rowMinLength );
  }
  return result;
}

void
BlackOilTables::buildAllTables( localIndex const ipOil,
                                localIndex const ipWater,
                                localIndex const ipGas,
                                array1d< integer > const & phaseTypes,
                                path_array const & tableFiles )
{

  // check if both oil and gas are defined
  auto lower = LvArray::sortedArrayManipulation::find( phaseTypes.begin(), phaseTypes.size(), int(ipOil) );
  bool const containsOil = lower != phaseTypes.size();
  GEOSX_ERROR_IF( !containsOil, "The oil phase must be defined for all PVT models" );

  // reading data from files
  const array1d< array1d< array1d< real64 > > > phaseTables = readAllTables( ipWater, phaseTypes, tableFiles );

  // finalize table construction
  for( localIndex i = 0; i < phaseTypes.size(); ++i )
  {
    array1d< array1d< real64 > > const & phaseTable = phaseTables[i];

    // oil
    if( phaseTypes[i] == ipOil )
    {
      m_oilTable.resize( phaseTable.size());
      for( localIndex ii = 0; ii < phaseTable.size(); ++ii )
      {
        m_oilTable[ii].resize( phaseTable[ii].size());
        for( localIndex jj = 0; jj < phaseTable[ii].size(); ++jj )
        {
          m_oilTable[ii][jj] = phaseTable[ii][jj];
        }
      }
    }
    // gas
    else if( phaseTypes[i] == ipGas )
    {
      m_gasTable.resize( phaseTable.size());
      for( localIndex ii = 0; ii < phaseTable.size(); ++ii )
      {
        m_gasTable[ii].resize( phaseTable[ii].size());
        for( localIndex jj = 0; jj < phaseTable[ii].size(); ++jj )
        {
          m_gasTable[ii][jj] = phaseTable[ii][jj];
        }
      }
    }
    // water
    else if( phaseTypes[i] == ipWater )
    {
      m_waterTable.resize( phaseTable[0].size());
      for( localIndex ii = 0; ii < phaseTable[0].size(); ++ii )
      {
        m_waterTable[ii] = phaseTable[0][ii];
      }
    }
    else
    {
      GEOSX_ERROR( "Phase type not supported for Black Oil model" );
    }
  }
}

} // namespace PVTProps

} // namespace constitutive

} // namespace geosx
