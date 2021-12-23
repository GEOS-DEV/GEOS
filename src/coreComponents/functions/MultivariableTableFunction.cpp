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
 * @file MultivariableTableFunction.cpp
 */

#include "MultivariableTableFunction.hpp"

#include "common/DataTypes.hpp"
#include <algorithm>

namespace geosx
{

using namespace dataRepository;

MultivariableTableFunction::MultivariableTableFunction( const string & name,
                                                        Group * const parent ):
  FunctionBase( name, parent )
{}

template< typename T >
void MultivariableTableFunction::parseFile( string const & filename, array1d< T > & target )
{
  std::ifstream inputStream( filename.c_str() );
  GEOSX_THROW_IF( !inputStream, catalogName() << " " << getName() << ": could not read input file " << filename, InputError );

  // Read the file
  // TODO: Update this to handle large parallel jobs
  string lineString;
  while( std::getline( inputStream, lineString ) )
  {
    std::istringstream ss( lineString );
    while( ss.peek() == ',' || ss.peek() == ' ' )
    {
      ss.ignore();
    }
    T value;
    while( ss >> value )
    {
      target.emplace_back( value );
      while( ss.peek() == ',' || ss.peek() == ' ' )
      {
        ss.ignore();
      }
    }
  }

  inputStream.close();
}


void MultivariableTableFunction::setTableCoordinates( array1d< real64_array > const & coordinates )
{
  m_coordinates.resize( 0 );
  for( localIndex i = 0; i < coordinates.size(); ++i )
  {
    for( localIndex j = 1; j < coordinates[i].size(); ++j )
    {
      GEOSX_THROW_IF( coordinates[i][j] - coordinates[i][j-1] <= 0,
                      catalogName() << " " << getName() << ": coordinates must be strictly increasing, but axis " << i << " is not",
                      InputError );
    }
    m_coordinates.appendArray( coordinates[i].begin(), coordinates[i].end() );
  }
  reInitializeFunction();
}

void MultivariableTableFunction::setTableValues( real64_array values )
{
  m_values = std::move( values );
  reInitializeFunction();
}

void MultivariableTableFunction::initializeFunction()
{
  // ND Table
  // parseFile( m_voxelFile, m_values );
  // array1d< real64 > tmp;
  // for( localIndex ii = 0; ii < m_coordinateFiles.size(); ++ii )
  // {
  //   tmp.clear();
  //   parseFile( m_coordinateFiles[ii], tmp );
  //   m_coordinates.appendArray( tmp.begin(), tmp.end() );
  // }

  reInitializeFunction();
}

void MultivariableTableFunction::reInitializeFunction()
{
  // Setup index increment (assume data is in Fortran array order)
  localIndex increment = 1;
  for( localIndex ii = 0; ii < m_coordinates.size(); ++ii )
  {
    increment *= m_coordinates.sizeOfArray( ii );
  }
  if( m_coordinates.size() > 0 && !m_values.empty() ) // coordinates and values have been set
  {
    GEOSX_THROW_IF_NE_MSG( increment, m_values.size(),
                           catalogName() << " " << getName() << ": number of values does not match total number of table coordinates",
                           InputError );
  }
}

REGISTER_CATALOG_ENTRY( FunctionBase, MultivariableTableFunction, string const &, Group * const )

} /* namespace ANST */
