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
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseKernels.hpp"

#include "common/DataTypes.hpp"
#include <algorithm>

namespace geosx
{

using namespace dataRepository;

MultivariableTableFunction::MultivariableTableFunction( const string & name,
                                                        Group * const parent ):
  FunctionBase( name, parent )
{
  registerWrapper( viewKeyStruct::coordinatesString(), &m_tableCoordinates1D ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coordinates inputs for 1D tables" );

  registerWrapper( viewKeyStruct::valuesString(), &m_values ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Values for 1D tables" );

  registerWrapper( viewKeyStruct::coordinateFilesString(), &m_coordinateFiles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of coordinate file names for ND Table" );

  registerWrapper( viewKeyStruct::voxelFileString(), &m_voxelFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Voxel file name for ND Table" );
}

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
  // Read in data
  if( m_coordinates.size() > 0 )
  {
    // This function appears to be already initialized
    // Apparently, this can be called multiple times during unit tests?
  }
  else if( m_coordinateFiles.empty() )
  {
    // 1D Table
    m_coordinates.appendArray( m_tableCoordinates1D.begin(), m_tableCoordinates1D.end() );
    GEOSX_THROW_IF_NE_MSG( m_tableCoordinates1D.size(), m_values.size(),
                           catalogName() << " " << getName() << ": 1D table function coordinates and values must have the same length",
                           InputError );
  }
  else
  {
    // ND Table
    parseFile( m_voxelFile, m_values );
    array1d< real64 > tmp;
    for( localIndex ii = 0; ii < m_coordinateFiles.size(); ++ii )
    {
      tmp.clear();
      parseFile( m_coordinateFiles[ii], tmp );
      m_coordinates.appendArray( tmp.begin(), tmp.end() );
    }
  }

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

real64 MultivariableTableFunction::evaluate( real64 const * const input ) const
{
  return 0;
}

real64 MultivariableTableFunction::evaluate( arrayView1d< real64 const > const & input,
                                             arrayView1d< real64 > const & output ) const
{
  createAndLaunch< parallelDevicePolicy<> >( m_coordinates.toViewConst(), m_values.toView(), input, output, output );
  return 2;
}

template< typename POLICY >
void
MultivariableTableFunction::createAndLaunch( ArrayOfArraysView< real64 const > const & coordinates,
                                             arrayView1d< real64 > const & values,
                                             arrayView1d< real64 const > const & input,
                                             arrayView1d< real64 > const & output,
                                             arrayView1d< real64 > const & output_derivatives )
{
  if( 0 )
  {
    CompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( coordinates.size(), [&] ( auto ND )
    {
      integer constexpr NUM_DIMS = ND();
      MultivariableTableFunction::KernelWrapper< NUM_DIMS, 2 > kernel( coordinates, values, input, output, output_derivatives );
      MultivariableTableFunction::KernelWrapper< NUM_DIMS, 2 >::template launch< POLICY >( input.size() / NUM_DIMS, kernel );
    } );
  }
  else if( 1 )
  {
    CompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( coordinates.size(), [&] ( auto ND )
    {
      integer constexpr NUM_DIMS = ND();
      MultivariableTableFunction::KernelWrapper< NUM_DIMS, 3 > kernel( coordinates, values, input, output, output_derivatives );
      MultivariableTableFunction::KernelWrapper< NUM_DIMS, 3 >::template launch< POLICY >( input.size() / NUM_DIMS, kernel );
    } );
  }
}



REGISTER_CATALOG_ENTRY( FunctionBase, MultivariableTableFunction, string const &, Group * const )

} /* namespace ANST */
