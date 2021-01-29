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
 * @file TableFunction.cpp
 */

#include "TableFunction.hpp"
#include "common/DataTypes.hpp"
#include <algorithm>

namespace geosx
{

namespace dataRepository
{
namespace keys
{
std::string const tableCoordinates = "coordinates";
std::string const tableValues = "values";
std::string const tableInterpolation = "interpolation";
std::string const coordinateFiles = "coordinateFiles";
std::string const voxelFile = "voxelFile";
std::string const valueType = "valueType";
}
}

using namespace dataRepository;



TableFunction::TableFunction( const std::string & name,
                              Group * const parent ):
  FunctionBase( name, parent ),
  m_tableCoordinates1D(),
  m_coordinateFiles(),
  m_voxelFile(),
  m_interpolationMethod( InterpolationType::Linear ),
  m_coordinates(),
  m_values(),
  m_dimensions( 0 ),
  m_size(),
  m_indexIncrement(),
  m_corners(),
  m_numCorners( 0 )
{
  registerWrapper( keys::tableCoordinates, &m_tableCoordinates1D )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Coordinates inputs for 1D tables" );

  registerWrapper( keys::tableValues, &m_values )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Values for 1D tables" );

  registerWrapper( keys::coordinateFiles, &m_coordinateFiles )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "List of coordinate file names for ND Table" );

  registerWrapper( keys::voxelFile, &m_voxelFile )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Voxel file name for ND Table" );

  registerWrapper( keys::tableInterpolation, &m_interpolationMethod )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Interpolation method. Valid options:\n* " + EnumStrings< InterpolationType >::concat( "\n* " ) )->
    setApplyDefaultValue( m_interpolationMethod );
}

TableFunction::~TableFunction()
{}


template< typename T >
void TableFunction::parseFile( array1d< T > & target, string const & filename, char delimiter )
{
  std::ifstream inputStream( filename.c_str());
  std::string lineString;
  T value;

  GEOSX_ERROR_IF( !inputStream, "Could not read input file: " << filename );

  // Read the file
  // TODO: Update this to handle large parallel jobs
  while( !inputStream.eof())
  {
    std::getline( inputStream, lineString );
    std::istringstream ss( lineString );

    while( ss.peek() == delimiter || ss.peek() == ' ' )
    {
      ss.ignore();
    }
    while( ss>>value )
    {
      target.emplace_back( value );
      while( ss.peek() == delimiter || ss.peek() == ' ' )
      {
        ss.ignore();
      }
    }
  }

  inputStream.close();
}

void TableFunction::setInterpolationMethod( InterpolationType const method )
{
  m_interpolationMethod = method;
  reInitializeFunction();
}

void TableFunction::setTableCoordinates( array1d< real64_array > coordinates )
{
  m_coordinates.resize( 0 );
  for( localIndex i = 0; i < coordinates.size(); ++i )
  {
    for( localIndex j = 1; j < coordinates[i].size(); ++j )
    {
      GEOSX_ERROR_IF( coordinates[i][j] - coordinates[i][j-1] <= 0,
                      "In the table, the coordinates must be strictly increasing, but axis " << i << "is not" );
    }
    m_coordinates.appendArray( coordinates[i].begin(), coordinates[i].end() );
  }
  reInitializeFunction();
}

void TableFunction::setTableValues( real64_array values )
{
  m_values = std::move( values );
  reInitializeFunction();
}

void TableFunction::initializeFunction()
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
    m_dimensions = 1;
    m_coordinates.appendArray( m_tableCoordinates1D.begin(), m_tableCoordinates1D.end() );
    m_size.emplace_back( m_tableCoordinates1D.size());

    // Check to make sure that the table dimensions match
    GEOSX_ERROR_IF( m_size[0] != m_values.size(), "1D Table function coordinates and values must have the same length." );
  }
  else
  {
    // ND Table
    m_dimensions = LvArray::integerConversion< localIndex >( m_coordinateFiles.size());

    parseFile( m_values, m_voxelFile, ',' );
    array1d< real64 > tmp;
    for( localIndex ii=0; ii<m_dimensions; ++ii )
    {
      tmp.resize( 0 );
      parseFile( tmp, m_coordinateFiles[ii], ',' );
      m_coordinates.appendArray( tmp.begin(), tmp.end() );
      m_size.emplace_back( m_coordinates[ii].size());
    }
  }

  reInitializeFunction();
}

void TableFunction::reInitializeFunction()
{
  m_dimensions = LvArray::integerConversion< localIndex >( m_coordinates.size());
  m_size.resize( m_dimensions );

  // Setup index increment (assume data is in Fortran array order)
  localIndex increment = 1;
  m_indexIncrement.resize( m_dimensions );
  for( localIndex ii=0; ii<m_dimensions; ++ii )
  {
    m_size[ii] = LvArray::integerConversion< localIndex >( m_coordinates[ii].size() );
    m_indexIncrement[ii] = increment;
    increment *= m_size[ii];
  }

  if( m_coordinates.size() > 0 && m_values.size() > 0 ) // coordinates and values have been set
  {
    GEOSX_ERROR_IF( increment != m_values.size(), "Table dimensions do not match!" );
  }

  // Build a quick map to help with linear interpolation
  m_numCorners = static_cast< localIndex >(pow( 2, m_dimensions ));
  for( localIndex ii=0; ii<m_numCorners; ++ii )
  {
    for( localIndex jj=0; jj<m_dimensions; ++jj )
    {
      m_corners[jj][ii] = int(ii / pow( 2, jj )) % 2;
    }
  }

  // Create the kernel wrapper
  m_kernelWrapper.create( m_interpolationMethod,
                          m_coordinates.toViewConst(),
                          m_values.toViewConst(),
                          m_dimensions,
                          m_size.toViewConst(),
                          m_indexIncrement.toViewConst(),
                          m_corners,
                          m_numCorners );

}

TableFunction::KernelWrapper TableFunction::createKernelWrapper() const
{
  return TableFunction::KernelWrapper( m_interpolationMethod,
                                       m_coordinates.toViewConst(),
                                       m_values.toViewConst(),
                                       m_dimensions,
                                       m_size.toViewConst(),
                                       m_indexIncrement.toViewConst(),
                                       m_corners,
                                       m_numCorners );
}

real64 TableFunction::evaluate( real64 const * const input ) const
{
  real64 scalarValue = 0;
  stackArray1d< real64, maxDimensions > derivativesArray( m_dimensions );

  // interpolate in table, return scalar value (derivatives are discarded)
  m_kernelWrapper.compute( input, scalarValue, derivativesArray );
  return scalarValue;
}

TableFunction::KernelWrapper::KernelWrapper( TableFunction::InterpolationType interpolationMethod,
                                             ArrayOfArraysView< real64 const > const & coordinates,
                                             arrayView1d< real64 const > const & values,
                                             localIndex dimensions,
                                             arrayView1d< localIndex const > const & size,
                                             arrayView1d< localIndex const > const & indexIncrement,
                                             localIndex const (&corners)[TableFunction::maxDimensions][16],
                                             localIndex const numCorners )
  :
  m_interpolationMethod( interpolationMethod ),
  m_coordinates( coordinates ),
  m_values( values ),
  m_dimensions( dimensions ),
  m_size( size ),
  m_indexIncrement( indexIncrement ),
  m_numCorners( numCorners )
{
  LvArray::tensorOps::copy< TableFunction::maxDimensions, 16 >( m_corners, corners );
}

TableFunction::KernelWrapper::KernelWrapper()
  :
  m_interpolationMethod(),
  m_coordinates(),
  m_values(),
  m_dimensions( 0 ),
  m_size(),
  m_indexIncrement(),
  m_numCorners( 0 )
{}

void TableFunction::KernelWrapper::create( TableFunction::InterpolationType interpolationMethod,
                                           ArrayOfArraysView< real64 const > const & coordinates,
                                           arrayView1d< real64 const > const & values,
                                           localIndex dimensions,
                                           arrayView1d< localIndex const > const & size,
                                           arrayView1d< localIndex const > const & indexIncrement,
                                           localIndex const (&corners)[TableFunction::maxDimensions][16],
                                           localIndex const numCorners )
{
  m_interpolationMethod = interpolationMethod;
  m_coordinates = coordinates;
  m_values = values;
  m_dimensions = dimensions;
  m_size = size;
  m_indexIncrement = indexIncrement;
  m_numCorners = numCorners;

  LvArray::tensorOps::copy< TableFunction::maxDimensions, 16 >( m_corners, corners );
}


REGISTER_CATALOG_ENTRY( FunctionBase, TableFunction, std::string const &, Group * const )

} /* namespace ANST */
