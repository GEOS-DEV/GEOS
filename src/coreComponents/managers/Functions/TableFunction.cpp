/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file TableFunction.cpp
 */

#include "TableFunction.hpp"
#include "common/DataTypes.hpp"
#include "codingUtilities/IOUtilities.hpp"
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


TableFunction::TableFunction( const std::string& name,
                              Group * const parent ):
  FunctionBase( name, parent ),
  m_coordinates(),
  m_values(),
  m_dimensions(0),
  m_size(),
  m_indexIncrement(),
  m_corners(),
  m_numCorners(0)
{
  registerWrapper<real64_array>(keys::tableCoordinates)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Table coordinates inputs for 1D tables");

  registerWrapper<real64_array>(keys::tableValues)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Table Values for 1D tables");

  registerWrapper<string_array>(keys::coordinateFiles)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("List of coordinate file names");

  registerWrapper<string>(keys::voxelFile)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Voxel file name");

  registerWrapper<string>(keys::tableInterpolation)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Interpolation method");
}

TableFunction::~TableFunction()
{}

void TableFunction::InitializeFunction()
{
  // Read in data
  string_array const & coordinateFiles = getReference<string_array>(keys::coordinateFiles);
  if (coordinateFiles.empty())
  {
    // 1D Table
    m_dimensions = 1;
    real64_array const & coordinates = getReference<real64_array>(keys::tableCoordinates);
    real64_array const & tmpValues = getReference<real64_array>(keys::tableValues);
    localIndex tableSize = tmpValues.size();

    if (coordinates.size() == tableSize)
    {
      m_coordinates.push_back(coordinates);
      for (localIndex ii=0; ii<tmpValues.size(); ii++)
      {
        m_values.push_back(tmpValues[ii]);
      }
      m_size.push_back(tableSize);
    }
  }
  else
  {
    m_dimensions = integer_conversion<localIndex>(coordinateFiles.size());
    m_coordinates.resize(m_dimensions);

    // TODO: Read these files on rank 0, then broadcast
    string const& voxelFile = getReference<string>(keys::voxelFile);
    IOUtilities::parse_file( m_values, voxelFile, ',' );
    for (localIndex ii=0 ; ii<m_dimensions ; ++ii)
    {
      IOUtilities::parse_file( m_coordinates[ii], coordinateFiles[ii], ',' );
      m_size.push_back(m_coordinates[ii].size());
    }
  }

  reInitializeFunction();
}

void TableFunction::reInitializeFunction()
{

  // Setup index increment (assume data is in Fortran array order)
  localIndex increment = 1;
  m_indexIncrement.resize(m_dimensions);
  for (localIndex ii=0 ; ii<m_dimensions ; ++ii)
  {
    m_size[ii] = m_coordinates[ii].size();
    m_indexIncrement[ii] = increment;
    increment *= m_size[ii];
  }

  // Error checking
  GEOS_ERROR_IF( increment != m_values.size(), "Table dimensions do not match!");

  // Build a quick map to help with linear interpolation
  m_numCorners = static_cast<localIndex>(pow(2, m_dimensions));
  for (localIndex ii=0 ; ii<m_numCorners ; ++ii)
  {
    for (localIndex jj=0 ; jj<m_dimensions ; ++jj)
    {
      m_corners[jj][ii] = int(ii / pow(2, jj)) % 2;
    }
  }
}


real64 TableFunction::Evaluate( real64 const * const input ) const
{
  localIndex bounds[m_maxDimensions][2];
  real64 weights[m_maxDimensions][2];

  // Determine position, weights
  for (localIndex ii=0 ; ii<m_dimensions ; ++ii)
  {
    if (input[ii] <= m_coordinates[ii][0])
    {
      // Coordinate is to the left of this axis
      bounds[ii][0] = 0;
      bounds[ii][1] = 0;
      weights[ii][0] = 0;
      weights[ii][1] = 1;
    }
    else if (input[ii] >= m_coordinates[ii][m_size[ii] - 1])
    {
      // Coordinate is to the right of this axis
      bounds[ii][0] = m_size[ii] - 1;
      bounds[ii][1] = bounds[ii][0];
      weights[ii][0] = 1;
      weights[ii][1] = 0;
    }
    else
    {
      // Find the coordinate index
      ///TODO make this fast
      // Note: lower_bound uses a binary search...  If we assume coordinates are
      // evenly spaced, we can speed things up considerably
      auto lower = std::lower_bound(m_coordinates[ii].begin(), m_coordinates[ii].end(), input[ii]);
      bounds[ii][1] = integer_conversion<localIndex>(std::distance(m_coordinates[ii].begin(), lower));
      bounds[ii][0] = bounds[ii][1] - 1;

      real64 dx = m_coordinates[ii][bounds[ii][1]] - m_coordinates[ii][bounds[ii][0]];
      weights[ii][0] = 1.0 - (input[ii] - m_coordinates[ii][bounds[ii][0]]) / dx;
      weights[ii][1] = 1.0 - weights[ii][0];
    }
  }

  // Linear interpolation
  real64 weightedValue = 0.0;
  for (localIndex ii=0 ; ii<m_numCorners ; ++ii)
  {
    // Find array index
    localIndex tableIndex = 0;
    for (localIndex jj=0 ; jj<m_dimensions ; ++jj)
    {
      tableIndex += bounds[jj][m_corners[jj][ii]] * m_indexIncrement[jj];
    }

    // Determine weighted value
    real64 cornerValue = m_values[tableIndex];
    for (localIndex jj=0 ; jj<m_dimensions ; ++jj)
    {
      cornerValue *= weights[jj][m_corners[jj][ii]];
    }
    weightedValue += cornerValue;
  }

  return weightedValue;
}

REGISTER_CATALOG_ENTRY( FunctionBase, TableFunction, std::string const &, Group * const )

} /* namespace ANST */
