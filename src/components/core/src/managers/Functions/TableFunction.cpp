/*
 * TableFunction.cpp
 *
 *  Created on: June 16, 2017
 *      Author: sherman
 */

#include "TableFunction.hpp"
#include "common/DataTypes.hpp"
//#include "mathpresso.h"
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
                              ManagedGroup * const parent ):
  FunctionBase( name, parent ),
  m_coordinates(),
  m_values(),
  m_dimensions(0),
  m_size(),
  m_indexIncrement(),
  m_corners(),
  m_numCorners(0)
{}

TableFunction::~TableFunction()
{
  // TODO Auto-generated destructor stub
}


void TableFunction::FillDocumentationNode()
{
  FunctionBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Table function");

  docNode->AllocateChildNode( keys::tableCoordinates,
                              keys::tableCoordinates,
                              -1,
                              "real64_array",
                              "real64_array",
                              "Table coordinates",
                              "Table coordinates inputs for 1D tables",
                              "",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::tableValues,
                              keys::tableValues,
                              -1,
                              "real64_array",
                              "real64_array",
                              "Table Values",
                              "Table Values for 1D tables",
                              "",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::coordinateFiles,
                              keys::coordinateFiles,
                              -1,
                              "string_array",
                              "string_array",
                              "List of coordinate file names",
                              "List of coordinate file names",
                              "",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::voxelFile,
                              keys::voxelFile,
                              -1,
                              "string",
                              "string",
                              "Voxel file name",
                              "Voxel file name",
                              "none",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::tableInterpolation,
                              keys::tableInterpolation,
                              -1,
                              "string",
                              "string",
                              "Interpolation method",
                              "Interpolation method",
                              "linear",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::valueType,
                              keys::valueType,
                              -1,
                              "string",
                              "string",
                              "Value Type",
                              "Value Type",
                              "real64",
                              "",
                              1,
                              1,
                              0 );
}

void TableFunction::BuildDataStructure( ManagedGroup * const domain )
{}

void TableFunction::InitializeFunction()
{
  // Read in data
  string_array const & coordinateFiles = getReference<string_array>(keys::coordinateFiles);
  if (coordinateFiles[0].empty())
  {
    // 1D Table
    m_dimensions = 1;
    real64_array const & coordinates = getReference<real64_array>(keys::tableCoordinates);
    real64_array const & tmpValues = getReference<real64_array>(keys::tableValues);
    localIndex tableSize = tmpValues.size();

    if (coordinates.size() == tableSize)
    {
      m_coordinates.push_back(coordinates);
      m_values.insert(std::end(m_values), std::begin(tmpValues), std::end(tmpValues));
      m_size.push_back(tableSize);
    }
  }
  else
  {
    m_dimensions = coordinateFiles.size();
    m_coordinates.resize(m_dimensions);

    // TODO: Read these files on rank 0, then broadcast
    view_rtype<string> voxelFile = getData<string>(keys::voxelFile);
    IOUtilities::parse_file( m_values, voxelFile, ',' );
    for (localIndex ii=0 ; ii<m_dimensions ; ++ii)
    {
      IOUtilities::parse_file( m_coordinates[ii], coordinateFiles[ii], ',' );
      m_size.push_back(m_coordinates[ii].size());
    }
  }

  // Setup index increment (assume data is in Fortran array order)
  localIndex increment = 1;
  for (localIndex ii=0 ; ii<m_dimensions ; ++ii)
  {
    m_indexIncrement.push_back(increment);
    increment *= m_size[ii];
  }

  // Error checking
  if (increment != m_values.size())
  {
    throw std::invalid_argument("Table dimensions do not match!");
  }

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

void TableFunction::Evaluate( dataRepository::ManagedGroup const * const group,
                              real64 const time,
                              lSet const & set,
                              real64_array & result ) const
{
  FunctionBase::EvaluateT<TableFunction>( group, time, set, result );
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
      bounds[ii][1] = multidimensionalArray::integer_conversion<localIndex>(std::distance(m_coordinates[ii].begin(), lower));
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

REGISTER_CATALOG_ENTRY( FunctionBase, TableFunction, std::string const &, ManagedGroup * const )

} /* namespace ANST */
