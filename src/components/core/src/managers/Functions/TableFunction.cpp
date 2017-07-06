/*
 * TableFunction.cpp
 *
 *  Created on: June 16, 2017
 *      Author: sherman
 */

#include "TableFunction.hpp"
#include "common/DataTypes.hpp"
#include <mathpresso/mathpresso.h>
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
                              ManagedGroup * const parent ) :
  FunctionBase( name, parent ),
  m_dimensions(0),
  m_coordinates(),
  m_values()
{}

TableFunction::~TableFunction()
{
  // TODO Auto-generated destructor stub
}


void TableFunction::FillDocumentationNode( dataRepository::ManagedGroup * const domain )
{
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
{
  RegisterDocumentationNodes();
}

void TableFunction::InitializeFunction()
{
  if (true)
  {
    // 1D Table
    m_dimensions = 1;
    view_rtype<real64_array> coordinates = getData<real64_array>(keys::tableCoordinates);
    view_rtype<real64_array> tmpValues = getData<real64_array>(keys::tableValues);
    localIndex tableSize = tmpValues.size();

    if (coordinates.size() == tableSize)
    {
      m_coordinates.push_back(coordinates);
      m_values.insert(std::end(m_values), std::begin(tmpValues), std::end(tmpValues));
      
    }
    else
    {
      throw std::invalid_argument("Table dimensions do not match!");
    }
  }
  else
  {
    // TODO: Port existing method from geos
  }
}

double TableFunction::Evaluate(double* input)
{
  // Simple Lookup
  // TODO: Use existing table functions
  localIndex lowerBounds[m_dimensions];
  for (localIndex ii=0; ii<m_dimensions; ++ii)
  {
    auto lower = std::lower_bound(m_coordinates[ii].begin(), m_coordinates[ii].end(), input[ii]);
    lowerBounds[ii] = std::distance(m_coordinates[ii].begin(), lower);

    if (lowerBounds[ii] >= m_coordinates[ii].size())
    {
      lowerBounds[ii] -= 1;
    }
  }

  return m_values[lowerBounds[0]];
}

REGISTER_CATALOG_ENTRY( FunctionBase, TableFunction, std::string const &, ManagedGroup * const )

} /* namespace ANST */
