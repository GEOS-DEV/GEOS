/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VTKHybridPartitioning.cpp
 */

#include "VTKHybridPartitioning.hpp"

#include "common/TimingMacros.hpp"
#include "common/format/Format.hpp"
#include "mesh/generators/VTKSuperCellPartitioning.hpp"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace geos
{
namespace vtk
{
namespace
{

using Clock = std::chrono::steady_clock;
constexpr localIndex NUM_CONSTRAINTS = 3;
constexpr real64 INTEGER_WEIGHT_RESOLUTION = 1000.0;

struct CellModel
{
  unsigned char vtkType;
  localIndex numNodes;
  localIndex numFaces;
  localIndex quadratureProxy;
  localIndex stateScalarProxy;
  char const * name;
};

/**
 * Conservative first-order load proxies. They are balancing models, not
 * operation counts. Keeping this table centralized makes the assumptions easy
 * to calibrate and test without changing topology code.
 */
CellModel cellModel( unsigned char const vtkType )
{
  switch( vtkType )
  {
    case VTK_TETRA: return { VTK_TETRA, 4, 4, 1, 4, "tetrahedron" };
    case VTK_PYRAMID: return { VTK_PYRAMID, 5, 5, 5, 5, "pyramid" };
    case VTK_WEDGE: return { VTK_WEDGE, 6, 5, 6, 6, "wedge" };
    case VTK_HEXAHEDRON: return { VTK_HEXAHEDRON, 8, 6, 8, 8, "hexahedron" };
    default: return { vtkType, 0, 0, 0, 0, "unsupported" };
  }
}

template< typename FUNC >
void forEachCellFace( unsigned char const vtkType, FUNC && func )
{
  switch( vtkType )
  {
    case VTK_TETRA:
    {
      constexpr localIndex faces[4][4] = {
        { 0, 1, 2, -1 }, { 0, 3, 1, -1 }, { 1, 3, 2, -1 }, { 2, 3, 0, -1 }
      };
      for( localIndex face = 0; face < 4; ++face )
      {
        func( 3, faces[face] );
      }
      return;
    }
    case VTK_PYRAMID:
    {
      constexpr localIndex faces[5][4] = {
        { 0, 3, 2, 1 }, { 0, 1, 4, -1 }, { 1, 2, 4, -1 },
        { 2, 3, 4, -1 }, { 3, 0, 4, -1 }
      };
      func( 4, faces[0] );
      for( localIndex face = 1; face < 5; ++face )
      {
        func( 3, faces[face] );
      }
      return;
    }
    case VTK_WEDGE:
    {
      constexpr localIndex faces[5][4] = {
        { 0, 2, 1, -1 }, { 3, 4, 5, -1 }, { 0, 1, 4, 3 },
        { 1, 2, 5, 4 }, { 2, 0, 3, 5 }
      };
      func( 3, faces[0] );
      func( 3, faces[1] );
      for( localIndex face = 2; face < 5; ++face )
      {
        func( 4, faces[face] );
      }
      return;
    }
    case VTK_HEXAHEDRON:
    {
      constexpr localIndex faces[6][4] = {
        { 0, 4, 7, 3 }, { 1, 2, 6, 5 }, { 0, 1, 5, 4 },
        { 3, 7, 6, 2 }, { 0, 3, 2, 1 }, { 4, 5, 6, 7 }
      };
      for( localIndex face = 0; face < 6; ++face )
      {
        func( 4, faces[face] );
      }
      return;
    }
    default:
      GEOS_ERROR( GEOS_FMT( "Unsupported VTK cell type {} in hybrid face extraction",
                            static_cast< integer >( vtkType ) ) );
  }
}

struct MeshArrays
{
  vtkIdType const * offsets = nullptr;
  vtkIdType const * connectivity = nullptr;
  unsigned char const * cellTypes = nullptr;
  vtkIdType const * globalPointIds = nullptr;
  vtkIdType const * globalCellIds = nullptr;
};

MeshArrays getMeshArrays( vtkUnstructuredGrid const & mesh )
{
  vtkUnstructuredGrid & mutableMesh = const_cast< vtkUnstructuredGrid & >( mesh );
  // Conversion only changes the integral storage width, not mesh semantics.
  vtkCellArray * const cells = mutableMesh.GetCells();
  GEOS_ERROR_IF( cells == nullptr, "Root-local VTK mesh has no cell connectivity" );
  cells->ConvertToDefaultStorage();

  vtkIdTypeArray * const offsets = vtkIdTypeArray::FastDownCast( cells->GetOffsetsArray() );
  vtkIdTypeArray * const connectivity = vtkIdTypeArray::FastDownCast( cells->GetConnectivityArray() );
  vtkUnsignedCharArray * const cellTypes = mutableMesh.GetCellTypesArray();
  vtkIdTypeArray * const pointIds = vtkIdTypeArray::FastDownCast(
    mutableMesh.GetPointData()->GetGlobalIds() );
  vtkIdTypeArray * const cellIds = vtkIdTypeArray::FastDownCast(
    mutableMesh.GetCellData()->GetGlobalIds() );

  GEOS_ERROR_IF( offsets == nullptr || connectivity == nullptr,
                 "Hybrid partitioning requires vtkIdType cell connectivity" );
  GEOS_ERROR_IF( cellTypes == nullptr, "Hybrid partitioning requires a VTK cell-type array" );
  GEOS_ERROR_IF( pointIds == nullptr || cellIds == nullptr,
                 "Hybrid partitioning requires vtkIdType global point and cell IDs" );
  GEOS_ERROR_IF_NE_MSG( pointIds->GetNumberOfTuples(), mutableMesh.GetNumberOfPoints(),
                        "Global point ID tuple count does not match the VTK mesh" );
  GEOS_ERROR_IF_NE_MSG( cellIds->GetNumberOfTuples(), mutableMesh.GetNumberOfCells(),
                        "Global cell ID tuple count does not match the VTK mesh" );

  return { offsets->GetPointer( 0 ),
           connectivity->GetPointer( 0 ),
           cellTypes->GetPointer( 0 ),
           pointIds->GetPointer( 0 ),
           cellIds->GetPointer( 0 ) };
}

struct FaceRecord
{
  std::array< globalIndex, 4 > points = {};
  localIndex cell = -1;
  uint8_t size = 0;
};

bool faceKeyLess( FaceRecord const & lhs, FaceRecord const & rhs )
{
  if( lhs.size != rhs.size )
  {
    return lhs.size < rhs.size;
  }
  for( uint8_t i = 0; i < lhs.size; ++i )
  {
    if( lhs.points[i] != rhs.points[i] )
    {
      return lhs.points[i] < rhs.points[i];
    }
  }
  return false;
}

bool sameFaceKey( FaceRecord const & lhs, FaceRecord const & rhs )
{
  if( lhs.size != rhs.size )
  {
    return false;
  }
  for( uint8_t i = 0; i < lhs.size; ++i )
  {
    if( lhs.points[i] != rhs.points[i] )
    {
      return false;
    }
  }
  return true;
}

struct PointVertexRecord
{
  globalIndex point = -1;
  localIndex vertex = -1;
};

struct SuperCellRecord
{
  globalIndex superCell = -1;
  localIndex cell = -1;
};

struct EdgeContribution
{
  localIndex first = -1;
  localIndex second = -1;
  real64 weight = 0.0;
};

struct AggregatedEdge
{
  localIndex first = -1;
  localIndex second = -1;
  pmet_idx_t weight = 0;
};

uint64_t rankPairKey( int first, int second )
{
  if( first > second )
  {
    std::swap( first, second );
  }
  return (static_cast< uint64_t >( static_cast< uint32_t >( first ) ) << 32) |
         static_cast< uint32_t >( second );
}

void appendRankPair( std::unordered_map< uint64_t, int64_t > & counts,
                     int const first,
                     int const second,
                     int64_t const increment )
{
  if( first == second || increment == 0 )
  {
    return;
  }
  uint64_t const key = rankPairKey( first, second );
  int64_t & value = counts[key];
  value += increment;
  GEOS_ERROR_IF_LT_MSG( value, 0, "Hybrid rank-neighbor multiplicity became negative" );
  if( value == 0 )
  {
    counts.erase( key );
  }
}

stdVector< int > uniqueEntityParts( arraySlice1d< localIndex const > const vertices,
                                    arrayView1d< int64_t const > const parts,
                                    localIndex const movedVertex = -1,
                                    int const destination = -1 )
{
  stdVector< int > result;
  result.reserve( vertices.size() );
  for( localIndex const vertex : vertices )
  {
    int const part = vertex == movedVertex ? destination : LvArray::integerConversion< int >( parts[vertex] );
    result.push_back( part );
  }
  std::sort( result.begin(), result.end() );
  result.erase( std::unique( result.begin(), result.end() ), result.end() );
  return result;
}

void appendPartPairs( stdVector< int > const & parts,
                      std::unordered_map< uint64_t, int64_t > & counts,
                      int64_t const increment )
{
  for( std::size_t i = 0; i < parts.size(); ++i )
  {
    for( std::size_t j = i + 1; j < parts.size(); ++j )
    {
      appendRankPair( counts, parts[i], parts[j], increment );
    }
  }
}

void appendPartPairDeltas( stdVector< int > const & oldParts,
                           stdVector< int > const & newParts,
                           std::map< uint64_t, int64_t > & deltas )
{
  for( std::size_t i = 0; i < oldParts.size(); ++i )
  {
    for( std::size_t j = i + 1; j < oldParts.size(); ++j )
    {
      --deltas[rankPairKey( oldParts[i], oldParts[j] )];
    }
  }
  for( std::size_t i = 0; i < newParts.size(); ++i )
  {
    for( std::size_t j = i + 1; j < newParts.size(); ++j )
    {
      ++deltas[rankPairKey( newParts[i], newParts[j] )];
    }
  }
}

vtkDataArray * validateWeightField( vtkUnstructuredGrid const & mesh,
                                    string const & fieldName )
{
  if( fieldName.empty() )
  {
    return nullptr;
  }

  vtkUnstructuredGrid & mutableMesh = const_cast< vtkUnstructuredGrid & >( mesh );
  vtkDataArray * const field = mutableMesh.GetCellData()->GetArray( fieldName.c_str() );
  GEOS_ERROR_IF( field == nullptr,
                 GEOS_FMT( "Hybrid partition weight field '{}' was not found in VTK cell data", fieldName ) );
  GEOS_ERROR_IF_NE_MSG( field->GetNumberOfComponents(), 1,
                        GEOS_FMT( "Hybrid partition weight field '{}' must be scalar", fieldName ) );
  GEOS_ERROR_IF_NE_MSG( field->GetNumberOfTuples(), mutableMesh.GetNumberOfCells(),
                        GEOS_FMT( "Hybrid partition weight field '{}' tuple count does not match the 3D mesh",
                                  fieldName ) );
  return field;
}

void buildPartitionVertices( vtkUnstructuredGrid const & mesh,
                             MeshArrays const & arrays,
                             HybridPartitionTopology & topology )
{
  vtkUnstructuredGrid & mutableMesh = const_cast< vtkUnstructuredGrid & >( mesh );
  localIndex const numCells = LvArray::integerConversion< localIndex >( mutableMesh.GetNumberOfCells() );
  topology.cellToVertex.resize( numCells );

  vtkIdTypeArray * const superCellIds = vtkIdTypeArray::FastDownCast(
    mutableMesh.GetCellData()->GetArray( "SuperCellId" ) );

  if( superCellIds == nullptr )
  {
    array1d< localIndex > capacities( numCells );
    capacities.setValues< serialPolicy >( 1 );
    topology.vertexToCells.resizeFromCapacities< serialPolicy >( numCells, capacities.data() );
    topology.representativeGlobalCellIds.resize( numCells );
    for( localIndex cell = 0; cell < numCells; ++cell )
    {
      topology.cellToVertex[cell] = cell;
      topology.vertexToCells.emplaceBack( cell, cell );
      topology.representativeGlobalCellIds[cell] = arrays.globalCellIds[cell];
    }
    return;
  }

  GEOS_ERROR_IF_NE_MSG( superCellIds->GetNumberOfTuples(), mutableMesh.GetNumberOfCells(),
                        "SuperCellId tuple count does not match the 3D mesh" );

  stdVector< SuperCellRecord > records;
  records.reserve( numCells );
  for( localIndex cell = 0; cell < numCells; ++cell )
  {
    records.push_back( { superCellIds->GetValue( cell ), cell } );
  }
  std::sort( records.begin(), records.end(), []( SuperCellRecord const & lhs, SuperCellRecord const & rhs )
  {
    return lhs.superCell < rhs.superCell ||
           (lhs.superCell == rhs.superCell && lhs.cell < rhs.cell);
  } );

  stdVector< localIndex > capacities;
  for( std::size_t begin = 0; begin < records.size(); )
  {
    std::size_t end = begin + 1;
    while( end < records.size() && records[end].superCell == records[begin].superCell )
    {
      ++end;
    }
    capacities.push_back( LvArray::integerConversion< localIndex >( end - begin ) );
    begin = end;
  }

  localIndex const numVertices = LvArray::integerConversion< localIndex >( capacities.size() );
  topology.vertexToCells.resizeFromCapacities< serialPolicy >( numVertices, capacities.data() );
  topology.representativeGlobalCellIds.resize( numVertices );

  localIndex vertex = 0;
  for( std::size_t begin = 0; begin < records.size(); ++vertex )
  {
    std::size_t end = begin + 1;
    globalIndex representative = arrays.globalCellIds[records[begin].cell];
    while( end < records.size() && records[end].superCell == records[begin].superCell )
    {
      ++end;
    }
    for( std::size_t record = begin; record < end; ++record )
    {
      localIndex const cell = records[record].cell;
      topology.cellToVertex[cell] = vertex;
      topology.vertexToCells.emplaceBack( vertex, cell );
      representative = std::min( representative, static_cast< globalIndex >( arrays.globalCellIds[cell] ) );
    }
    topology.representativeGlobalCellIds[vertex] = representative;
    begin = end;
  }
}

void buildExactTopology( vtkUnstructuredGrid const & mesh,
                         MeshArrays const & arrays,
                         HybridPartitionTopology & topology )
{
  vtkUnstructuredGrid & mutableMesh = const_cast< vtkUnstructuredGrid & >( mesh );
  localIndex const numCells = LvArray::integerConversion< localIndex >( mutableMesh.GetNumberOfCells() );
  int64_t totalFaces = 0;
  int64_t totalIncidences = 0;
  for( localIndex cell = 0; cell < numCells; ++cell )
  {
    CellModel const model = cellModel( arrays.cellTypes[cell] );
    totalFaces += model.numFaces;
    totalIncidences += model.numNodes;
  }

  stdVector< FaceRecord > faces;
  stdVector< PointVertexRecord > pointVertices;
  faces.reserve( LvArray::integerConversion< std::size_t >( totalFaces ) );
  pointVertices.reserve( LvArray::integerConversion< std::size_t >( totalIncidences ) );

  for( localIndex cell = 0; cell < numCells; ++cell )
  {
    CellModel const model = cellModel( arrays.cellTypes[cell] );
    vtkIdType const begin = arrays.offsets[cell];
    vtkIdType const end = arrays.offsets[cell + 1];
    GEOS_ERROR_IF_NE_MSG( end - begin, model.numNodes,
                          GEOS_FMT( "Linear {} cell {} has {} points; expected {}",
                                    model.name, arrays.globalCellIds[cell], end - begin, model.numNodes ) );
    vtkIdType const * const cellPoints = arrays.connectivity + begin;

    std::array< vtkIdType, 8 > localPoints = {};
    for( localIndex node = 0; node < model.numNodes; ++node )
    {
      vtkIdType const localPoint = cellPoints[node];
      GEOS_ERROR_IF( localPoint < 0 || localPoint >= mutableMesh.GetNumberOfPoints(),
                     GEOS_FMT( "Cell {} references invalid local point {}",
                               arrays.globalCellIds[cell], localPoint ) );
      localPoints[node] = localPoint;
      pointVertices.push_back( { arrays.globalPointIds[localPoint], topology.cellToVertex[cell] } );
    }
    std::sort( localPoints.begin(), localPoints.begin() + model.numNodes );
    GEOS_ERROR_IF( std::adjacent_find( localPoints.begin(), localPoints.begin() + model.numNodes ) !=
                   localPoints.begin() + model.numNodes,
                   GEOS_FMT( "Cell {} contains duplicate point references", arrays.globalCellIds[cell] ) );

    forEachCellFace( arrays.cellTypes[cell], [&]( localIndex const faceSize, localIndex const * const faceNodes )
    {
      FaceRecord record;
      record.cell = cell;
      record.size = LvArray::integerConversion< uint8_t >( faceSize );
      for( localIndex node = 0; node < faceSize; ++node )
      {
        record.points[node] = arrays.globalPointIds[cellPoints[faceNodes[node]]];
      }
      std::sort( record.points.begin(), record.points.begin() + faceSize );
      GEOS_ERROR_IF( std::adjacent_find( record.points.begin(), record.points.begin() + faceSize ) !=
                     record.points.begin() + faceSize,
                     GEOS_FMT( "Cell {} has a face with duplicate global point IDs", arrays.globalCellIds[cell] ) );
      faces.push_back( record );
    } );
  }

  std::sort( faces.begin(), faces.end(), faceKeyLess );
  stdVector< std::array< localIndex, 2 > > interiorFaces;
  interiorFaces.reserve( faces.size() / 2 );
  for( std::size_t begin = 0; begin < faces.size(); )
  {
    std::size_t end = begin + 1;
    while( end < faces.size() && sameFaceKey( faces[begin], faces[end] ) )
    {
      ++end;
    }
    std::size_t const incidence = end - begin;
    if( incidence == 2 )
    {
      GEOS_ERROR_IF_EQ_MSG( faces[begin].cell, faces[begin + 1].cell,
                            "A VTK cell generated the same canonical face more than once" );
      interiorFaces.push_back( { topology.cellToVertex[faces[begin].cell],
                                 topology.cellToVertex[faces[begin + 1].cell] } );
    }
    else if( incidence > 2 )
    {
      string pointList;
      for( uint8_t node = 0; node < faces[begin].size; ++node )
      {
        pointList += node == 0 ? GEOS_FMT( "{}", faces[begin].points[node] )
                              : GEOS_FMT( ",{}", faces[begin].points[node] );
      }
      GEOS_ERROR( GEOS_FMT( "Non-manifold face [{}] is shared by {} 3D cells; "
                            "hybrid partitioning requires at most two",
                            pointList, incidence ) );
    }
    begin = end;
  }
  stdVector< FaceRecord >().swap( faces );

  topology.faceToVertices.resize( interiorFaces.size(), 2 );
  topology.exactFaceWeights.resize( interiorFaces.size() );
  topology.exactFaceWeights.setValues< serialPolicy >( 1.0 );
  for( localIndex face = 0; face < topology.faceToVertices.size( 0 ); ++face )
  {
    topology.faceToVertices( face, 0 ) = interiorFaces[face][0];
    topology.faceToVertices( face, 1 ) = interiorFaces[face][1];
  }

  std::sort( pointVertices.begin(), pointVertices.end(), []( PointVertexRecord const & lhs,
                                                             PointVertexRecord const & rhs )
  {
    return lhs.point < rhs.point || (lhs.point == rhs.point && lhs.vertex < rhs.vertex);
  } );
  pointVertices.erase( std::unique( pointVertices.begin(), pointVertices.end(),
                                   []( PointVertexRecord const & lhs, PointVertexRecord const & rhs )
  {
    return lhs.point == rhs.point && lhs.vertex == rhs.vertex;
  } ), pointVertices.end() );

  stdVector< localIndex > pointCapacities;
  stdVector< globalIndex > pointIds;
  pointCapacities.reserve( mutableMesh.GetNumberOfPoints() );
  pointIds.reserve( mutableMesh.GetNumberOfPoints() );
  for( std::size_t begin = 0; begin < pointVertices.size(); )
  {
    std::size_t end = begin + 1;
    while( end < pointVertices.size() && pointVertices[end].point == pointVertices[begin].point )
    {
      ++end;
    }
    pointIds.push_back( pointVertices[begin].point );
    pointCapacities.push_back( LvArray::integerConversion< localIndex >( end - begin ) );
    begin = end;
  }

  localIndex const numUsedPoints = LvArray::integerConversion< localIndex >( pointCapacities.size() );
  topology.pointToVertices.resizeFromCapacities< serialPolicy >( numUsedPoints, pointCapacities.data() );
  topology.pointGlobalIds.resize( numUsedPoints );
  topology.exactPointWeights.resize( numUsedPoints );
  topology.exactPointWeights.setValues< serialPolicy >( 1.0 );

  localIndex point = 0;
  for( std::size_t begin = 0; begin < pointVertices.size(); ++point )
  {
    std::size_t end = begin + 1;
    while( end < pointVertices.size() && pointVertices[end].point == pointVertices[begin].point )
    {
      ++end;
    }
    topology.pointGlobalIds[point] = pointVertices[begin].point;
    for( std::size_t record = begin; record < end; ++record )
    {
      topology.pointToVertices.emplaceBack( point, pointVertices[record].vertex );
    }
    begin = end;
  }
}

void buildVertexWeights( vtkUnstructuredGrid const & mesh,
                         MeshArrays const & arrays,
                         HybridPartitionOptions const & options,
                         HybridPartitionTopology & topology )
{
  localIndex const numCells = topology.cellToVertex.size();
  localIndex const numVertices = topology.vertexToCells.size();
  vtkDataArray * const fvmField = validateWeightField( mesh, options.fvmWeightField );
  vtkDataArray * const femField = validateWeightField( mesh, options.femWeightField );
  vtkDataArray * const memoryField = validateWeightField( mesh, options.memoryWeightField );

  stdVector< std::array< long double, NUM_CONSTRAINTS > > rawWeights( numVertices );
  bool anyPositive[NUM_CONSTRAINTS] = { false, false, false };
  for( localIndex cell = 0; cell < numCells; ++cell )
  {
    CellModel const model = cellModel( arrays.cellTypes[cell] );
    real64 values[NUM_CONSTRAINTS] = {
      1.0 + model.numFaces + model.numNodes,
      1.0 + model.quadratureProxy + model.numNodes + model.numNodes * model.numNodes,
      1.0 + model.numNodes + model.numFaces + model.stateScalarProxy
    };
    vtkDataArray * const fields[NUM_CONSTRAINTS] = { fvmField, femField, memoryField };
    for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
    {
      if( fields[constraint] != nullptr )
      {
        values[constraint] = fields[constraint]->GetTuple1( cell );
      }
      GEOS_ERROR_IF( !std::isfinite( values[constraint] ) || values[constraint] < 0.0,
                     GEOS_FMT( "Hybrid partition constraint {} has invalid weight {} on global cell {}",
                               constraint, values[constraint], arrays.globalCellIds[cell] ) );
      anyPositive[constraint] = anyPositive[constraint] || values[constraint] > 0.0;
      rawWeights[topology.cellToVertex[cell]][constraint] += values[constraint];
    }
  }
  for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
  {
    GEOS_ERROR_IF( !anyPositive[constraint],
                   GEOS_FMT( "Hybrid partition constraint {} has no positive cell weights", constraint ) );
  }

  if( options.fractureWeight > 0 )
  {
    for( localIndex vertex = 0; vertex < numVertices; ++vertex )
    {
      if( topology.vertexToCells.sizeOfArray( vertex ) > 1 )
      {
        for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
        {
          rawWeights[vertex][constraint] += options.fractureWeight;
        }
      }
    }
  }

  topology.vertexWeights.resize( numVertices, NUM_CONSTRAINTS );
  for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
  {
    long double sum = 0.0;
    int64_t positive = 0;
    for( localIndex vertex = 0; vertex < numVertices; ++vertex )
    {
      if( rawWeights[vertex][constraint] > 0.0 )
      {
        sum += rawWeights[vertex][constraint];
        ++positive;
      }
    }
    GEOS_ERROR_IF( positive == 0 || !(sum > 0.0),
                   GEOS_FMT( "Hybrid partition constraint {} cannot be normalized", constraint ) );
    long double const mean = sum / positive;

    for( localIndex vertex = 0; vertex < numVertices; ++vertex )
    {
      long double const raw = rawWeights[vertex][constraint];
      if( raw == 0.0 )
      {
        topology.vertexWeights( vertex, constraint ) = 0;
        continue;
      }
      long double scaled = std::round( raw / mean * INTEGER_WEIGHT_RESOLUTION );
      scaled = std::max< long double >( scaled, 1.0 );
      GEOS_ERROR_IF( scaled > std::numeric_limits< pmet_idx_t >::max(),
                     GEOS_FMT( "Hybrid vertex weight overflow in constraint {}", constraint ) );
      topology.vertexWeights( vertex, constraint ) = static_cast< pmet_idx_t >( scaled );
    }
  }
}

void buildReverseIncidence( HybridPartitionTopology & topology )
{
  localIndex const numVertices = topology.vertexToCells.size();
  array1d< localIndex > faceCapacities( numVertices );
  for( localIndex face = 0; face < topology.faceToVertices.size( 0 ); ++face )
  {
    localIndex const first = topology.faceToVertices( face, 0 );
    localIndex const second = topology.faceToVertices( face, 1 );
    ++faceCapacities[first];
    if( second != first )
    {
      ++faceCapacities[second];
    }
  }
  topology.vertexToFaces.resizeFromCapacities< serialPolicy >( numVertices, faceCapacities.data() );
  for( localIndex face = 0; face < topology.faceToVertices.size( 0 ); ++face )
  {
    localIndex const first = topology.faceToVertices( face, 0 );
    localIndex const second = topology.faceToVertices( face, 1 );
    topology.vertexToFaces.emplaceBack( first, face );
    if( second != first )
    {
      topology.vertexToFaces.emplaceBack( second, face );
    }
  }

  array1d< localIndex > pointCapacities( numVertices );
  for( localIndex point = 0; point < topology.pointToVertices.size(); ++point )
  {
    for( localIndex const vertex : topology.pointToVertices[point] )
    {
      ++pointCapacities[vertex];
    }
  }
  topology.vertexToPoints.resizeFromCapacities< serialPolicy >( numVertices, pointCapacities.data() );
  for( localIndex point = 0; point < topology.pointToVertices.size(); ++point )
  {
    for( localIndex const vertex : topology.pointToVertices[point] )
    {
      topology.vertexToPoints.emplaceBack( vertex, point );
    }
  }
}

void buildWeightedGraph( HybridPartitionOptions const & options,
                         HybridPartitionTopology & topology )
{
  localIndex const numVertices = topology.vertexToCells.size();
  stdVector< EdgeContribution > contributions;
  contributions.reserve( topology.faceToVertices.size( 0 ) +
                         topology.pointToVertices.toViewConst().getOffsets()[topology.pointToVertices.size()] );

  auto const addContribution = [&contributions]( localIndex first, localIndex second, real64 const weight )
  {
    if( first == second || weight <= 0.0 )
    {
      return;
    }
    if( first > second )
    {
      std::swap( first, second );
    }
    contributions.push_back( { first, second, weight } );
  };

  for( localIndex face = 0; face < topology.faceToVertices.size( 0 ); ++face )
  {
    addContribution( topology.faceToVertices( face, 0 ),
                     topology.faceToVertices( face, 1 ),
                     options.fvmCommunicationWeight * topology.exactFaceWeights[face] );
  }

  for( localIndex point = 0; point < topology.pointToVertices.size(); ++point )
  {
    arraySlice1d< localIndex const > const vertices = topology.pointToVertices[point];
    if( vertices.size() < 2 )
    {
      continue;
    }
    localIndex anchor = vertices[0];
    for( localIndex const vertex : vertices )
    {
      if( topology.representativeGlobalCellIds[vertex] < topology.representativeGlobalCellIds[anchor] ||
          (topology.representativeGlobalCellIds[vertex] == topology.representativeGlobalCellIds[anchor] &&
           vertex < anchor) )
      {
        anchor = vertex;
      }
    }
    for( localIndex const vertex : vertices )
    {
      if( vertex != anchor )
      {
        addContribution( anchor, vertex,
                         options.femCommunicationWeight * topology.exactPointWeights[point] );
      }
    }
  }

  std::sort( contributions.begin(), contributions.end(), []( EdgeContribution const & lhs,
                                                              EdgeContribution const & rhs )
  {
    return lhs.first < rhs.first || (lhs.first == rhs.first && lhs.second < rhs.second);
  } );

  stdVector< AggregatedEdge > edges;
  edges.reserve( contributions.size() );
  for( std::size_t begin = 0; begin < contributions.size(); )
  {
    std::size_t end = begin + 1;
    long double weight = contributions[begin].weight;
    while( end < contributions.size() &&
           contributions[end].first == contributions[begin].first &&
           contributions[end].second == contributions[begin].second )
    {
      weight += contributions[end].weight;
      ++end;
    }
    long double scaled = std::round( weight * INTEGER_WEIGHT_RESOLUTION );
    scaled = std::max< long double >( scaled, 1.0 );
    GEOS_ERROR_IF( scaled > std::numeric_limits< pmet_idx_t >::max(),
                   GEOS_FMT( "Hybrid graph edge ({},{}) weight overflow",
                             contributions[begin].first, contributions[begin].second ) );
    edges.push_back( { contributions[begin].first,
                       contributions[begin].second,
                       static_cast< pmet_idx_t >( scaled ) } );
    begin = end;
  }
  stdVector< EdgeContribution >().swap( contributions );

  array1d< pmet_idx_t > capacities( numVertices );
  for( AggregatedEdge const & edge : edges )
  {
    ++capacities[edge.first];
    ++capacities[edge.second];
  }
  topology.graph.resizeFromCapacities< serialPolicy >( numVertices, capacities.data() );
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > weightRows;
  weightRows.resizeFromCapacities< serialPolicy >( numVertices, capacities.data() );
  for( AggregatedEdge const & edge : edges )
  {
    topology.graph.emplaceBack( edge.first, edge.second );
    weightRows.emplaceBack( edge.first, edge.weight );
    topology.graph.emplaceBack( edge.second, edge.first );
    weightRows.emplaceBack( edge.second, edge.weight );
  }

  pmet_idx_t const numAdjacencies = topology.graph.toViewConst().getOffsets()[numVertices];
  topology.edgeWeights.resize( numAdjacencies );
  std::copy( weightRows.toViewConst().getValues(),
             weightRows.toViewConst().getValues() + numAdjacencies,
             topology.edgeWeights.begin() );

  // Validate sorted, duplicate-free rows and exact weighted symmetry.
  auto const graph = topology.graph.toViewConst();
  for( localIndex vertex = 0; vertex < numVertices; ++vertex )
  {
    pmet_idx_t previous = -1;
    pmet_idx_t const rowBegin = graph.getOffsets()[vertex];
    pmet_idx_t const rowEnd = graph.getOffsets()[vertex + 1];
    for( pmet_idx_t entry = rowBegin; entry < rowEnd; ++entry )
    {
      pmet_idx_t const neighbor = graph.getValues()[entry];
      GEOS_ERROR_IF( neighbor < 0 || neighbor >= numVertices || neighbor == vertex,
                     GEOS_FMT( "Invalid hybrid graph edge {} -> {}", vertex, neighbor ) );
      GEOS_ERROR_IF( neighbor <= previous,
                     GEOS_FMT( "Hybrid graph row {} is not sorted and duplicate-free", vertex ) );
      GEOS_ERROR_IF( topology.edgeWeights[entry] <= 0,
                     GEOS_FMT( "Hybrid graph edge {} -> {} has nonpositive weight", vertex, neighbor ) );
      previous = neighbor;

      pmet_idx_t const reverseBegin = graph.getOffsets()[neighbor];
      pmet_idx_t const reverseEnd = graph.getOffsets()[neighbor + 1];
      pmet_idx_t const * const reverse = std::lower_bound( graph.getValues() + reverseBegin,
                                                           graph.getValues() + reverseEnd,
                                                           static_cast< pmet_idx_t >( vertex ) );
      GEOS_ERROR_IF( reverse == graph.getValues() + reverseEnd || *reverse != vertex,
                     GEOS_FMT( "Hybrid graph edge {} -> {} is not symmetric", vertex, neighbor ) );
      pmet_idx_t const reverseEntry = reverse - graph.getValues();
      GEOS_ERROR_IF_NE_MSG( topology.edgeWeights[entry], topology.edgeWeights[reverseEntry],
                            "Hybrid graph has asymmetric edge weights" );
    }
  }
}

struct ExactStats
{
  int64_t cutFaces = 0;
  real64 weightedFaceCut = 0.0;
  int64_t sharedPoints = 0;
  int64_t pointReplication = 0;
  real64 weightedPointReplication = 0.0;
  real64 objective = 0.0;
  stdVector< stdVector< int > > neighbors;
};

ExactStats exactStats( HybridPartitionTopology const & topology,
                       arrayView1d< int64_t const > const parts,
                       int const numParts,
                       HybridPartitionOptions const & options )
{
  ExactStats result;
  std::unordered_set< uint64_t > rankPairs;
  rankPairs.reserve( static_cast< std::size_t >( numParts ) * 8 );
  stdVector< int > pointParts;

  for( localIndex face = 0; face < topology.faceToVertices.size( 0 ); ++face )
  {
    int const first = LvArray::integerConversion< int >( parts[topology.faceToVertices( face, 0 )] );
    int const second = LvArray::integerConversion< int >( parts[topology.faceToVertices( face, 1 )] );
    if( first != second )
    {
      ++result.cutFaces;
      result.weightedFaceCut += topology.exactFaceWeights[face];
      rankPairs.insert( rankPairKey( first, second ) );
    }
  }

  for( localIndex point = 0; point < topology.pointToVertices.size(); ++point )
  {
    pointParts.clear();
    pointParts.reserve( topology.pointToVertices.sizeOfArray( point ) );
    for( localIndex const vertex : topology.pointToVertices[point] )
    {
      pointParts.push_back( LvArray::integerConversion< int >( parts[vertex] ) );
    }
    std::sort( pointParts.begin(), pointParts.end() );
    pointParts.erase( std::unique( pointParts.begin(), pointParts.end() ), pointParts.end() );
    if( pointParts.size() > 1 )
    {
      ++result.sharedPoints;
      int64_t const replication = LvArray::integerConversion< int64_t >( pointParts.size() - 1 );
      result.pointReplication += replication;
      result.weightedPointReplication += replication * topology.exactPointWeights[point];
      for( std::size_t i = 0; i < pointParts.size(); ++i )
      {
        for( std::size_t j = i + 1; j < pointParts.size(); ++j )
        {
          rankPairs.insert( rankPairKey( pointParts[i], pointParts[j] ) );
        }
      }
    }
  }

  result.objective = options.fvmCommunicationWeight * result.weightedFaceCut +
                     options.femCommunicationWeight * result.weightedPointReplication +
                     options.neighborPenalty * rankPairs.size();

  result.neighbors.resize( numParts );
  for( uint64_t const pair : rankPairs )
  {
    int const first = static_cast< int >( pair >> 32 );
    int const second = static_cast< int >( pair & 0xFFFFFFFFu );
    result.neighbors[first].push_back( second );
    result.neighbors[second].push_back( first );
  }
  for( stdVector< int > & neighbors : result.neighbors )
  {
    std::sort( neighbors.begin(), neighbors.end() );
  }
  return result;
}

class RefinementState
{
public:
  RefinementState( HybridPartitionTopology const & topology,
                   HybridPartitionOptions const & options,
                   int const numParts,
                   array1d< int64_t > & parts )
    : m_topology( topology ),
      m_options( options ),
      m_numParts( numParts ),
      m_parts( parts ),
      m_loads( numParts * NUM_CONSTRAINTS, 0.0L ),
      m_limits( NUM_CONSTRAINTS, 0.0L ),
      m_partCounts( numParts, 0 )
  {
    std::array< long double, NUM_CONSTRAINTS > totals = {};
    for( localIndex vertex = 0; vertex < m_topology.vertexToCells.size(); ++vertex )
    {
      int const part = LvArray::integerConversion< int >( m_parts[vertex] );
      GEOS_ERROR_IF( part < 0 || part >= m_numParts,
                     GEOS_FMT( "Hybrid partition contains invalid part {} for vertex {}", part, vertex ) );
      ++m_partCounts[part];
      for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
      {
        long double const weight = m_topology.vertexWeights( vertex, constraint );
        load( part, constraint ) += weight;
        totals[constraint] += weight;
      }
    }
    for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
    {
      m_limits[constraint] = (1.0L + m_options.imbalance[constraint]) *
                             totals[constraint] / m_numParts;
    }

    if( m_options.neighborPenalty > 0.0 )
    {
      for( localIndex face = 0; face < m_topology.faceToVertices.size( 0 ); ++face )
      {
        int const first = LvArray::integerConversion< int >( m_parts[m_topology.faceToVertices( face, 0 )] );
        int const second = LvArray::integerConversion< int >( m_parts[m_topology.faceToVertices( face, 1 )] );
        appendRankPair( m_rankPairMultiplicities, first, second, 1 );
      }
      for( localIndex point = 0; point < m_topology.pointToVertices.size(); ++point )
      {
        appendPartPairs( uniqueEntityParts( m_topology.pointToVertices[point], m_parts.toViewConst() ),
                         m_rankPairMultiplicities, 1 );
      }
    }
  }

  bool partIsOverloaded( int const part ) const
  {
    for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
    {
      if( load( part, constraint ) > m_limits[constraint] + 1.0e-9L )
      {
        return true;
      }
    }
    return false;
  }

  bool anyPartIsOverloaded() const
  {
    for( int part = 0; part < m_numParts; ++part )
    {
      if( partIsOverloaded( part ) )
      {
        return true;
      }
    }
    return false;
  }

  bool canMove( localIndex const vertex, int const destination ) const
  {
    int const source = LvArray::integerConversion< int >( m_parts[vertex] );
    if( source == destination || destination < 0 || destination >= m_numParts || m_partCounts[source] <= 1 )
    {
      return false;
    }
    for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
    {
      if( load( destination, constraint ) + m_topology.vertexWeights( vertex, constraint ) >
          m_limits[constraint] + 1.0e-9L )
      {
        return false;
      }
    }
    return true;
  }

  real64 gain( localIndex const vertex, int const destination ) const
  {
    int const source = LvArray::integerConversion< int >( m_parts[vertex] );
    if( source == destination )
    {
      return 0.0;
    }
    long double result = 0.0;
    for( localIndex const face : m_topology.vertexToFaces[vertex] )
    {
      localIndex const firstVertex = m_topology.faceToVertices( face, 0 );
      localIndex const secondVertex = m_topology.faceToVertices( face, 1 );
      if( firstVertex == secondVertex )
      {
        continue;
      }
      int const firstOld = LvArray::integerConversion< int >( m_parts[firstVertex] );
      int const secondOld = LvArray::integerConversion< int >( m_parts[secondVertex] );
      int const firstNew = firstVertex == vertex ? destination : firstOld;
      int const secondNew = secondVertex == vertex ? destination : secondOld;
      result += m_options.fvmCommunicationWeight * m_topology.exactFaceWeights[face] *
                (static_cast< int >( firstOld != secondOld ) - static_cast< int >( firstNew != secondNew ));
    }
    for( localIndex const point : m_topology.vertexToPoints[vertex] )
    {
      stdVector< int > const oldParts = uniqueEntityParts( m_topology.pointToVertices[point], m_parts.toViewConst() );
      stdVector< int > const newParts = uniqueEntityParts( m_topology.pointToVertices[point], m_parts.toViewConst(),
                                                          vertex, destination );
      result += m_options.femCommunicationWeight * m_topology.exactPointWeights[point] *
                (static_cast< int64_t >( oldParts.size() ) - static_cast< int64_t >( newParts.size() ));
    }

    if( m_options.neighborPenalty > 0.0 )
    {
      std::map< uint64_t, int64_t > deltas;
      pairDeltas( vertex, destination, deltas );
      for( auto const & [key, delta] : deltas )
      {
        auto const oldIt = m_rankPairMultiplicities.find( key );
        int64_t const oldCount = oldIt == m_rankPairMultiplicities.end() ? 0 : oldIt->second;
        int64_t const newCount = oldCount + delta;
        GEOS_ERROR_IF_LT_MSG( newCount, 0, "Hybrid move produced a negative neighbor multiplicity" );
        result += m_options.neighborPenalty *
                  (static_cast< int >( oldCount > 0 ) - static_cast< int >( newCount > 0 ));
      }
    }
    return static_cast< real64 >( result );
  }

  void apply( localIndex const vertex, int const destination )
  {
    int const source = LvArray::integerConversion< int >( m_parts[vertex] );
    if( m_options.neighborPenalty > 0.0 )
    {
      std::map< uint64_t, int64_t > deltas;
      pairDeltas( vertex, destination, deltas );
      for( auto const & [key, delta] : deltas )
      {
        int64_t & count = m_rankPairMultiplicities[key];
        count += delta;
        GEOS_ERROR_IF_LT_MSG( count, 0, "Hybrid neighbor multiplicity became negative" );
        if( count == 0 )
        {
          m_rankPairMultiplicities.erase( key );
        }
      }
    }
    for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
    {
      long double const weight = m_topology.vertexWeights( vertex, constraint );
      load( source, constraint ) -= weight;
      load( destination, constraint ) += weight;
    }
    --m_partCounts[source];
    ++m_partCounts[destination];
    m_parts[vertex] = destination;
  }

  long double normalizedDestinationLoad( localIndex const vertex, int const destination ) const
  {
    long double maximum = 0.0;
    for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
    {
      if( m_limits[constraint] > 0.0 )
      {
        maximum = std::max( maximum,
                            (load( destination, constraint ) +
                             m_topology.vertexWeights( vertex, constraint )) / m_limits[constraint] );
      }
    }
    return maximum;
  }

private:
  long double & load( int const part, localIndex const constraint )
  {
    return m_loads[part * NUM_CONSTRAINTS + constraint];
  }

  long double load( int const part, localIndex const constraint ) const
  {
    return m_loads[part * NUM_CONSTRAINTS + constraint];
  }

  void pairDeltas( localIndex const vertex,
                   int const destination,
                   std::map< uint64_t, int64_t > & deltas ) const
  {
    for( localIndex const face : m_topology.vertexToFaces[vertex] )
    {
      localIndex const firstVertex = m_topology.faceToVertices( face, 0 );
      localIndex const secondVertex = m_topology.faceToVertices( face, 1 );
      if( firstVertex == secondVertex )
      {
        continue;
      }
      int const firstOld = LvArray::integerConversion< int >( m_parts[firstVertex] );
      int const secondOld = LvArray::integerConversion< int >( m_parts[secondVertex] );
      int const firstNew = firstVertex == vertex ? destination : firstOld;
      int const secondNew = secondVertex == vertex ? destination : secondOld;
      if( firstOld != secondOld )
      {
        --deltas[rankPairKey( firstOld, secondOld )];
      }
      if( firstNew != secondNew )
      {
        ++deltas[rankPairKey( firstNew, secondNew )];
      }
    }
    for( localIndex const point : m_topology.vertexToPoints[vertex] )
    {
      appendPartPairDeltas( uniqueEntityParts( m_topology.pointToVertices[point], m_parts.toViewConst() ),
                            uniqueEntityParts( m_topology.pointToVertices[point], m_parts.toViewConst(),
                                               vertex, destination ),
                            deltas );
    }
  }

  HybridPartitionTopology const & m_topology;
  HybridPartitionOptions const & m_options;
  int m_numParts;
  array1d< int64_t > & m_parts;
  stdVector< long double > m_loads;
  stdVector< long double > m_limits;
  stdVector< int64_t > m_partCounts;
  std::unordered_map< uint64_t, int64_t > m_rankPairMultiplicities;
};

stdVector< localIndex > orderedVertices( HybridPartitionTopology const & topology )
{
  stdVector< localIndex > order( topology.vertexToCells.size() );
  std::iota( order.begin(), order.end(), 0 );
  std::sort( order.begin(), order.end(), [&topology]( localIndex const lhs, localIndex const rhs )
  {
    return topology.representativeGlobalCellIds[lhs] < topology.representativeGlobalCellIds[rhs] ||
           (topology.representativeGlobalCellIds[lhs] == topology.representativeGlobalCellIds[rhs] && lhs < rhs);
  } );
  return order;
}

void ensureNonemptyParts( HybridPartitionTopology const & topology,
                          int const numParts,
                          stdVector< localIndex > const & order,
                          array1d< int64_t > & parts )
{
  stdVector< int64_t > counts( numParts, 0 );
  for( int64_t const part : parts )
  {
    GEOS_ERROR_IF( part < 0 || part >= numParts,
                   GEOS_FMT( "METIS returned invalid partition {}", part ) );
    ++counts[part];
  }

  for( int emptyPart = 0; emptyPart < numParts; ++emptyPart )
  {
    if( counts[emptyPart] != 0 )
    {
      continue;
    }
    int donor = -1;
    for( int part = 0; part < numParts; ++part )
    {
      if( counts[part] > 1 &&
          (donor < 0 || counts[part] > counts[donor] ||
           (counts[part] == counts[donor] && part < donor)) )
      {
        donor = part;
      }
    }
    GEOS_ERROR_IF( donor < 0,
                   GEOS_FMT( "Cannot populate empty hybrid partition {} without emptying another partition",
                             emptyPart ) );

    localIndex selected = -1;
    pmet_idx_t selectedWeight = std::numeric_limits< pmet_idx_t >::max();
    for( localIndex const vertex : order )
    {
      if( parts[vertex] != donor )
      {
        continue;
      }
      pmet_idx_t weight = 0;
      for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
      {
        weight += topology.vertexWeights( vertex, constraint );
      }
      if( selected < 0 || weight < selectedWeight )
      {
        selected = vertex;
        selectedWeight = weight;
      }
    }
    GEOS_ERROR_IF_LT_MSG( selected, 0, "Failed to select a vertex for an empty hybrid partition" );
    parts[selected] = emptyPart;
    --counts[donor];
    ++counts[emptyPart];
  }
}

int64_t repairBalance( HybridPartitionTopology const & topology,
                       HybridPartitionOptions const & options,
                       int const numParts,
                       stdVector< localIndex > const & order,
                       array1d< int64_t > & parts )
{
  RefinementState state( topology, options, numParts, parts );
  int64_t moves = 0;
  for( integer pass = 0; pass < 4 && state.anyPartIsOverloaded(); ++pass )
  {
    bool changed = false;
    for( localIndex const vertex : order )
    {
      int const source = LvArray::integerConversion< int >( parts[vertex] );
      if( !state.partIsOverloaded( source ) )
      {
        continue;
      }
      int bestDestination = -1;
      long double bestLoad = std::numeric_limits< long double >::max();
      for( int destination = 0; destination < numParts; ++destination )
      {
        if( state.canMove( vertex, destination ) )
        {
          long double const load = state.normalizedDestinationLoad( vertex, destination );
          if( load < bestLoad ||
              (std::abs( load - bestLoad ) < 1.0e-18L && destination < bestDestination) )
          {
            bestLoad = load;
            bestDestination = destination;
          }
        }
      }
      if( bestDestination >= 0 )
      {
        state.apply( vertex, bestDestination );
        ++moves;
        changed = true;
      }
    }
    if( !changed )
    {
      break;
    }
  }
  GEOS_WARNING_IF( state.anyPartIsOverloaded(),
                   "Hybrid balance repair could not satisfy every tolerance; an indivisible vertex or "
                   "conflicting multi-constraint loads may be responsible" );
  return moves;
}

int64_t refineExactObjective( HybridPartitionTopology const & topology,
                              HybridPartitionOptions const & options,
                              int const numParts,
                              stdVector< localIndex > const & order,
                              array1d< int64_t > & parts )
{
  if( options.refinementPasses <= 0 || numParts == 1 )
  {
    return 0;
  }

  array1d< int64_t > const baselineParts = parts;
  real64 const baselineObjective = evaluateHybridObjective( topology, parts.toViewConst(), numParts, options );
  RefinementState state( topology, options, numParts, parts );
  stdVector< int > markers( numParts, -1 );
  stdVector< int > candidates;
  candidates.reserve( std::min( numParts, 32 ) );
  int marker = 0;
  int64_t moves = 0;

  for( integer pass = 0; pass < options.refinementPasses; ++pass )
  {
    int64_t passMoves = 0;
    for( localIndex const vertex : order )
    {
      ++marker;
      if( marker == std::numeric_limits< int >::max() )
      {
        std::fill( markers.begin(), markers.end(), -1 );
        marker = 0;
      }
      candidates.clear();
      int const source = LvArray::integerConversion< int >( parts[vertex] );
      auto addCandidate = [&]( int const destination )
      {
        if( destination != source && markers[destination] != marker )
        {
          markers[destination] = marker;
          candidates.push_back( destination );
        }
      };

      for( localIndex const face : topology.vertexToFaces[vertex] )
      {
        addCandidate( LvArray::integerConversion< int >( parts[topology.faceToVertices( face, 0 )] ) );
        addCandidate( LvArray::integerConversion< int >( parts[topology.faceToVertices( face, 1 )] ) );
      }
      for( localIndex const point : topology.vertexToPoints[vertex] )
      {
        for( localIndex const neighbor : topology.pointToVertices[point] )
        {
          addCandidate( LvArray::integerConversion< int >( parts[neighbor] ) );
        }
      }
      std::sort( candidates.begin(), candidates.end() );

      int bestDestination = -1;
      real64 bestGain = 0.0;
      for( int const destination : candidates )
      {
        if( !state.canMove( vertex, destination ) )
        {
          continue;
        }
        real64 const gain = state.gain( vertex, destination );
        real64 constexpr epsilon = 1.0e-12;
        if( gain > bestGain + epsilon ||
            (std::abs( gain - bestGain ) <= epsilon && gain > epsilon && destination < bestDestination) )
        {
          bestGain = gain;
          bestDestination = destination;
        }
      }
      if( bestDestination >= 0 )
      {
        state.apply( vertex, bestDestination );
        ++passMoves;
        ++moves;
      }
    }
    if( passMoves == 0 )
    {
      break;
    }
  }

  real64 const finalObjective = evaluateHybridObjective( topology, parts.toViewConst(), numParts, options );
  real64 const tolerance = 1.0e-10 * std::max( 1.0, std::abs( baselineObjective ) );
  if( finalObjective > baselineObjective + tolerance )
  {
    GEOS_WARNING( GEOS_FMT( "Hybrid refinement objective unexpectedly increased from {} to {}; reverting",
                            baselineObjective, finalObjective ) );
    parts = baselineParts;
    return 0;
  }
  return moves;
}

localIndex findRoot( stdVector< localIndex > & parent, localIndex vertex )
{
  localIndex root = vertex;
  while( parent[root] != root )
  {
    root = parent[root];
  }
  while( parent[vertex] != vertex )
  {
    localIndex const next = parent[vertex];
    parent[vertex] = root;
    vertex = next;
  }
  return root;
}

int64_t maxConnectedComponents( HybridPartitionTopology const & topology,
                                arrayView1d< int64_t const > const parts,
                                int const numParts )
{
  localIndex const numVertices = topology.vertexToCells.size();
  stdVector< localIndex > parent( numVertices );
  std::iota( parent.begin(), parent.end(), 0 );
  for( localIndex face = 0; face < topology.faceToVertices.size( 0 ); ++face )
  {
    localIndex const first = topology.faceToVertices( face, 0 );
    localIndex const second = topology.faceToVertices( face, 1 );
    if( first == second || parts[first] != parts[second] )
    {
      continue;
    }
    localIndex const firstRoot = findRoot( parent, first );
    localIndex const secondRoot = findRoot( parent, second );
    if( firstRoot != secondRoot )
    {
      parent[secondRoot] = firstRoot;
    }
  }

  stdVector< std::unordered_set< localIndex > > roots( numParts );
  for( localIndex vertex = 0; vertex < numVertices; ++vertex )
  {
    roots[parts[vertex]].insert( findRoot( parent, vertex ) );
  }
  int64_t maximum = 0;
  for( auto const & partRoots : roots )
  {
    maximum = std::max( maximum, LvArray::integerConversion< int64_t >( partRoots.size() ) );
  }
  return maximum;
}

HybridPartitionMetrics buildMetrics( HybridPartitionTopology const & topology,
                                     arrayView1d< int64_t const > const parts,
                                     int const numParts,
                                     HybridPartitionOptions const & options,
                                     ExactStats const & initialStats,
                                     ExactStats const & finalStats )
{
  HybridPartitionMetrics metrics;
  metrics.numCells = topology.cellToVertex.size();
  metrics.numPartitionVertices = topology.vertexToCells.size();
  metrics.numGraphEdges = topology.graph.toViewConst().getOffsets()[topology.graph.size()] / 2;
  metrics.numInteriorFaces = topology.faceToVertices.size( 0 );
  metrics.numCutFaces = finalStats.cutFaces;
  metrics.weightedFaceCut = finalStats.weightedFaceCut;
  metrics.numSharedPoints = finalStats.sharedPoints;
  metrics.pointReplication = finalStats.pointReplication;
  metrics.weightedPointReplication = finalStats.weightedPointReplication;
  metrics.initialObjective = initialStats.objective;
  metrics.finalObjective = finalStats.objective;
  metrics.initialCutFaces = initialStats.cutFaces;
  metrics.initialPointReplication = initialStats.pointReplication;
  metrics.estimatedRootGraphPeakBytes = topology.estimatedPeakBytes;
  metrics.maxConnectedComponents = maxConnectedComponents( topology, parts, numParts );

  int64_t totalNeighbors = 0;
  for( stdVector< int > const & neighbors : finalStats.neighbors )
  {
    int64_t const count = LvArray::integerConversion< int64_t >( neighbors.size() );
    metrics.maxRankNeighbors = std::max( metrics.maxRankNeighbors, count );
    totalNeighbors += count;
  }
  metrics.averageRankNeighbors = static_cast< real64 >( totalNeighbors ) / numParts;

  metrics.imbalanceByConstraint.resize( NUM_CONSTRAINTS );
  stdVector< long double > loads( numParts * NUM_CONSTRAINTS, 0.0L );
  for( localIndex vertex = 0; vertex < topology.vertexToCells.size(); ++vertex )
  {
    for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
    {
      loads[parts[vertex] * NUM_CONSTRAINTS + constraint] += topology.vertexWeights( vertex, constraint );
    }
  }
  for( localIndex constraint = 0; constraint < NUM_CONSTRAINTS; ++constraint )
  {
    long double total = 0.0;
    long double maximum = 0.0;
    for( int part = 0; part < numParts; ++part )
    {
      long double const value = loads[part * NUM_CONSTRAINTS + constraint];
      total += value;
      maximum = std::max( maximum, value );
    }
    long double const average = total / numParts;
    metrics.imbalanceByConstraint[constraint] =
      average > 0.0 ? static_cast< real64 >( maximum / average - 1.0 ) : 0.0;
  }

  if( options.diagnostics > 0 )
  {
    for( int part = 0; part < numParts; ++part )
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "Hybrid rank {:4}: loads [{:10.0f}, {:10.0f}, {:10.0f}], neighbors {}",
                                 part,
                                 static_cast< real64 >( loads[part * NUM_CONSTRAINTS] ),
                                 static_cast< real64 >( loads[part * NUM_CONSTRAINTS + 1] ),
                                 static_cast< real64 >( loads[part * NUM_CONSTRAINTS + 2] ),
                                 finalStats.neighbors[part].size() ) );
    }
  }
  return metrics;
}

} // namespace

bool isHybridPartitioningSupported( vtkUnstructuredGrid const & mesh,
                                    string & reason )
{
  vtkUnstructuredGrid & mutableMesh = const_cast< vtkUnstructuredGrid & >( mesh );
  if( mutableMesh.GetNumberOfCells() == 0 )
  {
    reason = "the root 3D mesh is empty";
    return false;
  }
  if( mutableMesh.GetNumberOfCells() > std::numeric_limits< localIndex >::max() )
  {
    reason = "the 3D cell count exceeds GEOS localIndex capacity";
    return false;
  }
  for( vtkIdType cell = 0; cell < mutableMesh.GetNumberOfCells(); ++cell )
  {
    unsigned char const type = static_cast< unsigned char >( mutableMesh.GetCellType( cell ) );
    CellModel const model = cellModel( type );
    if( model.numNodes == 0 )
    {
      reason = GEOS_FMT( "cell {} has unsupported VTK type {} (hybrid supports linear tetrahedra, "
                         "pyramids, wedges, and hexahedra)", cell, static_cast< integer >( type ) );
      return false;
    }
    if( mutableMesh.GetCellSize( cell ) != model.numNodes )
    {
      reason = GEOS_FMT( "cell {} is a non-linear or malformed {} with {} points instead of {}",
                         cell, model.name, mutableMesh.GetCellSize( cell ), model.numNodes );
      return false;
    }
  }
  reason.clear();
  return true;
}

int64_t estimateHybridPartitionMemory( vtkUnstructuredGrid const & mesh )
{
  vtkUnstructuredGrid & mutableMesh = const_cast< vtkUnstructuredGrid & >( mesh );
  long double totalFaces = 0.0;
  long double totalIncidences = 0.0;
  for( vtkIdType cell = 0; cell < mutableMesh.GetNumberOfCells(); ++cell )
  {
    CellModel const model = cellModel( static_cast< unsigned char >( mutableMesh.GetCellType( cell ) ) );
    if( model.numNodes == 0 )
    {
      return std::numeric_limits< int64_t >::max();
    }
    totalFaces += model.numFaces;
    totalIncidences += model.numNodes;
  }
  long double const numCells = mutableMesh.GetNumberOfCells();
  long double const numPoints = mutableMesh.GetNumberOfPoints();

  // Peak includes sort records, edge contributions, retained exact incidence,
  // CSR/weights, METIS work estimates, and refinement state with a 20% margin.
  long double const bytes = 1.2L * (
    totalFaces * sizeof( FaceRecord ) +
    totalIncidences * sizeof( PointVertexRecord ) +
    (totalFaces / 2.0L + totalIncidences) * sizeof( EdgeContribution ) +
    (totalFaces + totalIncidences) * (sizeof( localIndex ) + sizeof( real64 )) +
    numCells * (8.0L * sizeof( pmet_idx_t ) + 5.0L * sizeof( localIndex )) +
    numPoints * (sizeof( globalIndex ) + 3.0L * sizeof( localIndex )) );
  if( bytes >= std::numeric_limits< int64_t >::max() )
  {
    return std::numeric_limits< int64_t >::max();
  }
  return static_cast< int64_t >( std::ceil( bytes ) );
}

HybridPartitionTopology buildHybridPartitionTopology(
  vtkUnstructuredGrid const & mesh,
  SuperCellInfo const * const superCells,
  HybridPartitionOptions const & options )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( superCells );
  string reason;
  GEOS_ERROR_IF( !isHybridPartitioningSupported( mesh, reason ),
                 GEOS_FMT( "Cannot build hybrid partition topology: {}", reason ) );
  GEOS_ERROR_IF( options.imbalance.size() != NUM_CONSTRAINTS,
                 GEOS_FMT( "partitionImbalance must contain exactly {} values", NUM_CONSTRAINTS ) );
  GEOS_ERROR_IF( options.fvmCommunicationWeight < 0.0 || options.femCommunicationWeight < 0.0 ||
                 options.neighborPenalty < 0.0,
                 "Hybrid communication and neighbor weights must be nonnegative" );
  GEOS_ERROR_IF( options.fvmCommunicationWeight == 0.0 && options.femCommunicationWeight == 0.0,
                 "At least one hybrid communication weight must be positive" );
  for( real64 const tolerance : options.imbalance )
  {
    GEOS_ERROR_IF( !std::isfinite( tolerance ) || tolerance < 0.0,
                   "Hybrid imbalance tolerances must be finite and nonnegative" );
  }

  MeshArrays const arrays = getMeshArrays( mesh );
  HybridPartitionTopology topology;
  topology.estimatedPeakBytes = estimateHybridPartitionMemory( mesh );
  buildPartitionVertices( mesh, arrays, topology );
  buildExactTopology( mesh, arrays, topology );
  buildVertexWeights( mesh, arrays, options, topology );
  buildReverseIncidence( topology );
  buildWeightedGraph( options, topology );
  return topology;
}

real64 evaluateHybridObjective( HybridPartitionTopology const & topology,
                                arrayView1d< int64_t const > const vertexParts,
                                int const numParts,
                                HybridPartitionOptions const & options )
{
  GEOS_ERROR_IF_NE_MSG( vertexParts.size(), topology.vertexToCells.size(),
                        "Hybrid objective partition size does not match the partition-vertex count" );
  GEOS_ERROR_IF( numParts <= 0, "Hybrid objective requires at least one partition" );
  return exactStats( topology, vertexParts, numParts, options ).objective;
}

HybridPartitionResult partitionHybridMeshOnRoot(
  vtkUnstructuredGrid const & cells3D,
  SuperCellInfo const * const superCells,
  HybridPartitionOptions const & options,
  int const numParts )
{
  GEOS_MARK_FUNCTION;
  GEOS_ERROR_IF( numParts <= 0, "Hybrid partitioning requires at least one target rank" );

  Clock::time_point const buildBegin = Clock::now();
  HybridPartitionTopology topology = buildHybridPartitionTopology( cells3D, superCells, options );
  Clock::time_point const buildEnd = Clock::now();
  localIndex const numVertices = topology.vertexToCells.size();
  GEOS_ERROR_IF( numParts > numVertices,
                 GEOS_FMT( "Cannot partition {} atomic vertices over {} nonempty ranks", numVertices, numParts ) );

  Clock::time_point const partitionBegin = Clock::now();
  array1d< pmet_idx_t > const metisParts = metis::partitionWeighted(
    topology.graph.toViewConst(),
    topology.edgeWeights.toViewConst(),
    topology.vertexWeights.toViewConst(),
    numParts,
    options.imbalance.toViewConst(),
    options.seed );
  array1d< int64_t > vertexParts( numVertices );
  for( localIndex vertex = 0; vertex < numVertices; ++vertex )
  {
    vertexParts[vertex] = metisParts[vertex];
  }
  Clock::time_point const partitionEnd = Clock::now();

  stdVector< localIndex > const order = orderedVertices( topology );
  ensureNonemptyParts( topology, numParts, order, vertexParts );
  repairBalance( topology, options, numParts, order, vertexParts );
  ExactStats const initialStats = exactStats( topology, vertexParts.toViewConst(), numParts, options );

  Clock::time_point const refinementBegin = Clock::now();
  int64_t const refinementMoves = refineExactObjective( topology, options, numParts, order, vertexParts );
  Clock::time_point const refinementEnd = Clock::now();
  ExactStats const finalStats = exactStats( topology, vertexParts.toViewConst(), numParts, options );

  HybridPartitionResult result;
  result.cellParts.resize( topology.cellToVertex.size() );
  for( localIndex cell = 0; cell < topology.cellToVertex.size(); ++cell )
  {
    result.cellParts[cell] = vertexParts[topology.cellToVertex[cell]];
  }
  result.metrics = buildMetrics( topology, vertexParts.toViewConst(), numParts, options,
                                 initialStats, finalStats );
  result.metrics.refinementMoves = refinementMoves;
  result.metrics.graphBuildSeconds = std::chrono::duration< real64 >( buildEnd - buildBegin ).count();
  result.metrics.initialPartitionSeconds = std::chrono::duration< real64 >( partitionEnd - partitionBegin ).count();
  result.metrics.exactRefinementSeconds = std::chrono::duration< real64 >( refinementEnd - refinementBegin ).count();
  result.rankNeighbors = finalStats.neighbors;

  GEOS_LOG_RANK_0( GEOS_FMT(
    "Hybrid partition: {} cells, {} atomic vertices, {} graph edges; objective {:.3f} -> {:.3f}, "
    "cut faces {} -> {}, point replication {} -> {}, max neighbors {}, imbalance [{:.3f}%, {:.3f}%, {:.3f}%], "
    "build {:.3f}s, METIS {:.3f}s, refine {:.3f}s ({} moves), estimated peak {:.1f} MiB",
    result.metrics.numCells,
    result.metrics.numPartitionVertices,
    result.metrics.numGraphEdges,
    result.metrics.initialObjective,
    result.metrics.finalObjective,
    result.metrics.initialCutFaces,
    result.metrics.numCutFaces,
    result.metrics.initialPointReplication,
    result.metrics.pointReplication,
    result.metrics.maxRankNeighbors,
    100.0 * result.metrics.imbalanceByConstraint[0],
    100.0 * result.metrics.imbalanceByConstraint[1],
    100.0 * result.metrics.imbalanceByConstraint[2],
    result.metrics.graphBuildSeconds,
    result.metrics.initialPartitionSeconds,
    result.metrics.exactRefinementSeconds,
    result.metrics.refinementMoves,
    result.metrics.estimatedRootGraphPeakBytes / (1024.0 * 1024.0) ) );

  return result;
}

} // namespace vtk
} // namespace geos
