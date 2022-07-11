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

// Source includes
#include "VTKPolyDataWriterInterface.hpp"

#include "common/Logger.hpp"
#include "common/TypeDispatch.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"

// TPL includes
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

// System includes
#include <unordered_set>

namespace geosx
{

using namespace dataRepository;

namespace vtk
{

VTKPolyDataWriterInterface::VTKPolyDataWriterInterface( string name ):
  m_outputDir( "." ),
  m_outputName( std::move( name ) ),
  m_pvd( m_outputName + ".pvd" ),
  m_plotLevel( PlotLevel::LEVEL_1 ),
  m_previousCycle( -1 ),
  m_outputMode( VTKOutputMode::BINARY ),
  m_outputRegionType( VTKRegionTypes::ALL )
{}

static int toVTKCellType( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Vertex:        return VTK_VERTEX;
    case ElementType::Line:          return VTK_LINE;
    case ElementType::Triangle:      return VTK_TRIANGLE;
    case ElementType::Quadrilateral: return VTK_QUAD;
    case ElementType::Polygon:       return VTK_POLYGON;
    case ElementType::Tetrahedron:   return VTK_TETRA;
    case ElementType::Pyramid:       return VTK_PYRAMID;
    case ElementType::Wedge:         return VTK_WEDGE;
    case ElementType::Hexahedron:    return VTK_HEXAHEDRON;
    case ElementType::Prism5:        return VTK_PENTAGONAL_PRISM;
    case ElementType::Prism6:        return VTK_HEXAGONAL_PRISM;
    case ElementType::Polyhedron:    return VTK_POLYHEDRON;
  }
  return VTK_EMPTY_CELL;
}

static std::vector< int > getVtkToGeosxNodeOrdering( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Vertex:        return { 0 };
    case ElementType::Line:          return { 0, 1 };
    case ElementType::Triangle:      return { 0, 1, 2 };
    case ElementType::Quadrilateral: return { 0, 1, 2, 3 }; // TODO check
    case ElementType::Polygon:       return { 0, 1, 2, 3, 4, 5, 6, 7, 8 }; // TODO
    case ElementType::Tetrahedron:   return { 0, 1, 2, 3 };
    case ElementType::Pyramid:       return { 0, 1, 3, 2, 4 };
    case ElementType::Wedge:         return { 0, 4, 2, 1, 5, 3 };
    case ElementType::Hexahedron:    return { 0, 1, 3, 2, 4, 5, 7, 6 };
    case ElementType::Prism5:        return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    case ElementType::Prism6:        return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    case ElementType::Polyhedron:    return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }; // TODO
  }
  return {};
}

/**
 * @brief Gets the vertices coordinates as a VTK Object for @p nodeManager
 * @param[in] nodeManager the NodeManager associated with the domain being written
 * @return a VTK object storing all nodes of the mesh
 */
vtkSmartPointer< vtkPoints >
getVtkPoints( NodeManager const & nodeManager,
              arrayView1d< localIndex const > const & nodeIndices )
{
  localIndex const numNodes = LvArray::integerConversion< localIndex >( nodeIndices.size() );
  auto points = vtkSmartPointer< vtkPoints >::New();
  points->SetNumberOfPoints( numNodes );
  auto const coord = nodeManager.referencePosition().toViewConst();
  forAll< parallelHostPolicy >( numNodes, [=, pts = points.GetPointer()]( localIndex const k )
  {
    localIndex const v = nodeIndices[k];
    pts->SetPoint( k, coord[v][0], coord[v][1], coord[v][2] );
  } );
  return points;
}

/**
 * @brief Gets the cell connectivities and the vertices coordinates as VTK objects for a specific WellElementSubRegion.
 * @param[in] subRegion the WellElementSubRegion to be output
 * @param[in] nodeManager the NodeManager associated with the DomainPartition being written.
 * @return a pair containing a VTKPoints (with the information on the vertices and their coordinates)
 * and a VTKCellArray (with the cell connectivities).
 */
std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
getWell( WellElementSubRegion const & subRegion,
         NodeManager const & nodeManager )
{
  // some notes about WellElementSubRegion:
  // - if the well represented by this subRegion is not on this rank, esr.size() = 0
  // - otherwise, esr.size() is equal to the number of well elements of the well on this rank
  // Each well element has two nodes, shared with the previous and next well elements, respectively
  auto points = vtkSmartPointer< vtkPoints >::New();
  // if esr.size() == 0, we set the number of points and cells to zero
  // if not, we set the number of points to esr.size()+1 and the number of cells to esr.size()
  localIndex const numPoints = subRegion.size() > 0 ? subRegion.size() + 1 : 0;
  points->SetNumberOfPoints( numPoints );
  auto cellsArray = vtkSmartPointer< vtkCellArray >::New();
  cellsArray->SetNumberOfCells( subRegion.size() );
  localIndex const numberOfNodesPerElement = subRegion.numNodesPerElement();
  GEOSX_ERROR_IF_NE( numberOfNodesPerElement, 2 );
  std::vector< vtkIdType > connectivity( numberOfNodesPerElement );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const referencePosition = nodeManager.referencePosition();

  // note that if esr.size() == 0, we don't have any point or cell to add below and we just return

  for( localIndex edge = 0; edge < subRegion.size(); edge++ )
  {
    localIndex const firstPoint = subRegion.nodeList()[edge][0];
    auto point = referencePosition[firstPoint];
    points->SetPoint( edge, point[0], point[1], point[2] );
    connectivity[0] = edge;
    connectivity[1] = edge + 1; // We can do that because of the pattern in which the wells are stored
    cellsArray->InsertNextCell( numberOfNodesPerElement, connectivity.data() );
  }

  if( subRegion.size() > 0 )
  {
    localIndex const lastPoint = subRegion.nodeList()[subRegion.size() - 1][1];
    auto point = referencePosition[lastPoint];
    points->SetPoint( subRegion.size(), point[0], point[1], point[2] );
  }

  return std::make_pair( points, cellsArray );
}

/**
 * @brief Gets the cell connectivities and the vertices coordinates as VTK objects for a specific FaceElementSubRegion.
 * @param[in] subRegion the FaceElementSubRegion to be output
 * @param[in] nodeManager the NodeManager associated with the DomainPartition being written.
 * @return a pair containing a VTKPoints (with the information on the vertices and their coordinates)
 * and a VTKCellArray (with the cell connectivities).
 */
std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
getSurface( FaceElementSubRegion const & subRegion,
            NodeManager const & nodeManager )
{
  // Get unique node set composing the surface
  auto & nodeListPerElement = subRegion.nodeList();
  auto cellsArray = vtkSmartPointer< vtkCellArray >::New();
  cellsArray->SetNumberOfCells( subRegion.size() );
  std::unordered_map< localIndex, localIndex > geosx2VTKIndexing;
  geosx2VTKIndexing.reserve( subRegion.size() * subRegion.numNodesPerElement() );
  localIndex nodeIndexInVTK = 0;
  std::vector< vtkIdType > connectivity( subRegion.numNodesPerElement() );
  std::vector< int > const vtkOrdering = getVtkToGeosxNodeOrdering( subRegion.getElementType() );

  for( localIndex ei = 0; ei < subRegion.size(); ei++ )
  {
    auto const & elem = nodeListPerElement[ei];
    for( localIndex i = 0; i < elem.size(); ++i )
    {
      auto const & VTKIndexPos = geosx2VTKIndexing.find( elem[vtkOrdering[i]] );
      if( VTKIndexPos == geosx2VTKIndexing.end() )
      {
        connectivity[i] = geosx2VTKIndexing[elem[vtkOrdering[i]]] = nodeIndexInVTK++;
      }
      else
      {
        connectivity[i] = VTKIndexPos->second;
      }
    }
    cellsArray->InsertNextCell( elem.size(), connectivity.data() );
  }

  auto points = vtkSmartPointer< vtkPoints >::New();
  points->SetNumberOfPoints( geosx2VTKIndexing.size() );
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const referencePosition = nodeManager.referencePosition();

  for( auto nodeIndex: geosx2VTKIndexing )
  {
    auto point = referencePosition[nodeIndex.first];
    points->SetPoint( nodeIndex.second, point[0], point[1], point[2] );
  }

  return std::make_pair( points, cellsArray );
}

/**
 * @brief Gets the cell connectivities and the vertices coordinates as VTK objects for a specific EmbeddedSurafaceSubRegion.
 * @param[in] subRegion the EmbeddedSurfaceSubRegion to be output
 * @param[in] elemManager the elemManager associated with the DomainPartition being written.
 * @param[in] nodeManager the NodeManager associated with the DomainPartition being written.
 * @param[in] edgeManager the edgeManager associated with the DomainPartition being written.
 * @return a pair containing a VTKPoints (with the information on the vertices and their coordinates)
 * and a VTKCellArray (with the cell connectivities).
 */
std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
getEmbeddedSurface( EmbeddedSurfaceSubRegion const & subRegion,
                    EmbeddedSurfaceNodeManager const & nodeManager )
{
  auto cellsArray = vtkSmartPointer< vtkCellArray >::New();
  auto points = vtkSmartPointer< vtkPoints >::New();

  localIndex const numNodes = nodeManager.size();
  auto const intersectionPoints = nodeManager.referencePosition();

  points->SetNumberOfPoints( numNodes );
  for( localIndex pointIndex = 0; pointIndex < numNodes; ++pointIndex )
  {
    auto const pointCoords = intersectionPoints[pointIndex];
    points->SetPoint( pointIndex, pointCoords[0], pointCoords[1], pointCoords[2] );
  }

  auto const toNodesMap = subRegion.nodeList().toViewConst();
  std::vector< vtkIdType > connectivity( 10 );
  for( localIndex cellIndex = 0; cellIndex < subRegion.size(); ++cellIndex )
  {
    auto const nodes = toNodesMap[cellIndex];
    connectivity.resize( nodes.size() );
    for( localIndex i = 0; i < nodes.size(); ++i )
    {
      connectivity[i] = nodes[i];
    }
    cellsArray->InsertNextCell( nodes.size(), connectivity.data() );
  }

  return std::make_pair( points, cellsArray );
}

struct CellData
{
  std::vector< int > cellTypes;
  vtkSmartPointer< vtkCellArray > cells;
  array1d< localIndex > nodes;
};

/**
 * @brief Gets the cell connectivities as a VTK object for the CellElementRegion @p region
 * @param[in] region the CellElementRegion to be written
 * @param[in] numNodes number of local nodes
 * @return a struct consisting of:
 *         - a list of types for each cell,
 *         - a VTK object containing the connectivity information
 *         - a list of relevant node indices in order in which they must be stored
 */
CellData getVtkCells( CellElementRegion const & region, localIndex const numNodes )
{
  localIndex const numElems = region.getNumberOfElements< CellElementSubRegion >();
  if( numElems == 0 )
  {
    return { {}, vtkSmartPointer< vtkCellArray >::New(), {} };
  }

  // 1. Mark (in parallel) relevant nodes
  std::vector< localIndex > newNodeIndices( numNodes ); // temporary use as a marker array
  region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    auto const nodeList = subRegion.nodeList().toViewConst();
    forAll< parallelHostPolicy >( subRegion.size(), [&, nodeList]( localIndex const cellIdx )
    {
      auto const nodes = nodeList[cellIdx];
      for( localIndex i = 0; i < nodes.size(); ++i )
      {
        // use atomic write to avoid technical UB
        RAJA::atomicExchange< parallelHostAtomic >( &newNodeIndices[nodes[i]], 1 );
      }
    } );
  } );

  // 2. Assign new node IDs (serial step)
  array1d< localIndex > relevantNodes;
  relevantNodes.reserve( numNodes );
  localIndex newNodeIdx = 0;
  for( localIndex nodeIdx = 0; nodeIdx < numNodes; ++nodeIdx )
  {
    if( newNodeIndices[nodeIdx] > 0 )
    {
      relevantNodes.emplace_back( nodeIdx );
      newNodeIndices[nodeIdx] = newNodeIdx++;
    }
  }

  // 3. Write connectivity using new node IDs
  localIndex const numConns = [&]
  {
    localIndex numConn = 0;
    region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
    {
      numConn += subRegion.size() * subRegion.numNodesPerElement();
    } );
    return numConn;
  }();

  std::vector< int > cellTypes;
  cellTypes.reserve( numElems );

  auto const offsets = vtkSmartPointer< vtkIdTypeArray >::New();
  offsets->SetNumberOfTuples( numElems + 1 );

  auto const connectivity = vtkSmartPointer< vtkIdTypeArray >::New();
  connectivity->SetNumberOfTuples( numConns );

  // 4. Write connectivity using new node IDs
  localIndex elemOffset = 0;
  localIndex connOffset = 0;
  region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    cellTypes.insert( cellTypes.end(), subRegion.size(), toVTKCellType( subRegion.getElementType() ) );
    std::vector< int > const vtkOrdering = getVtkToGeosxNodeOrdering( subRegion.getElementType() );
    localIndex const nodesPerElem = subRegion.numNodesPerElement();
    auto const nodeList = subRegion.nodeList().toViewConst();

    forAll< parallelHostPolicy >( subRegion.size(), [=, &connectivity, &offsets]( localIndex const c )
    {
      localIndex const elemConnOffset = connOffset + c * nodesPerElem;
      auto const nodes = nodeList[c];
      for( localIndex i = 0; i < nodesPerElem; ++i )
      {
        connectivity->SetTypedComponent( elemConnOffset + i, 0, newNodeIndices[nodes[vtkOrdering[i]]] );
      }
      offsets->SetTypedComponent( elemOffset + c, 0, elemConnOffset );
    } );

    elemOffset += subRegion.size();
    connOffset += subRegion.size() * nodesPerElem;
  } );
  offsets->SetTypedComponent( elemOffset, 0, connOffset );

  auto cellsArray = vtkSmartPointer< vtkCellArray >::New();
  cellsArray->SetData( offsets, connectivity );

  return { std::move( cellTypes ), cellsArray, std::move( relevantNodes ) };
}

/**
 * @brief Writes timestamp information required by VisIt
 * @param[in] ug the VTK unstructured grid.
 * @param[in] time the current time-step
 */
void writeTimestamp( vtkUnstructuredGrid * ug,
                     real64 const time )
{
  auto t = vtkSmartPointer< vtkDoubleArray >::New();
  t->SetName( "TIME" );
  t->SetNumberOfTuples( 1 );
  t->SetTuple1( 0, time );
  ug->GetFieldData()->AddArray( t );
}

/**
 * @brief Writes a field from @p wrapper.
 * @param[in] wrapper a wrapper around the field to be written
 * @param[in] offset the cell index offset at which to start writing data (in case of multiple subregions)
 * @param[in,out] data a VTK data container, must be a vtkAOSDataArrayTemplate of the correct value type
 */
void writeField( WrapperBase const & wrapper,
                 localIndex const offset,
                 vtkDataArray * data )
{
  types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    using T = typename ArrayType::ValueType;
    vtkAOSDataArrayTemplate< T > * typedData = vtkAOSDataArrayTemplate< T >::FastDownCast( data );
    auto const sourceArray = Wrapper< ArrayType >::cast( wrapper ).reference().toViewConst();

    forAll< parallelHostPolicy >( sourceArray.size( 0 ), [sourceArray, offset, typedData]( localIndex const i )
    {
      LvArray::forValuesInSlice( sourceArray[i], [&, compIndex = 0]( T const & value ) mutable
      {
        typedData->SetTypedComponent( offset + i, compIndex++, value );
      } );
    } );
  } );
}

/**
 * @brief Writes a field from @p wrapper using an index list.
 * @param[in] wrapper a wrapper around the field to be written
 * @param[in] indices a list of indices into @p wrapper that will be written
 * @param[in] offset the cell index offset at which to start writing data (in case of multiple subregions)
 * @param[in,out] data a VTK data container, must be a vtkAOSDataArrayTemplate of the correct value type
 */
void writeField( WrapperBase const & wrapper,
                 arrayView1d< localIndex const > const & indices,
                 localIndex const offset,
                 vtkDataArray * data )
{
  types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    using T = typename ArrayType::ValueType;
    vtkAOSDataArrayTemplate< T > * typedData = vtkAOSDataArrayTemplate< T >::FastDownCast( data );
    auto const sourceArray = Wrapper< ArrayType >::cast( wrapper ).reference().toViewConst();

    forAll< parallelHostPolicy >( indices.size(), [=]( localIndex const i )
    {
      LvArray::forValuesInSlice( sourceArray[indices[i]], [&, compIndex = 0]( T const & value ) mutable
      {
        typedData->SetTypedComponent( offset + i, compIndex++, value );
      } );
    } );
  } );
}

/**
 * @brief Build/expand a list of strings used as default component labels on demand.
 * @param size number of labels requested
 * @return a span over range of strings (stored permanently in memory)
 */
Span< string const > getDefaultLabels( localIndex const size )
{
  static std::vector< string > labels;
  localIndex oldSize = LvArray::integerConversion< localIndex >( labels.size() );
  std::generate_n( std::back_inserter( labels ), size - oldSize, [&] { return std::to_string( oldSize++ ); } );
  return { labels.begin(), labels.begin() + size };
}

/**
 * @brief Checks consistency of user-provided per-dimension labels (must match the size of array).
 * @param wrapper the array wrapper
 * @param dim dimension index to check
 */
template< typename T, int NDIM, typename PERM >
void checkLabels( Wrapper< Array< T, NDIM, PERM > > const & wrapper, int const dim )
{
  GEOSX_ERROR_IF_NE_MSG( LvArray::integerConversion< localIndex >( wrapper.getDimLabels( dim ).size() ),
                         wrapper.reference().size( dim ),
                         "VTK writer: component names are set, but don't match the array size.\n"
                         "This is likely a bug in physics module (solver or constitutive model)." );
}

/**
 * @brief Get a list of component labels for a particular dimension of an array.
 * @param wrapper array wrapper
 * @param dim dimension index
 * @return a span over range of strings representing labels
 */
template< typename T, int NDIM, typename PERM >
Span< string const > getDimLabels( Wrapper< Array< T, NDIM, PERM > > const & wrapper,
                                   int const dim )
{
  Span< string const > const labels = wrapper.getDimLabels( dim );
  if( labels.empty() )
  {
    return getDefaultLabels( wrapper.reference().size( dim ) );
  }
  checkLabels( wrapper, dim );
  return labels;
}

/**
 * @brief Build a multidimensional component name out of distinct dimension-wise labels.
 * @tparam Ts types of indices
 * @tparam Is dummy template argument required for positional expansion of the pack
 * @param dimLabels per-dimension component labels
 * @param indices per-dimension component indices
 * @return combined component name
 */
template< typename ... Ts, integer ... Is >
string makeComponentName( Span< string const >(&dimLabels)[sizeof...(Ts)],
                          std::integer_sequence< integer, Is... >,
                          Ts const & ... indices )
{
  return stringutilities::concat( '/', dimLabels[Is][indices] ... );
}

/**
 * @brief Specialized component metadata handler for 1D arrays.
 * @param wrapper GEOSX typed wrapper over source array
 * @param data VTK typed data array
 */
template< typename T, typename PERM >
void setComponentMetadata( Wrapper< Array< T, 1, PERM > > const & GEOSX_UNUSED_PARAM( wrapper ),
                           vtkAOSDataArrayTemplate< T > * data )
{
  data->SetNumberOfComponents( 1 );
}

/**
 * @brief Specialized component metadata handler for 2D arrays.
 * @param wrapper GEOSX typed wrapper over source array
 * @param data VTK typed data array
 *
 * This exists because we want to keep default VTK handling for unlabeled components
 * (i.e. X/Y/Z for 1-3 components, numeric indices for higher) for the time being.
 * This function can be removed if we force each physics package to always set its labels.
 */
template< typename T, typename PERM >
void setComponentMetadata( Wrapper< Array< T, 2, PERM > > const & wrapper,
                           vtkAOSDataArrayTemplate< T > * data )
{
  auto const view = wrapper.referenceAsView();
  data->SetNumberOfComponents( view.size( 1 ) );

  Span< string const > const labels = wrapper.getDimLabels( 1 );
  if( !labels.empty() )
  {
    checkLabels( wrapper, 1 );
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      data->SetComponentName( i, labels[i].c_str() );
    }
  }
}

/**
 * @brief Produces a temporary array slice from a view that can be looped over.
 * @param view the source view
 * @return a fake slice that does not point to real data but has correct dims/strides.
 * @note The slice is only valid as long as the @p view is in scope.
 *       Values in the slice may be uninitialized and should not be used.
 */
template< typename T, int NDIM, int USD >
ArraySlice< T const, NDIM - 1, USD - 1 >
makeTemporarySlice( ArrayView< T const, NDIM, USD > const & view )
{
  // The following works in all compilers, but technically invokes undefined behavior:
  // return ArraySlice< T, NDIM - 1, USD - 1 >( nullptr, view.dims() + 1, view.strides() + 1 );
  localIndex const numComp = LvArray::indexing::multiplyAll< NDIM - 1 >( view.dims() + 1 );
  static array1d< T > arr;
  arr.template resizeWithoutInitializationOrDestruction( numComp );
  return ArraySlice< T const, NDIM - 1, USD - 1 >( arr.data(), view.dims() + 1, view.strides() + 1 );
}

/**
 * @brief Generic component metadata handler for multidimensional arrays.
 * @param wrapper GEOSX typed wrapper over source array
 * @param data VTK typed data array
 */
template< typename T, int NDIM, typename PERM >
void setComponentMetadata( Wrapper< Array< T, NDIM, PERM > > const & wrapper,
                           vtkAOSDataArrayTemplate< T > * data )
{
  data->SetNumberOfComponents( wrapper.numArrayComp() );

  Span< string const > labels[NDIM-1];
  for( integer dim = 1; dim < NDIM; ++dim )
  {
    labels[dim-1] = getDimLabels( wrapper, dim );
  }

  auto const view = wrapper.referenceAsView();
  auto const slice = view.size( 0 ) > 0 ? view[0] : makeTemporarySlice( view );

  integer compIndex = 0;
  LvArray::forValuesInSliceWithIndices( slice, [&]( T const &, auto const ... indices )
  {
    using idx_seq = std::make_integer_sequence< integer, sizeof...(indices) >;
    data->SetComponentName( compIndex++, makeComponentName( labels, idx_seq{}, indices ... ).c_str() );
  } );
}

template< class SUBREGION = Group >
void writeElementField( Group const & subRegions,
                        string const & field,
                        vtkCellData * cellData )
{
  // instantiate vtk array of the correct type
  vtkSmartPointer< vtkDataArray > data;
  localIndex numElements = 0;
  bool first = true;
  int numDims = 0;
  subRegions.forSubGroups< SUBREGION >( [&]( SUBREGION const & subRegion )
  {
    numElements += subRegion.size();
    WrapperBase const & wrapper = subRegion.getWrapperBase( field );
    if( first )
    {
      types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
      {
        using ArrayType = decltype( array );
        using T = typename ArrayType::ValueType;
        auto typedData = vtkAOSDataArrayTemplate< T >::New();
        data.TakeReference( typedData );
        setComponentMetadata( Wrapper< ArrayType >::cast( wrapper ), typedData );
      } );
      first = false;
      numDims = wrapper.numArrayDims();
    }
    else
    {
      // Sanity check
      GEOSX_ERROR_IF_NE_MSG( wrapper.numArrayDims(), numDims,
                             "VTK writer: sanity check failed for " << field << " (inconsistent array dimensions)" );
      GEOSX_ERROR_IF_NE_MSG( wrapper.numArrayComp(), data->GetNumberOfComponents(),
                             "VTK writer: sanity check failed for " << field << " (inconsistent array sizes)" );
    }
  } );

  data->SetNumberOfTuples( numElements );
  data->SetName( field.c_str() );

  // write each subregion in turn, keeping track of element offset
  localIndex offset = 0;
  subRegions.forSubGroups< SUBREGION >( [&]( SUBREGION const & subRegion )
  {
    WrapperBase const & wrapper = subRegion.getWrapperBase( field );
    writeField( wrapper, offset, data.GetPointer() );
    offset += subRegion.size();
  } );
  cellData->AddArray( data );
}

void VTKPolyDataWriterInterface::writeNodeFields( NodeManager const & nodeManager,
                                                  arrayView1d< localIndex const > const & nodeIndices,
                                                  vtkPointData * pointData ) const
{
  for( auto const & wrapperIter : nodeManager.wrappers() )
  {
    auto const & wrapper = *wrapperIter.second;
    if( isFieldPlotEnabled( wrapper ) )
    {
      vtkSmartPointer< vtkDataArray > data;
      types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
      {
        using ArrayType = decltype( array );
        using T = typename ArrayType::ValueType;
        auto typedData = vtkAOSDataArrayTemplate< T >::New();
        data.TakeReference( typedData );
        setComponentMetadata( Wrapper< ArrayType >::cast( wrapper ), typedData );
      } );

      data->SetNumberOfTuples( nodeIndices.size() );
      data->SetName( wrapper.getName().c_str() );

      writeField( wrapper, nodeIndices, 0, data.GetPointer() );
      pointData->AddArray( data );
    }
  }
}

template< class SUBREGION >
void VTKPolyDataWriterInterface::writeElementFields( ElementRegionBase const & region,
                                                     vtkCellData * cellData ) const
{
  std::unordered_set< string > materialFields;
  conduit::Node fakeRoot;
  Group materialData( "averagedMaterialData", fakeRoot );
  region.forElementSubRegions< SUBREGION >( [&]( SUBREGION const & subRegion )
  {
    // Register a dummy group for each subregion
    Group & subReg = materialData.registerGroup( subRegion.getName() );
    subReg.resize( subRegion.size() );

    // Collect a list of plotted constitutive fields and create wrappers containing averaged data
    subRegion.getConstitutiveModels().forSubGroups( [&]( Group const & material )
    {
      material.forWrappers( [&]( WrapperBase const & wrapper )
      {
        string const fieldName = constitutive::ConstitutiveBase::makeFieldName( material.getName(), wrapper.getName() );
        if( isFieldPlotEnabled( wrapper.getPlotLevel(), fieldName ) )
        {
          subReg.registerWrapper( wrapper.averageOverSecondDim( fieldName, subReg ) );
          materialFields.insert( fieldName );
        }
      } );
    } );
  } );

  // Write averaged material data
  for( string const & field : materialFields )
  {
    writeElementField( materialData, field, cellData );
  }

  // Collect a list of regular fields (filter out material field wrappers)
  // TODO: this can be removed if we stop hanging constitutive wrappers on the mesh
  std::unordered_set< string > regularFields;
  region.forElementSubRegions< SUBREGION >( [&]( ElementSubRegionBase const & subRegion )
  {
    for( auto const & wrapperIter : subRegion.wrappers() )
    {
      if( isFieldPlotEnabled( *wrapperIter.second ) && materialFields.count( wrapperIter.first ) == 0 )
      {
        regularFields.insert( wrapperIter.first );
      }
    }
  } );

  // Write regular fields
  for( string const & field : regularFields )
  {
    writeElementField< SUBREGION >( region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() ), field, cellData );
  }
}

void VTKPolyDataWriterInterface::writeCellElementRegions( real64 const time,
                                                          integer const cycle,
                                                          ElementRegionManager const & elemManager,
                                                          NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion const & region )
  {
    CellData VTKCells = getVtkCells( region, nodeManager.size() );
    vtkSmartPointer< vtkPoints > const VTKPoints = getVtkPoints( nodeManager, VTKCells.nodes );

    auto const ug = vtkSmartPointer< vtkUnstructuredGrid >::New();
    ug->SetCells( VTKCells.cellTypes.data(), VTKCells.cells );
    ug->SetPoints( VTKPoints );

    writeTimestamp( ug.GetPointer(), time );
    writeElementFields< CellElementSubRegion >( region, ug->GetCellData() );
    writeNodeFields( nodeManager, VTKCells.nodes, ug->GetPointData() );
    writeUnstructuredGrid( cycle, region.getName(), ug.GetPointer() );
  } );
}

void VTKPolyDataWriterInterface::writeWellElementRegions( real64 const time,
                                                          integer const cycle,
                                                          ElementRegionManager const & elemManager,
                                                          NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< WellElementRegion >( [&]( WellElementRegion const & region )
  {
    auto const & subRegion = region.getSubRegion< WellElementSubRegion >( 0 );
    auto const ug = vtkSmartPointer< vtkUnstructuredGrid >::New();
    auto const VTKWell = getWell( subRegion, nodeManager );
    ug->SetPoints( VTKWell.first );
    ug->SetCells( VTK_LINE, VTKWell.second );
    writeTimestamp( ug.GetPointer(), time );
    writeElementFields< WellElementSubRegion >( region, ug->GetCellData() );
    writeUnstructuredGrid( cycle, region.getName(), ug.GetPointer() );
  } );
}

void VTKPolyDataWriterInterface::writeSurfaceElementRegions( real64 const time,
                                                             integer const cycle,
                                                             ElementRegionManager const & elemManager,
                                                             NodeManager const & nodeManager,
                                                             EmbeddedSurfaceNodeManager const & embSurfNodeManager ) const
{
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & region )
  {
    auto const ug = vtkSmartPointer< vtkUnstructuredGrid >::New();
    if( region.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::embeddedElement )
    {
      auto const & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

      auto const VTKSurface = getEmbeddedSurface( subRegion, embSurfNodeManager );
      ug->SetPoints( VTKSurface.first );
      ug->SetCells( VTK_POLYGON, VTKSurface.second );

      writeElementFields< EmbeddedSurfaceSubRegion >( region, ug->GetCellData() );
    }
    else if( region.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::faceElement )
    {
      auto const & subRegion = region.getSubRegion< FaceElementSubRegion >( 0 );

      auto const VTKSurface = getSurface( subRegion, nodeManager );
      ug->SetPoints( VTKSurface.first );

      if( subRegion.numNodesPerElement() == 8 )
      {
        ug->SetCells( VTK_HEXAHEDRON, VTKSurface.second );
      }
      else if( subRegion.numNodesPerElement() == 6 )
      {
        ug->SetCells( VTK_WEDGE, VTKSurface.second );
      }
      else
      {
        GEOSX_ERROR( "Elements with " << subRegion.numNodesPerElement() << " nodes can't be output "
                                      << "in the FaceElementRegion " << region.getName() );
      }
      writeElementFields< FaceElementSubRegion >( region, ug->GetCellData() );
    }
    writeTimestamp( ug.GetPointer(), time );
    writeUnstructuredGrid( cycle, region.getName(), ug.GetPointer() );
  } );
}

static string getCycleSubFolder( integer const cycle )
{
  return GEOSX_FMT( "{:06d}", cycle );
}

static string getRegionFileName( integer const rank, string const & regionName )
{
  int const width = static_cast< int >( std::log10( MpiWrapper::commSize() ) ) + 1;
  return GEOSX_FMT( "{:>0{}}_{}.vtu", rank, width, regionName );
}

void VTKPolyDataWriterInterface::writeVtmFile( integer const cycle,
                                               ElementRegionManager const & elemManager,
                                               VTKVTMWriter const & vtmWriter ) const
{
  GEOSX_ASSERT_EQ_MSG( MpiWrapper::commRank(), 0, "Must only be called on rank 0" );

  int const mpiSize = MpiWrapper::commSize();
  auto addRegion = [&]( ElementRegionBase const & region )
  {
    if( !vtmWriter.hasBlock( region.getCatalogName() ) )
    {
      vtmWriter.addBlock( region.getCatalogName() );
    }
    vtmWriter.addSubBlock( region.getCatalogName(), region.getName() );
    for( int i = 0; i < mpiSize; i++ )
    {
      string const dataSetFile = joinPath( getCycleSubFolder( cycle ), getRegionFileName( i, region.getName() ) );
      vtmWriter.addDataToSubBlock( region.getCatalogName(), region.getName(), dataSetFile, i );
    }
  };

  // Output each of the region types
  if( ( m_outputRegionType == VTKRegionTypes::CELL ) || ( m_outputRegionType == VTKRegionTypes::ALL ) )
  {
    elemManager.forElementRegions< CellElementRegion >( addRegion );
  }

  if( ( m_outputRegionType == VTKRegionTypes::WELL ) || ( m_outputRegionType == VTKRegionTypes::ALL ) )
  {
    elemManager.forElementRegions< WellElementRegion >( addRegion );
  }

  if( ( m_outputRegionType == VTKRegionTypes::SURFACE ) || ( m_outputRegionType == VTKRegionTypes::ALL ) )
  {
    elemManager.forElementRegions< SurfaceElementRegion >( addRegion );
  }

  vtmWriter.save();
}

void VTKPolyDataWriterInterface::writeUnstructuredGrid( integer const cycle,
                                                        string const & name,
                                                        vtkUnstructuredGrid * ug ) const
{
  auto const vtuWriter = vtkSmartPointer< vtkXMLUnstructuredGridWriter >::New();
  vtuWriter->SetInputData( ug );
  string const vtuFilePath = joinPath( m_outputDir,
                                       m_outputName,
                                       getCycleSubFolder( cycle ),
                                       getRegionFileName( MpiWrapper::commRank(), name ) );
  vtuWriter->SetFileName( vtuFilePath.c_str() );
  if( m_outputMode == VTKOutputMode::BINARY )
  {
    vtuWriter->SetDataModeToBinary();
  }
  else if( m_outputMode == VTKOutputMode::ASCII )
  {
    vtuWriter->SetDataModeToAscii();
  }
  vtuWriter->Write();
}

void VTKPolyDataWriterInterface::write( real64 const time,
                                        integer const cycle,
                                        DomainPartition const & domain )
{
  // This guard prevents crashes observed on MacOS due to a floating point exception
  // triggered inside VTK by a progress indicator
#if defined(__APPLE__) && defined(__MACH__)
  LvArray::system::FloatingPointExceptionGuard guard;
#endif

  string const stepSubFolder = joinPath( m_outputName, getCycleSubFolder( cycle ) );
  int const rank = MpiWrapper::commRank();
  if( rank == 0 )
  {
    if( m_previousCycle == -1 )
    {
      makeDirsForPath( joinPath( m_outputDir, m_outputName ) );
    }
    makeDirectory( joinPath( m_outputDir, stepSubFolder ) );
  }
  MpiWrapper::barrier( MPI_COMM_GEOSX );

  ElementRegionManager const & elemManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getElemManager();
  NodeManager const & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getNodeManager();
  EmbeddedSurfaceNodeManager const & embSurfNodeManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getEmbSurfNodeManager();
  writeCellElementRegions( time, cycle, elemManager, nodeManager );
  writeWellElementRegions( time, cycle, elemManager, nodeManager );
  writeSurfaceElementRegions( time, cycle, elemManager, nodeManager, embSurfNodeManager );

  if( rank == 0 )
  {
    string const vtmName = stepSubFolder + ".vtm";
    VTKVTMWriter vtmWriter( joinPath( m_outputDir, vtmName ) );
    writeVtmFile( cycle, elemManager, vtmWriter );

    if( cycle != m_previousCycle )
    {
      m_pvd.addData( time, vtmName );
      m_pvd.save();
    }
  }

  m_previousCycle = cycle;
}

bool VTKPolyDataWriterInterface::isFieldPlotEnabled( PlotLevel const wrapperPlotLevel,
                                                     string const & wrapperName ) const
{
  // check if the logLevel is sufficient high for plotting
  bool plotEnabled = ( wrapperPlotLevel <= m_plotLevel );

  // override the logLevel if the fieldNames list was provided
  if( !m_fieldNames.empty() )
  {
    auto search = m_fieldNames.find( wrapperName );
    plotEnabled = ( search != m_fieldNames.end() );
  }
  return plotEnabled;
}


} // namespace vtk
} // namespace geosx
