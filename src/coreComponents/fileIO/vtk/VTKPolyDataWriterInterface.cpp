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

static string paddedRank( MPI_Comm const & comm, int const rank = -1 )
{
  int const width = LvArray::integerConversion< int >( std::to_string( MpiWrapper::commSize( comm ) ).size() );
  return stringutilities::padValue( rank >= 0 ? rank : MpiWrapper::commRank( comm ), width );
}

static int toVTKCellType( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Line:          return VTK_LINE;
    case ElementType::Triangle:      return VTK_TRIANGLE;
    case ElementType::Quadrilateral: return VTK_QUAD;
    case ElementType::Polygon:       return VTK_POLYGON;
    case ElementType::Tetrahedron:    return VTK_TETRA;
    case ElementType::Pyramid:       return VTK_PYRAMID;
    case ElementType::Prism:         return VTK_WEDGE;
    case ElementType::Hexahedron:    return VTK_HEXAHEDRON;
    case ElementType::Polyhedron:    return VTK_POLYHEDRON;
  }
  return VTK_EMPTY_CELL;
}

static std::vector< int > getVTKNodeOrdering( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Line:          return { 0, 1 };
    case ElementType::Triangle:      return { 0, 1, 2 };
    case ElementType::Quadrilateral: return { 0, 1, 2, 3 }; // TODO check
    case ElementType::Polygon:       return { 0, 1, 2, 3, 4, 5, 6, 7, 8 }; // TODO
    case ElementType::Tetrahedron:    return { 1, 0, 2, 3 };
    case ElementType::Pyramid:       return { 0, 3, 2, 1, 4, 0, 0, 0 };
    case ElementType::Prism:         return { 0, 4, 2, 1, 5, 3, 0, 0 };
    case ElementType::Hexahedron:    return { 0, 1, 3, 2, 4, 5, 7, 6 };
    case ElementType::Polyhedron:    return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }; // TODO
  }
  return {};
}

/**
 * @brief Ask rank @p rank for the number of elements in its ElementRegionBase @p er.
 * @param[in] region the element region for which we want to know the number of elements
 * @param[out] nbElemsInRegion output array
 * @return the number of elements in the region for the asked rank
 */
std::vector< localIndex >
gatherNbElementsInRegion( ElementRegionBase const & region,
                          MPI_Comm const & comm = MPI_COMM_GEOSX )
{
  localIndex const nbElems = region.getNumberOfElements();
  std::vector< localIndex > nbElemsInRegion( MpiWrapper::commSize( comm ) );
  MpiWrapper::gather( &nbElems, 1, nbElemsInRegion.data(), 1, 0, comm );
  return nbElemsInRegion;
}

/**
 * @brief Gets the vertices coordinates as a VTK Object for @p nodeManager
 * @param[in] nodeManager the NodeManager associated with the domain being written
 * @return a VTK object storing all nodes of the mesh
 */
vtkSmartPointer< vtkPoints >
getVtkPoints( NodeManager const & nodeManager )
{
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  points->SetNumberOfPoints( nodeManager.size() );
  auto const coord = nodeManager.referencePosition();
  for( localIndex v = 0; v < nodeManager.size(); v++ )
  {
    points->SetPoint( v, coord[v][0], coord[v][1], coord[v][2] );
  }
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
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  // if esr.size() == 0, we set the number of points and cells to zero
  // if not, we set the number of points to esr.size()+1 and the number of cells to esr.size()
  localIndex const numPoints = subRegion.size() > 0 ? subRegion.size() + 1 : 0;
  points->SetNumberOfPoints( numPoints );
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
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
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( subRegion.size() );
  std::unordered_map< localIndex, localIndex > geosx2VTKIndexing;
  geosx2VTKIndexing.reserve( subRegion.size() * subRegion.numNodesPerElement() );
  localIndex nodeIndexInVTK = 0;
  std::vector< vtkIdType > connectivity( subRegion.numNodesPerElement() );
  std::vector< int > vtkOrdering = getVTKNodeOrdering( subRegion.getElementType() );

  for( localIndex ei = 0; ei < subRegion.size(); ei++ )
  {
    auto const & elem = nodeListPerElement[ei];
    for( localIndex i = 0; i < elem.size(); i++ )
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

  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
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
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & intersectionPoints = nodeManager.referencePosition();

  points->SetNumberOfPoints( intersectionPoints.size( 0 ) );
  for( localIndex pointIndex = 0; pointIndex < intersectionPoints.size( 0 ); pointIndex++ )
  {
    points->SetPoint( pointIndex, intersectionPoints[pointIndex][0], intersectionPoints[pointIndex][1], intersectionPoints[pointIndex][2] );
  }

  EmbeddedSurfaceSubRegion::NodeMapType const & toNodesMap = subRegion.nodeList();
  array1d< vtkIdType > connectivity( 10 );
  for( localIndex cellIndex = 0; cellIndex < subRegion.size(); cellIndex++ )
  {
    connectivity.resize( toNodesMap.sizeOfArray( cellIndex ) );
    for( localIndex i = 0; i < connectivity.size(); ++i )
    {
      connectivity[i] = subRegion.nodeList( cellIndex, i );
    }
    cellsArray->InsertNextCell( connectivity.size(), connectivity.data() );
  }

  return std::make_pair( points, cellsArray );
}

/**
 * @brief Gets the cell connectivities as a VTK object for the CellElementRegion @p er
 * @param[in] region the CellElementRegion to be written
 * @return a pair, first value is a table with the same size than the total number of element in the CellElementRegion
 * contaning the type of the cells, the second value is a VTK object containing the connectivity information
 */
std::pair< std::vector< int >, vtkSmartPointer< vtkCellArray > >
getVtkCells( CellElementRegion const & region )
{
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( region.getNumberOfElements< CellElementRegion >() );
  std::vector< int > cellType;
  cellType.reserve( region.getNumberOfElements< CellElementRegion >() );
  region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    std::vector< vtkIdType > connectivity( subRegion.numNodesPerElement() );
    std::vector< int > vtkOrdering = getVTKNodeOrdering( subRegion.getElementType() );
    int vtkCellType = toVTKCellType( subRegion.getElementType() );
    for( localIndex c = 0; c < subRegion.size(); c++ )
    {
      for( std::size_t i = 0; i < connectivity.size(); i++ )
      {
        connectivity[i] = subRegion.nodeList( c, vtkOrdering[i] );
      }
      cellType.push_back( vtkCellType );
      cellsArray->InsertNextCell( subRegion.numNodesPerElement(), connectivity.data() );
    }
  } );
  return std::make_pair( cellType, cellsArray );
}

/**
 * @brief Writes timestamp information required by VisIt
 * @param[in] ug the VTK unstructured grid.
 * @param[in] time the current time-step
 */
void writeTimestamp( vtkUnstructuredGrid & ug,
                     real64 const time )
{
  vtkDoubleArray * t = vtkDoubleArray::New();
  t->SetName( "TIME" );
  t->SetNumberOfTuples( 1 );
  t->SetTuple1( 0, time );
  ug.GetFieldData()->AddArray( t );
}

/**
 * @brief Writes a field from \p wrapperBase
 * @details Sets the number of components, the number of value and fill the VTK data structure using
 * a wrapper around a field.
 * @param[in] wrapperBase a wrapper around the field to be written
 * @param[in] offset the cell index offset at which to start writing data (in case of multiple subregions)
 * @param[in,out] data a VTK data container, must be a vtkAOSDataArrayTemplate of the correct value type
 */
void writeField( WrapperBase const & wrapper,
                 localIndex const offset,
                 vtkDataArray & data )
{
  types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    using T = typename ArrayType::ValueType;
    vtkAOSDataArrayTemplate< T > & typedData = *vtkAOSDataArrayTemplate< T >::FastDownCast( &data );
    auto const sourceArray = Wrapper< ArrayType >::cast( wrapper ).reference().toViewConst();

    // TODO: check if parallel host policy is faster/slower
    forAll< serialPolicy >( sourceArray.size( 0 ), [sourceArray, offset, &typedData]( localIndex const i )
    {
      int compIndex = 0;
      LvArray::forValuesInSlice( sourceArray[i], [&]( T const & value )
      {
        typedData.SetTypedComponent( offset + i, compIndex++, value );
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
                           vtkAOSDataArrayTemplate< T > & data )
{
  data.SetNumberOfComponents( 1 );
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
                           vtkAOSDataArrayTemplate< T > & data )
{
  auto const view = wrapper.referenceAsView();
  data.SetNumberOfComponents( view.size( 1 ) );

  Span< string const > const labels = wrapper.getDimLabels( 1 );
  if( !labels.empty() )
  {
    checkLabels( wrapper, 1 );
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      data.SetComponentName( i, labels[i].c_str() );
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
                           vtkAOSDataArrayTemplate< T > & data )
{
  data.SetNumberOfComponents( wrapper.numArrayComp() );

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
    data.SetComponentName( compIndex++, makeComponentName( labels, idx_seq{}, indices ... ).c_str() );
  } );
}

template< class SUBREGION = Group >
void writeElementField( Group const & subRegions,
                        string const & field,
                        vtkCellData & cellData )
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
        data = typedData;
        setComponentMetadata( Wrapper< ArrayType >::cast( wrapper ), *typedData );
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
    writeField( wrapper, offset, *data );
    offset += subRegion.size();
  } );
  cellData.AddArray( data );
}

void VTKPolyDataWriterInterface::writeNodeFields( NodeManager const & nodeManager,
                                                  vtkPointData & pointData ) const
{
  for( auto const & wrapperIter : nodeManager.wrappers() )
  {
    auto const & wrapper = *wrapperIter.second;
    if( wrapper.getPlotLevel() <= m_plotLevel )
    {
      vtkSmartPointer< vtkDataArray > data;
      types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
      {
        using ArrayType = decltype( array );
        using T = typename ArrayType::ValueType;
        auto typedData = vtkAOSDataArrayTemplate< T >::New();
        data = typedData;
        setComponentMetadata( Wrapper< ArrayType >::cast( wrapper ), *typedData );
      } );

      data->SetNumberOfTuples( nodeManager.size() );
      data->SetName( wrapper.getName().c_str() );

      writeField( wrapper, 0, *data );
      pointData.AddArray( data );
    }
  }
}

template< class SUBREGION >
void VTKPolyDataWriterInterface::writeElementFields( ElementRegionBase const & region,
                                                     vtkCellData & cellData ) const
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
        if( wrapper.getPlotLevel() <= m_plotLevel )
        {
          string const fieldName = constitutive::ConstitutiveBase::makeFieldName( material.getName(), wrapper.getName() );
          subReg.registerWrapper( fieldName, wrapper.averageOverSecondDim( fieldName, subReg ) );
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
      if( wrapperIter.second->getPlotLevel() <= m_plotLevel && materialFields.count( wrapperIter.first ) == 0 )
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
    if( region.getNumberOfElements< CellElementSubRegion >() != 0 )
    {
      vtkSmartPointer< vtkUnstructuredGrid > const ug = vtkUnstructuredGrid::New();
      auto VTKPoints = getVtkPoints( nodeManager );
      ug->SetPoints( VTKPoints );
      auto VTKCells = getVtkCells( region );
      ug->SetCells( VTKCells.first.data(), VTKCells.second );
      writeTimestamp( *ug, time );
      writeElementFields< CellElementSubRegion >( region, *ug->GetCellData() );
      writeNodeFields( nodeManager, *ug->GetPointData() );
      writeUnstructuredGrid( cycle, region.getName(), *ug );
    }
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
    vtkSmartPointer< vtkUnstructuredGrid > const ug = vtkUnstructuredGrid::New();
    auto const VTKWell = getWell( subRegion, nodeManager );
    ug->SetPoints( VTKWell.first );
    ug->SetCells( VTK_LINE, VTKWell.second );
    writeTimestamp( *ug, time );
    writeElementFields< WellElementSubRegion >( region, *ug->GetCellData() );
    writeUnstructuredGrid( cycle, region.getName(), *ug );
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
    vtkSmartPointer< vtkUnstructuredGrid > const ug = vtkUnstructuredGrid::New();
    if( region.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::embeddedElement )
    {
      auto const & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

      auto const VTKSurface = getEmbeddedSurface( subRegion, embSurfNodeManager );
      ug->SetPoints( VTKSurface.first );
      ug->SetCells( VTK_POLYGON, VTKSurface.second );

      writeElementFields< EmbeddedSurfaceSubRegion >( region, *ug->GetCellData() );
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
      writeElementFields< FaceElementSubRegion >( region, *ug->GetCellData() );
    }
    writeTimestamp( *ug, time );
    writeUnstructuredGrid( cycle, region.getName(), *ug );
  } );
}


void VTKPolyDataWriterInterface::writeVtmFile( integer const cycle,
                                               ElementRegionManager const & elemManager,
                                               VTKVTMWriter const & vtmWriter ) const
{
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  auto addRegion = [&]( ElementRegionBase const & region )
  {
    if( mpiRank == 0 )
    {
      if ( !vtmWriter.hasBlock( region.getCatalogName() ) )
      {
        vtmWriter.addBlock( region.getCatalogName() );
      }
    }

    std::vector< localIndex > const nbElemsInRegion = gatherNbElementsInRegion( region, MPI_COMM_GEOSX );
    vtmWriter.addSubBlock( region.getCatalogName(), region.getName() );
    for( int i = 0; i < mpiSize; i++ )
    {
      if( mpiRank == 0 )
      {
        string const dataSetFile = GEOSX_FMT( "{:06d}/{}_{}.vtu", cycle, paddedRank( MPI_COMM_GEOSX, i ), region.getName() );
        vtmWriter.addDataToSubBlock( region.getCatalogName(), region.getName(), dataSetFile, i );
      }
    }
  };

  // Output each of the region types
  if ( ( m_outputRegionType == VTKRegionTypes::CELL ) || ( m_outputRegionType == VTKRegionTypes::ALL ) )
  {
    elemManager.forElementRegions< CellElementRegion >( addRegion );
  }

  if ( ( m_outputRegionType == VTKRegionTypes::WELL ) || ( m_outputRegionType == VTKRegionTypes::ALL ) )
  {
    elemManager.forElementRegions< WellElementRegion >( addRegion );
  }

  if ( ( m_outputRegionType == VTKRegionTypes::SURFACE ) || ( m_outputRegionType == VTKRegionTypes::ALL ) )
  {
    elemManager.forElementRegions< SurfaceElementRegion >( addRegion );
  }

  if( mpiRank == 0 )
  {
    vtmWriter.save();
  }
}

void VTKPolyDataWriterInterface::writeUnstructuredGrid( integer const cycle,
                                                        string const & name,
                                                        vtkUnstructuredGrid & ug ) const
{
  string const cycleSubFolder = joinPath( m_outputDir, getCycleSubFolder( cycle ) );
  vtkSmartPointer< vtkXMLUnstructuredGridWriter > const vtuWriter = vtkXMLUnstructuredGridWriter::New();
  vtuWriter->SetInputData( &ug );
  string const vtuFileName = paddedRank( MPI_COMM_GEOSX ) + "_" + name + ".vtu";
  string const vtuFilePath = joinPath( cycleSubFolder, vtuFileName );
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

string VTKPolyDataWriterInterface::getCycleSubFolder( integer const cycle ) const
{
  return joinPath( m_outputName, GEOSX_FMT( "{:06d}", cycle ) );
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

  string const stepSubFolder = getCycleSubFolder( cycle );
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
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

  string const vtmName = stepSubFolder + ".vtm";
  VTKVTMWriter vtmWriter( joinPath( m_outputDir, vtmName ) );
  writeVtmFile( cycle, elemManager, vtmWriter );

  if( cycle != m_previousCycle )
  {
    m_pvd.addData( time, vtmName );
    m_pvd.save();
  }
  m_previousCycle = cycle;
}

} // namespace vtk
} // namespace geosx
