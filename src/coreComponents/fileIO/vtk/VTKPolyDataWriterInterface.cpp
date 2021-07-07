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

// Source includes
#include "VTKPolyDataWriterInterface.hpp"

#include "dataRepository/Group.hpp"

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
  m_outputMode( VTKOutputMode::BINARY )
{}

string paddedRank( MPI_Comm const & comm, int const rank = -1 )
{
  int const width = LvArray::integerConversion< int >( std::to_string( MpiWrapper::commSize( comm ) ).size() );
  return stringutilities::padValue( rank >= 0 ? rank : MpiWrapper::commRank( comm ), width );
}

/**
 * @brief Gets the VTK cell identifier
 * @param[in] elementType the type of the element (using the abaqus nomenclature)
 * @return the VTK cell identifier
 */
int toVTKCellType( const string & elementType )
{
  static const std::map< string, int > geosx2VTKCellTypes =
  {
    { "C3D4", VTK_TETRA },
    { "C3D8", VTK_HEXAHEDRON },
    { "C3D6", VTK_WEDGE },
    { "C3D5", VTK_PYRAMID }
  };

  GEOSX_THROW_IF( geosx2VTKCellTypes.count( elementType ) == 0,
                  "Element type not recognized for VTK output: " << elementType,
                  std::runtime_error );
  return geosx2VTKCellTypes.at( elementType );
}

/**
 * @brief Ask rank @p rank for the number of elements in its ElementRegionBase @p er.
 * @param[in] er the element region for which we want to know the number of elements
 * @param[out] nbElemsInRegion output array
 * @return the number of elements in the region for the asked rank
 */
std::vector< localIndex >
gatherNbElementsInRegion( ElementRegionBase const & er,
                          MPI_Comm const & comm = MPI_COMM_GEOSX )
{
  localIndex const nbElems = er.getNumberOfElements();
  std::vector< localIndex > nbElemsInRegion( MpiWrapper::commSize( comm ) );
  MpiWrapper::gather( &nbElems, 1, nbElemsInRegion.data(), 1, 0, comm );
  return nbElemsInRegion;
}

vtkSmartPointer< vtkPoints >
getVtkPoints( NodeManager const & nodeManager )
{
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  points->SetNumberOfPoints( nodeManager.size() );
  auto connectivity = nodeManager.referencePosition();
  for( localIndex v = 0; v < nodeManager.size(); v++ )
  {
    points->SetPoint( v, connectivity[v][0], connectivity[v][1], connectivity[v][2] );
  }
  return points;
}

std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
getWell( WellElementSubRegion const & esr,
         NodeManager const & nodeManager )
{
  // some notes about WellElementSubRegion:
  // - if the well represented by this subRegion is not on this rank, esr.size() = 0
  // - otherwise, esr.size() is equal to the number of well elements of the well on this rank
  // Each well element has two nodes, shared with the previous and next well elements, respectively
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  // if esr.size() == 0, we set the number of points and cells to zero
  // if not, we set the number of points to esr.size()+1 and the number of cells to esr.size()
  localIndex const numPoints = esr.size() > 0 ? esr.size() + 1 : 0;
  points->SetNumberOfPoints( numPoints );
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( esr.size() );
  localIndex const numberOfNodesPerElement = esr.numNodesPerElement();
  GEOSX_ERROR_IF_NE( numberOfNodesPerElement, 2 );
  std::vector< vtkIdType > connectivity( numberOfNodesPerElement );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const referencePosition = nodeManager.referencePosition();

  // note that if esr.size() == 0, we don't have any point or cell to add below and we just return

  for( localIndex edge = 0; edge < esr.size(); edge++ )
  {
    localIndex const firstPoint = esr.nodeList()[edge][0];
    auto point = referencePosition[firstPoint];
    points->SetPoint( edge, point[0], point[1], point[2] );
    connectivity[0] = edge;
    connectivity[1] = edge + 1; // We can do that because of the pattern in which the wells are stored
    cellsArray->InsertNextCell( numberOfNodesPerElement, connectivity.data() );
  }

  if( esr.size() > 0 )
  {
    localIndex const lastPoint = esr.nodeList()[esr.size() - 1][1];
    auto point = referencePosition[lastPoint];
    points->SetPoint( esr.size(), point[0], point[1], point[2] );
  }

  return std::make_pair( points, cellsArray );
}

std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
getSurface( FaceElementSubRegion const & esr,
            NodeManager const & nodeManager )
{
  // Get unique node set composing the surface
  auto & nodeListPerElement = esr.nodeList();
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( esr.size() );
  std::unordered_map< localIndex, localIndex > geosx2VTKIndexing;
  geosx2VTKIndexing.reserve( esr.size() * esr.numNodesPerElement() );
  localIndex nodeIndexInVTK = 0;
  std::vector< vtkIdType > connectivity( esr.numNodesPerElement() );
  std::vector< int > vtkOrdering = esr.getVTKNodeOrdering();

  for( localIndex ei = 0; ei < esr.size(); ei++ )
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

std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
getEmbeddedSurface( EmbeddedSurfaceSubRegion const & esr,
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

  EmbeddedSurfaceSubRegion::NodeMapType const & toNodesMap = esr.nodeList();
  array1d< vtkIdType > connectivity( 10 );
  for( localIndex cellIndex = 0; cellIndex < esr.size(); cellIndex++ )
  {
    connectivity.resize( toNodesMap.sizeOfArray( cellIndex ) );
    for( localIndex i = 0; i < connectivity.size(); ++i )
    {
      connectivity[i] = esr.nodeList( cellIndex, i );
    }
    cellsArray->InsertNextCell( connectivity.size(), connectivity.data() );
  }

  return std::make_pair( points, cellsArray );
}

std::pair< std::vector< int >, vtkSmartPointer< vtkCellArray > >
getVtkCells( CellElementRegion const & er )
{
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( er.getNumberOfElements< CellElementRegion >() );
  std::vector< int > cellType;
  cellType.reserve( er.getNumberOfElements< CellElementRegion >() );
  er.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & esr )
  {
    std::vector< vtkIdType > connectivity( esr.numNodesPerElement() );
    std::vector< int > vtkOrdering = esr.getVTKNodeOrdering();
    int vtkCellType = toVTKCellType( esr.getElementTypeString() );
    for( localIndex c = 0; c < esr.size(); c++ )
    {
      for( std::size_t i = 0; i < connectivity.size(); i++ )
      {
        connectivity[i] = esr.nodeList( c, vtkOrdering[i] );
      }
      cellType.push_back( vtkCellType );
      cellsArray->InsertNextCell( esr.numNodesPerElement(), connectivity.data() );
    }
  } );
  return std::make_pair( cellType, cellsArray );
}

void writeTimestamp( vtkUnstructuredGrid & ug,
                     real64 const time )
{
  vtkDoubleArray * t = vtkDoubleArray::New();
  t->SetName( "TIME" );
  t->SetNumberOfTuples( 1 );
  t->SetTuple1( 0, time );
  ug.GetFieldData()->AddArray( t );
}

void writeField( WrapperBase const & wrapper,
                 localIndex const offset,
                 vtkDataArray & data )
{
  rtTypes::applyArrayTypeLambda2( rtTypes::typeID( wrapper.getTypeId() ), true, [&]( auto array, auto scalar )
  {
    using T = decltype( scalar );
    vtkAOSDataArrayTemplate< T > & typedData = *vtkAOSDataArrayTemplate< T >::FastDownCast( &data );
    auto const sourceArray = wrapper.cast< decltype( array ) >().reference().toViewConst();

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
  GEOSX_ERROR_IF_NE_MSG( wrapper.getDimLabels( dim ).size(), wrapper.reference().size( dim ),
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
 * @brief Generic component metadata handler for multidimensional arrays.
 * @param wrapper GEOSX typed wrapper over source array
 * @param data VTK typed data array
 */
template< typename T, int NDIM, typename PERM >
void setComponentMetadata( Wrapper< Array< T, NDIM, PERM > > const & wrapper,
                           vtkAOSDataArrayTemplate< T > & data )
{
  auto const view = wrapper.referenceAsView();
  data.SetNumberOfComponents( view.size() / view.size( 0 ) );

  Span< string const > labels[NDIM-1];
  for( integer dim = 1; dim < NDIM; ++dim )
  {
    labels[dim-1] = getDimLabels( wrapper, dim );
  }

  integer compIndex = 0;
  LvArray::forValuesInSliceWithIndices( view[0], [&]( T const &, auto const ... indices )
  {
    using idx_seq = std::make_integer_sequence< integer, sizeof...(indices) >;
    data.SetComponentName( compIndex++, makeComponentName( labels, idx_seq{}, indices ... ).c_str() );
  } );
}

template< class SUBREGION = Group >
void writeElementField( vtkCellData & cellData,
                        Group const & subRegions,
                        string const & field )
{
  // instantiate vtk array of the correct type
  vtkSmartPointer< vtkDataArray > data;
  localIndex numElements = 0;
  bool first = true;
  int numDims = 0;
  subRegions.forSubGroups< SUBREGION >( [&]( SUBREGION const & esr )
  {
    numElements += esr.size();
    WrapperBase const & wrapper = esr.getWrapperBase( field );
    if( first )
    {
      rtTypes::applyArrayTypeLambda2( rtTypes::typeID( wrapper.getTypeId() ), true, [&]( auto array, auto scalar )
      {
        auto typedData = vtkAOSDataArrayTemplate< decltype( scalar ) >::New();
        data = typedData;
        setComponentMetadata( wrapper.cast< decltype( array ) >(), *typedData );
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
  subRegions.forSubGroups< SUBREGION >( [&]( SUBREGION const & esr )
  {
    WrapperBase const & wrapper = esr.getWrapperBase( field );
    writeField( wrapper, offset, *data );
    offset += esr.size();
  } );
  cellData.AddArray( data );
}

void VTKPolyDataWriterInterface::writeNodeFields( vtkPointData & pointData,
                                                  NodeManager const & nodeManager ) const
{
  for( auto const & wrapperIter : nodeManager.wrappers() )
  {
    auto const & wrapper = *wrapperIter.second;
    if( wrapper.getPlotLevel() <= m_plotLevel )
    {
      vtkSmartPointer< vtkDataArray > data;
      rtTypes::applyArrayTypeLambda2( rtTypes::typeID( wrapper.getTypeId() ), true, [&]( auto array, auto scalar )
      {
        auto typedData = vtkAOSDataArrayTemplate< decltype( scalar ) >::New();
        data = typedData;
        setComponentMetadata( wrapper.cast< decltype( array ) >(), *typedData );
      } );

      data->SetNumberOfTuples( nodeManager.size() );
      data->SetName( wrapper.getName().c_str() );

      writeField( wrapper, 0, *data );
      pointData.AddArray( data );
    }
  }
}

template< class SUBREGION >
void VTKPolyDataWriterInterface::writeElementFields( vtkCellData & cellData,
                                                     ElementRegionBase const & er ) const
{
  std::unordered_set< string > regularFields;
  std::unordered_set< string > materialFields;
  er.forElementSubRegions< SUBREGION >( [&]( SUBREGION const & esr )
  {
    for( auto const & wrapperIter : esr.wrappers() )
    {
      string const & name = wrapperIter.first;
      auto const & wrapper = *wrapperIter.second;
      if( wrapper.getPlotLevel() <= m_plotLevel )
      {
        // Slightly cursed: split fields into regular and material fields based on whether wrapper name begins
        // with a known material name. I can't think of an easier way at the moment, but there ought to be one.
        auto const startsWithMaterial = [&]( string const & mat ){ return name.substr( 0, mat.size() ) == mat; };
        if( std::find_if( er.getMaterialList().begin(), er.getMaterialList().end(), startsWithMaterial ) != er.getMaterialList().end() )
        {
          materialFields.insert( name );
        }
        else
        {
          regularFields.insert( name );
        }
      }
    }
  } );

  // Just write regular fields normally
  for( string const & field : regularFields )
  {
    writeElementField< SUBREGION >( cellData, er.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() ), field );
  }

  // Very cursed: create a fake group to host averaged-over-gauss-point wrappers of material fields
  conduit::Node fakeRoot;
  Group materialData( "averagedMaterialData", fakeRoot );
  er.forElementSubRegions< SUBREGION >( [&]( SUBREGION const & esr )
  {
    // Can't register type SUBREGION because of SurfaceElementSubRegion's use of getParent()
    Group & subReg = materialData.registerGroup( esr.getName() );
    subReg.resize( esr.size() );
    for( string const & field : materialFields )
    {
      subReg.registerWrapper( field, esr.getWrapperBase( field ).averageOverSecondDim( field, subReg ) );
    }
  } );

  // Write averaged material data
  for( string const & field : materialFields )
  {
    writeElementField( cellData, materialData, field );
  }
}

void VTKPolyDataWriterInterface::writeCellElementRegions( real64 const time,
                                                          ElementRegionManager const & elemManager,
                                                          NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion const & er )
  {
    if( er.getNumberOfElements< CellElementSubRegion >() != 0 )
    {
      vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
      auto VTKPoints = getVtkPoints( nodeManager );
      ug->SetPoints( VTKPoints );
      auto VTKCells = getVtkCells( er );
      ug->SetCells( VTKCells.first.data(), VTKCells.second );
      writeTimestamp( *ug, time );
      writeElementFields< CellElementSubRegion >( *ug->GetCellData(), er );
      writeNodeFields( *ug->GetPointData(), nodeManager );
      writeUnstructuredGrid( *ug, time, er.getName() );
    }
  } );
}

void VTKPolyDataWriterInterface::writeWellElementRegions( real64 const time,
                                                          ElementRegionManager const & elemManager,
                                                          NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< WellElementRegion >( [&]( WellElementRegion const & er )
  {
    WellElementSubRegion const & esr = dynamicCast< WellElementSubRegion const & >( er.getSubRegion( 0 ) );
    vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
    auto VTKWell = getWell( esr, nodeManager );
    ug->SetPoints( VTKWell.first );
    ug->SetCells( VTK_LINE, VTKWell.second );
    writeElementFields< WellElementSubRegion >( *ug->GetCellData(), er );
    writeUnstructuredGrid( *ug, time, er.getName() );
  } );
}

void VTKPolyDataWriterInterface::writeSurfaceElementRegions( real64 const time,
                                                             ElementRegionManager const & elemManager,
                                                             NodeManager const & nodeManager,
                                                             EmbeddedSurfaceNodeManager const & embSurfNodeManager ) const
{
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & er )
  {
    vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
    if( er.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::embeddedElement )
    {
      EmbeddedSurfaceSubRegion const & esr = dynamicCast< EmbeddedSurfaceSubRegion const & >( er.getSubRegion( 0 ) );

      auto VTKSurface = getEmbeddedSurface( esr, embSurfNodeManager );
      ug->SetPoints( VTKSurface.first );
      ug->SetCells( VTK_POLYGON, VTKSurface.second );

      writeElementFields< EmbeddedSurfaceSubRegion >( *ug->GetCellData(), er );
    }
    else if( er.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::faceElement )
    {
      FaceElementSubRegion const & esr = dynamicCast< FaceElementSubRegion const & >( er.getSubRegion( 0 ) );

      auto VTKSurface = getSurface( esr, nodeManager );

      ug->SetPoints( VTKSurface.first );
      if( esr.numNodesPerElement() == 8 )
      {
        ug->SetCells( VTK_HEXAHEDRON, VTKSurface.second );
      }
      else if( esr.numNodesPerElement() == 6 )
      {
        ug->SetCells( VTK_WEDGE, VTKSurface.second );
      }
      else
      {
        GEOSX_ERROR( "Elements with " << esr.numNodesPerElement() << " nodes can't be output "
                                      << "in the FaceElementRegion " << er.getName() );
      }
      writeElementFields< FaceElementSubRegion >( *ug->GetCellData(), er );
    }
    writeUnstructuredGrid( *ug, time, er.getName() );
  } );
}

/**
 * @brief Writes a VTM file for the time-step \p time.
 * @details a VTM file is a VTK Multiblock file. It contains relative path to different files organized in blocks.
 * @param[in] time the time-step
 * @param[in] elemManager the ElementRegionManager containing all the regions to be output and referred to in the VTM file
 * @param[in] vtmWriter a writer specialized for the VTM file format
 */
void writeVtmFile( real64 const time,
                   ElementRegionManager const & elemManager,
                   VTKVTMWriter const & vtmWriter )
{
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  auto writeSubBlocks = [&]( ElementRegionBase const & er )
  {
    std::vector< localIndex > const nbElemsInRegion = gatherNbElementsInRegion( er, MPI_COMM_GEOSX );
    vtmWriter.addSubBlock( er.getCatalogName(), er.getName() );
    for( int i = 0; i < mpiSize; i++ )
    {
      if( mpiRank == 0 && nbElemsInRegion[i] > 0 )
      {
        string const dataSetFile = std::to_string( time ) + "/" + paddedRank( MPI_COMM_GEOSX, i ) + "_" + er.getName() + ".vtu";
        vtmWriter.addDataToSubBlock( er.getCatalogName(), er.getName(), dataSetFile, i );
      }
    }
  };
  if( mpiRank == 0 )
  {
    vtmWriter.addBlock( CellElementRegion::catalogName() );
    vtmWriter.addBlock( WellElementRegion::catalogName() );
    vtmWriter.addBlock( SurfaceElementRegion::catalogName() );
  }

  elemManager.forElementRegions< CellElementRegion >( writeSubBlocks );
  elemManager.forElementRegions< WellElementRegion >( writeSubBlocks );
  elemManager.forElementRegions< SurfaceElementRegion >( writeSubBlocks );

  if( mpiRank == 0 )
  {
    vtmWriter.save();
  }
}

void VTKPolyDataWriterInterface::writeUnstructuredGrid( vtkUnstructuredGrid & ug,
                                                        real64 const time,
                                                        string const & name ) const
{
  string const timeStepSubFolder = joinPath( m_outputDir, VTKPolyDataWriterInterface::getTimeStepSubFolder( time ) );
  vtkSmartPointer< vtkXMLUnstructuredGridWriter > vtuWriter = vtkXMLUnstructuredGridWriter::New();
  vtuWriter->SetInputData( &ug );
  string const vtuFileName = paddedRank( MPI_COMM_GEOSX ) + "_" + name + ".vtu";
  string const vtuFilePath = joinPath( timeStepSubFolder, vtuFileName );
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

string VTKPolyDataWriterInterface::getTimeStepSubFolder( real64 const time ) const
{
  return joinPath( m_outputName, std::to_string( time ) );
}

void VTKPolyDataWriterInterface::write( real64 const time,
                                        integer const cycle,
                                        DomainPartition const & domain )
{
  string const stepSubFolder = getTimeStepSubFolder( time );
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
  writeCellElementRegions( time, elemManager, nodeManager );
  writeWellElementRegions( time, elemManager, nodeManager );
  writeSurfaceElementRegions( time, elemManager, nodeManager, embSurfNodeManager );

  string const vtmName = stepSubFolder + ".vtm";
  VTKVTMWriter vtmWriter( joinPath( m_outputDir, vtmName ) );
  writeVtmFile( time, elemManager, vtmWriter );

  if( cycle != m_previousCycle )
  {
    m_pvd.addData( time, vtmName );
    m_pvd.save();
  }
  m_previousCycle = cycle;
}

} // namespace vtk
} // namespace geosx
