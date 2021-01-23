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
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkExtentTranslator.h>

// System includes
#include <unordered_set>
#include <sys/stat.h>



namespace geosx
{
namespace vtk
{


/*!
 * @brief Get the DataSet file path
 * @details the DataSet file path is the path to the .vtu mesh file
 * @param[in] er The ElementRegion
 * @param[in] time the time-step
 * @param[in] rank the rank to be written
 */
string GetDataSetFilePath( ElementRegionBase const & er, double time, int rank )
{
  return std::to_string( time ) + "/" + stringutilities::PadValue( rank, std::to_string( MpiWrapper::commSize() ).size() ) + "_" + er.getName() + ".vtu";
}

const std::map< string, int > geosx2VTKCellTypes =
{
  { "C3D4", VTK_TETRA },
  { "C3D8", VTK_HEXAHEDRON },
  { "C3D6", VTK_WEDGE },
  { "C3D5", VTK_PYRAMID }
};

/*!
 * @brief Gets the VTK cell identifier
 * @param[in] elementType the type of the element (using the abaqus nomenclature)
 * @return the VTK cell identifier
 */
int ToVTKCellType( const string & elementType )
{
  int vtkIdentifier = VTK_EMPTY_CELL;
  try
  {
    vtkIdentifier = geosx2VTKCellTypes.at( elementType );
  }
  catch( const std::out_of_range & outOfRange )
  {
    GEOSX_ERROR( "Element type " << elementType << " not recognized for VTK output " );
  }
  return vtkIdentifier;
}

void VTKPolyDataWriterInterface::gatherNbElementsInRegion( ElementRegionBase const & er,
                                                           array1d< localIndex > & nbElemsInRegion ) const
{
  localIndex nbElems = er.getNumberOfElements();
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  nbElemsInRegion.resize( mpiSize );
  MpiWrapper::gather( &nbElems, 1, nbElemsInRegion.data(), 1, 0, MPI_COMM_GEOSX );
}

VTKPolyDataWriterInterface::VTKPolyDataWriterInterface( string const & outputName ):
  m_outputFolder( outputName ),
  m_pvd( outputName + ".pvd" ),
  m_previousCycle( -1 )
{
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    mode_t mode = 0733;
    int errorCode = mkdir( outputName.c_str(), mode );
    if( errorCode == -1 )
    {
      if( errno == EEXIST )
      {
        GEOSX_WARNING( "Path " << outputName << " already exists from a previous simulation" );
      }
      else
      {
        GEOSX_ERROR( "Fail to create main directory " << outputName << " for the VTK output" );
      }
    }
  }
  MpiWrapper::barrier();
}

vtkSmartPointer< vtkPoints >  VTKPolyDataWriterInterface::getVtkPoints( NodeManager const & nodeManager ) const
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

std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >VTKPolyDataWriterInterface::getWell( WellElementSubRegion const & esr,
                                                                                                               NodeManager const & nodeManager ) const
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
    connectivity[1] = edge+1; // We can do that because of the pattern in which the wells are stored
    cellsArray->InsertNextCell( numberOfNodesPerElement, connectivity.data() );
  }

  if( esr.size() > 0 )
  {
    localIndex const lastPoint = esr.nodeList()[ esr.size() -1  ][1];
    auto point = referencePosition[lastPoint];
    points->SetPoint( esr.size(), point[0], point[1], point[2] );
  }

  return std::make_pair( points, cellsArray );
}
std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
VTKPolyDataWriterInterface::getSurface( FaceElementSubRegion const & esr,
                                        NodeManager const & nodeManager ) const
{
  // Get unique node set composing the surface
  auto & nodeListPerElement = esr.nodeList();
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( esr.size() );
  std::unordered_map< localIndex, localIndex > geosx2VTKIndexing;
  geosx2VTKIndexing.reserve( esr.size() * esr.numNodesPerElement() );
  localIndex nodeIndexInVTK = 0;
  std::vector< vtkIdType > connectivity( esr.numNodesPerElement() );
  std::vector< int >  vtkOrdering = esr.getVTKNodeOrdering();

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
VTKPolyDataWriterInterface::getEmbeddedSurface( EmbeddedSurfaceSubRegion const & esr,
                                                NodeManager const & nodeManager ) const
{
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & intersectionPoints = nodeManager.embSurfNodesPosition();

  points->SetNumberOfPoints( intersectionPoints.size( 0 ) );
  for( localIndex pointIndex = 0; pointIndex < intersectionPoints.size( 0 ); pointIndex++ )
  {
    points->SetPoint( pointIndex, intersectionPoints[pointIndex][0], intersectionPoints[pointIndex][1], intersectionPoints[pointIndex][2] );
  }

  EmbeddedSurfaceSubRegion::NodeMapType const & toNodesMap = esr.nodeList();
  std::vector< vtkIdType > connectivity( 10 );
  for( localIndex cellIndex = 0; cellIndex < esr.size(); cellIndex++ )
  {
    connectivity.resize( toNodesMap.sizeOfArray( cellIndex ) );
    for( std::size_t i = 0; i < connectivity.size(); ++i )
    {
      connectivity[i] = esr.nodeList( cellIndex, i );
    }
    cellsArray->InsertNextCell( connectivity.size(), connectivity.data() );
  }

  return std::make_pair( points, cellsArray );
}

std::pair< std::vector< int >, vtkSmartPointer< vtkCellArray > >
VTKPolyDataWriterInterface::getVtkCells( CellElementRegion const & er ) const
{
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( er.getNumberOfElements< CellElementRegion >() );
  std::vector< int > cellType;
  cellType.reserve( er.getNumberOfElements< CellElementRegion >() );
  er.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & esr )
  {
    std::vector< vtkIdType > connectivity( esr.numNodesPerElement() );
    std::vector< int > vtkOrdering = esr.getVTKNodeOrdering();
    int vtkCellType = ToVTKCellType( esr.getElementTypeString() );
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

void VTKPolyDataWriterInterface::writeField( WrapperBase const & wrapperBase,
                                             vtkSmartPointer< VTKGEOSXData > data,
                                             localIndex size, localIndex & count ) const
{
  std::type_info const & typeID = wrapperBase.getTypeId();
  if( typeID==typeid(r1_array) )
  {
    // We need a special case for the R1 array
    // Because those array are not stored the same way than the classical array2d which is
    // preferable to store a vector on each element.
    data->SetNumberOfComponents( 3 );
  }
  rtTypes::applyArrayTypeLambda2( rtTypes::typeID( typeID ),
                                  true,
                                  [&]( auto array, auto GEOSX_UNUSED_PARAM( Type ) )->void
  {
    typedef decltype( array ) arrayType;
    Wrapper< arrayType > const & wrapperT = Wrapper< arrayType >::cast( wrapperBase );
    traits::ViewTypeConst< arrayType > const sourceArray = wrapperT.reference().toViewConst();
    if( typeID!=typeid(r1_array) )
    {
      integer nbOfComponents = 1;
      for( localIndex i = 1; i < arrayType::NDIM; i++ )
      {
        nbOfComponents = nbOfComponents * sourceArray.size( i );
      }
      data->SetNumberOfComponents( nbOfComponents );
    }
    for( localIndex i = 0; i < size; i++ )
    {
      LvArray::forValuesInSlice( sourceArray[i], [&]( auto const & value )
      {
        data->customInsertValue( count++, value );
      } );
    }
  } );

}

void VTKPolyDataWriterInterface::writeNodeFields( vtkSmartPointer< vtkPointData > const pointdata,
                                                  NodeManager const & nodeManager ) const
{
  for( auto const & wrapperIter : nodeManager.wrappers() )
  {
    auto const & wrapper = *wrapperIter.second;
    if( wrapper.getPlotLevel() <= m_plotLevel )
    {
      vtkSmartPointer< VTKGEOSXData > data = VTKGEOSXData::New();
      data->SetNumberOfValues( nodeManager.size() );
      data->SetName( wrapper.getName().c_str() );;
      localIndex count = 0;
      writeField( wrapper, data, nodeManager.size(), count );
      pointdata->AddArray( data );
    }
  }
}

template< class SUBREGION >
void VTKPolyDataWriterInterface::writeElementFields( vtkSmartPointer< vtkCellData > const celldata,
                                                     ElementRegionBase const & er ) const
{
  std::unordered_set< string > allFields;
  er.forElementSubRegions< SUBREGION >( [&]( auto const & esr )
  {
    for( auto const & wrapperIter : esr.wrappers() )
    {
      auto const * const wrapper = wrapperIter.second;
      if( wrapper->getPlotLevel() <= m_plotLevel )
      {
        allFields.insert( wrapperIter.first );
      }
    }
  } );

  for( auto const & field : allFields )
  {
    vtkSmartPointer< VTKGEOSXData > data = VTKGEOSXData::New();
    data->SetNumberOfValues( er.getNumberOfElements< SUBREGION >() );
    data->SetName( field.c_str() );

    localIndex count = 0;
    er.forElementSubRegions< SUBREGION >( [&]( auto const & esr )
    {
      auto const & wrapper = *esr.getWrapperBase( field );
      writeField( wrapper, data, esr.size(), count );
    } );
    celldata->AddArray( data );
  }
}
void VTKPolyDataWriterInterface::writeCellElementRegions( real64 time,
                                                          ElementRegionManager const & elemManager,
                                                          NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion const & er )->void
  {
    if( er.getNumberOfElements< CellElementSubRegion >() != 0 )
    {
      vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
      auto VTKPoints = getVtkPoints( nodeManager );
      ug->SetPoints( VTKPoints );
      auto VTKCells = getVtkCells( er );
      ug->SetCells( VTKCells.first.data(), VTKCells.second );
      writeElementFields< CellElementSubRegion >( ug->GetCellData(), er );
      writeNodeFields( ug->GetPointData(), nodeManager );
      writeUnstructuredGrid( ug, time, er.getName() );
    }
  } );
}

void VTKPolyDataWriterInterface::writeWellElementRegions( real64 time, ElementRegionManager const & elemManager,
                                                          NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< WellElementRegion >( [&]( WellElementRegion const & er )->void
  {
    auto esr = er.getSubRegion( 0 )->groupCast< WellElementSubRegion const * >();
    vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
    auto VTKWell = getWell( *esr, nodeManager );
    ug->SetPoints( VTKWell.first );
    ug->SetCells( VTK_LINE, VTKWell.second );
    writeElementFields< WellElementSubRegion >( ug->GetCellData(), er );
    writeUnstructuredGrid( ug, time, er.getName() );
  } );
}

void VTKPolyDataWriterInterface::writeSurfaceElementRegions( real64 time,
                                                             ElementRegionManager const & elemManager,
                                                             NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & er )->void
  {
    vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
    if( er.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::embeddedElement )
    {
      auto esr = er.getSubRegion( 0 )->groupCast< EmbeddedSurfaceSubRegion const * >();

      auto VTKSurface = getEmbeddedSurface( *esr, nodeManager );
      ug->SetPoints( VTKSurface.first );
      ug->SetCells( VTK_POLYGON, VTKSurface.second );

      writeElementFields< EmbeddedSurfaceSubRegion >( ug->GetCellData(), er );
    }
    else if( er.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::faceElement )
    {
      auto esr = er.getSubRegion( 0 )->groupCast< FaceElementSubRegion const * >();

      auto VTKSurface = getSurface( *esr, nodeManager );

      ug->SetPoints( VTKSurface.first );
      if( esr->numNodesPerElement() == 8 )
      {
        ug->SetCells( VTK_HEXAHEDRON, VTKSurface.second );
      }
      else if( esr->numNodesPerElement() == 6 )
      {
        ug->SetCells( VTK_WEDGE, VTKSurface.second );
      }
      else
      {
        GEOSX_ERROR( "Elements with " << esr->numNodesPerElement() << " nodes can't be output "
                                      << "in the FaceElementRegion " << er.getName() );
      }
      writeElementFields< FaceElementSubRegion >( ug->GetCellData(), er );
    }
    writeUnstructuredGrid( ug, time, er.getName() );
  } );
}

void VTKPolyDataWriterInterface::writeVtmFile( real64 time,
                                               ElementRegionManager const & elemManager,
                                               VTKVTMWriter const & vtmWriter ) const
{
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  auto writeSubBlocks = [&]( ElementRegionBase const & er )
  {
    array1d< localIndex > nbElemsInRegion;
    gatherNbElementsInRegion( er, nbElemsInRegion );
    vtmWriter.addSubBlock( er.getCatalogName(), er.getName() );
    for( int i = 0; i < mpiSize; i++ )
    {
      string const dataSetFile = GetDataSetFilePath( er, time, i );
      if( nbElemsInRegion[i] > 0 )
      {
        if( mpiRank == 0 )
        {
          vtmWriter.addDataToSubBlock( er.getCatalogName(), er.getName(), dataSetFile, i );
        }
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

void VTKPolyDataWriterInterface::createTimeStepSubFolder( real64 time ) const
{
  int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    mode_t mode = 0773;
    string const timeStepSubFolder = getTimeStepSubFolder( time );
    int errorCode = mkdir( timeStepSubFolder.c_str(), mode );
    if( errorCode == -1 )
    {
      if( errno == EEXIST )
      {
        GEOSX_WARNING( "Path " << timeStepSubFolder << " already exists from a previous simulation" );
      }
      else
      {
        GEOSX_ERROR( "Fail to create the timestep directory " << timeStepSubFolder << " for the VTK output" );
      }
    }
  }
  MpiWrapper::barrier();
}

void VTKPolyDataWriterInterface::writeUnstructuredGrid( vtkSmartPointer< vtkUnstructuredGrid > ug,
                                                        double time,
                                                        string const & name ) const
{
  string timeStepSubFolder = VTKPolyDataWriterInterface::getTimeStepSubFolder( time );
  vtkSmartPointer< vtkXMLUnstructuredGridWriter > vtuWriter =vtkXMLUnstructuredGridWriter::New();
  vtuWriter->SetInputData( ug );
  string vtuFilePath = timeStepSubFolder + "/" +
                       stringutilities::PadValue( MpiWrapper::commRank(), std::to_string( MpiWrapper::commSize() ).size() ) +"_" + name + ".vtu";
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

string VTKPolyDataWriterInterface::getTimeStepSubFolder( real64 time ) const
{
  return m_outputFolder + "/" + std::to_string( time );
}

void VTKPolyDataWriterInterface::write( real64 time, integer cycle, DomainPartition const & domain )
{
  createTimeStepSubFolder( time );
  ElementRegionManager const & elemManager = *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  NodeManager const & nodeManager = *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getNodeManager();
  writeCellElementRegions( time, elemManager, nodeManager );
  writeWellElementRegions( time, elemManager, nodeManager );
  writeSurfaceElementRegions( time, elemManager, nodeManager );
  string vtmPath = getTimeStepSubFolder( time ) + ".vtm";
  VTKVTMWriter vtmWriter( vtmPath );
  writeVtmFile( time, elemManager, vtmWriter );
  if( cycle != m_previousCycle )
  {
    m_pvd.addData( time, vtmPath );
    m_pvd.save();
  }
  m_previousCycle = cycle;
}
}
}
