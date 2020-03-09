/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "VTKMultiBlockDataWriterInterface.hpp"

#include "dataRepository/WrapperBase.hpp"

#include "Array2DVTKDataArray.hpp"
#include "CellVTKDataArray.hpp"

#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>



namespace geosx
{
namespace vtk
{
  VTKPolyDataWriterInterface::VTKPolyDataWriterInterface( string const & outputName ) :
    m_outputFolder( outputName ),
    m_pvd( outputName )
  {
    mode_t mode = 0733;
    mkdir( outputName.c_str(), mode );
  }

  vtkSmartPointer< vtkPoints >  VTKPolyDataWriterInterface::GetVTKPoints( NodeManager  const * const nodeManager ) const
  {
    vtkSmartPointer< vtkPoints > points = vtkPoints::New();
    points->SetNumberOfPoints( nodeManager->size() );
    auto connectivity = nodeManager->referencePosition();
    for( localIndex v = 0 ; v < nodeManager->size(); v++)
    {
      points->SetPoint( v, connectivity[v][0], connectivity[v][1], connectivity[v][2] );
    }
    return points;
  }

  std::tuple< vtkSmartPointer< vtkPoints >,  vtkSmartPointer< vtkCellArray > >VTKPolyDataWriterInterface::GetWell( WellElementSubRegion  const * const esr, NodeManager const * const nodeManager ) const
  {
    vtkSmartPointer< vtkPoints > points = vtkPoints::New();
    points->SetNumberOfPoints( esr->size() + 1 );
    vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
    cellsArray->SetNumberOfCells( esr->size() );
    for( localIndex edge = 0 ; edge < esr->size(); edge++)
    {
      localIndex firstPoint = esr->nodeList()[edge][0];
      auto point = nodeManager->referencePosition()[firstPoint];
      points->SetPoint( edge, point[0], point[1], point[2] );
      std::vector< vtkIdType > connectivity(esr->numNodesPerElement() );
      connectivity[0] = edge;
      connectivity[1] = edge+1;
      cellsArray->InsertNextCell( 2, connectivity.data() );
    }
    if( esr->size() > 0 )
    {
      localIndex lastPoint = esr->nodeList()[ esr->size() -1  ][1];
      auto point = nodeManager->referencePosition()[lastPoint];
      points->SetPoint( esr->size(), point[0], point[1], point[2] );
    }
    return std::make_tuple( points, cellsArray );
  }

  std::tuple< vtkSmartPointer< vtkCellArray >,  std::vector< int> > VTKPolyDataWriterInterface::GetVTKCells( CellElementRegion const * const er ) const
  {
    vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
    cellsArray->SetNumberOfCells( er->getNumberOfElements< CellElementRegion >() );
    std::vector< int > cellType;
    cellType.reserve( er->getNumberOfElements< CellElementRegion >() );
    er->forElementSubRegions<CellElementSubRegion>([&](CellElementSubRegion const * const esr )->void
    {
      for( localIndex c = 0; c < esr->size(); c++ )
      {
        std::vector< vtkIdType > connectivity(esr->numNodesPerElement() );
        if( esr->GetElementTypeString() == "C3D8" )
        {
          connectivity[0] = esr->nodeList(c,0);
          connectivity[1] = esr->nodeList(c,1);
          connectivity[2] = esr->nodeList(c,3);
          connectivity[3] = esr->nodeList(c,2);
          connectivity[4] = esr->nodeList(c,4);
          connectivity[5] = esr->nodeList(c,5);
          connectivity[6] = esr->nodeList(c,7);
          connectivity[7] = esr->nodeList(c,6);
        }
        cellType.push_back( VTK_HEXAHEDRON );
        cellsArray->InsertNextCell( esr->numNodesPerElement(), connectivity.data() );
      }
    });
    return std::make_tuple(cellsArray, cellType);
  }

  void VTKPolyDataWriterInterface::SetCellFields( vtkSmartPointer< vtkCellData > & celldata, CellElementRegion const * const er ) const
  {
      for( auto const & wrapperIter : er->wrappers() )
      {
        auto const * const wrapper = wrapperIter.second;

        if( wrapper->getPlotLevel() < m_plotLevel )
        {
          // the field name is the key to the map
          string const fieldName = wrapper->getName();
          /*
          std::type_info const & typeID = wrapper->get_typeid();
          rtTypes::TypeIDs fieldType = rtTypes::typeID(wrapper->get_typeid());
          if( !geosxToVTKTypeMap.count( typeID ) )
            continue;
          int dimension = 0;
          if( fieldType == rtTypes::TypeIDs::r1_array_id )
          {
            dimension = 3;
          }
          else
          {
            dimension = 1;
          }
          cellFields.insert(std::make_tuple(fieldName, geosxToVTKTypeMap.at( typeID ), dimension, fieldType) );
          */
        }
      }
  }
  void VTKPolyDataWriterInterface::WriteMeshFiles( double time, DomainPartition const * const domain ) const
  {
    string timeStepSubFolder = VTKPolyDataWriterInterface::GetTimeStepSubFolder( time );
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    ElementRegionManager const * const  elemManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
    elemManager->forElementRegions<CellElementRegion>([&](CellElementRegion const * const er)->void
    {
      vtkSmartPointer< vtkUnstructuredGrid > vtkUg = vtkUnstructuredGrid::New();
      vtkSmartPointer< vtkXMLUnstructuredGridWriter > vtkUgWriter = vtkXMLUnstructuredGridWriter::New();
      auto VTKPoints = GetVTKPoints( domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager() );
      vtkUg->SetPoints( VTKPoints );
      auto VTKCells = GetVTKCells( er );
      vtkUg->SetCells( std::get<1>(VTKCells).data(), std::get<0>(VTKCells) );
      string vtuFilePath = timeStepSubFolder + "/" + std::to_string( mpiRank) +"_" + er->getName() + ".vtu";
      vtkUgWriter->SetFileName( vtuFilePath.c_str() );
      vtkUgWriter->SetInputData( vtkUg );
      vtkUgWriter->SetDataModeToAscii();
      vtkUgWriter->Write();
    });
  }

  void VTKPolyDataWriterInterface::WriteWellFiles( double time, DomainPartition const * const domain ) const
  {
    string timeStepSubFolder = VTKPolyDataWriterInterface::GetTimeStepSubFolder( time );
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    ElementRegionManager const * const  elemManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
    elemManager->forElementRegions<WellElementRegion>([&](WellElementRegion const * const er)->void
    {
      auto esr = er->GetSubRegion(0)->group_cast<WellElementSubRegion const *>();
      auto VTKWell = GetWell( esr, domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager() );
      vtkSmartPointer<vtkUnstructuredGrid> ug = vtkUnstructuredGrid::New();
      ug->SetPoints(std::get<0>(VTKWell));
      ug->SetCells(VTK_LINE, std::get<1>(VTKWell));
      vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriter =vtkXMLUnstructuredGridWriter::New();
      vtuWriter->SetInputData( ug );
      string vtuFilePath = timeStepSubFolder + "/" + std::to_string( mpiRank) +"_" + er->getName() + ".vtu";
      vtuWriter->SetFileName( vtuFilePath.c_str() );
      vtuWriter->SetDataModeToAscii();
      vtuWriter->Write();
    });
  }

  void VTKPolyDataWriterInterface::WriteVTMFile( double time, DomainPartition const * const domain, VTKVTMWriter const& vtmWriter ) const
  {
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    int const mpiSize = MpiWrapper::Comm_size(MPI_COMM_GEOSX);
    if( mpiRank == 0 )
    {
      // Cells
      vtmWriter.AddBlock( CellElementRegion::CatalogName() );
      ElementRegionManager const * const  elemManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
      elemManager->forElementRegions<CellElementRegion>([&](CellElementRegion const * const er)->void
      {

        vtmWriter.AddSubBlock( CellElementRegion::CatalogName(), er->getName() );
        for(int i = 0; i < mpiSize; i++ )
        {
          vtmWriter.AddDataToSubBlock(CellElementRegion::CatalogName(), er->getName(), std::to_string(time) + "/" + std::to_string( i ) + "_" + er->getName() + ".vtu", i );
        }
      });
      
      // Wells
      vtmWriter.AddBlock( WellElementRegion::CatalogName() );
      elemManager->forElementRegions<WellElementRegion>([&](WellElementRegion const * const er)->void
      {
        vtmWriter.AddSubBlock( WellElementRegion::CatalogName(), er->getName() );
        for(int i = 0; i < mpiSize; i++ )
        {
          vtmWriter.AddDataToSubBlock( WellElementRegion::CatalogName(), er->getName(), std::to_string(time) + "/" + std::to_string( i ) +"_" + er->getName() + ".vtu", i );
        }
      });

      vtmWriter.Save();
    }
  }

  void VTKPolyDataWriterInterface::CreateTimeStepSubFolder( double time ) const
  {
    int const mpiRank = MpiWrapper::Comm_rank(MPI_COMM_GEOSX);
    if( mpiRank == 0 )
    {
      mode_t mode = 0773;
      string timeStepSubFolder = m_outputFolder + "/" + std::to_string( time );
      mkdir( timeStepSubFolder.c_str(), mode );
    }
    MpiWrapper::Barrier();
  }

  string VTKPolyDataWriterInterface::GetTimeStepSubFolder( double time ) const
  {
    return  m_outputFolder + "/" + std::to_string( time );
  }

  void VTKPolyDataWriterInterface::Write( real64 time, DomainPartition const * const domain  ) const
  {
    CreateTimeStepSubFolder( time );
    string vtmPath = m_outputFolder + "/" + std::to_string( time ) + ".vtm";
    VTKVTMWriter vtmWriter( vtmPath );
    WriteMeshFiles( time, domain);
    WriteWellFiles( time, domain);
    WriteVTMFile( time, domain, vtmWriter );
    m_pvd.AddData( time, vtmPath );
    m_pvd.Save();
  }
}
}
