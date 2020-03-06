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

#include "Array2DVTKDataArray.hpp"
#include "CellVTKDataArray.hpp"
#include <vtkUnstructuredGrid.h>



namespace geosx
{
namespace vtk
{
  VTKPolyDataWriterInterface::VTKPolyDataWriterInterface( string const & outputName ) :
    m_writer( vtkXMLUnstructuredGridWriter::New() ),
    m_outputFolder( outputName )
  {
    mode_t mode = 0733;
    mkdir( outputName.c_str(), mode );
  }

  void VTKPolyDataWriterInterface::SetFileName( double time ) const
  {
    string fileName =  m_outputFolder + "/"+ m_outputFolder + "_"+ std::to_string( time ) + ".vtu";
    m_writer->SetFileName( fileName.c_str() );
  }

  void VTKPolyDataWriterInterface::LinkMesh( DomainPartition const * const domain ) const
  {
    auto VTKPoints = GetVTKPoints( domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager() );
    auto VTKCells = GetVTKCells( domain->getMeshBody(0)->getMeshLevel(0)->getElemManager()->GetRegion(0)->GetSubRegion(0)->group_cast< CellElementSubRegion const * const >() );
    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkUnstructuredGrid::New();
    ug->SetPoints( VTKPoints );
  //  ug->SetCells( VTK_HEXAHEDRON, VTKCells );
    m_writer->SetInputData( ug );
  }

  vtkSmartPointer< vtkPoints >  VTKPolyDataWriterInterface::GetVTKPoints( NodeManager  const * const nodeManager ) const
  {
    Array2DVTKDataArray<  const real64, nodes::REFERENCE_POSITION_USD > * 
      pointsArray = new Array2DVTKDataArray<  const real64, nodes::REFERENCE_POSITION_USD >( nodeManager->referencePosition() );
    vtkSmartPointer< vtkPoints > points = vtkPoints::New();
    points->SetData( pointsArray );
    return points;
  }

  vtkSmartPointer< vtkCellArray > VTKPolyDataWriterInterface::GetVTKCells( CellElementSubRegion const * const esr ) const
  {
    CellConnectivityVTKDataArray * connectivity = new CellConnectivityVTKDataArray( esr);
    CellOffsetVTKDataArray * offset = new CellOffsetVTKDataArray( esr);

    vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
    cellsArray->SetData( offset, connectivity );
    return cellsArray;
  }

  void VTKPolyDataWriterInterface::Write( real64 time, DomainPartition const * const domain  ) const
  {
    SetFileName( time );
    LinkMesh( domain );
    m_writer->SetCompressorTypeToNone();
    m_writer->SetDataModeToAscii();
    m_writer->Write();
  }
}
}
