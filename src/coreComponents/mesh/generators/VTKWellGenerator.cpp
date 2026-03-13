/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file VTKWellGenerator.cpp
 *
 */

#include "VTKWellGenerator.hpp"

#include "mesh/LogLevelsInfo.hpp"
#include "mesh/generators/VTKUtilities.hpp"
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPolyLine.h>
#include <vtkCellArrayIterator.h>

namespace geos
{
using namespace dataRepository;

VTKWellGenerator::VTKWellGenerator( string const & name, Group * const parent ):
  WellGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Path to the well file" );

  addLogLevel< logInfo::VTKSteps >();
}

void VTKWellGenerator::fillPolylineDataStructure( )
{
  GEOS_MARK_FUNCTION;

  vtkSmartPointer< vtkMultiProcessController > controller = vtk::getController();
  vtkMultiProcessController::SetGlobalController( controller );

  GEOS_LOG_RANK_0( GEOS_FMT( "{} '{}': reading well from {}", catalogName(), getName(), m_filePath ) );
  {
    GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, "  reading the dataset..." );
    vtk::AllMeshes allMeshes = vtk::loadAllMeshes( m_filePath, "main", string_array());
    vtkSmartPointer< vtkDataSet > loadedMesh = allMeshes.getMainMesh();
    controller->Broadcast( loadedMesh, 0 );

    vtkSmartPointer< vtkPolyData > polyData = vtkPolyData::SafeDownCast( loadedMesh );

    // load points
    vtkPoints * points = polyData->GetPoints();
    m_polyNodeCoords.resize( points->GetNumberOfPoints(), 3 );
    globalIndex ipoint = 0;
    for( vtkIdType c = 0; c < points->GetNumberOfPoints(); ++c, ++ipoint )
    {
      real64 point[3];
      points->GetPoint( c, point );
      LvArray::tensorOps::copy< 3 >( m_polyNodeCoords[ipoint], point );
    }

    GEOS_ERROR_IF( polyData->GetLines()->GetNumberOfCells() == 0,
                   GEOS_FMT( "Error! Your VTK file {} doesn't contain any well",
                             m_filePath ), this->getDataContext());

    // load edges
    vtkNew< vtkIdList > cellPts;
    auto iter = ::vtk::TakeSmartPointer( polyData->GetLines()->NewIterator() );
    iter->GetCellAtId( 0, cellPts );

    const globalIndex numCells = polyData->GetLines()->GetNumberOfCells();

    if( numCells > 1 && cellPts->GetNumberOfIds() == 2 )
    {
      // The well is stored as individual Line segments (2 pts each) rather than a single PolyLine.
      // Rebuild the full ordered point list by iterating all cells.
      GEOS_LOG_RANK_0( GEOS_FMT( "{}: VTK file {} contains {} individual Line segments. Concatenating into a single polyline.",
                                  this->getName(), m_filePath, numCells ) );
      vtkNew< vtkIdList > segPts;
      cellPts->Reset();
      vtkNew< vtkIdList > firstCell;
      iter->GetCellAtId( 0, firstCell );
      cellPts->InsertNextId( firstCell->GetId( 0 ) );
      for( vtkIdType cellId = 0; cellId < numCells; ++cellId )
      {
        iter->GetCellAtId( cellId, segPts );
        cellPts->InsertNextId( segPts->GetId( 1 ) );
      }
    }
    else
    {
      GEOS_LOG_RANK_0_IF( numCells > 1,
                          GEOS_FMT( "{}: Warning! Your VTK file {} contains {} PolyLine cells. Only the first one will be read.",
                                    this->getName(), m_filePath, numCells ) );
    }

    const globalIndex nbSegments = cellPts->GetNumberOfIds() - 1;
    m_segmentToPolyNodeMap.resizeDimension< 0 >( nbSegments );
    m_segmentToPolyNodeMap.resizeDimension< 1 >( m_numNodesPerElem );

    globalIndex iseg = 0;
    for( vtkIdType pointId = 0; pointId < nbSegments; ++pointId )
    {
      m_segmentToPolyNodeMap[iseg][0] = cellPts->GetId( pointId );
      m_segmentToPolyNodeMap[iseg][1] = cellPts->GetId( pointId + 1 );
      ++iseg;
    }
  }
}

REGISTER_CATALOG_ENTRY( MeshComponentBase, VTKWellGenerator, string const &, Group * const )
}
