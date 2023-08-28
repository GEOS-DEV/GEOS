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

/*
 * @file VTKWellGenerator.cpp
 *
 */

#include "VTKWellGenerator.hpp"

#include <vtkPolyData.h>
#include <vtkFieldData.h>
#include <vtkCellData.h>

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

}

void VTKWellGenerator::fillPolylineDataStructure( )
{
  GEOS_MARK_FUNCTION;

  vtkSmartPointer< vtkMultiProcessController > controller = vtk::getController();
  vtkMultiProcessController::SetGlobalController( controller );

  GEOS_LOG_RANK_0( GEOS_FMT( "{} '{}': reading well from {}", catalogName(), getName(), m_filePath ) );
  {
    GEOS_LOG_LEVEL_RANK_0( 2, "  reading the dataset..." );
    vtkSmartPointer< vtkDataSet > loadedMesh = vtk::loadMesh( m_filePath, "main" );
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

    // load edges
    polyData->GetLines()->InitTraversal();
    vtkNew< vtkIdList > idList;
    polyData->GetLines()->GetNextCell( idList );

    const globalIndex nbSegments = idList->GetNumberOfIds() - 1;
    m_segmentToPolyNodeMap.resizeDimension< 0 >( nbSegments );
    m_segmentToPolyNodeMap.resizeDimension< 1 >( m_numNodesPerElem );

    globalIndex iseg = 0;
    for( vtkIdType pointId = 0; pointId < nbSegments; ++pointId )
    {
      m_segmentToPolyNodeMap[iseg][0] = idList->GetId( pointId );
      m_segmentToPolyNodeMap[iseg][1] = idList->GetId( pointId + 1 );
      ++iseg;
    }
  }
}

REGISTER_CATALOG_ENTRY( WellGeneratorBase, VTKWellGenerator, string const &, Group * const )
}
