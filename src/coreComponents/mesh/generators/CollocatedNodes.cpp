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

#include "CollocatedNodes.hpp"

#include "common/MpiWrapper.hpp"


#include <vtkPointData.h>

namespace geos::vtk
{

CollocatedNodes::CollocatedNodes( string const & faceBlockName,
                                  vtkSmartPointer< vtkDataSet > faceMesh,
                                  bool performGlobalCheck )
{
  // The vtk field to the collocated nodes for fractures.
  string const COLLOCATED_NODES = "collocated_nodes";

  vtkIdTypeArray const * collocatedNodes = vtkIdTypeArray::FastDownCast( faceMesh->GetPointData()->GetArray( COLLOCATED_NODES.c_str() ) );

  // With performGlobalCheck, a NULL on the local rank is OK if any other rank found the field
  // (vtk does not provide fields for empty meshes). Otherwise we require it on every rank.
  // After MPI max a non-zero result means the field exists on at
  // least one rank, so the field is not globally missing.
  bool const isMissing = performGlobalCheck
                         ? ( MpiWrapper::max( reinterpret_cast< std::uintptr_t >( collocatedNodes ) ) == 0 )
                         : ( collocatedNodes == nullptr );

  if( isMissing )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "Available point data fields in '{}':", faceBlockName ) );
    for( int i = 0; i < faceMesh->GetPointData()->GetNumberOfArrays(); ++i )
    {
      GEOS_LOG_RANK_0( GEOS_FMT( " - {} of type '{}'",
                                 faceMesh->GetPointData()->GetArrayName( i ),
                                 faceMesh->GetPointData()->GetArray( i )->GetDataTypeAsString() ) );
    }
    GEOS_ERROR( GEOS_FMT( "Could not find valid field \"{}\" for fracture \"{}\"{}.",
                          COLLOCATED_NODES,
                          faceBlockName,
                          performGlobalCheck ? "" : " on this rank" ) );
  }

  if( collocatedNodes )
  {
    init( collocatedNodes );
  }
}

void CollocatedNodes::init( vtkIdTypeArray const * collocatedNodes )
{
  vtkIdType const numTuples = collocatedNodes->GetNumberOfTuples();
  int const numComponents = collocatedNodes->GetNumberOfComponents();
  m_collocatedNodes.resize( numTuples );
  for( vtkIdType i = 0; i < numTuples; ++i )
  {
    m_collocatedNodes[i].reserve( numComponents );
  }

  for( vtkIdType i = 0; i < numTuples; ++i )
  {
    for( int j = 0; j < numComponents; ++j )
    {
      vtkIdType const tmp = collocatedNodes->GetTypedComponent( i, j );
      if( tmp > -1 )
      {
        m_collocatedNodes[i].emplace_back( tmp );
      }
    }
  }
}


} // geos::vtk
