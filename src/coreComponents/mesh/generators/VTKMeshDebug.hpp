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

/**
 * @file VTKMeshDebug.hpp
 */

#ifndef GEOS_MESH_GENERATORS_VTKMESHDEBUG_HPP_
#define GEOS_MESH_GENERATORS_VTKMESHDEBUG_HPP_

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkPartitionedDataSet.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include <mpi.h>

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

namespace geos::vtk::meshDebug
{

inline bool enabled()
{
  static bool const value = ::getenv( "GEOS_MESH_DEBUG" ) != nullptr;
  return value;
}

inline double elapsedSeconds()
{
  using Clock = std::chrono::steady_clock;
  static Clock::time_point const start = Clock::now();
  return std::chrono::duration< double >( Clock::now() - start ).count();
}

inline void rankAndSize( MPI_Comm comm, int & rank, int & size )
{
  rank = 0;
  size = 1;
#ifdef GEOS_USE_MPI
  int initialized = 0;
  MPI_Initialized( &initialized );
  if( initialized )
  {
    MPI_Comm_rank( comm, &rank );
    MPI_Comm_size( comm, &size );
  }
#else
  (void) comm;
#endif
}

inline void log( MPI_Comm comm, char const * message )
{
  if( !enabled() )
  {
    return;
  }

  int rank = 0;
  int size = 1;
  rankAndSize( comm, rank, size );
  std::fprintf( stdout,
                "[GEOS_MESH_DEBUG t=%10.6f rank=%d/%d] %s\n",
                elapsedSeconds(), rank, size, message );
  std::fflush( stdout );
}

inline void logf( MPI_Comm comm, char const * format, ... )
{
  if( !enabled() )
  {
    return;
  }

  char message[1024];
  va_list args;
  va_start( args, format );
  std::vsnprintf( message, sizeof( message ), format, args );
  va_end( args );
  log( comm, message );
}

inline void logDataSet( MPI_Comm comm, char const * message, vtkDataSet * mesh )
{
  if( !enabled() )
  {
    return;
  }

  int rank = 0;
  int size = 1;
  rankAndSize( comm, rank, size );

  if( mesh == nullptr )
  {
    std::fprintf( stdout,
                  "[GEOS_MESH_DEBUG t=%10.6f rank=%d/%d] %s dataset=null\n",
                  elapsedSeconds(), rank, size, message );
    std::fflush( stdout );
    return;
  }

  double bounds[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  if( mesh->GetNumberOfCells() > 0 || mesh->GetNumberOfPoints() > 0 )
  {
    mesh->GetBounds( bounds );
  }

  int const cellArrays = mesh->GetCellData() == nullptr ? -1 : mesh->GetCellData()->GetNumberOfArrays();
  int const pointArrays = mesh->GetPointData() == nullptr ? -1 : mesh->GetPointData()->GetNumberOfArrays();
  std::fprintf( stdout,
                "[GEOS_MESH_DEBUG t=%10.6f rank=%d/%d] %s type=%s cells=%lld points=%lld cellArrays=%d pointArrays=%d bounds=[%g,%g,%g,%g,%g,%g]\n",
                elapsedSeconds(), rank, size, message,
                mesh->GetClassName(),
                static_cast< long long >( mesh->GetNumberOfCells() ),
                static_cast< long long >( mesh->GetNumberOfPoints() ),
                cellArrays, pointArrays,
                bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5] );
  std::fflush( stdout );
}

inline void logPartitionedDataSet( MPI_Comm comm, char const * message, vtkPartitionedDataSet * parts )
{
  if( !enabled() )
  {
    return;
  }

  int rank = 0;
  int size = 1;
  rankAndSize( comm, rank, size );

  if( parts == nullptr )
  {
    std::fprintf( stdout,
                  "[GEOS_MESH_DEBUG t=%10.6f rank=%d/%d] %s partitions=null\n",
                  elapsedSeconds(), rank, size, message );
    std::fflush( stdout );
    return;
  }

  long long totalCells = 0;
  long long totalPoints = 0;
  long long nonEmptyParts = 0;
  unsigned int const numParts = parts->GetNumberOfPartitions();
  for( unsigned int p = 0; p < numParts; ++p )
  {
    vtkDataSet * part = vtkDataSet::SafeDownCast( parts->GetPartition( p ) );
    if( part != nullptr )
    {
      long long const cells = static_cast< long long >( part->GetNumberOfCells() );
      totalCells += cells;
      totalPoints += static_cast< long long >( part->GetNumberOfPoints() );
      nonEmptyParts += cells > 0 ? 1 : 0;
    }
  }

  std::fprintf( stdout,
                "[GEOS_MESH_DEBUG t=%10.6f rank=%d/%d] %s partitions=%u nonEmpty=%lld totalCells=%lld totalPoints=%lld\n",
                elapsedSeconds(), rank, size, message,
                numParts, nonEmptyParts, totalCells, totalPoints );
  std::fflush( stdout );
}

} // namespace geos::vtk::meshDebug

#endif // GEOS_MESH_GENERATORS_VTKMESHDEBUG_HPP_
