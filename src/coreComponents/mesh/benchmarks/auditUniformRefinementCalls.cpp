/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/** @file auditUniformRefinementCalls.cpp
 * Linux preload instrumentation for the positive-refinement interval.
 * This benchmark/test library forwards every call; it does not change GEOS.
 */
#include "mesh/generators/VTKUniformRefinement.hpp"
#include "mesh/generators/VTKMeshScattering.hpp"
#include "mesh/generators/VTKUtilities.hpp"
#if defined(GEOS_USE_PARMETIS)
#include "mesh/generators/ParMETISInterface.hpp"
#endif
#if defined(GEOS_USE_SCOTCH)
#include "mesh/generators/PTScotchInterface.hpp"
#endif

#include <vtkAppendFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkCleanUnstructuredGrid.h>
#include <vtkMergePoints.h>
#if defined(GEOS_USE_MPI)
#include <vtkRedistributeDataSetFilter.h>
#endif
#include <vtkStaticCleanPolyData.h>
#include <vtkStaticCleanUnstructuredGrid.h>
#include <vtkUnstructuredGrid.h>

#include <cstdio>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <stdexcept>

namespace
{
struct Counters
{
  std::atomic< unsigned > gather{}, scatter{}, redistribute{}, append{}, clean{}, repartition{};
};
std::atomic< Counters * > active{ nullptr };

template< typename FUNCTION >
FUNCTION next( FUNCTION intercepted )
{
  static_assert( sizeof( FUNCTION ) == sizeof( void * ) );
  void * address;
  std::memcpy( &address, &intercepted, sizeof( address ) );
  Dl_info info{};
  if( !dladdr( address, &info ) || !info.dli_sname ) throw std::runtime_error( "Cannot identify audit symbol" );
  address = dlsym( RTLD_NEXT, info.dli_sname );
  if( !address ) throw std::runtime_error( std::string( "Cannot forward audit symbol " ) + info.dli_sname );
  FUNCTION original;
  std::memcpy( &original, &address, sizeof( original ) );
  return original;
}

struct Interval
{
  Counters counters;
  Counters * previous;
  int levels;
  MPI_Comm comm;
  Interval( int levelCount, MPI_Comm communicator ): levels( levelCount ), comm( communicator )
  {
    previous = active.exchange( &counters );
  }
  ~Interval()
  {
    active.store( previous );
    int rank = 0;
    int mpiEnabled = 0;
#if defined(GEOS_USE_MPI)
    mpiEnabled = 1;
    PMPI_Comm_rank( comm, &rank );
#endif
    // ProblemManager captures stderr for external-library error detection.
    std::fprintf( stdout, "uniform_refinement_audit rank=%d levels=%d mpi=%d gather=%u scatter=%u redistribute=%u append=%u clean=%u repartition=%u\n",
                  rank, levels, mpiEnabled, counters.gather.load(), counters.scatter.load(), counters.redistribute.load(),
                  counters.append.load(), counters.clean.load(), counters.repartition.load() );
    std::fflush( stdout );
  }
};
} // namespace

#define AUDIT_NEW( CLASS, COUNTER ) \
  CLASS * CLASS::New() \
  { \
    if( auto * counters = active.load() ) ++counters->COUNTER; \
    static auto original = next( &CLASS::New ); \
    return original(); \
  }

AUDIT_NEW( vtkAppendFilter, append )
AUDIT_NEW( vtkCleanPolyData, clean )
AUDIT_NEW( vtkCleanUnstructuredGrid, clean )
AUDIT_NEW( vtkStaticCleanPolyData, clean )
AUDIT_NEW( vtkStaticCleanUnstructuredGrid, clean )
AUDIT_NEW( vtkMergePoints, clean )
#if defined(GEOS_USE_MPI)
AUDIT_NEW( vtkRedistributeDataSetFilter, repartition )
#endif
#undef AUDIT_NEW

#if defined(GEOS_USE_MPI)
extern "C" int MPI_Gather( void const * send, int sendCount, MPI_Datatype sendType,
                            void * receive, int receiveCount, MPI_Datatype receiveType,
                            int root, MPI_Comm comm )
{
  if( auto * counters = active.load() ) ++counters->gather;
  return PMPI_Gather( send, sendCount, sendType, receive, receiveCount, receiveType, root, comm );
}

extern "C" int MPI_Gatherv( void const * send, int sendCount, MPI_Datatype sendType,
                             void * receive, int const * counts, int const * displacements,
                             MPI_Datatype receiveType, int root, MPI_Comm comm )
{
  if( auto * counters = active.load() ) ++counters->gather;
  return PMPI_Gatherv( send, sendCount, sendType, receive, counts, displacements, receiveType, root, comm );
}
#endif

namespace geos::vtk
{
UniformRefinementResult refineUniformly( AllMeshes & meshes, int levels, UniformRefinementOptions const & options, MPI_Comm comm )
{
  static auto original = next( &refineUniformly );
  if( levels <= 0 ) return original( meshes, levels, options, comm );
  Interval interval( levels, comm );
  if( std::getenv( "GEOS_REFINEMENT_AUDIT_PROBE" ) )
  {
    // Deliberate harmless calls prove the instrumentation detects forbidden work.
    vtkAppendFilter::New()->Delete();
    vtkCleanPolyData::New()->Delete();
#if defined(GEOS_USE_MPI)
    vtkRedistributeDataSetFilter::New()->Delete();
    MPI_Gather( nullptr, 0, MPI_BYTE, nullptr, 0, MPI_BYTE, 0, comm );
#endif
  }
  return original( meshes, levels, options, comm );
}

vtkSmartPointer< vtkDataSet > scatterMesh( ScatterMethod method, vtkDataSet & mesh,
                                          arrayView1d< integer const > partitions, MPI_Comm comm )
{
  if( auto * counters = active.load() ) ++counters->scatter;
  static auto original = next( &scatterMesh );
  return original( method, mesh, partitions, comm );
}

vtkSmartPointer< vtkUnstructuredGrid > scatterByRankAssignment( vtkUnstructuredGrid * mesh,
                                                               stdVector< integer > assignment, MPI_Comm comm )
{
  if( auto * counters = active.load() ) ++counters->scatter;
  static auto original = next( &scatterByRankAssignment );
  return original( mesh, std::move( assignment ), comm );
}

AllMeshes redistributeMeshes( integer logLevel, vtkSmartPointer< vtkDataSet > mesh,
                              stdMap< string, vtkSmartPointer< vtkDataSet > > & fractures, MPI_Comm comm,
                              ScatterMethod scatter, arrayView1d< int const > partitions, PartitionMethod method,
                              int partitionRefinement, int fractureWeight, int useIds, string const & indexName )
{
  if( auto * counters = active.load() ) ++counters->redistribute;
  static auto original = next( &redistributeMeshes );
  return original( logLevel, std::move( mesh ), fractures, comm, scatter, partitions, method,
                   partitionRefinement, fractureWeight, useIds, indexName );
}
} // namespace geos::vtk

#if defined(GEOS_USE_PARMETIS)
namespace geos::parmetis
{
array1d< pmet_idx_t > partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                arrayView1d< pmet_idx_t const > const & distribution,
                                pmet_idx_t parts, MPI_Comm comm, int refinements )
{
  if( auto * counters = active.load() ) ++counters->repartition;
  static auto original = next( &partition );
  return original( graph, distribution, parts, comm, refinements );
}
array1d< pmet_idx_t > partitionWeighted( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                        arrayView1d< pmet_idx_t const > const & weights,
                                        arrayView1d< pmet_idx_t const > const & distribution,
                                        pmet_idx_t parts, MPI_Comm comm, int refinements )
{
  if( auto * counters = active.load() ) ++counters->repartition;
  static auto original = next( &partitionWeighted );
  return original( graph, weights, distribution, parts, comm, refinements );
}
} // namespace geos::parmetis
#endif

#if defined(GEOS_USE_SCOTCH)
namespace geos::ptscotch
{
array1d< int64_t > partition( ArrayOfArraysView< int64_t const, int64_t > const & graph, int64_t parts, MPI_Comm comm )
{
  if( auto * counters = active.load() ) ++counters->repartition;
  static auto original = next( &partition );
  return original( graph, parts, comm );
}
} // namespace geos::ptscotch
#endif
