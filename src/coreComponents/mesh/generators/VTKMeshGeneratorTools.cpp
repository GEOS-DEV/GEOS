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

/**
 * @file VTKMeshGeneratorTools.hpp
 */

#include "VTKMeshGeneratorTools.hpp"

#include <vtkAppendFilter.h>
#include <vtkDIYUtilities.h>

// NOTE: do NOT include anything from GEOSX here.
// See full explanation in VTKMeshGeneratorTools.hpp.

namespace geosx
{
namespace vtk
{

vtkSmartPointer< vtkUnstructuredGrid >
redistribute( vtkPartitionedDataSet & localParts,
              MPI_Comm mpiComm )
{
  // The code below is modified from vtkDIYKdTreeUtilities::Exchange():
  // https://gitlab.kitware.com/vtk/vtk/-/blob/7037a148605bf9628710d8b729c22f27dd0ede93/Filters/ParallelDIY2/vtkDIYKdTreeUtilities.cxx#L289
  // We cannot call that function directly because vtkDIYKdTreeUtilities.hpp is a private header in VTK.
  // We also make simplifying assumptions about the nature of input (e.g. one partition per target rank).

  diy::mpi::communicator comm( mpiComm );
  assert( static_cast< int >( localParts.GetNumberOfPartitions() ) == comm.size() );

  using VectorOfUG = std::vector< vtkSmartPointer< vtkUnstructuredGrid > >;

  diy::Master master( comm, 1, -1,
                      []() { return static_cast< void * >( new VectorOfUG() ); },
                      []( void * b ) { delete static_cast< VectorOfUG * >( b ); } );

  diy::ContiguousAssigner const assigner( comm.size(), comm.size() );
  diy::RegularDecomposer< diy::DiscreteBounds > decomposer( 1, diy::interval( 0, comm.size() - 1 ), comm.size() );
  decomposer.decompose( comm.rank(), assigner, master );
  assert( master.size() == 1 );

  int const myRank = comm.rank();
  diy::all_to_all( master, assigner, [myRank, &localParts]( VectorOfUG * block, diy::ReduceProxy const & rp )
  {
    if( rp.in_link().size() == 0 )
    {
      // enqueue blocks to send.
      block->resize( localParts.GetNumberOfPartitions() );
      for( unsigned int partId = 0; partId < localParts.GetNumberOfPartitions(); ++partId )
      {
        if( auto part = vtkUnstructuredGrid::SafeDownCast( localParts.GetPartition( partId ) ) )
        {
          int const targetRank = static_cast< int >( partId );
          if( targetRank == myRank )
          {
            // short-circuit messages to self.
            block->push_back( part );
          }
          else
          {
            rp.enqueue< vtkDataSet * >( rp.out_link().target( targetRank ), part );
          }
        }
      }
    }
    else
    {
      for( int i = 0; i < rp.in_link().size(); ++i )
      {
        int const gid = rp.in_link().target( i ).gid;
        while( rp.incoming( gid ) )
        {
          vtkDataSet * ptr = nullptr;
          rp.dequeue< vtkDataSet * >( rp.in_link().target( i ), ptr );

          vtkSmartPointer< vtkUnstructuredGrid > sptr;
          sptr.TakeReference( vtkUnstructuredGrid::SafeDownCast( ptr ) );
          block->push_back( sptr );
        }
      }
    }
  } );

  vtkNew< vtkAppendFilter > appender;
  appender->MergePointsOn();
  for( auto & ug : *master.block< VectorOfUG >( 0 ) )
  {
    appender->AddInputDataObject( ug );
  }
  appender->Update();

  return vtkUnstructuredGrid::SafeDownCast( appender->GetOutputDataObject( 0 ) );
}

} // namespace vtk
} // namespace geosx
