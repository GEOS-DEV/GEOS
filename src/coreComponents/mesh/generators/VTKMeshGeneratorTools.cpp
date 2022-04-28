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
#include <vtkDIYGhostUtilities.h>

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

  using BlockType = std::vector< vtkSmartPointer< vtkUnstructuredGrid > >;

  diy::Master master( comm, 1, -1,
                      [] { return static_cast< void * >( new BlockType() ); },
                      []( void * b ) { delete static_cast< BlockType * >( b ); } );

  diy::ContiguousAssigner const assigner( comm.size(), comm.size() );
  diy::RegularDecomposer< diy::DiscreteBounds > decomposer( 1, diy::interval( 0, comm.size() - 1 ), comm.size() );
  decomposer.decompose( comm.rank(), assigner, master );
  assert( master.size() == 1 );

  int const myRank = comm.rank();
  diy::all_to_all( master, assigner, [myRank, &localParts]( BlockType * block, diy::ReduceProxy const & rp )
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
  for( auto & ug : *master.block< BlockType >( 0 ) )
  {
    appender->AddInputDataObject( ug );
  }
  appender->Update();

  return vtkUnstructuredGrid::SafeDownCast( appender->GetOutputDataObject( 0 ) );
}

std::vector< vtkBoundingBox >
exchangeBoundingBoxes( vtkDataSet & dataSet, MPI_Comm mpiComm )
{
  // The code below is modified from vtkDIYGhostUtilities::ExchangeBoundingBoxes():
  // https://gitlab.kitware.com/vtk/vtk/-/blob/1f0e4b2d0be7cd328795131642b5bf7984f681c1/Parallel/DIY/vtkDIYGhostUtilities.txx#L300
  // It makes some simplifications (e.g. just one input dataset per rank).

  using BlockType = std::map< int, vtkBoundingBox >;

  diy::mpi::communicator comm( mpiComm );
  diy::Master master( comm, 1, -1,
                      [] { return static_cast< void * >( new BlockType() ); },
                      []( void * b ) { delete static_cast< BlockType * >( b ); } );

  diy::ContiguousAssigner const assigner( comm.size(), comm.size() );
  diy::RegularDecomposer< diy::DiscreteBounds > decomposer( 1, diy::interval( 0, comm.size() - 1 ), comm.size() );
  decomposer.decompose( comm.rank(), assigner, master );
  assert( master.size() == 1 );

  diy::all_to_all( master, assigner, [&dataSet]( BlockType * block, diy::ReduceProxy const & rp )
  {
    int myBlockId = rp.gid();
    if( rp.round() == 0 )
    {
      vtkBoundingBox bb( dataSet.GetBounds() );
      for( int i = 0; i < rp.out_link().size(); ++i )
      {
        diy::BlockID const blockId = rp.out_link().target( i );
        if( blockId.gid != myBlockId )
        {
          rp.enqueue( blockId, bb.GetMinPoint(), 3 );
          rp.enqueue( blockId, bb.GetMaxPoint(), 3 );
        }
      }
    }
    else
    {
      double minPoint[3], maxPoint[3];
      for( int i = 0; i < static_cast< int >( rp.in_link().size() ); ++i )
      {
        diy::BlockID const blockId = rp.in_link().target( i );
        if( blockId.gid != myBlockId )
        {
          rp.dequeue( blockId, minPoint, 3 );
          rp.dequeue( blockId, maxPoint, 3 );

          block->emplace( blockId.gid, vtkBoundingBox( minPoint[0], maxPoint[0], minPoint[1],
                                                       maxPoint[1], minPoint[2], maxPoint[2] ) );
        }
      }
    }
  } );

  BlockType & boxMap = *master.block< BlockType >( 0 );
  boxMap.emplace( comm.rank(), vtkBoundingBox( dataSet.GetBounds() ) );
  assert( static_cast< int >( boxMap.size() ) == comm.size() );

  std::vector< vtkBoundingBox > boxes;
  boxes.reserve( boxMap.size() );
  for( auto const & rankBox : boxMap )
  {
    boxes.push_back( rankBox.second );
  }
  return boxes;
}

} // namespace vtk
} // namespace geosx
