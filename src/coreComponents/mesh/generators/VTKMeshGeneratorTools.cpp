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
 * @file VTKMeshGeneratorTools.hpp
 */

#include "VTKMeshGeneratorTools.hpp"

#include "LvArray/src/system.hpp"

#include <vtkBoundingBox.h>
#include <vtkAppendFilter.h>
#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkVariant.h>
#include <vtkVersionMacros.h>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#ifdef GEOS_USE_MPI
#include <vtkDIYGhostUtilities.h>
#include <vtkDIYUtilities.h>
#endif

// Do not include GEOS headers that transitively include Format.hpp here.
// See full explanation in VTKMeshGeneratorTools.hpp.

namespace geos::vtk
{

vtkSmartPointer< vtkUnstructuredGrid >
appendMeshParts( stdVector< vtkUnstructuredGrid * > const & meshes )
{
  stdVector< vtkSmartPointer< vtkUnstructuredGrid > > inputs;
  for( auto * mesh : meshes )
  {
    if( !mesh ) continue;
    vtkSmartPointer< vtkUnstructuredGrid > input = mesh;
#if VTK_VERSION_NUMBER == VTK_VERSION_CHECK( 9, 7, 0 )
    input = vtkSmartPointer< vtkUnstructuredGrid >::New();
    input->ShallowCopy( mesh );
    // DIY deserialization can change an ID array's concrete class. VTK's
    // append filter recognizes only vtkIdTypeArray for topological merging.
    for( vtkDataSetAttributes * attributes : { static_cast< vtkDataSetAttributes * >( input->GetPointData() ),
                                             static_cast< vtkDataSetAttributes * >( input->GetCellData() ) } )
    {
      auto * original = attributes->GetGlobalIds();
      if( !original || vtkIdTypeArray::SafeDownCast( original ) ) continue;
      auto ids = vtkSmartPointer< vtkIdTypeArray >::New();
      ids->SetName( original->GetName() );
      ids->SetComponentName( 0, original->GetComponentName( 0 ) );
      ids->SetNumberOfValues( original->GetNumberOfTuples() );
      for( vtkIdType i = 0; i < original->GetNumberOfTuples(); ++i ) ids->SetValue( i, original->GetVariantValue( i ).ToLongLong() );
      attributes->SetGlobalIds( ids );
    }
#endif
    inputs.emplace_back( std::move( input ) );
  }
  vtkNew< vtkAppendFilter > appender;
  appender->MergePointsOn();
  for( auto const & input : inputs ) appender->AddInputDataObject( input );
  appender->Update();
  vtkSmartPointer< vtkUnstructuredGrid > result = appender->GetOutput();
#if VTK_VERSION_NUMBER == VTK_VERSION_CHECK( 9, 7, 0 )
  // A single nonempty grid is shallow-copied by VTK without tuple conversion.
  // Keep that path free of additional ID maps and writes to shared arrays.
  if( std::count_if( inputs.begin(), inputs.end(), []( auto const & input )
      { return input->GetNumberOfPoints() > 0 || input->GetNumberOfCells() > 0; } ) <= 1 ) return result;
  // VTK 9.7's CopyTuple fallback copies vtkIdTypeArray through double when
  // appending multiple inputs. Recopy these arrays with typed access, using
  // the filter's first-occurrence point order and concatenated cell order.
  bool const allPointIds = std::all_of( inputs.begin(), inputs.end(), []( auto const & input )
  { return input->GetNumberOfPoints() == 0 || vtkIdTypeArray::SafeDownCast( input->GetPointData()->GetGlobalIds() ); } );
  std::unordered_map< vtkIdType, vtkIdType > pointIndices;
  if( allPointIds )
  {
    for( auto const & input : inputs )
    {
      auto * ids = vtkIdTypeArray::SafeDownCast( input->GetPointData()->GetGlobalIds() );
      for( vtkIdType p = 0; p < input->GetNumberOfPoints(); ++p ) pointIndices.emplace( ids->GetValue( p ), pointIndices.size() );
    }
    if( pointIndices.size() != static_cast< std::size_t >( result->GetNumberOfPoints() ) )
      throw std::runtime_error( "VTK append did not preserve the expected global-ID point topology" );
  }
  for( bool const points : { true, false } )
  {
    if( points && !allPointIds ) continue;
    vtkDataSetAttributes * output = points ? static_cast< vtkDataSetAttributes * >( result->GetPointData() )
                                          : static_cast< vtkDataSetAttributes * >( result->GetCellData() );
    for( int a = 0; a < output->GetNumberOfArrays(); ++a )
    {
      auto * target = vtkIdTypeArray::SafeDownCast( output->GetAbstractArray( a ) );
      if( !target ) continue;
      vtkIdType offset = 0;
      for( auto const & input : inputs )
      {
        vtkDataSetAttributes * source = points ? static_cast< vtkDataSetAttributes * >( input->GetPointData() )
                                              : static_cast< vtkDataSetAttributes * >( input->GetCellData() );
        auto * array = vtkIdTypeArray::SafeDownCast( target == output->GetGlobalIds() ? source->GetGlobalIds()
                                                     : target->GetName() ? source->GetAbstractArray( target->GetName() ) : nullptr );
        vtkIdType const count = points ? input->GetNumberOfPoints() : input->GetNumberOfCells();
        if( array )
        {
          auto * ids = vtkIdTypeArray::SafeDownCast( input->GetPointData()->GetGlobalIds() );
          for( vtkIdType i = 0; i < count; ++i )
          {
            vtkIdType const destination = points ? pointIndices.at( ids->GetValue( i ) ) : offset + i;
            for( int c = 0; c < target->GetNumberOfComponents(); ++c ) target->SetTypedComponent( destination, c, array->GetTypedComponent( i, c ) );
          }
        }
        offset += count;
      }
      target->Modified();
    }
  }
#endif
  return result;
}

#ifdef GEOS_USE_MPI

vtkSmartPointer< vtkUnstructuredGrid >
redistribute( vtkPartitionedDataSet & localParts,
              MPI_Comm mpiComm )
{
  // VTK's XML writer, used by DIY to serialize the exchanged grids, can
  // evaluate empty arrays while reporting progress. Do not let those
  // internal operations trip GEOS' floating-point exception handler.
  LvArray::system::FloatingPointExceptionGuard guard;

  // The code below is modified from vtkDIYKdTreeUtilities::Exchange():
  // https://gitlab.kitware.com/vtk/vtk/-/blob/7037a148605bf9628710d8b729c22f27dd0ede93/Filters/ParallelDIY2/vtkDIYKdTreeUtilities.cxx#L289
  // We cannot call that function directly because vtkDIYKdTreeUtilities.hpp is a private header in VTK.
  // We also make simplifying assumptions about the nature of input (e.g. one partition per target rank).

  diy::mpi::communicator comm( mpiComm );
  assert( static_cast< int >( localParts.GetNumberOfPartitions() ) == comm.size() );

  using BlockType = stdVector< vtkSmartPointer< vtkUnstructuredGrid > >;

  diy::Master master( comm, 1, -1,
                      [] { return static_cast< void * >( new BlockType() ); },
                      []( void * b ) { delete static_cast< BlockType * >( b ); } );

  diy::ContiguousAssigner const assigner( comm.size(), comm.size() );
  diy::RegularDecomposer< diy::DiscreteBounds > decomposer( 1, diy::interval( 0, comm.size() - 1 ), comm.size() );
  decomposer.decompose( comm.rank(), assigner, master );
  assert( master.size() == 1 );

  int const myRank = comm.rank();
  diy::all_to_all( master, assigner, [myRank, &localParts]( BlockType * block, diy::ReduceProxy const & reduceProxy )
  {
    if( reduceProxy.in_link().size() == 0 )
    {
      // enqueue blocks to send.
      block->reserve( localParts.GetNumberOfPartitions() );
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
            reduceProxy.enqueue< vtkDataSet * >( reduceProxy.out_link().target( targetRank ), part );
          }
        }
      }
    }
    else
    {
      for( int i = 0; i < reduceProxy.in_link().size(); ++i )
      {
        int const gid = reduceProxy.in_link().target( i ).gid;
        while( reduceProxy.incoming( gid ) )
        {
          vtkDataSet * ptr = nullptr;
          reduceProxy.dequeue< vtkDataSet * >( reduceProxy.in_link().target( i ), ptr );

          vtkSmartPointer< vtkUnstructuredGrid > sptr;
          sptr.TakeReference( vtkUnstructuredGrid::SafeDownCast( ptr ) );
          block->push_back( sptr );
        }
      }
    }
  } );

  // At this point of the process, it is legitimate to have ranks with no cells for the cases with fractures.
  // But this leaves us with a technical problem since `vtkAppendFilter`
  // (that will be used to merge the different pieces of the meshes)
  // discards the empty the data sets it merges.
  // The definition of "empty" in its context is having no points nor cells...
  //
  // However, some other information which was defined in the discarded data sets gets lost too!
  // In particular, the cell, points and field data were _defined_, but _legitimately_ _empty_.
  // After the `vtkAppendFilter` processing, the definition of those fields is no more available.
  //
  // This leaves us with a specific case when importing fields from vtk,
  // since we need to take extra care of the empty data sets, while we should not have to do this.
  // To circumvent this issue, we gather the point, cell and field data by hand,
  // before registering them by hand again into the final `vtkUnstructuredGrid`.

  // This little structure stores the information we'll need to register back into
  // the cell, points and field data into the final vtkUnstructuredGrid.
  // We should not need it outside of this function, but if we needed, make it a little more solid.
  struct FieldMetaInfo
  {
    enum Location
    {
      CELL,
      POINT,
      FIELD
    };
    std::string name;
    int numComponents;
    int dataType;
    Location location;

    bool operator<( FieldMetaInfo const & other ) const
    {
      return std::tie( name, numComponents, dataType, location ) < std::tie( other.name, other.numComponents, other.dataType, other.location );
    }
  };
  // First step is to gather all the field information
  std::set< FieldMetaInfo > fieldMetaInfo;
  for( unsigned int i = 0; i < master.size(); ++i )
  {
    for( vtkUnstructuredGrid * ug: *master.block< BlockType >( i ) )
    {
      if( !ug )
      {
        break;
      }
      // vtkFieldData::GetArray() returns nullptr for vtkStringArray and other
      // non-vtkDataArray objects. Empty-rank reconstruction uses CreateArray
      // on the stored VTK type, so the scan must use GetAbstractArray().
      for( int c = 0; c < ug->GetCellData()->GetNumberOfArrays(); ++c )
      {
        vtkAbstractArray * array = ug->GetCellData()->GetAbstractArray( c );
        fieldMetaInfo.insert( { array->GetName(), array->GetNumberOfComponents(), array->GetDataType(), FieldMetaInfo::Location::CELL } );
      }
      for( int c = 0; c < ug->GetPointData()->GetNumberOfArrays(); ++c )
      {
        vtkAbstractArray * array = ug->GetPointData()->GetAbstractArray( c );
        fieldMetaInfo.insert( { array->GetName(), array->GetNumberOfComponents(), array->GetDataType(), FieldMetaInfo::Location::POINT } );
      }
      for( int c = 0; c < ug->GetFieldData()->GetNumberOfArrays(); ++c )
      {
        vtkAbstractArray * array = ug->GetFieldData()->GetAbstractArray( c );
        fieldMetaInfo.insert( { array->GetName(), array->GetNumberOfComponents(), array->GetDataType(), FieldMetaInfo::Location::FIELD } );
      }
    }
  }

  stdVector< vtkUnstructuredGrid * > meshes;
  for( unsigned int i = 0; i < master.size(); ++i )
  {
    for( vtkUnstructuredGrid * ug: *master.block< BlockType >( i ) )
    {
      meshes.emplace_back( ug );
    }
  }
  auto result = appendMeshParts( meshes );
  // Now we register back the field info.
  if( result->GetNumberOfCells() == 0 )
  {
    for( FieldMetaInfo const & info: fieldMetaInfo )
    {
      vtkAbstractArray * array = vtkAbstractArray::CreateArray( info.dataType );
      array->SetNumberOfComponents( info.numComponents );
      array->SetNumberOfTuples( 0 );
      array->SetName( info.name.c_str() );
      if( info.location == 0 )
      {
        result->GetCellData()->AddArray( array );
      }
      if( info.location == 1 )
      {
        result->GetPointData()->AddArray( array );
      }
      if( info.location == 2 )
      {
        result->GetFieldData()->AddArray( array );
      }
    }
  }

  return result;
}

stdVector< vtkBoundingBox >
exchangeBoundingBoxes( vtkDataSet & dataSet, MPI_Comm mpiComm )
{
  // The code below is modified from vtkDIYGhostUtilities::ExchangeBoundingBoxes():
  // https://gitlab.kitware.com/vtk/vtk/-/blob/1f0e4b2d0be7cd328795131642b5bf7984f681c1/Parallel/DIY/vtkDIYGhostUtilities.txx#L300
  // It makes some simplifications (e.g. just one input dataset per rank).

  using BlockType = stdMap< int, vtkBoundingBox >;

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

  stdVector< vtkBoundingBox > boxes;
  boxes.reserve( boxMap.size() );
  for( auto const & rankBox : boxMap )
  {
    boxes.push_back( rankBox.second );
  }
  return boxes;
}

#else

vtkSmartPointer< vtkUnstructuredGrid >
redistribute( vtkPartitionedDataSet & localParts, MPI_Comm mpiComm )
{
  static_cast< void >( mpiComm );
  assert( localParts.GetNumberOfPartitions() == 1 );
  auto result = vtkSmartPointer< vtkUnstructuredGrid >::New();
  if( auto * part = vtkUnstructuredGrid::SafeDownCast( localParts.GetPartition( 0 ) ) ) result->ShallowCopy( part );
  return result;
}

stdVector< vtkBoundingBox >
exchangeBoundingBoxes( vtkDataSet & dataSet, MPI_Comm mpiComm )
{
  static_cast< void >( mpiComm );
  return { vtkBoundingBox( dataSet.GetBounds() ) };
}

#endif

} // namespace geos::vtk
