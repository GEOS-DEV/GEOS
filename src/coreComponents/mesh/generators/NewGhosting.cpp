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

#include "NewGhosting.hpp"

#include "NewGlobalNumbering.hpp"

#include "Pods.hpp"

#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"

#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_VectorOut.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_RowMatrixTransposer.h>
#include <Epetra_Vector.h>

#include "Indices.hpp"

#include <vtkCellData.h>
#include <vtkPointData.h>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

#include <algorithm>
#include <utility>

namespace geos::ghosting
{

struct MaxGlbIdcs
{
  NodeGlbIdx nodes;
  EdgeGlbIdx edges;
  FaceGlbIdx faces;
  CellGlbIdx cells;
};

void to_json( json & j,
              const MaxGlbIdcs & mo )
{
  j = json{ { "nodes", mo.nodes },
            { "edges", mo.edges },
            { "faces", mo.faces },
            { "cells", mo.cells } };
}

void g( void * in,
        void * inout,
        int * len,
        MPI_Datatype * dataType )
{
  GEOS_ASSERT_EQ( *len, 1 );

  MaxGlbIdcs const * i = reinterpret_cast<MaxGlbIdcs const *>(in);
  MaxGlbIdcs * io = reinterpret_cast<MaxGlbIdcs *>(inout);

  io->nodes = std::max( i->nodes, io->nodes );
  io->edges = std::max( i->edges, io->edges );
  io->faces = std::max( i->faces, io->faces );
  io->cells = std::max( i->cells, io->cells );
}

MaxGlbIdcs gatherOffset( vtkSmartPointer< vtkDataSet > mesh,
                         EdgeGlbIdx const & maxEdgeId,
                         FaceGlbIdx const & maxFaceId )
{
  MaxGlbIdcs offsets{ NodeGlbIdx{ 0 }, maxEdgeId, maxFaceId, CellGlbIdx{ 0 } };

  auto const extract = []( vtkDataArray * globalIds ) -> vtkIdType
  {
    vtkIdTypeArray * gids = vtkIdTypeArray::FastDownCast( globalIds );
    Span< vtkIdType > const s( (vtkIdType *) gids->GetPointer( 0 ), gids->GetNumberOfTuples() );
    return *std::max_element( s.begin(), s.end() );
  };

  offsets.nodes = NodeGlbIdx{ extract( mesh->GetPointData()->GetGlobalIds() ) };
  offsets.cells = CellGlbIdx{ extract( mesh->GetCellData()->GetGlobalIds() ) };

  // Otherwise, use `MPI_Type_create_struct`.
  static_assert( std::is_same_v< NodeGlbIdx::UnderlyingType, EdgeGlbIdx::UnderlyingType > );
  static_assert( std::is_same_v< NodeGlbIdx::UnderlyingType, FaceGlbIdx::UnderlyingType > );
  static_assert( std::is_same_v< NodeGlbIdx::UnderlyingType, CellGlbIdx::UnderlyingType > );

//  MPI_Datatype const underlying = internal::getMpiType< NodeGlbIdx::UnderlyingType >();
  MPI_Datatype t;
//  MPI_Type_contiguous( sizeof( MatrixOffsets ) / sizeof( underlying ), underlying, &t );
  MPI_Type_contiguous( 4, MPI_LONG_LONG_INT, &t );
  MPI_Type_commit( &t );

  MPI_Op op;
  MPI_Op_create( g, true, &op );

  MaxGlbIdcs result( offsets );

  MPI_Allreduce( &offsets, &result, 1, t, op, MPI_COMM_WORLD );

  return result;
}

struct MeshGraph
{
  std::map< CellGlbIdx, std::set< FaceGlbIdx > > c2f;  // TODO What about the metadata (e.g. flip the face)
  std::map< FaceGlbIdx, std::set< EdgeGlbIdx > > f2e;
  std::map< EdgeGlbIdx, std::tuple< NodeGlbIdx, NodeGlbIdx > > e2n; // TODO use Edge here?
  std::set< NodeGlbIdx > n;
};

void to_json( json & j,
              const MeshGraph & v )  // For display
{
  j = json{ { "c2f", v.f2e },
            { "f2e", v.f2e },
            { "e2n", v.e2n },
            { "n", v.n } };
}

/**
 * @brief Builds the graph information for the owned elements only.
 * @param mesh
 * @param buckets
 * @param offsets
 * @param curRank
 * @return The tuple of first the owned geometrical quantities (ie keys are owned)
 * and second the geometrical quantities that are present on the rank but not owned.
 */
std::tuple< MeshGraph, MeshGraph > buildMeshGraph( vtkSmartPointer< vtkDataSet > mesh,  // TODO give a sub-mesh?
                                                   Buckets const & buckets,
                                                   BucketOffsets const & offsets,
                                                   MpiRank curRank )
{
  MeshGraph owned, present;

  auto const isCurrentRankOwning = [&curRank]( std::set< MpiRank > const & ranks ) -> bool
  {
    return curRank == *std::min_element( std::cbegin( ranks ), std::cend( ranks ) );
  };

  for( auto const & [ranks, nodes]: buckets.nodes )
  {
    std::set< NodeGlbIdx > & n = isCurrentRankOwning( ranks ) ? owned.n : present.n;
    n.insert( std::cbegin( nodes ), std::cend( nodes ) );
  }

  for( auto const & [ranks, edges]: buckets.edges )
  {
    auto & e2n = isCurrentRankOwning( ranks ) ? owned.e2n : present.e2n;
    EdgeGlbIdx egi = offsets.edges.at( ranks );  // TODO hack
    for( Edge const & edge: edges )
    {
      e2n[egi] = edge;
      ++egi;
    }
  }
  // The `e2n` is a mapping for all the geometrical entities, not only the one owned like `result.e2n`.
  // TODO check that it is really useful.
  std::map< EdgeGlbIdx, std::tuple< NodeGlbIdx, NodeGlbIdx > > e2n;
  for( auto const & m: { owned.e2n, present.e2n } )
  {
    e2n.insert( std::cbegin( m ), std::cend( m ) );
  }

  // Simple inversion
  std::map< std::tuple< NodeGlbIdx, NodeGlbIdx >, EdgeGlbIdx > n2e;
  for( auto const & [e, n]: e2n )
  {
    n2e[n] = e;  // TODO what about ownership?
  }

  for( auto const & [ranks, faces]: buckets.faces )
  {
    auto & f2e = isCurrentRankOwning( ranks ) ? owned.f2e : present.f2e;

    FaceGlbIdx fgi = offsets.faces.at( ranks );
    for( Face face: faces )  // Intentional copy for the future `emplace_back`.
    {
      face.emplace_back( face.front() );  // Trick to build the edges.
      for( std::size_t i = 0; i < face.size() - 1; ++i )
      {
        NodeGlbIdx const & n0 = face[i], & n1 = face[i + 1];
        std::pair< NodeGlbIdx, NodeGlbIdx > const p0 = std::make_pair( n0, n1 );
        std::pair< NodeGlbIdx, NodeGlbIdx > const p1 = std::minmax( n0, n1 );
        f2e[fgi].insert( n2e.at( p1 ) );
        bool const flipped = p0 != p1;  // TODO store somewhere.
      }
      ++fgi;
    }
  }

  std::map< std::vector< NodeGlbIdx >, FaceGlbIdx > n2f;
  for( auto const & [ranks, faces]: buckets.faces )
  {
    FaceGlbIdx fgi = offsets.faces.at( ranks );
    for( Face const & face: faces )  // TODO hack
    {
      n2f[face] = fgi;
      ++fgi;
    }
  }

  vtkIdTypeArray const * globalPtIds = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );
  vtkIdTypeArray const * globalCellIds = vtkIdTypeArray::FastDownCast( mesh->GetCellData()->GetGlobalIds() );  // TODO do the mapping beforehand
  for( vtkIdType c = 0; c < mesh->GetNumberOfCells(); ++c )
  {
    vtkCell * cell = mesh->GetCell( c );
    CellGlbIdx const gci{ globalCellIds->GetValue( c ) };
    // TODO copy paste?
    for( auto f = 0; f < cell->GetNumberOfFaces(); ++f )
    {
      vtkCell * face = cell->GetFace( f );
      vtkIdList * pids = face->GetPointIds();
      std::vector< NodeGlbIdx > faceNodes( pids->GetNumberOfIds() );
      for( std::size_t i = 0; i < faceNodes.size(); ++i )
      {
        vtkIdType const nli = face->GetPointId( i );
        vtkIdType const ngi = globalPtIds->GetValue( nli );
        faceNodes[i] = NodeGlbIdx{ ngi };
      }
      std::vector< NodeGlbIdx > const reorderedFaceNodes = reorderFaceNodes( faceNodes );
      owned.c2f[gci].insert( n2f.at( reorderedFaceNodes ) );
      // TODO... bool const flipped = ... compare nodes and reorderedNodes. Or ask `reorderFaceNodes` to tell
    }
  }

  return { std::move( owned ), std::move( present ) };
}

std::unique_ptr< Epetra_CrsMatrix > makeTranspose( Epetra_CrsMatrix & input,
                                                   bool makeDataContiguous = true )
{
  Epetra_RowMatrixTransposer transposer( &input );  // This process does not modify the original `ghosted` matrix.
  Epetra_CrsMatrix * tr = nullptr;  // The transposer returns a pointer we must handle ourselves.
  transposer.CreateTranspose( makeDataContiguous, tr );

  std::unique_ptr< Epetra_CrsMatrix > ptr;
  ptr.reset( tr );
  return ptr;
}

struct Ownerships
{
  std::map< NodeGlbIdx, MpiRank > nodes;
  std::map< EdgeGlbIdx, MpiRank > edges;
  std::map< FaceGlbIdx, MpiRank > faces;
  std::map< CellGlbIdx, MpiRank > cells;
  std::set< MpiRank > neighbors;
};

Epetra_CrsMatrix multiply( int commSize,
                           Epetra_CrsMatrix const & indicator,
                           Epetra_CrsMatrix & upward )
{
  Epetra_Map const & ownedMap = upward.RowMap();
  Epetra_Map const & mpiMap = indicator.RangeMap();

  // Upward (n -> e -> f -> c)

  Epetra_CrsMatrix result_u0_0( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( upward, false, indicator, true, result_u0_0, false );
  result_u0_0.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-0.mat", result_u0_0 );

  Epetra_CrsMatrix result_u0_1( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( upward, false, result_u0_0, false, result_u0_1, false );
  result_u0_1.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-1.mat", result_u0_1 );

  Epetra_CrsMatrix result_u0_2( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( upward, false, result_u0_1, false, result_u0_2, false );
  result_u0_2.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-2.mat", result_u0_2 );

  // Downward (c -> f -> e -> n)
  auto tDownward = makeTranspose( upward );  // TODO check the algorithm to understand what's more relevant.

  Epetra_CrsMatrix result_d0_0( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( *tDownward, false, result_u0_2, false, result_d0_0, false );
  result_d0_0.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-4.mat", result_d0_0 );

  Epetra_CrsMatrix result_d0_1( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( *tDownward, false, result_d0_0, false, result_d0_1, false );
  result_d0_1.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-5.mat", result_d0_1 );

  Epetra_CrsMatrix result_d0_2( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( *tDownward, false, result_d0_1, false, result_d0_2, false );
  result_d0_2.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-6.mat", result_d0_2 );

  return result_d0_2;
}

class FindGeometricalType
{
public:
  enum Geom
  {
    NODE,
    EDGE,
    FACE,
    CELL
  };

  FindGeometricalType( std::size_t const edgeOffset,
                       std::size_t const faceOffset,
                       std::size_t const cellOffset )
    : m_edgeOffset( intConv< int >( edgeOffset ) ),
      m_faceOffset( intConv< int >( faceOffset ) ),
      m_cellOffset( intConv< int >( cellOffset ) )
  { }

  Geom operator()( int const & index ) const
  {
    if( index < m_edgeOffset )
    {
      return Geom::NODE;
    }
    else if( index < m_faceOffset )
    {
      return Geom::EDGE;
    }
    else if( index < m_cellOffset )
    {
      return Geom::FACE;
    }
    else
    {
      return Geom::CELL;
    }
  }

  template< typename GLOBAL_INDEX >
  [[nodiscard]] GLOBAL_INDEX as( int const & index ) const
  {
    if constexpr( std::is_same_v< GLOBAL_INDEX, NodeGlbIdx > )
    {
      return NodeGlbIdx{ intConv< NodeGlbIdx::UnderlyingType >( index ) };
    }
    else if constexpr( std::is_same_v< GLOBAL_INDEX, EdgeGlbIdx > )
    {
      return EdgeGlbIdx{ intConv< EdgeGlbIdx::UnderlyingType >( index - m_edgeOffset ) };
    }
    else if constexpr( std::is_same_v< GLOBAL_INDEX, FaceGlbIdx > )
    {
      return FaceGlbIdx{ intConv< FaceGlbIdx::UnderlyingType >( index - m_faceOffset ) };
    }
    else if constexpr( std::is_same_v< GLOBAL_INDEX, CellGlbIdx > )
    {
      return CellGlbIdx{ intConv< CellGlbIdx::UnderlyingType >( index - m_cellOffset ) };
    }
    else
    {
      GEOS_ERROR( "Internal Error." );
    }
  }

private:
  int const m_edgeOffset;
  int const m_faceOffset;
  int const m_cellOffset;
};

std::tuple< MeshGraph, Ownerships > assembleAdjacencyMatrix( MeshGraph const & owned,
                                                             MeshGraph const & present,
                                                             MaxGlbIdcs const & gis,
                                                             MpiRank curRank )
{
  std::size_t const edgeOffset = gis.nodes.get() + 1;
  std::size_t const faceOffset = edgeOffset + gis.edges.get() + 1;
  std::size_t const cellOffset = faceOffset + gis.faces.get() + 1;

  std::size_t const n = cellOffset + gis.cells.get() + 1;  // Total number of entries in the graph.

  std::size_t const numOwnedNodes = std::size( owned.n );
  std::size_t const numOwnedEdges = std::size( owned.e2n );
  std::size_t const numOwnedFaces = std::size( owned.f2e );
  std::size_t const numOwnedCells = std::size( owned.c2f );
  std::size_t const numOwned = numOwnedNodes + numOwnedEdges + numOwnedFaces + numOwnedCells;

  std::size_t const numOtherNodes = std::size( present.n );
  std::size_t const numOtherEdges = std::size( present.e2n );
  std::size_t const numOtherFaces = std::size( present.f2e );
  std::size_t const numOther = numOtherNodes + numOtherEdges + numOtherFaces;

  std::vector< int > ownedGlbIdcs, numEntriesPerRow;  // TODO I couldn't use a vector of `std::size_t`
  std::vector< std::vector< int > > indices;
  ownedGlbIdcs.reserve( numOwned );
  numEntriesPerRow.reserve( numOwned );
  indices.reserve( numOwned );

  std::vector< int > otherGlbIdcs;  // TODO I couldn't use a vector of `std::size_t`
  otherGlbIdcs.reserve( numOther );
  for( NodeGlbIdx const & ngi: present.n )
  {
    otherGlbIdcs.emplace_back( ngi.get() );
  }
  for( auto const & [egi, _]: present.e2n )
  {
    otherGlbIdcs.emplace_back( egi.get() + edgeOffset );
  }
  for( auto const & [fgi, _]: present.f2e )
  {
    otherGlbIdcs.emplace_back( fgi.get() + faceOffset );
  }
  GEOS_ASSERT_EQ( numOther, std::size( otherGlbIdcs ) );

  for( NodeGlbIdx const & ngi: owned.n )
  {
    auto const i = ngi.get();
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( 1 );
    std::vector< int > const tmp( 1, ownedGlbIdcs.back() );
    indices.emplace_back( tmp );
  }
  for( auto const & [egi, nodes]: owned.e2n )
  {
    auto const i = egi.get() + edgeOffset;
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( std::tuple_size_v< decltype( nodes ) > + 1 );  // `+1` comes from the identity
    std::vector< int > const tmp{ int( std::get< 0 >( nodes ).get() ), int( std::get< 1 >( nodes ).get() ), ownedGlbIdcs.back() };
    indices.emplace_back( tmp );
  }
  for( auto const & [fgi, edges]: owned.f2e )
  {
    auto const i = fgi.get() + faceOffset;
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( std::size( edges ) + 1 );  // `+1` comes from the identity
    std::vector< int > tmp;
    tmp.reserve( numEntriesPerRow.back() );
    for( EdgeGlbIdx const & egi: edges )
    {
      tmp.emplace_back( egi.get() + edgeOffset );
    }
    tmp.emplace_back( ownedGlbIdcs.back() );
    indices.emplace_back( tmp );
  }
  for( auto const & [cgi, faces]: owned.c2f )
  {
    auto const i = cgi.get() + cellOffset;
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( std::size( faces ) + 1 );  // `+1` comes from the identity
    std::vector< int > tmp;
    tmp.reserve( numEntriesPerRow.back() );
    for( FaceGlbIdx const & fgi: faces )
    {
      tmp.emplace_back( fgi.get() + faceOffset );
    }
    tmp.emplace_back( ownedGlbIdcs.back() );
    indices.emplace_back( tmp );
  }

  GEOS_ASSERT_EQ( numOwned, std::size( ownedGlbIdcs ) );
  GEOS_ASSERT_EQ( numOwned, std::size( numEntriesPerRow ) );
  GEOS_ASSERT_EQ( numOwned, std::size( indices ) );
  for( std::size_t i = 0; i < numOwned; ++i )
  {
    GEOS_ASSERT_EQ( indices[i].size(), std::size_t( numEntriesPerRow[i] ) );
  }

  Epetra_MpiComm const & comm = Epetra_MpiComm( MPI_COMM_GEOSX );
  Epetra_Map const ownedMap( n, numOwned, ownedGlbIdcs.data(), 0, comm );

  Epetra_CrsMatrix upward( Epetra_DataAccess::Copy, ownedMap, numEntriesPerRow.data(), true );

  for( std::size_t i = 0; i < numOwned; ++i )
  {
    std::vector< int > const & rowIndices = indices[i];
    std::vector< double > const rowValues( std::size( rowIndices ), 1. );
    GEOS_ASSERT_EQ( std::size( rowIndices ), std::size_t( numEntriesPerRow[i] ) );
    upward.InsertGlobalValues( ownedGlbIdcs[i], std::size( rowIndices ), rowValues.data(), rowIndices.data() );
  }

//  upward.FillComplete();
  upward.FillComplete( ownedMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/adj.mat", upward );

  int const commSize( MpiWrapper::commSize() );

  // Now let's build the domain indicator matrix.
  // It's rectangular, one dimension being the number of MPI ranks, the other the number of nodes in the mesh graph.
  Epetra_Map const mpiMap( commSize, 0, comm );  // Let the current rank get the appropriate index in this map.
  Epetra_CrsMatrix indicator( Epetra_DataAccess::Copy, mpiMap, numOwned + numOther, true );

  std::vector< double > const ones( n, 1. );
  indicator.InsertGlobalValues( curRank.get(), numOwned, ones.data(), ownedGlbIdcs.data() );
  indicator.InsertGlobalValues( curRank.get(), numOther, ones.data(), otherGlbIdcs.data() );
  indicator.FillComplete( ownedMap, mpiMap );

  // TODO Could we use an Epetra_Vector as a diagonal matrix?
  std::vector< double > myRank( 1, curRank.get() );
  Epetra_CrsMatrix ownership( Epetra_DataAccess::Copy, ownedMap, 1, true );
  for( auto const & i: ownedGlbIdcs )
  {
    ownership.InsertGlobalValues( i, 1, myRank.data(), &i );
  }
  ownership.FillComplete( ownedMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/ownership.mat", ownership );

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "indicator.NumGlobalCols() = " << indicator.NumGlobalCols() );
    GEOS_LOG_RANK( "indicator.NumGlobalRows() = " << indicator.NumGlobalRows() );
    GEOS_LOG_RANK( "ownership.NumGlobalCols() = " << ownership.NumGlobalCols() );
    GEOS_LOG_RANK( "ownership.NumGlobalRows() = " << ownership.NumGlobalRows() );
    GEOS_LOG_RANK( "ownership diag = " << std::boolalpha << ownership.LowerTriangular() and ownership.UpperTriangular() );
  }
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/indicator.mat", indicator );

  Epetra_CrsMatrix ghosted( multiply( commSize, indicator, upward ) );
  ghosted.PutScalar( 1. );  // This can be done after the `FillComplete`.
  ghosted.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/ghosted.mat", ghosted );

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "ghosted->NumGlobalCols() = " << ghosted.NumGlobalCols() );
    GEOS_LOG_RANK( "ghosted->NumGlobalRows() = " << ghosted.NumGlobalRows() );
  }
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/ghosted-after.mat", ghosted );

  Epetra_CrsMatrix ghostInfo( Epetra_DataAccess::Copy, mpiMap, 1, false );
  EpetraExt::MatrixMatrix::Multiply( ghosted, true, ownership, false, ghostInfo, false );
  ghostInfo.FillComplete( ownedMap, mpiMap );
  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "ghostInfo->NumGlobalCols() = " << ghostInfo.NumGlobalCols() );
    GEOS_LOG_RANK( "ghostInfo->NumGlobalRows() = " << ghostInfo.NumGlobalRows() );
  }

  auto tDownward = makeTranspose( upward );  // TODO give it to multiply!
//  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/tDownward-before.mat", *tDownward );
  Epetra_Vector const zeros( ownedMap, true );
  tDownward->ReplaceDiagonalValues( zeros );  // TODO directly build zeros in the function call.
//  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/tDownward-after.mat", *tDownward );

  int extracted = 0;
  std::vector< double > extractedValues( n );
  std::vector< int > extractedIndices( n );
  ghostInfo.ExtractGlobalRowCopy( curRank.get(), int( n ), extracted, extractedValues.data(), extractedIndices.data() );
  extractedValues.resize( extracted );
  extractedIndices.resize( extracted );
  GEOS_LOG_RANK( "extracted = " << extracted );

  std::vector< int > missingIndices;
  missingIndices.reserve( n );

//  int extractedEdges = 0;
//  std::vector< double > extractedEdgesValues( n );
//  std::vector< int > extractedEdgesIndices( n );
  FindGeometricalType const getGeomType( edgeOffset, faceOffset, cellOffset );
  using Geom = FindGeometricalType::Geom;

  Ownerships ownerships;
  for( int i = 0; i < extracted; ++i )  // TODO Can we `zip` instead?
  {
    int const & index = extractedIndices[i];

    MpiRank const owner = MpiRank{ int( extractedValues[i] ) };
    if( owner != curRank )
    {
      ownerships.neighbors.insert( owner );
    }

    switch( getGeomType( index ) )
    {
      case Geom::NODE:
      {
        NodeGlbIdx const ngi = getGeomType.as< NodeGlbIdx >( index );
        ownerships.nodes.emplace( ngi, owner );
        break;
      }
      case Geom::EDGE:
      {
        EdgeGlbIdx const egi = getGeomType.as< EdgeGlbIdx >( index );
        ownerships.edges.emplace( egi, owner );
        // TODO make all the following check in on time with sets comparison instead of "point-wise" comparisons.
        // TODO same for the faces and cells...
        if( owned.e2n.find( egi ) == owned.e2n.cend() and present.e2n.find( egi ) == present.e2n.cend() )
        {
          missingIndices.emplace_back( index );
        }
        break;
      }
      case Geom::FACE:
      {
        FaceGlbIdx const fgi = getGeomType.as< FaceGlbIdx >( index );
        ownerships.faces.emplace( fgi, owner );
        if( owned.f2e.find( fgi ) == owned.f2e.cend() and present.f2e.find( fgi ) == present.f2e.cend() )
        {
          missingIndices.emplace_back( index );
        }
        break;
      }
      case Geom::CELL:
      {
        CellGlbIdx const cgi = getGeomType.as< CellGlbIdx >( index );
        ownerships.cells.emplace( cgi, owner );
        if( owned.c2f.find( cgi ) == owned.c2f.cend() )
        {
          missingIndices.emplace_back( index );
        }
        break;
      }
      default:
      {
        GEOS_ERROR( "Internal error" );
      }
    }
  }

  GEOS_LOG( "Gathering the missing mappings." );

  std::size_t const numMissingIndices = std::size( missingIndices );
  std::size_t const numGlobalMissingIndices = MpiWrapper::sum( numMissingIndices );

//  Epetra_Map const missingIndicesMap( intConv< long long int >( numGlobalMissingIndices ), intConv< int >( numMissingIndices ), 0, comm );
  Epetra_Map const missingIndicesMap( intConv< int >( numGlobalMissingIndices ), intConv< int >( numMissingIndices ), 0, comm );

//  std::size_t const offset = MpiWrapper::prefixSum< std::size_t >( numMissingIndices );
  int const offset = *missingIndicesMap.MyGlobalElements();
//  GEOS_LOG_RANK( "numMissingIndices, numGlobalMissingIndices, offset = " << numMissingIndices << ", " << numGlobalMissingIndices << ", " << offset );

  Epetra_CrsMatrix missing( Epetra_DataAccess::Copy, missingIndicesMap, 1, true );
  for( std::size_t i = 0; i < numMissingIndices; ++i )
  {
    missing.InsertGlobalValues( offset + i, 1, ones.data(), &missingIndices[i] );
  }
  missing.FillComplete( ownedMap, missingIndicesMap );

  // TODO Don't put ones in downward: reassemble with the ordering (e.g. edges point to a first and last node)
  Epetra_CrsMatrix missingMappings( Epetra_DataAccess::Copy, missingIndicesMap, 1, false );
  EpetraExt::MatrixMatrix::Multiply( missing, false, *tDownward, true, missingMappings, false );
  missingMappings.FillComplete( ownedMap, missingIndicesMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/missingMappings.mat", missingMappings );

  MeshGraph ghosts;

  int ext = 0;
  std::vector< double > extValues( n );
  std::vector< int > extIndices( n );
//  double * extValues = nullptr;
//  int * extIndices = nullptr;

  for( int i = 0; i < int( numMissingIndices ); ++i )
  {
    missingMappings.ExtractGlobalRowCopy( offset + i, int( n ), ext, extValues.data(), extIndices.data() );
//    missingMappings.ExtractGlobalRowView( intConv< int >( offset + i ), ext, extValues );
//    missingMappings.ExtractGlobalRowView( intConv< int >( offset + i ), ext, extValues, extIndices );
//    Span< double > const s0( extValues, ext );
//    Span< int > const s1( extIndices, ext );
//    GEOS_LOG( "ext 0 = " << std::size( s0 ) );
//    GEOS_LOG( "ext 1 = " << std::size( s1 ) );
//    GEOS_LOG_RANK( "ext = " << ext );
    int const index = missingIndices[i];
    switch( getGeomType( index ) )
    {
      case Geom::EDGE:
      {
        GEOS_ASSERT_EQ( ext, 2 );  // TODO check that val != 0?
        EdgeGlbIdx const egi = getGeomType.as< EdgeGlbIdx >( index );
        NodeGlbIdx const ngi0 = getGeomType.as< NodeGlbIdx >( extIndices[0] );
        NodeGlbIdx const ngi1 = getGeomType.as< NodeGlbIdx >( extIndices[1] );
        ghosts.e2n[egi] = std::minmax( { ngi0, ngi1 } );
        break;
      }
      case Geom::FACE:
      {
        GEOS_ASSERT_EQ( ext, 4 );  // TODO temporary check for the dev
        FaceGlbIdx const fgi = getGeomType.as< FaceGlbIdx >( index );
        std::set< EdgeGlbIdx > & tmp = ghosts.f2e[fgi];
        for( int ii = 0; ii < ext; ++ii )
        {
          tmp.insert( EdgeGlbIdx{ intConv< globalIndex >( extIndices[ii] - edgeOffset ) } );
        }
        break;
      }
      case Geom::CELL:
      {
        GEOS_ASSERT_EQ( ext, 6 );  // TODO temporary check for the dev
        CellGlbIdx const cgi = getGeomType.as< CellGlbIdx >( index );
        std::set< FaceGlbIdx > & tmp = ghosts.c2f[cgi];
        for( int ii = 0; ii < ext; ++ii )
        {
          tmp.insert( FaceGlbIdx{ intConv< globalIndex >( extIndices[ii] - faceOffset ) } );
        }
        break;
      }
      default:
      {
        GEOS_ERROR( "Internal error." );
      }
    }
  }

  if( curRank == 2_mpi )
  {
    GEOS_LOG_RANK( "neighboringConnexions = " << json( ghosts ) );
    std::map< MpiRank, std::set< CellGlbIdx > > tmp;
    for( auto const & [geom, rk]: ownerships.cells )
    {
      tmp[rk].insert( geom );
    }
    for( auto const & [rk, geom]: tmp )
    {
      GEOS_LOG_RANK( "ghost geom = " << rk << ": " << json( geom ) );
    }
  }
  GEOS_LOG_RANK( "my final neighbors are " << json( ownerships.neighbors ) );

  return { std::move( ghosts ), std::move( ownerships ) };


  // TODO to build the edges -> nodes map containing the all the ghost (i.e. the ghosted ones as well),
  // Let's create a vector (or matrix?) full of ones where we have edges and multiply using the adjacency matrix.
}

/**
 * @brief
 * @tparam GI A global index type
 * @param m
 * @return
 */
template< class GI >
std::tuple< std::vector< MpiRank >, std::vector< GI > >
buildGhostRankAndL2G( std::map< GI, MpiRank > const & m )
{
  std::size_t const size = std::size( m );

  std::vector< MpiRank > ghostRank;
  ghostRank.reserve( size );
  std::vector< GI > l2g;
  l2g.reserve( size );
  for( auto const & [t, rank]: m )
  {
    ghostRank.emplace_back( rank );
    l2g.emplace_back( t );
  }

  return std::make_tuple( ghostRank, l2g );
}

void buildPods( MeshGraph const & owned,
                MeshGraph const & present,
                MeshGraph const & ghosts,
                Ownerships const & ownerships )
{
  std::size_t const numNodes = std::size( ownerships.nodes );
  std::size_t const numEdges = std::size( ownerships.edges );
  std::size_t const numFaces = std::size( ownerships.faces );

  auto [ghostRank, l2g] = buildGhostRankAndL2G( ownerships.edges );

  NodeMgrImpl const nodeMgr( NodeLocIdx{ intConv< localIndex >( numNodes ) } );
  EdgeMgrImpl const edgeMgr( EdgeLocIdx{ intConv< localIndex >( numEdges ) }, std::move( ghostRank ), std::move( l2g ) );
  FaceMgrImpl const faceMgr( FaceLocIdx{ intConv< localIndex >( numFaces ) } );
}

std::unique_ptr< generators::MeshMappings > doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                                                              std::set< MpiRank > const & neighbors )
{
  auto const [buckets, offsets] = doTheNewGlobalNumbering( mesh, neighbors );

  // Now we exchange the data with our neighbors.
  MpiRank const curRank{ MpiWrapper::commRank() };

  GEOS_LOG_RANK( "offsets on rank " << curRank << " -> " << json( offsets ) );

  MpiRank const nextRank = curRank + 1_mpi;
  MaxGlbIdcs const matrixOffsets = gatherOffset( mesh, offsets.edges.at( { nextRank } ) - 1_egi, offsets.faces.at( { nextRank } ) - 1_fgi );
  std::cout << "matrixOffsets on rank " << curRank << " -> " << json( matrixOffsets ) << std::endl;

  auto const [owned, present] = buildMeshGraph( mesh, buckets, offsets, curRank );  // TODO change into buildOwnedMeshGraph?
//  if( curRank == 1_mpi )
//  {
//    GEOS_LOG_RANK( "My owned is " << json( owned ) );
//    GEOS_LOG_RANK( "My present is " << json( present ) );
//  }
//  MpiWrapper::barrier();

  auto const [ghosts, ownerships] = assembleAdjacencyMatrix( owned, present, matrixOffsets, curRank );

  buildPods( owned, present, ghosts, ownerships );

  return {};
}

std::unique_ptr< generators::MeshMappings > doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                                                              std::set< int > const & neighbors )
{
  std::set< MpiRank > neighbors_;
  for( int const & rank: neighbors )
  {
    neighbors_.insert( MpiRank{ rank } );
  }
  GEOS_LOG_RANK( "my initial neighbors are " << json( neighbors_ ) );

  return doTheNewGhosting( mesh, neighbors_ );
}

}  // end of namespace geos::ghosting
