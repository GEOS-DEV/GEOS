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

#include "BuildPods.hpp"

#include "Pods.hpp"

#include "Indices.hpp"

#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_Comm.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>

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

void to_json( json & j,
              const EdgeInfo & v )  // For display
{
  j = json{ { "index", v.index },
            { "start", v.start } };
}

void to_json( json & j,
              const FaceInfo & v )  // For display
{
  j = json{ { "index",     v.index },
            { "isFlipped", v.isFlipped },
            { "start",     v.start } };
}

void to_json( json & j,
              const MeshGraph & v )  // For display
{
  j = json{ { "c2f", v.f2e },
            { "f2e", v.f2e },
            { "e2n", v.e2n },
            { "n2pos", v.n2pos } };
}


std::map< NodeGlbIdx, vtkIdType > buildNgiToVtk( vtkSmartPointer< vtkDataSet > mesh )
{
  std::map< NodeGlbIdx, vtkIdType > res;

  vtkIdTypeArray const * globalPtIds = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );

  for( vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i )
  {
    res[NodeGlbIdx{ globalPtIds->GetValue( i ) }] = i;
  }

  return res;
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

  std::map< NodeGlbIdx, vtkIdType > const n2vtk = buildNgiToVtk( mesh );
  for( auto const & [ranks, nodes]: buckets.nodes )
  {
    std::map< NodeGlbIdx, std::array< double, 3 > > & n2pos = isCurrentRankOwning( ranks ) ? owned.n2pos : present.n2pos;
    for( NodeGlbIdx const & ngi: nodes )
    {
      double const * pos = mesh->GetPoint( n2vtk.at( ngi ) );
      n2pos[ngi] = { pos[0], pos[1], pos[2] };
    }
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
        EdgeInfo const info = { n2e.at( p1 ), std::uint8_t{ p0 != p1 } };
        f2e[fgi].emplace_back( info );
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
      bool isFlipped;
      std::uint8_t start;
      std::vector< NodeGlbIdx > const reorderedFaceNodes = reorderFaceNodes( faceNodes, isFlipped, start );
      owned.c2f[gci].emplace_back( FaceInfo{ n2f.at( reorderedFaceNodes ), isFlipped, start } );
    }
  }

  return { std::move( owned ), std::move( present ) };
}

//using TriLocIdx = fluent::NamedType< localIndex, struct TagTriLocIdx, fluent::Arithmetic, fluent::ImplicitlyConvertibleTo >;
//using TriGlbIdx = fluent::NamedType< globalIndex, struct TagTriGlbIdx, fluent::Arithmetic, fluent::ImplicitlyConvertibleTo >;
//template<>
//struct ::std::is_signed< TriLocIdx > : ::std::is_signed< typename TriLocIdx::UnderlyingType >
//{
//};

using TriLocIdx = localIndex;
using TriGlbIdx = globalIndex;

using TriMap = Tpetra::Map< TriLocIdx, TriGlbIdx >;
// We use matrices of integers for some of our usages
// But it seems that Tpetra wants *signed* integers (at least the way we install it).
// Consider using *unsigned* integers if you manage to make it work/link.
using TriScalarInt = std::int32_t;
using TriCrsMatrix = Tpetra::CrsMatrix< TriScalarInt , TriLocIdx, TriGlbIdx >;
using TriDblCrsMatrix = Tpetra::CrsMatrix< double, TriLocIdx, TriGlbIdx >;
using TriComm = Teuchos::Comm< int >;

template< typename T, typename... ARGS >
Teuchos::RCP< T > make_rcp( ARGS && ... args )
{
  return Teuchos::rcp( new T( std::forward< ARGS >( args )... ) );
}

void to_json( json & j,
              const GhostSend & v )
{
  j = json{ { "nodes", v.nodes },
            { "edges", v.edges },
            { "faces", v.faces },
            { "cells", v.cells } };
}

void to_json( json & j,
              const GhostRecv & v )
{
  j = json{ { "nodes", v.nodes },
            { "edges", v.edges },
            { "faces", v.faces },
            { "cells", v.cells } };
}


TriCrsMatrix multiply( int commSize,
                       TriCrsMatrix const & indicator,
                       TriCrsMatrix const & upward )
{
  Teuchos::RCP< TriMap const > ownedMap = upward.getRowMap();
  Teuchos::RCP< TriMap const > mpiMap = indicator.getRangeMap();

  // Upward (n -> e -> f -> c)

  TriCrsMatrix result_u0_0( ownedMap, commSize );
  Tpetra::MatrixMatrix::Multiply( upward, false, indicator, true, result_u0_0, false );
  result_u0_0.fillComplete( mpiMap, ownedMap );
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/result-0.mat", result_u0_0 );

  TriCrsMatrix result_u0_1( ownedMap, commSize );
  Tpetra::MatrixMatrix::Multiply( upward, false, result_u0_0, false, result_u0_1, false );
  result_u0_1.fillComplete( mpiMap, ownedMap );
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/result-1.mat", result_u0_1 );

  TriCrsMatrix result_u0_2( ownedMap, commSize );
  Tpetra::MatrixMatrix::Multiply( upward, false, result_u0_1, false, result_u0_2, false );
  result_u0_2.fillComplete( mpiMap, ownedMap );
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/result-2.mat", result_u0_2 );

  // Downward (c -> f -> e -> n)

  TriCrsMatrix result_d0_0( ownedMap, commSize );
  Tpetra::MatrixMatrix::Multiply( upward, true, result_u0_2, false, result_d0_0, false );
  result_d0_0.fillComplete( mpiMap, ownedMap );
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/result-4.mat", result_d0_0 );

  TriCrsMatrix result_d0_1( ownedMap, commSize );
  Tpetra::MatrixMatrix::Multiply( upward, true, result_d0_0, false, result_d0_1, false );
  result_d0_1.fillComplete( mpiMap, ownedMap );
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/result-5.mat", result_d0_1 );

  TriCrsMatrix result_d0_2( ownedMap, commSize );
  Tpetra::MatrixMatrix::Multiply( upward, true, result_d0_1, false, result_d0_2, false );
  result_d0_2.fillComplete( mpiMap, ownedMap );
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/result-6.mat", result_d0_2 );

  MpiWrapper::barrier();

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

  FindGeometricalType( NodeGlbIdx const & maxNodeGlbIdx,
                       EdgeGlbIdx const & maxEdgeGlbIdx,
                       FaceGlbIdx const & maxFaceGlbIdx,
                       CellGlbIdx const & maxCellGlbIdx )
    : m_edgeOffset( intConv< int >( maxNodeGlbIdx.get() + 1 ) ),
      m_faceOffset( intConv< int >( m_edgeOffset + TriGlbIdx( maxEdgeGlbIdx.get() ) + TriGlbIdx( 1 ) ) ),
      m_cellOffset( intConv< int >( m_faceOffset + TriGlbIdx( maxFaceGlbIdx.get() ) + TriGlbIdx( 1 ) ) ),
      m_numEntries( intConv< int >( m_cellOffset + TriGlbIdx( maxCellGlbIdx.get() ) + TriGlbIdx( 1 ) ) )
  { }

  [[nodiscard]] TriGlbIdx numEntries() const
  {
    return m_numEntries;
  }

  [[nodiscard]] Geom getGeometricalType( TriGlbIdx const & index ) const
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

  [[nodiscard]] NodeGlbIdx toNodeGlbIdx( TriGlbIdx const & index ) const
  {
    return NodeGlbIdx{ intConv< NodeGlbIdx::UnderlyingType >( index ) };
  }

  [[nodiscard]] EdgeGlbIdx toEdgeGlbIdx( TriGlbIdx const & index ) const
  {
    return EdgeGlbIdx{ intConv< EdgeGlbIdx::UnderlyingType >( index - m_edgeOffset ) };
  }

  [[nodiscard]] FaceGlbIdx toFaceGlbIdx( TriGlbIdx const & index ) const
  {
    return FaceGlbIdx{ intConv< FaceGlbIdx::UnderlyingType >( index - m_faceOffset ) };
  }

  [[nodiscard]] CellGlbIdx toCellGlbIdx( TriGlbIdx const & index ) const
  {
    return CellGlbIdx{ intConv< CellGlbIdx::UnderlyingType >( index - m_cellOffset ) };
  }

  [[nodiscard]] TriGlbIdx fromNodeGlbIdx( NodeGlbIdx const & ngi ) const
  {
    return TriGlbIdx( ngi.get() );
  }

  [[nodiscard]] TriGlbIdx fromEdgeGlbIdx( EdgeGlbIdx const & egi ) const
  {
    return TriGlbIdx( egi.get() ) + m_edgeOffset;
  }

  [[nodiscard]] TriGlbIdx fromFaceGlbIdx( FaceGlbIdx const & fgi ) const
  {
    return TriGlbIdx( fgi.get() ) + m_faceOffset;
  }

  [[nodiscard]] TriGlbIdx fromCellGlbIdx( CellGlbIdx const & cgi ) const
  {
    return TriGlbIdx( cgi.get() ) + m_cellOffset;
  }

private:
  TriGlbIdx const m_edgeOffset;
  TriGlbIdx const m_faceOffset;
  TriGlbIdx const m_cellOffset;
  TriGlbIdx const m_numEntries;
};

template< std::uint8_t N >
std::size_t encode( std::size_t const & basis,
                    std::array< std::size_t, N > const & array )
{
  std::size_t result = 0;
  for( auto i = 0; i < N; ++i )
  {
    result *= basis;
    result += array[i];
  }
  return result;
}

template< std::uint8_t N >
std::array< std::size_t, N > decode( std::size_t const & basis,
                                     std::size_t const & input )
{
  std::array< std::size_t, N > result, pows;

  for( std::size_t i = 0, j = N - 1, pow = 1; i < N; ++i, --j )
  {
    pows[j] = pow;
    pow *= basis;
  }

  std::size_t dec = input;
  for( auto i = 0; i < N; ++i )
  {
    std::lldiv_t const div = std::lldiv( dec, pows[i] );
    result[i] = div.quot;
    dec = div.rem;
  }
  return result;
}

struct Adjacency
{
  std::vector< TriGlbIdx > ownedNodesIdcs;
  std::vector< TriGlbIdx > ownedGlbIdcs;
  std::vector< TriGlbIdx > otherGlbIdcs;
  Teuchos::Array< std::size_t > numEntriesPerRow;
  std::vector< std::vector< TriGlbIdx > > indices;
  std::vector< std::vector< TriScalarInt > > values;
};

Adjacency buildAdjacency( MeshGraph const & owned,
                          MeshGraph const & present,
                          FindGeometricalType const & convert )
{
  Adjacency adjacency;

  std::size_t const numOwned = std::size( owned.n2pos ) + std::size( owned.e2n ) + std::size( owned.f2e ) + std::size( owned.c2f );
  std::size_t const numOther = std::size( present.n2pos ) + std::size( present.e2n ) + std::size( present.f2e );

  // Aliases
  std::vector< TriGlbIdx > & ownedNodesIdcs = adjacency.ownedNodesIdcs;
  std::vector< TriGlbIdx > & ownedGlbIdcs = adjacency.ownedGlbIdcs;
  std::vector< TriGlbIdx > & otherGlbIdcs = adjacency.otherGlbIdcs;
  Teuchos::Array< std::size_t > & numEntriesPerRow = adjacency.numEntriesPerRow;
  std::vector< std::vector< TriGlbIdx > > & indices = adjacency.indices;
  std::vector< std::vector< TriScalarInt > > & values = adjacency.values;

  ownedNodesIdcs.reserve( std::size( owned.n2pos ) );
  ownedGlbIdcs.reserve( numOwned );
  otherGlbIdcs.reserve( numOther );
  numEntriesPerRow.reserve( numOwned );
  indices.reserve( numOwned );
  values.reserve( numOwned );

  for( auto const & [ngi, _]: present.n2pos )
  {
    otherGlbIdcs.emplace_back( convert.fromNodeGlbIdx( ngi ) );
  }
  for( auto const & [egi, _]: present.e2n )
  {
    otherGlbIdcs.emplace_back( convert.fromEdgeGlbIdx( egi ) );
  }
  for( auto const & [fgi, _]: present.f2e )
  {
    otherGlbIdcs.emplace_back( convert.fromFaceGlbIdx( fgi ) );
  }
  std::sort( std::begin( otherGlbIdcs ), std::end( otherGlbIdcs ) );
  GEOS_ASSERT_EQ( numOther, std::size( otherGlbIdcs ) );

  for( auto const & [ngi, _]: owned.n2pos )
  {
    // Nodes depend on no other geometrical entity,
    // so we only have one entry `1` in the diagonal of the matrix,
    // because we add the identity to the adjacency matrix.
    TriGlbIdx const i = convert.fromNodeGlbIdx( ngi );
    ownedNodesIdcs.emplace_back( i );
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.push_back( 0 + 1 );  // `+1` comes from the diagonal
    indices.emplace_back( 1, ownedGlbIdcs.back() );
    values.emplace_back( 1, ownedGlbIdcs.back() );
  }
  for( auto const & [egi, nodes]: owned.e2n )
  {
    // Edges always depend on exactly 2 nodes, so this value can be hard coded.
    // Also, edges have two different direction (starting from one point of the other).
    // To keep the symmetry with the faces (see below),
    // - we store the local index (0 or 1) to express this information.
    // - we store the number of underlying nodes (2) on the diagonal.
    size_t constexpr numNodes = std::tuple_size_v< decltype( nodes ) >;

    ownedGlbIdcs.emplace_back( convert.fromEdgeGlbIdx( egi ) );

    numEntriesPerRow.push_back( numNodes + 1 );  // `+1` comes from the diagonal
    indices.emplace_back( std::vector< TriGlbIdx >{ convert.fromNodeGlbIdx( std::get< 0 >( nodes ) ),
                                                    convert.fromNodeGlbIdx( std::get< 1 >( nodes ) ),
                                                    ownedGlbIdcs.back() } );
    // Note that when storing a value that can be `0`, we always add `1`,
    // (for edges but also later for faces and cells),
    // to be sure that there will always be some a noticeable figure where we need one.
    values.emplace_back( std::vector< TriScalarInt >{ 0 + 1, 1 + 1, numNodes } );
  }
  for( auto const & [fgi, edges]: owned.f2e )
  {
    // Faces point to multiple edges, but their edges can be flipped w.r.t. the canonical edges
    // (ie minimal node global index comes first).
    // In order not to lose this information (held by an `EdgeInfo` instance), we serialise
    // - the `EdgeInfo` instance,
    // - plus the order in which the edges are appearing,
    // as an integer and store it as an entry in the matrix.
    // On the diagonal, we'll store the number of edges of the face.
    std::size_t const numEdges = std::size( edges );

    ownedGlbIdcs.emplace_back( convert.fromFaceGlbIdx( fgi ) );

    numEntriesPerRow.push_back( numEdges + 1 );  // `+1` comes from the diagonal
    std::vector< TriGlbIdx > & ind = indices.emplace_back( numEntriesPerRow.back() );
    std::vector< TriScalarInt > & val = values.emplace_back( numEntriesPerRow.back() );
    for( std::size_t i = 0; i < numEdges; ++i )
    {
      EdgeInfo const & edgeInfo = edges[i];
      ind[i] = convert.fromEdgeGlbIdx( edgeInfo.index );
      std::size_t const v = 1 + encode< 2 >( numEdges, { edgeInfo.start, i } );
      val[i] = intConv< TriScalarInt >( v );
    }
    ind.back() = ownedGlbIdcs.back();
    val.back() = intConv< TriScalarInt >( numEdges );
  }
  for( auto const & [cgi, faces]: owned.c2f )
  {
    // The same comment as for faces and edges applies for cells and faces (see above).
    // The main differences being that
    // - the `FaceInfo` as a little more information,
    // - the diagonal stores the type of the cell which conveys the number of face, nodes, ordering...
    std::size_t const numFaces = std::size( faces );

    ownedGlbIdcs.emplace_back( convert.fromCellGlbIdx( cgi ) );

    numEntriesPerRow.push_back( numFaces + 1 );  // `+1` comes from the diagonal
    std::vector< TriGlbIdx > & ind = indices.emplace_back( numEntriesPerRow.back() );
    std::vector< TriScalarInt > & val = values.emplace_back( numEntriesPerRow.back() );
    for( std::size_t i = 0; i < numFaces; ++i )
    {
      FaceInfo const & faceInfo = faces[i];
      ind[i] = convert.fromFaceGlbIdx( faceInfo.index );
      std::size_t const v = 1 + encode< 3 >( numFaces, { faceInfo.isFlipped, faceInfo.start, i } );
      val[i] = intConv< TriScalarInt >( v );
    }
    ind.back() = ownedGlbIdcs.back();
    val.back() = intConv< TriScalarInt >( numFaces );  // TODO This should be Hex and the not the number of faces...
  }
  std::sort( std::begin( ownedGlbIdcs ), std::end( ownedGlbIdcs ) );

  GEOS_ASSERT_EQ( numOwned, std::size( ownedGlbIdcs ) );
  GEOS_ASSERT_EQ( numOwned, intConv< std::size_t >( numEntriesPerRow.size() ) );
  GEOS_ASSERT_EQ( numOwned, std::size( indices ) );
  GEOS_ASSERT_EQ( numOwned, std::size( values ) );
  for( std::size_t i = 0; i < numOwned; ++i )
  {
    GEOS_ASSERT_EQ( indices[i].size(), std::size_t( numEntriesPerRow[i] ) );
  }

  return adjacency;
}

std::tuple< MeshGraph, GhostRecv, GhostSend > performGhosting( MeshGraph const & owned,
                                                               MeshGraph const & present,
                                                               MaxGlbIdcs const & gis,
                                                               MpiRank curRank )
{
  using Teuchos::RCP;

  FindGeometricalType const convert( gis.nodes, gis.edges, gis.faces, gis.cells );
  using Geom = FindGeometricalType::Geom;

  std::size_t const n = convert.numEntries();  // Total number of entries in the graph.

  Adjacency const adjacency = buildAdjacency( owned, present, convert );
  std::size_t const numOwned = std::size( adjacency.ownedGlbIdcs );
  std::size_t const numOther = std::size( adjacency.otherGlbIdcs );

  // Aliases
  std::vector< TriGlbIdx > const & ownedNodesIdcs = adjacency.ownedNodesIdcs;
  std::vector< TriGlbIdx > const & ownedGlbIdcs = adjacency.ownedGlbIdcs;
  std::vector< TriGlbIdx > const & otherGlbIdcs = adjacency.otherGlbIdcs;
  Teuchos::Array< std::size_t > const & numEntriesPerRow = adjacency.numEntriesPerRow;
  std::vector< std::vector< TriGlbIdx > > const & indices = adjacency.indices;
  std::vector< std::vector< TriScalarInt > > const  & values = adjacency.values;

  RCP< TriComm const > const comm = make_rcp< Teuchos::MpiComm< int > const >( MPI_COMM_GEOSX );
  auto const ownedMap = make_rcp< TriMap const >( Tpetra::global_size_t( n ), ownedGlbIdcs.data(), TriLocIdx( numOwned ), TriGlbIdx{ 0 }, comm );

  // The `upward` matrix offers a representation of the graph connections.
  // The matrix is square. Each size being the number of geometrical entities.
  TriCrsMatrix upward( ownedMap, numEntriesPerRow() );

  for( std::size_t i = 0; i < numOwned; ++i )
  {
    std::vector< TriGlbIdx > const & rowIndices = indices[i];
    std::vector< TriScalarInt > const & rowValues = values[i];
    GEOS_ASSERT_EQ( std::size( rowIndices ), numEntriesPerRow[i] );
    GEOS_ASSERT_EQ( std::size( rowValues ), numEntriesPerRow[i] );
    upward.insertGlobalValues( ownedGlbIdcs[i], TriLocIdx( std::size( rowIndices ) ), rowValues.data(), rowIndices.data() );
  }

  upward.fillComplete( ownedMap, ownedMap );
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/upward.mat", upward );

  int const commSize( MpiWrapper::commSize() );

  // Now let's build the domain indicator matrix.
  // It's rectangular, number of rows being the number of MPI ranks,
  // number of columns being the number of nodes in the mesh graph,
  // ie the number of geometrical entities in the mesh.
  // It contains one for every time the geometrical entity is present on the rank.
  // By present, we do not meant that the ranks owns it: just that it's already available on the rank.
  auto const mpiMap = make_rcp< TriMap const >( Tpetra::global_size_t( commSize ), TriGlbIdx{ 0 }, comm );  // Let the current rank get the appropriate index in this map.
  TriCrsMatrix indicator( mpiMap, numOwned + numOther );

  std::vector< TriScalarInt > const ones( std::max( { numOwned, numOther, std::size_t( 1 ) } ), 1 );
  indicator.insertGlobalValues( TriGlbIdx( curRank.get() ), TriLocIdx( numOwned ), ones.data(), ownedGlbIdcs.data() );
  indicator.insertGlobalValues( TriGlbIdx( curRank.get() ), TriLocIdx( numOther ), ones.data(), otherGlbIdcs.data() );
  indicator.fillComplete( ownedMap, mpiMap );

  // The `ownership` matrix is a diagonal square matrix.
  // Each size being the number of geometrical entities.
  // The value of the diagonal term will be the owning rank.
  // By means of matrices multiplications, it will be possible to exchange the ownerships information across the ranks.
  // TODO Could we use an Epetra_Vector as a diagonal matrix?
  std::vector< TriScalarInt > myRank( 1, curRank.get() );
  TriCrsMatrix ownership( ownedMap, 1 );
  for( TriGlbIdx const & i: ownedGlbIdcs )
  {
    ownership.insertGlobalValues( i, TriLocIdx{ 1 }, myRank.data(), &i );
  }
  ownership.fillComplete( ownedMap, ownedMap );
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/ownership.mat", ownership );

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "indicator.NumGlobalCols() = " << indicator.getGlobalNumCols() );
    GEOS_LOG_RANK( "indicator.NumGlobalRows() = " << indicator.getGlobalNumRows() );
    GEOS_LOG_RANK( "ownership.NumGlobalCols() = " << ownership.getGlobalNumCols() );
    GEOS_LOG_RANK( "ownership.NumGlobalRows() = " << ownership.getGlobalNumRows() );
  }
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/indicator.mat", indicator );

  // The `ghostingFootprint` matrix is rectangular,
  // the number of columns being the number of MPI ranks,
  // the number of rows being the number of nodes in the mesh graph.
  //
  // For any MPI rank (ie column), a non-zero entry means
  // that the corresponding geometrical entry has to be eventually present onto the rank
  // for the ghosting to be effective.
  //
  // From `ghostingFootprint` we can extract where the current rank has to send any owned graph node.
  TriCrsMatrix ghostingFootprint( multiply( commSize, indicator, upward ) );
  ghostingFootprint.resumeFill();
  ghostingFootprint.setAllToScalar( 1. );
  ghostingFootprint.fillComplete( mpiMap, ownedMap );
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/ghostingFootprint.mat", ghostingFootprint );
  // We put `1` everywhere there's a non-zero entry, so we'll be able to compose with the `ownership` matrix.

  // FIXME TODO WARNING From `ghostingFootprint` extract where I have to send
  // Then, perform the ghostingFootprint * ownership matrix multiplication,
  // so I get to know from whom I'll receive.

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "ghosted->NumGlobalCols() = " << ghostingFootprint.getGlobalNumCols() );
    GEOS_LOG_RANK( "ghosted->NumGlobalRows() = " << ghostingFootprint.getGlobalNumRows() );
  }

  // The `ghostExchange` matrix is rectangular,
  // the number of columns being the number of nodes in the mesh graph,
  // the number of rows being the number of MPI ranks.
  //
  // As the result of the multiplication between the `ghostingFootprint` matrix and the `ownership` matrix,
  // for each row owned (ie at the current MPI rank index),
  // the value of the `ghostExchange` matrix term will provide the actual owning rank for all the .
  //
  // From `ghostExchange` we can extract which other rank will send to the current rank any graph node.
  TriCrsMatrix ghostExchange( mpiMap, 10000 ); // TODO having `0` there should be working!
  Tpetra::MatrixMatrix::Multiply( ghostingFootprint, true, ownership, false, ghostExchange, false );
  ghostExchange.fillComplete( ownedMap, mpiMap );

  // TODO Do I have to work with `ghostingFootprint` if I already have `ghostExchange` which may convey more information?
  // TODO Maybe because of the ownership of the ranks: one is also the "scaled" transposed of the other.

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "ghostExchange->NumGlobalCols() = " << ghostExchange.getGlobalNumCols() );
    GEOS_LOG_RANK( "ghostExchange->NumGlobalRows() = " << ghostExchange.getGlobalNumRows() );
  }
  Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/ghostInfo.mat", ghostExchange );

  std::size_t extracted = 0;
  Teuchos::Array< TriScalarInt > extractedValues;
  Teuchos::Array< TriGlbIdx > extractedIndices;

  GhostSend send;
  for( TriGlbIdx const & index: ownedGlbIdcs )
  {
    std::size_t const length = ghostingFootprint.getNumEntriesInGlobalRow( index );
    extractedValues.resize( length );
    extractedIndices.resize( length );
    ghostingFootprint.getGlobalRowCopy( index, extractedIndices(), extractedValues(), extracted );
    GEOS_ASSERT_EQ( extracted, length );

    std::set< MpiRank > neighbors;
    for( std::size_t i = 0; i < extracted; ++i )
    {
      MpiRank const rank( extractedIndices[i] );
      if( rank != curRank )
      {
        neighbors.insert( rank );
      }
    }

    if( std::empty( neighbors ) )
    {
      continue;
    }

    switch( convert.getGeometricalType( index ) )
    {
      case Geom::NODE:
      {
        send.nodes.emplace( convert.toNodeGlbIdx( index ), std::move( neighbors ) );
        break;
      }
      case Geom::EDGE:
      {
        send.edges.emplace( convert.toEdgeGlbIdx( index ), std::move( neighbors ) );
        break;
      }
      case Geom::FACE:
      {
        send.faces.emplace( convert.toFaceGlbIdx( index ), std::move( neighbors ) );
        break;
      }
      case Geom::CELL:
      {
        send.cells.emplace( convert.toCellGlbIdx( index ), std::move( neighbors ) );
        break;
      }
      default:
      {
        GEOS_ERROR( "Internal error" );
      }
    }
  }

  std::size_t const numNeededIndices = ghostExchange.getNumEntriesInGlobalRow( TriGlbIdx( curRank.get() ) );
  extractedValues.resize( numNeededIndices );
  extractedIndices.resize( numNeededIndices );
  ghostExchange.getGlobalRowCopy( TriGlbIdx( curRank.get() ), extractedIndices(), extractedValues(), extracted );
  GEOS_ASSERT_EQ( extracted, numNeededIndices );

  std::set< TriGlbIdx > const allNeededIndices( std::cbegin( extractedIndices ), std::cend( extractedIndices ) );
  std::set< TriGlbIdx > receivedIndices;  // The graph nodes that my neighbors will send me.
  std::set_difference( std::cbegin( allNeededIndices ), std::cend( allNeededIndices ),
                       std::cbegin( ownedGlbIdcs ), std::cend( ownedGlbIdcs ),
                       std::inserter( receivedIndices, std::end( receivedIndices ) ) );
  std::vector< TriGlbIdx > notPresentIndices;  // The graphs nodes that are nor owned neither present by/on the current rank.
  std::set_difference( std::cbegin( receivedIndices ), std::cend( receivedIndices ),
                       std::cbegin( otherGlbIdcs ), std::cend( otherGlbIdcs ),
                       std::back_inserter( notPresentIndices ) );
  std::vector< TriGlbIdx > notPresentNodes;
  std::copy_if( std::cbegin( notPresentIndices ), std::cend( notPresentIndices ),
                std::back_inserter( notPresentNodes ), [&]( TriGlbIdx const & i )
                {
                  return convert.getGeometricalType( i ) == Geom::NODE;
                } );
  GEOS_ASSERT_EQ( std::size( allNeededIndices ), numNeededIndices );

  GhostRecv recv;
  for( std::size_t i = 0; i < extracted; ++i )
  {
    TriGlbIdx const & index = extractedIndices[i];
    if( receivedIndices.find( index ) == std::cend( receivedIndices ) )  // TODO make a map `receivedIndices -> mpi rank`
    {
      continue;
    }

    MpiRank const sender( extractedValues[i] );
    switch( convert.getGeometricalType( index ) )
    {
      case Geom::NODE:
      {
        recv.nodes[convert.toNodeGlbIdx( index )] = sender;
        break;
      }
      case Geom::EDGE:
      {
        recv.edges[convert.toEdgeGlbIdx( index )] = sender;
        break;
      }
      case Geom::FACE:
      {
        recv.faces[convert.toFaceGlbIdx( index )] = sender;
        break;
      }
      case Geom::CELL:
      {
        recv.cells[convert.toCellGlbIdx( index )] = sender;
        break;
      }
      default:
      {
        GEOS_ERROR( "Internal error" );
      }
    }
  }

//  if( curRank == 1_mpi or curRank == 2_mpi )
//  {
//    GEOS_LOG_RANK( "recv.edges = " << json( recv.edges ) );
//    GEOS_LOG_RANK( "send.edges = " << json( send.edges ) );
//  }

  // At this point, each rank knows what it has to send to and what it is going to receive from the other ranks.
  //
  // The remaining action is about
  // - retrieving the additional graph information
  //   for the new geometrical quantities that will be sent by the neighbors.
  // - retrieving the positions of the ghosted nodes:
  //   knowing the index of the nodes is not enough.
  //
  // In order to do that, we build the `missingIndicator` matrix, which is rectangular:
  // - The number of columns is the number of graph nodes in the mesh graph.
  // - The number of rows is the total number of graph nodes that are missing on ranks.
  //   Note that the same quantity can be missing on multiple ranks, and that's handled.
  // - The non-zero terms equal `1`, meaning that a given quantity (column) is missing on a given rank (row).
  //
  // Combining `missingIndicator` with the adjacency matrix which conveys a lot of connections information,
  // we'll be able to create the final `missingIndices` matrix,
  // with `range` and `domain` maps (MPI ranks ownerships and offsets) are appropriately defined
  // such that the rows will be available to any ranks that need them (nothing more, nothing less).
  //
  // To get the `missingNodePos` matrix which will convey the ghosted nodes positions that are missing on the rank,
  // we first create a specific `missingNodesIndicator` matrix, which is alike `missingIndicator` but only for the nodes.
  // We could have used `missingIndicator` and ditched the superfluous information,
  // but the node positions are reals, where the connections are integers.
  // As long as we're using `Epetra`, where the type of matrix term is `double`, this makes no critical difference.
  // But when we switch to `Tpetra`, where we can select the type of the matrix term (and therefore use integers where we need integers),
  // this may have become an issue.
  // So the work has been done to separate `missingIndicator` and `missingNodesIndicator`.
  //
  // `missingNodesIndicator` is a rectangular like `missingIndicator`:
  // - The number of columns is the number of graph nodes in the mesh graph
  //   that actually are physical nodes in the mesh.
  // - The number of rows is the total number of graph nodes (that actually are physical nodes in the mesh) that are missing on ranks.
  // - The non-zero terms equal `1`, meaning that a given quantity (column) is missing on a given rank (row).
  std::size_t const numLocMissingCompound[2] = { std::size( notPresentIndices ), std::size( notPresentNodes ) };
  std::size_t numGlbMissingCompound[2] = { 0, 0 };
  MpiWrapper::allReduce( numLocMissingCompound, numGlbMissingCompound, 2, MPI_SUM );
  std::size_t const & numLocMissingIndices = numLocMissingCompound[0];
  std::size_t const & numLocMissingNodes = numLocMissingCompound[1];
  Tpetra::global_size_t const numGlbMissingIndices = numGlbMissingCompound[0];
  Tpetra::global_size_t const numGlbMissingNodes = numGlbMissingCompound[1];

  std::size_t const numGlbNodes = gis.nodes.get() + 1;  // TODO Use strongly typed integers
  std::size_t const numOwnedNodes = std::size( ownedNodesIdcs );  // TODO Use strongly typed integers

  auto const missingIndicesMap = make_rcp< TriMap const >( numGlbMissingIndices, numGlbMissingIndices, TriGlbIdx{ 0 }, comm );
  auto const missingNodesMap = make_rcp< TriMap const >( numGlbMissingNodes, numLocMissingNodes, TriGlbIdx{ 0 }, comm );

  // Following information is needed to compute the global row.
  TriGlbIdx const missingIndicesOffset = missingIndicesMap->getGlobalElement( TriLocIdx{ 0 } );  // TODO Hocus Pocus
  TriGlbIdx const missingNodesOffset = missingNodesMap->getGlobalElement( TriLocIdx{ 0 } );  // TODO Hocus Pocus

  // Indicator matrix for the all the missing quantities (nodes, edges, faces and cells).
  TriCrsMatrix missingIndicator( missingIndicesMap, 1 );
  for( std::size_t i = 0; i < numLocMissingIndices; ++i )
  {
    missingIndicator.insertGlobalValues( missingIndicesOffset + TriGlbIdx( i ), TriLocIdx( 1 ), ones.data(), &notPresentIndices[i] );
  }
  missingIndicator.fillComplete( ownedMap, missingIndicesMap );

  // `missingIndices` will contain the missing connectivity information.
  TriCrsMatrix missingIndices( missingIndicesMap, 1000 ); // TODO having `0` there should be working!
  Tpetra::MatrixMatrix::Multiply( missingIndicator, false, upward, false, missingIndices, false );
  missingIndices.fillComplete( ownedMap, missingIndicesMap );

  // Indicator matrix only for the nodes.
  // Note that it should be an integer matrix full of ones,
  // but it's a double matrix because it's going to be multiplied a double matrix.
  TriDblCrsMatrix missingNodesIndicator( missingNodesMap, 1 );
  for( std::size_t i = 0; i < numLocMissingNodes; ++i )
  {
    double const o = 1.;
    missingNodesIndicator.insertGlobalValues( missingNodesOffset + TriGlbIdx( i ), TriLocIdx( 1 ), &o, &notPresentNodes[i] );
  }
  auto const ownedNodesMap = make_rcp< TriMap const >( Tpetra::global_size_t( numGlbNodes ), ownedNodesIdcs.data(), TriLocIdx( numOwnedNodes ), TriGlbIdx{ 0 }, comm );
  missingNodesIndicator.fillComplete( ownedNodesMap, missingNodesMap );

  // The `nodePositions` matrix is rectangular.
  // - Its number of rows is the total number of nodes in the mesh.
  // - Its number of columns is 3: the x, y, and z coordinates of the nodes.
  TriDblCrsMatrix nodePositions( ownedNodesMap, 3 );
  std::vector< TriGlbIdx > const zot{ TriGlbIdx( 0 ), TriGlbIdx( 1 ), TriGlbIdx( 2 ) };  // zot: zero, one, two.
  for( auto const & [ngi, pos]: owned.n2pos )
  {
    nodePositions.insertGlobalValues( convert.fromNodeGlbIdx( ngi ), TriLocIdx( 3 ), pos.data(), zot.data() );
  }
  auto const threeMap = make_rcp< TriMap const >( Tpetra::global_size_t( 3 ), TriGlbIdx( 0 ), comm );
  nodePositions.fillComplete( threeMap, ownedNodesMap );

  // `missingNodePos` will contain the missing node positions.
  TriDblCrsMatrix missingNodePos( missingNodesMap, 1000 ); // TODO having `0` there should be working!
  Tpetra::MatrixMatrix::Multiply( missingNodesIndicator, false, nodePositions, false, missingNodePos, false );
  missingNodePos.fillComplete( threeMap, missingNodesMap );

  MeshGraph ghosts;
  Teuchos::Array< double > extractedNodePos( 3 );
  for( std::size_t i = 0; i < numLocMissingIndices; ++i )
  {
    TriGlbIdx const index = notPresentIndices[i];
    Geom const geometricalType = convert.getGeometricalType( index );

    std::size_t const length = missingIndices.getNumEntriesInGlobalRow( missingIndicesOffset + TriGlbIdx( i ) );
    extractedValues.resize( length );
    extractedIndices.resize( length );
    missingIndices.getGlobalRowCopy( missingIndicesOffset + TriGlbIdx( i ), extractedIndices(), extractedValues(), extracted );
    GEOS_ASSERT_EQ( extracted, length );
    if( geometricalType == Geom::NODE )
    {
      // The case of nodes is a bit different from the other cases,
      // because nodes do not rely on other geometrical quantities,
      // but we need to extract the position of the node instead.
      // In order to extract these positions, we use the other matrix `missingNodePos`.
      GEOS_ASSERT_EQ( length, 1 );
      std::size_t const lengthPos = missingNodePos.getNumEntriesInGlobalRow( missingNodesOffset + TriGlbIdx( i ) );
      GEOS_ASSERT_EQ( lengthPos, 3 );
      extractedIndices.resize( lengthPos );
      missingNodePos.getGlobalRowCopy( missingNodesOffset + TriGlbIdx( i ), extractedIndices(), extractedNodePos(), extracted );
      GEOS_ASSERT_EQ( extracted, lengthPos );
      std::array< double, 3 > & pos = ghosts.n2pos[convert.toNodeGlbIdx( index )];
      for( auto dim = 0; dim < 3; ++dim )
      {
        pos[extractedIndices[dim]] = extractedNodePos[dim];
      }
      continue;
    }

    auto const cit = std::find( std::cbegin( extractedIndices ), std::cend( extractedIndices ), index );
    std::size_t const numGeomQuantitiesIdx = intConv< std::size_t >( std::distance( std::cbegin( extractedIndices ), cit ) );
    TriScalarInt const & numGeomQuantities = extractedValues[numGeomQuantitiesIdx];
    GEOS_ASSERT_EQ( extracted, intConv< std::size_t >( numGeomQuantities + 1 ) );

    switch( geometricalType )
    {
      case Geom::EDGE:
      {
        TriScalarInt const & numNodes = numGeomQuantities;  // Alias
        GEOS_ASSERT_EQ( numNodes, 2 );
        std::array< NodeGlbIdx, 2 > order{};

        for( std::size_t ii = 0; ii < extracted; ++ii )
        {
          if( ii == numGeomQuantitiesIdx )
          {
            continue;
          }

          NodeGlbIdx const ngi = convert.toNodeGlbIdx( extractedIndices[ii] );
          integer const ord = integer( extractedValues[ii] - 1 );
          GEOS_ASSERT( ord == 0 or ord == 1 );
          order[ord] = ngi;
        }
        GEOS_ASSERT_EQ( std::size( order ), intConv< std::size_t >( numNodes ) );
        EdgeGlbIdx const egi = convert.toEdgeGlbIdx( index );
        std::tuple< NodeGlbIdx, NodeGlbIdx > & tmp = ghosts.e2n[egi];
        std::get< 0 >( tmp ) = order[0];
        std::get< 1 >( tmp ) = order[1];
        break;
      }
      case Geom::FACE:
      {
        TriScalarInt const & numEdges = numGeomQuantities;  // Alias
        std::map< integer, EdgeInfo > order;
        for( std::size_t ii = 0; ii < extracted; ++ii )
        {
          if( ii == numGeomQuantitiesIdx )
          {
            continue;
          }

          EdgeGlbIdx const egi = convert.toEdgeGlbIdx( extractedIndices[ii] );
          std::array< std::size_t, 2 > const decoded = decode< 2 >( numEdges, std::size_t( extractedValues[ii] - 1 ) );
          order[decoded[1]] = { egi, intConv< std::uint8_t >( decoded[0] ) };
          GEOS_ASSERT( decoded[0] == 0 or decoded[0] == 1 );
        }
        GEOS_ASSERT_EQ( std::size( order ), intConv< std::size_t >( numEdges ) );
        FaceGlbIdx const fgi = convert.toFaceGlbIdx( index );
        std::vector< EdgeInfo > & tmp = ghosts.f2e[fgi];
        tmp.resize( numEdges );
        for( auto const & [ord, edgeInfo]: order )
        {
          tmp[ord] = edgeInfo;
        }
//        GEOS_LOG_RANK( "faces ; index, extIndices, extValues, order = " << index << " | " << json( extIndices ) << " | " << json( extValues ) << " | " << json( order ) );
        break;
      }
      case Geom::CELL:
      {
        TriScalarInt const & numFaces = numGeomQuantities;  // Alias // TODO This should receive the cell type instead.
        std::map< integer, FaceInfo > order;
        for( std::size_t ii = 0; ii < extracted; ++ii )
        {
          if( ii == numGeomQuantitiesIdx )
          {
            continue;
          }

          FaceGlbIdx const fgi = convert.toFaceGlbIdx( extractedIndices[ii] );
          std::array< std::size_t, 3 > const decoded = decode< 3 >( numFaces, std::size_t( extractedValues[ii] - 1 ) );
          order[decoded[2]] = { fgi, intConv< bool >( decoded[0] ), intConv< std::uint8_t >( decoded[1] ) };
        }
        GEOS_ASSERT_EQ( std::size( order ), intConv< std::size_t >( numFaces ) );
        CellGlbIdx const cgi = convert.toCellGlbIdx( index );
        std::vector< FaceInfo > & tmp = ghosts.c2f[cgi];
        tmp.resize( numFaces );
        for( auto const & [ord, faceInfo]: order )
        {
          tmp[ord] = faceInfo;
        }
        break;
      }
      default:
      {
        GEOS_ERROR( "Internal error." );
      }
    }
  }

//  if( curRank == 2_mpi )
//  {
//    GEOS_LOG_RANK( "neighboringConnexions = " << json( ghosts ) );
//    std::map< MpiRank, std::set< CellGlbIdx > > tmp;
//    for( auto const & [geom, rk]: ownerships.cells )
//    {
//      tmp[rk].insert( geom );
//    }
//    for( auto const & [rk, geom]: tmp )
//    {
//      GEOS_LOG_RANK( "ghost geom = " << rk << ": " << json( geom ) );
//    }
//  }
//  GEOS_LOG_RANK( "my final neighbors are " << json( ownerships.neighbors ) );

//  if( curRank == 1_mpi )
//  {
//    GEOS_LOG_RANK( "ghosts_n2ps = " << json( ghosts.n2pos ) );
//  }

  return { std::move( ghosts ), std::move( recv ), std::move( send ) };
}

void doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                       std::set< MpiRank > const & neighbors,
                       MeshMappingImpl & meshMappings )
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

  auto const [ghosts, recv, send] = performGhosting( owned, present, matrixOffsets, curRank );

  buildPods( owned, present, ghosts, recv, send, meshMappings );
}

void doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                       std::set< int > const & neighbors,
                       MeshMappingImpl & meshMappings )
{
  std::set< MpiRank > neighbors_;
  for( int const & rank: neighbors )
  {
    neighbors_.insert( MpiRank{ rank } );
  }
  GEOS_LOG_RANK( "my initial neighbors are " << json( neighbors_ ) );

  return doTheNewGhosting( mesh, neighbors_, meshMappings );
}

}  // end of namespace geos::ghosting
