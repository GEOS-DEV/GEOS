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

/**
 * @brief Complete set of information when a face refers to an edge.
 * @details When a face refers to an edge, w.r.t. the canonical ordering of the edge,
 * there's one node that comes before the other.
 * The @c start parameter (which will be equal to either 0 or 1) will hold this information.
 */
struct EdgeInfo
{
  /// The global index of the edge the face is referring to.
  EdgeGlbIdx index;
  /// Which nodes comes first (and therefore which comes second) in the original edge numbering in the cell.
  std::uint8_t start; // Where the initial index in the canonical face ordering was in the initial face.
};

void to_json( json & j,
              const EdgeInfo & v )  // For display
{
  j = json{ { "index", v.index },
            { "start", v.start } };
}

/**
 * @brief Complete set of information when a cell refers to a face.
 * @details When a cell refers to a face, w.r.t. the canonical ordering of the face,
 * the nodes of the face (in the cell) can be shifted,
 * and travelled in a different direction (clockwise or counter-clockwise).
 * The @c isFlipped parameter informs about the travel discrepancy.
 * The @c start parameter informs about the shift.
 * @note This class does not refer to how multiple faces are ordered by a cell,
 * but to the precise information when refering to one given face.
 */
struct FaceInfo
{
  /// The global index of the face the cell is referring to.
  FaceGlbIdx index;
  /// Is the face travelled in the same direction as the canonical face.
  bool isFlipped;
  /// Where to start iterating over the nodes of the canonical faces to get back to the original face numbering in the cell.
  std::uint8_t start;
};

void to_json( json & j,
              const FaceInfo & v )  // For display
{
  j = json{ { "index",     v.index },
            { "isFlipped", v.isFlipped },
            { "start",     v.start } };
}

struct MeshGraph
{
  std::map< CellGlbIdx, std::vector< FaceInfo > > c2f;  // TODO What about the metadata (e.g. flip the face)
  std::map< FaceGlbIdx, std::vector< EdgeInfo > > f2e;
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
        EdgeInfo const info = { n2e.at( p1 ), std::uint8_t{ p0 == p1 } };
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

struct GhostSend
{
  std::map< NodeGlbIdx, std::set< MpiRank > > nodes;
  std::map< EdgeGlbIdx, std::set< MpiRank > > edges;
  std::map< FaceGlbIdx, std::set< MpiRank > > faces;
  std::map< CellGlbIdx, std::set< MpiRank > > cells;
};

void to_json( json & j,
              const GhostSend & v )
{
  j = json{ { "nodes", v.nodes },
            { "edges", v.edges },
            { "faces", v.faces },
            { "cells", v.cells } };
}

struct GhostRecv
{
  std::map< NodeGlbIdx, MpiRank > nodes;
  std::map< EdgeGlbIdx, MpiRank > edges;
  std::map< FaceGlbIdx, MpiRank > faces;
  std::map< CellGlbIdx, MpiRank > cells;
};

void to_json( json & j,
              const GhostRecv & v )
{
  j = json{ { "nodes", v.nodes },
            { "edges", v.edges },
            { "faces", v.faces },
            { "cells", v.cells } };
}


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

  FindGeometricalType( NodeGlbIdx const & maxNodeGlbIdx,
                       EdgeGlbIdx const & maxEdgeGlbIdx,
                       FaceGlbIdx const & maxFaceGlbIdx,
                       CellGlbIdx const & maxCellGlbIdx )
    : m_edgeOffset( intConv< int >( maxNodeGlbIdx.get() + 1 ) ),
      m_faceOffset( intConv< int >( m_edgeOffset + maxEdgeGlbIdx.get() + 1 ) ),
      m_cellOffset( intConv< int >( m_faceOffset + maxFaceGlbIdx.get() + 1 ) ),
      m_numEntries( intConv< int >( m_cellOffset + maxCellGlbIdx.get() + 1 ) )
  { }

  [[nodiscard]] int numEntries() const
  {
    return m_numEntries;
  }

  [[nodiscard]] Geom getGeometricalType( int const & index ) const
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

  [[nodiscard]] NodeGlbIdx toNodeGlbIdx( int const & index ) const
  {
    return NodeGlbIdx{ intConv< NodeGlbIdx::UnderlyingType >( index ) };
  }

  [[nodiscard]] EdgeGlbIdx toEdgeGlbIdx( int const & index ) const
  {
    return EdgeGlbIdx{ intConv< EdgeGlbIdx::UnderlyingType >( index - m_edgeOffset ) };
  }

  [[nodiscard]] FaceGlbIdx toFaceGlbIdx( int const & index ) const
  {
    return FaceGlbIdx{ intConv< FaceGlbIdx::UnderlyingType >( index - m_faceOffset ) };
  }

  [[nodiscard]] CellGlbIdx toCellGlbIdx( int const & index ) const
  {
    return CellGlbIdx{ intConv< CellGlbIdx::UnderlyingType >( index - m_cellOffset ) };
  }

  [[nodiscard]] int fromNodeGlbIdx( NodeGlbIdx const & ngi ) const
  {
    return intConv< int >( ngi.get() );
  }

  [[nodiscard]] int fromEdgeGlbIdx( EdgeGlbIdx const & egi ) const
  {
    return intConv< int >( egi.get() + m_edgeOffset );
  }

  [[nodiscard]] int fromFaceGlbIdx( FaceGlbIdx const & fgi ) const
  {
    return intConv< int >( fgi.get() + m_faceOffset );
  }

  [[nodiscard]] int fromCellGlbIdx( CellGlbIdx const & cgi ) const
  {
    return intConv< int >( cgi.get() + m_cellOffset );
  }

private:
  int const m_edgeOffset;
  int const m_faceOffset;
  int const m_cellOffset;
  int const m_numEntries;
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


std::tuple< MeshGraph, GhostRecv, GhostSend > assembleAdjacencyMatrix( MeshGraph const & owned,
                                                                       MeshGraph const & present,
                                                                       MaxGlbIdcs const & gis,
                                                                       MpiRank curRank )
{
  FindGeometricalType const convert( gis.nodes, gis.edges, gis.faces, gis.cells );
  using Geom = FindGeometricalType::Geom;

  std::size_t const n = convert.numEntries();  // Total number of entries in the graph.

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
  std::vector< std::vector< double > > values;
  ownedGlbIdcs.reserve( numOwned );
  numEntriesPerRow.reserve( numOwned );
  indices.reserve( numOwned );
  indices.reserve( numOwned );

  std::vector< int > otherGlbIdcs;  // TODO I couldn't use a vector of `std::size_t`
  otherGlbIdcs.reserve( numOther );
  for( NodeGlbIdx const & ngi: present.n )
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

  for( NodeGlbIdx const & ngi: owned.n )
  {
    // Nodes depend on no other geometrical entity,
    // so we only have one entry `1` in the diagonal of the matrix,
    // because we add the identity to the adjacency matrix.
    ownedGlbIdcs.emplace_back( convert.fromNodeGlbIdx( ngi ) );
    numEntriesPerRow.emplace_back( 0 + 1 );  // `+1` comes from the diagonal
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

    numEntriesPerRow.emplace_back( numNodes + 1 );  // `+1` comes from the diagonal
    indices.emplace_back( std::vector< int >{ convert.fromNodeGlbIdx( std::get< 0 >( nodes ) ),
                                              convert.fromNodeGlbIdx( std::get< 1 >( nodes ) ),
                                              ownedGlbIdcs.back() } );
    // Note that when storing a value that can be `0`, we always add `1`,
    // (for edges but also later for faces and cells),
    // to be sure that there will always be some a noticeable figure where we need one.
    values.emplace_back( std::vector< double >{ 0 + 1., 1 + 1., numNodes } );
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

    numEntriesPerRow.emplace_back( numEdges + 1 );  // `+1` comes from the diagonal
    std::vector< int > & ind = indices.emplace_back( numEntriesPerRow.back() );
    std::vector< double > & val = values.emplace_back( numEntriesPerRow.back() );
    for( std::size_t i = 0; i < numEdges; ++i )
    {
      EdgeInfo const & edgeInfo = edges[i];
      ind[i] = convert.fromEdgeGlbIdx( edgeInfo.index );
      std::size_t const v = 1 + encode< 2 >( numEdges, { edgeInfo.start, i } );
      val[i] = double( v );
    }
    ind.back() = ownedGlbIdcs.back();
    val.back() = double( numEdges );
  }
  for( auto const & [cgi, faces]: owned.c2f )
  {
    // The same comment as for faces and edges applies for cells and faces (see above).
    // The main differences being that
    // - the `FaceInfo` as a little more information,
    // - the diagonal stores the type of the cell which conveys the number of face, nodes, ordering...
    std::size_t const numFaces = std::size( faces );

    ownedGlbIdcs.emplace_back( convert.fromCellGlbIdx( cgi ) );

    numEntriesPerRow.emplace_back( numFaces + 1 );  // `+1` comes from the diagonal
    std::vector< int > & ind = indices.emplace_back( numEntriesPerRow.back() );
    std::vector< double > & val = values.emplace_back( numEntriesPerRow.back() );
    for( std::size_t i = 0; i < numFaces; ++i )
    {
      FaceInfo const & faceInfo = faces[i];
      ind[i] = convert.fromFaceGlbIdx( faceInfo.index );
      std::size_t const v = 1 + encode< 3 >( numFaces, { faceInfo.isFlipped, faceInfo.start, i } );
      val[i] = double( v );
    }
    ind.back() = ownedGlbIdcs.back();
    val.back() = double( numFaces );  // TODO This should be Hex and the not the number of faces...
  }
  std::sort( std::begin( ownedGlbIdcs ), std::end( ownedGlbIdcs ) );

  GEOS_ASSERT_EQ( numOwned, std::size( ownedGlbIdcs ) );
  GEOS_ASSERT_EQ( numOwned, std::size( numEntriesPerRow ) );
  GEOS_ASSERT_EQ( numOwned, std::size( indices ) );
  GEOS_ASSERT_EQ( numOwned, std::size( values ) );
  for( std::size_t i = 0; i < numOwned; ++i )
  {
    GEOS_ASSERT_EQ( indices[i].size(), std::size_t( numEntriesPerRow[i] ) );
  }

  Epetra_MpiComm const & comm = Epetra_MpiComm( MPI_COMM_GEOSX );
  Epetra_Map const ownedMap( n, numOwned, ownedGlbIdcs.data(), 0, comm );

  // The `upward` matrix offers a representation of the graph connections.
  // The matrix is square. Each size being the number of geometrical entities.
  Epetra_CrsMatrix upward( Epetra_DataAccess::Copy, ownedMap, numEntriesPerRow.data(), true );

  for( std::size_t i = 0; i < numOwned; ++i )
  {
    std::vector< int > const & rowIndices = indices[i];
    std::vector< double > const & rowValues = values[i];
    GEOS_ASSERT_EQ( std::size( rowIndices ), std::size_t( numEntriesPerRow[i] ) );
    GEOS_ASSERT_EQ( std::size( rowValues ), std::size_t( numEntriesPerRow[i] ) );
    upward.InsertGlobalValues( ownedGlbIdcs[i], std::size( rowIndices ), rowValues.data(), rowIndices.data() );
  }

  upward.FillComplete( ownedMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/adj.mat", upward );

  int const commSize( MpiWrapper::commSize() );

  // Now let's build the domain indicator matrix.
  // It's rectangular, number of rows being the number of MPI ranks,
  // number of columns being the number of nodes in the mesh graph,
  // ie the number of geometrical entities in the mesh.
  // It contains one for every time the geometrical entity is present on the rank.
  // By present, we do not meant that the ranks owns it: just that it's already available on the rank.
  Epetra_Map const mpiMap( commSize, 0, comm );  // Let the current rank get the appropriate index in this map.
  Epetra_CrsMatrix indicator( Epetra_DataAccess::Copy, mpiMap, numOwned + numOther, true );

  std::vector< double > const ones( n, 1. );
  indicator.InsertGlobalValues( curRank.get(), numOwned, ones.data(), ownedGlbIdcs.data() );
  indicator.InsertGlobalValues( curRank.get(), numOther, ones.data(), otherGlbIdcs.data() );
  indicator.FillComplete( ownedMap, mpiMap );

  // The `ownership` matrix is a diagonal square matrix.
  // Each size being the number of geometrical entities.
  // The value of the diagonal term will be the owning rank.
  // By means of matrices multiplications, it will be possible to exchange the ownerships information across the ranks.
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

  // The `ghostingFootprint` matrix is rectangular,
  // the number of columns being the number of MPI ranks,
  // the number of rows being the number of nodes in the mesh graph.
  //
  // For any MPI rank (ie column), a non-zero entry means
  // that the corresponding geometrical entry has to be eventually present onto the rank
  // for the ghosting to be effective.
  //
  // From `ghostingFootprint` we can extract where the current rank has to send any owned graph node.
  Epetra_CrsMatrix ghostingFootprint( multiply( commSize, indicator, upward ) );
  ghostingFootprint.PutScalar( 1. );  // This can be done after the `FillComplete`.
  ghostingFootprint.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/ghostingFootprint.mat", ghostingFootprint );
  // We put `1` everywhere there's a non-zero entry, so we'll be able to compose with the `ownership` matrix.

  // FIXME TODO WARNING From `ghostingFootprint` extract where I have to send
  // Then, perform the ghostingFootprint * ownership matrix multiplication,
  // so I get to know from whom I'll receive.

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "ghosted->NumGlobalCols() = " << ghostingFootprint.NumGlobalCols() );
    GEOS_LOG_RANK( "ghosted->NumGlobalRows() = " << ghostingFootprint.NumGlobalRows() );
  }

  // The `ghostExchange` matrix is rectangular,
  // the number of columns being the number of nodes in the mesh graph,
  // the number of rows being the number of MPI ranks.
  //
  // As the result of the multiplication of the `ghostingFootprint` matrix and the `ownership`,
  // for each row owned (ie at the current MPI rank index),
  // the value of the matrix term will provide the actual owning rank for all the .
  //
  // From `ghostExchange` we can extract which other rank will send to the current rank any graph node.
  Epetra_CrsMatrix ghostExchange( Epetra_DataAccess::Copy, mpiMap, 1, false );
  EpetraExt::MatrixMatrix::Multiply( ghostingFootprint, true, ownership, false, ghostExchange, false );
  ghostExchange.FillComplete( ownedMap, mpiMap );
  // TODO Do I have to work with `ghostingFootprint` if I already have `ghostExchange` which may convey more information?
  // TODO Maybe because of the ownership of the ranks: one is also the "scaled" transposed of the other.

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "ghostExchange->NumGlobalCols() = " << ghostExchange.NumGlobalCols() );
    GEOS_LOG_RANK( "ghostExchange->NumGlobalRows() = " << ghostExchange.NumGlobalRows() );
  }
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/ghostInfo.mat", ghostingFootprint );

  int extracted2 = 0;
  std::vector< double > extractedValues2( commSize );  // TODO improve with resize...
  std::vector< int > extractedIndices2( commSize );

  GhostSend send;
  for( int const & index: ownedGlbIdcs )
  {
    ghostingFootprint.ExtractGlobalRowCopy( index, commSize, extracted2, extractedValues2.data(), extractedIndices2.data() );
    std::set< MpiRank > * sendingTo = nullptr;
    switch( convert.getGeometricalType( index ) )
    {
      case Geom::NODE:
      {
        sendingTo = &send.nodes[convert.toNodeGlbIdx( index )];
        break;
      }
      case Geom::EDGE:
      {
        sendingTo = &send.edges[convert.toEdgeGlbIdx( index )];
        break;
      }
      case Geom::FACE:
      {
        sendingTo = &send.faces[convert.toFaceGlbIdx( index )];
        break;
      }
      case Geom::CELL:
      {
        sendingTo = &send.cells[convert.toCellGlbIdx( index )];
        break;
      }
    }

    for( int ii = 0; ii < extracted2; ++ii )
    {
      MpiRank const rank{ extractedIndices2[ii] };
      if( rank != curRank )
      {
        sendingTo->insert( rank );
      }
    }
  }

  int extracted = 0;
  int const length = ghostExchange.NumGlobalEntries( curRank.get() );
  std::vector< double > extractedValues( length );
  std::vector< int > extractedIndices( length );
  ghostExchange.ExtractGlobalRowCopy( curRank.get(), length, extracted, extractedValues.data(), extractedIndices.data() );
  GEOS_ASSERT_EQ( extracted, length );

  std::set< int > const allNeededIndices( std::cbegin( extractedIndices ), std::cend( extractedIndices ) );
  std::set< int > receivedIndices;  // The graph nodes that my neighbors will send me.
  std::set_difference( std::cbegin( allNeededIndices ), std::cend( allNeededIndices ),
                       std::cbegin( ownedGlbIdcs ), std::cend( ownedGlbIdcs ),
                       std::inserter( receivedIndices, std::end( receivedIndices ) ) );
  std::vector< int > notPresentIndices;  // The graphs nodes that are nor owned neither present by/on the current rank.
  std::set_difference( std::cbegin( receivedIndices ), std::cend( receivedIndices ),
                       std::cbegin( otherGlbIdcs ), std::cend( otherGlbIdcs ),
                       std::back_inserter( notPresentIndices ) );

  GhostRecv recv;
  for( int i = 0; i < extracted; ++i )
  {
    int const & index = extractedIndices[i];
    if( receivedIndices.find( index ) == std::cend( receivedIndices ) )  // TODO make a map `receivedIndices -> mpi rank`
    {
      continue;
    }

    MpiRank const sender{ int( extractedValues[i] ) };
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
    }
  }

//  if( curRank == 1_mpi or curRank == 2_mpi )
//  {
//    GEOS_LOG_RANK( "recv.edges = " << json( recv.edges ) );
//    GEOS_LOG_RANK( "send.edges = " << json( send.edges ) );
//  }

  // At this point, each rank knows what it has to send to and what it is going to receive from the other ranks.
  // The remaining action is about retrieving the additional graph information
  // for the new geometrical quantities that will be sent by the neighbors.
  //
  // In order to do that, we build the `missing` matrix, which is rectangular,
  // - the number of columns being the number of nodes in the mesh graph,
  // - the number of rows being the total number graph nodes that's missing on ranks.
  // - all its terms will be `1`.
  // Combined with the adjacency matrix which conveya a lot of connections information,
  // we'll be able to create the final `missingMappings`.
  // We also create `maps` with specific MPI ranks ownerships and offsets,
  // so we'll be able to retrieve the information on the ranks that need it (nothing more, nothing less).
  std::size_t const numMissingIndices = std::size( notPresentIndices );
  std::size_t const numGlobalMissingIndices = MpiWrapper::sum( numMissingIndices );  // A missing node can be missing on multiple ranks and that's OK.

//  Epetra_Map const missingIndicesMap( intConv< long long int >( numGlobalMissingIndices ), intConv< int >( numMissingIndices ), 0, comm );
  Epetra_Map const missingIndicesMap( intConv< int >( numGlobalMissingIndices ), intConv< int >( numMissingIndices ), 0, comm );

//  std::size_t const offset = MpiWrapper::prefixSum< std::size_t >( numMissingIndices );
  int const offset = *missingIndicesMap.MyGlobalElements();  // Needed to compute the global row
//  GEOS_LOG_RANK( "numMissingIndices, numGlobalMissingIndices, offset = " << numMissingIndices << ", " << numGlobalMissingIndices << ", " << offset );

  Epetra_CrsMatrix missing( Epetra_DataAccess::Copy, missingIndicesMap, 1, true );
  for( std::size_t i = 0; i < numMissingIndices; ++i )
  {
    missing.InsertGlobalValues( offset + i, 1, ones.data(), &notPresentIndices[i] );
  }
  missing.FillComplete( ownedMap, missingIndicesMap );

  auto tDownward = makeTranspose( upward );  // TODO give it to multiply!
  Epetra_CrsMatrix missingMappings( Epetra_DataAccess::Copy, missingIndicesMap, 1, false );
  EpetraExt::MatrixMatrix::Multiply( missing, false, *tDownward, true, missingMappings, false );
  missingMappings.FillComplete( ownedMap, missingIndicesMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/missingMappings.mat", missingMappings );

  MeshGraph ghosts;  // TODO most likely this will be wrong.

  int ext = 0;
  std::vector< double > extValues;
  std::vector< int > extIndices;

  for( int i = 0; i < int( numMissingIndices ); ++i )
  {
    int const length2 = missingMappings.NumGlobalEntries( offset + i );
    extValues.resize( length2 );
    extIndices.resize( length2 );
    missingMappings.ExtractGlobalRowCopy( offset + i, length2, ext, extValues.data(), extIndices.data() );
    GEOS_ASSERT_EQ( ext, length2 );
    int const index = notPresentIndices[i];
    switch( convert.getGeometricalType( index ) )
    {
      case Geom::NODE:
      {
        // Nothing to be done for nodes since they do not rely on any underlying information.
        break;
      }
      case Geom::EDGE:
      {
        auto const cit = std::find( std::cbegin( extIndices ), std::cend( extIndices ), index );
        std::ptrdiff_t const numNodesIdx = std::distance( std::cbegin( extIndices ), cit );
        int const numNodes = int( extValues[numNodesIdx] );
        GEOS_ASSERT_EQ( ext, numNodes + 1 );
        GEOS_ASSERT_EQ( numNodes, 2 );
        std::array< NodeGlbIdx, 2 > order{};

        for( int ii = 0; ii < ext; ++ii )
        {
          if( ii == numNodesIdx )
          {
            continue;
          }

          NodeGlbIdx const ngi = convert.toNodeGlbIdx( extIndices[ii] );
          integer const ord = integer( extValues[ii] - 1 );
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
        auto const cit = std::find( std::cbegin( extIndices ), std::cend( extIndices ), index );
        std::ptrdiff_t const numEdgesIdx = std::distance( std::cbegin( extIndices ), cit );
        int const numEdges = int( extValues[numEdgesIdx] );
        GEOS_ASSERT_EQ( ext, numEdges + 1 );
        std::map< integer, EdgeInfo > order;
        for( int ii = 0; ii < ext; ++ii )
        {
          if( ii == numEdgesIdx )
          {
            continue;
          }

          EdgeGlbIdx const egi = convert.toEdgeGlbIdx( extIndices[ii] );
          std::array< std::size_t, 2 > const decoded = decode< 2 >( numEdges, std::size_t( extValues[ii] - 1 ) );
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
        auto const cit = std::find( std::cbegin( extIndices ), std::cend( extIndices ), index );
        std::ptrdiff_t const numFacesIdx = std::distance( std::cbegin( extIndices ), cit );
        int const numFaces = int( extValues[numFacesIdx] );  // TODO This should receive the cell type instead.
        GEOS_ASSERT_EQ( ext, numFaces + 1 );
        std::map< integer, FaceInfo > order;
        for( int ii = 0; ii < ext; ++ii )
        {
          if( ii == numFacesIdx )
          {
            continue;
          }

          FaceGlbIdx const fgi = convert.toFaceGlbIdx( extIndices[ii] );
          std::array< std::size_t, 3 > const decoded = decode< 3 >( numFaces, std::size_t( extValues[ii] - 1 ) );
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
//        GEOS_LOG_RANK( "cells ; index, extIndices, extValues, order = " << index << " | " << json( extIndices ) << " | " << json( extValues ) << " | " << json( order ) );
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

  return { std::move( ghosts ), std::move( recv ), std::move( send ) };
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
                GhostRecv const & recv,
                GhostSend const & send )
{
//  std::size_t const numNodes = std::size( ownerships.nodes );
//  std::size_t const numEdges = std::size( ownerships.edges );
//  std::size_t const numFaces = std::size( ownerships.faces );
//
//  auto [ghostRank, l2g] = buildGhostRankAndL2G( ownerships.edges );
//
//  NodeMgrImpl const nodeMgr( NodeLocIdx{ intConv< localIndex >( numNodes ) } );
//  EdgeMgrImpl const edgeMgr( EdgeLocIdx{ intConv< localIndex >( numEdges ) }, std::move( ghostRank ), std::move( l2g ) );
//  FaceMgrImpl const faceMgr( FaceLocIdx{ intConv< localIndex >( numFaces ) } );
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

  auto const [ghosts, recv, send] = assembleAdjacencyMatrix( owned, present, matrixOffsets, curRank );

  buildPods( owned, present, ghosts, recv, send );

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
