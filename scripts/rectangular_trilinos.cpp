/*
 * Compile: mpicxx -std=c++20 -I /tmp/tmp.EtfQgEUXLg/src/cmake-build-main/_deps/namedtype-src/include -I /opt/GEOS/GEOS_TPL-263-252-ba785a2/trilinos/include/ /tmp/tmp.EtfQgEUXLg/scripts/rectangular_trilinos.cpp -L /opt/GEOS/GEOS_TPL-263-252-ba785a2/trilinos/lib -lepetra -lepetraext -ltpetra -lteuchoscomm -lteuchoscore -lkokkoscore -Wl,-rpath /opt/GEOS/GEOS_TPL-263-252-ba785a2/trilinos/lib -o /tmp/rectangular_trilinos
 * Launch: mpirun -n 3 /tmp/rectangular_trilinos
 */

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_Comm.hpp>

#include <MatrixMarket_Tpetra.hpp>
//#include <Kokkos_DefaultNode.hpp>

#include <mpi.h>

#include <NamedType/named_type.hpp>

#include <iostream>
#include <vector>

//using TriMap = Tpetra::Map< int, long long, Kokkos::DefaultNode::DefaultNodeType >;
//using TriCrsMatrix = Tpetra::CrsMatrix< double, int, long long, Kokkos::DefaultNode::DefaultNodeType >;

//using TriLocIdx = fluent::NamedType< int, struct TagTriLocIdx, fluent::Arithmetic, fluent::ImplicitlyConvertibleTo >;
//using TriGlbIdx = fluent::NamedType< long long int, struct TagTriGlbIdx, fluent::Arithmetic, fluent::ImplicitlyConvertibleTo >;
//template<>
//struct std::is_signed< TriLocIdx > : std::is_signed< typename TriLocIdx::UnderlyingType >
//{
//};

using TriLocIdx = int;
using TriGlbIdx = long long int;

using TriMap = Tpetra::Map< TriLocIdx, TriGlbIdx >;
using TriCrsMatrix = Tpetra::CrsMatrix< std::uint32_t, TriLocIdx, TriGlbIdx >;
using TriComm = Teuchos::Comm< int >;

using Teuchos::RCP;

template< typename T, typename... ARGS >
Teuchos::RCP< T > make_rcp( ARGS && ... args )
{
  return Teuchos::rcp( new T( std::forward< ARGS >( args )... ) );
}

int mainTpetra( int argc,
                char ** argv )
{
  int rank, size;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  {  // Scope for MPI race destructors
    Tpetra::global_size_t const numGlbCols{ 4 };
    Tpetra::global_size_t const numGlbRows( size );

    RCP< TriComm const > comm = make_rcp< Teuchos::MpiComm< int > const >( MPI_COMM_WORLD );
    std::vector< TriGlbIdx > const ownedGlbRows{ TriGlbIdx( rank ) };

    RCP< TriMap const > rowMap = make_rcp< TriMap const >( numGlbRows,
                                                           ownedGlbRows.data(),
                                                           TriLocIdx( std::size( ownedGlbRows ) ),
                                                           TriGlbIdx{ 0 },
                                                           comm );
    RCP< TriCrsMatrix > m = make_rcp< TriCrsMatrix >( rowMap, numGlbCols );

    std::uint32_t const r = rank;
    std::vector< std::uint32_t > const value{ 1 + 10 * r, 2 + 10 * r, 3 + 10 * r, 4 + 10 * r };
    std::vector< TriGlbIdx > bulk{ TriGlbIdx{ 3 }, TriGlbIdx{ 2 }, TriGlbIdx{ 1 }, TriGlbIdx{ 0 } };
    m->insertGlobalValues( TriGlbIdx( rank ), TriLocIdx{ 4 }, value.data(), bulk.data() );

    RCP< TriMap const > domainMap = make_rcp< TriMap const >( numGlbCols, TriGlbIdx{ 0 }, comm );
    RCP< TriMap const > rangeMap = make_rcp< TriMap const >( numGlbCols, TriGlbIdx{ 0 }, comm );

    m->fillComplete( domainMap, rangeMap );

    std::cout << m->getDomainMap() << std::endl;

    std::cout << "glb(rows, cols) = " << m->getGlobalNumRows() << " " << m->getGlobalNumCols() << std::endl;
    std::cout << "loc(rows, cols) = " << m->getLocalNumRows() << " " << m->getLocalNumCols() << std::endl;

//    Tpetra::MatrixMarket::Writer< TriCrsMatrix >::writeSparseFile( "/tmp/matrices/poc-tpetra.mat", m );
  }

  MPI_Finalize();

  return 0;
}

int mainEpetra( int argc,
                char ** argv )
{
  int rank, size;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  const int numGlbCols = 4;
  const int numGlbRows = size;

  Epetra_MpiComm const & comm = Epetra_MpiComm( MPI_COMM_WORLD );
  std::vector< int > const ownedGlbRows{ rank };

  Epetra_Map const rowMap( numGlbRows, 1, ownedGlbRows.data(), 0, comm );
  Epetra_CrsMatrix m( Epetra_DataAccess::Copy, rowMap, numGlbCols, true );

  std::vector< double > const value{ 1 + 10. * rank, 2 + 10. * rank, 3 + 10. * rank, 4 + 10. * rank };
  m.InsertGlobalValues( rank, numGlbCols, value.data(), std::vector< int >{ 3, 2, 1, 0 }.data() );

  Epetra_Map const domainMap( numGlbCols, 0, comm );
  Epetra_Map const rangeMap( numGlbRows, 0, comm );

//  m.FillComplete();
  m.FillComplete( domainMap, rangeMap );

//  m.ColMap().Print(std::cout);

  m.DomainMap().Print( std::cout );

  std::cout << "glb(rows, cols) = " << m.NumGlobalRows64() << " " << m.NumGlobalCols64() << std::endl;
  std::cout << "loc(rows, cols) = " << m.NumMyRows() << " " << m.NumMyCols() << std::endl;

//  m.Print( std::cout );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/poc-epetra.mat", m );

  MPI_Finalize();

  return 0;
}

int main( int argc,
          char ** argv )
{
  return mainTpetra( argc, argv );
//  return mainEpetra( argc, argv );
}