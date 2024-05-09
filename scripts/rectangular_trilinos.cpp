/*
 * Compile: mpicxx -std=c++20 -I /opt/GEOS/GEOS_TPL-263-252-ba785a2/trilinos/include/ rectangular_trilinos.cpp -L /opt/GEOS/GEOS_TPL-263-252-ba785a2/trilinos/lib -lepetra -Wl,-rpath /opt/GEOS/GEOS_TPL-263-252-ba785a2/trilinos/lib -o /tmp/rectangular_trilinos
 * Launch: mpirun -n 3 /tmp/rectangular_trilinos
 */

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>

#include <mpi.h>

#include <iostream>

int main( int argc,
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

  MPI_Finalize();

  return 0;
}