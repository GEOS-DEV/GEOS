/*
 * Compile: mpicxx -std=c++20 mpi_global_offset.cpp -o /tmp/mpi_global_offset
 * Launch: mpirun -n 3 /tmp/mpi_global_offset
 */

#include <mpi.h>

//#include <cstddef>
//#include <array>
#include <iostream>

struct S
{
  std::int32_t n;  // Compute offsets and not values... S.edgeOffset, S.faceOffset, S.cellOffset
  std::int32_t e;
  std::int32_t f;
};

void f_impl( S const * in,
             S * inout )
{
  inout->n = std::max( in->n, inout->n );
  inout->e = std::max( in->e, inout->e );
  inout->f = std::max( in->f, inout->f );
}

void f( void * in,
        void * inout,
        int * len,
        MPI_Datatype * dptr )
{
  f_impl( (S const *) in, (S *) inout );
}

int main( int argc,
          char ** argv )
{
  int rank, size;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  std::cout << "rank " << rank << std::endl;

//  constexpr int count = 3;
//  std::array< int, count > const lengths{ 1, 1, 1 };
//  std::array< MPI_Aint, count > const offsets{ offsetof( S, n ), offsetof( S, e ), offsetof( S, f ) };
//  std::array< MPI_Datatype, count > const types{ MPI_INT32_T, MPI_INT32_T, MPI_INT32_T };
//  MPI_Datatype t;
//  MPI_Type_create_struct( count, lengths.data(), offsets.data(), types.data(), &t );
//  MPI_Type_commit( &t );

  MPI_Datatype t;
  MPI_Type_contiguous( 3, MPI_INT32_T, &t );
  MPI_Type_commit( &t );

  MPI_Op myOp;
  MPI_Op_create( f, true, &myOp );

  S s{ rank, 2 * rank, - rank };
  S answer;

  MPI_Allreduce( &s, &answer, 1, t, myOp, MPI_COMM_WORLD );

//  std::cout << "inp on rank " << rank << " is [" << s.n << ", " << s.e << ", " << s.f << "]" << std::endl;
  std::cout << "ans on rank " << rank << " is [" << answer.n << ", " << answer.e << ", " << answer.f << "]" << std::endl;

  // Just be a little in order
  MPI_Barrier( MPI_COMM_WORLD );

  MPI_Finalize();

  return 0;
}
