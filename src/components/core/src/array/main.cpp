/*
 * main.cpp
 *
 *  Created on: Jul 30, 2016
 *      Author: rrsettgast
 */


#include <sys/time.h>
#include <stdint.h>
#include <string>
#include <math.h>

#include "MultidimensionalArray.hpp"

using namespace multidimensionalArray;
uint64_t GetTimeMs64()
{
  struct timeval tv;

  gettimeofday( &tv, NULL );

  uint64_t ret = tv.tv_usec;
  /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
  ret /= 1000;

  /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
  ret += ( tv.tv_sec * 1000 );

  return ret;

}




int main( int argc, char* argv[] )
{
  unsigned int seed = time( NULL );

  const int num_i = std::stoi( argv[1] );
  const int num_k = std::stoi( argv[2] );
  const int num_j = std::stoi( argv[3] );
  const int ITERATIONS = std::stoi( argv[4] );



  double A[num_i][num_k];
  double B[num_k][num_j];
  double C1[num_i][num_j];
  double C2a[num_i][num_j];
  double C2b[num_i][num_j];

  srand( seed );

  for( int i = 0 ; i < num_i ; ++i )
    for( int k = 0 ; k < num_k ; ++k )
      A[i][k] = rand();

  for( int k = 0 ; k < num_k ; ++k )
    for( int j = 0 ; j < num_j ; ++j )
    {
      B[k][j] = rand();
    }

  for( int i = 0 ; i < num_i ; ++i )
    for( int j = 0 ; j < num_j ; ++j )
    {
      C1[i][j] = 0.0;
      C2a[i][j] = 0.0;
      C2b[i][j] = 0.0;
    }

  double __restrict * const A1d = &(A[0][0]);
  double __restrict * const B1d = &(B[0][0]);
  double __restrict * const C1d = &(C1[0][0]);

  uint64_t startTime = GetTimeMs64();
  for( int iter = 0 ; iter < ITERATIONS ; ++iter )
  {
    for( int i = 0 ; i < num_i ; ++i )
    {
      for( int j = 0 ; j < num_j ; ++j )
      {
        for( int k = 0 ; k < num_k ; ++k )
        {
          C1d[ i*num_j+j ] += A1d[ i*num_k+k ] * B1d[ k*num_j+j ];
        }
      }
    }
  }
  uint64_t endTime = GetTimeMs64();
  double runTime1 = ( endTime - startTime ) / 1000.0;
//  printf( "[Native 1D Array]   Elapsed time: %6.3f seconds\n", ( endTime - startTime ) / 1000.0 );


  startTime = GetTimeMs64();
  for( int iter = 0 ; iter < ITERATIONS ; ++iter )
  {
    for( int i = 0 ; i < num_i ; ++i )
    {
      for( int j = 0 ; j < num_j ; ++j )
      {
        for( int k = 0 ; k < num_k ; ++k )
        {
          C2a[i][j] += A[i][k] * B[k][j];
        }
      }
    }
  }
  endTime = GetTimeMs64();
  double runTime2a = ( endTime - startTime ) / 1000.0;

//  printf( "[Native 2D Array]   Elapsed time: %6.3f seconds\n", ( endTime - startTime ) / 1000.0 );


  int64 lengthsA[] = { num_i , num_k };
  int64 lengthsB[] = { num_k , num_j };
  int64 lengthsC[] = { num_i , num_j };
  ArrayAccessor<double,2> arrayA( &(A[0][0]), lengthsA );
  ArrayAccessor<double,2> arrayB( &(B[0][0]), lengthsB );
  ArrayAccessor<double,2> arrayC( &(C2b[0][0]), lengthsC );

  startTime = GetTimeMs64();
  for( int iter = 0 ; iter < ITERATIONS ; ++iter )
  {
    for( int i = 0 ; i < num_i ; ++i )
    {
      for( int j = 0 ; j < num_j ; ++j )
      {
        for( int k = 0 ; k < num_k ; ++k )
        {
          arrayC[i][j] += arrayA[i][k] * arrayB[k][j];
        }
      }
    }
  }
  endTime = GetTimeMs64();
  double runTime2b = ( endTime - startTime ) / 1000.0;

//  printf( "[Randy's 2D Array]   Elapsed time: %6.3f seconds\n", ( endTime - startTime ) / 1000.0 );


  double error12a = 0.0;
  double error12b = 0.0;
  double error2a2b = 0.0;
  double error2a2c = 0.0;

  for( int i = 0 ; i < num_i ; ++i )
  {
    for( int j = 0 ; j < num_j ; ++j )
    {
      error12a  += pow( C1[i][j] - C2a[i][j] , 2 ) ;
      error12b  += pow( C1[i][j] - C2b[i][j] , 2 ) ;
      error2a2b += pow( C2a[i][j] - C2b[i][j] , 2 ) ;
    }
  }
  printf( "1d, 2d_native, 2db: %6.3f %6.3f %6.3f\n", runTime1, runTime2a, runTime2b );
  std::cout<<"error12a = "<<error12a<<std::endl;
  std::cout<<"error12b = "<<error12b<<std::endl;
  std::cout<<"error2a2b = "<<error2a2b<<std::endl;
  std::cout<<"error2a2c = "<<error2a2c<<std::endl;
  return 0;
}
