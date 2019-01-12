/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


#include "TensorT.h"
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
#include <vector>

/// returns the amount of cpu time use for this process
realT getcputime(void);



void function1( realT* __restrict__ a, realT* __restrict__ b, realT* __restrict__  c, realT* __restrict__  d, realT* __restrict__  e, const int n );

void function2( R2TensorT<3>& A, R2TensorT<3>& B, R2TensorT<3>& C, R2TensorT<3>& D, R2TensorT<3>& E, const int n );

int main(int argc, char* argv[] )
{
  realT iRANDMAX = 1.0 / RAND_MAX;

  std::srand(1234);

//   int num_iter = 1000000;
  int num_iter = atoi(argv[1]);


  realT a[9] = { realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand())  };
  realT b[9] = { realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand())  };
  realT c[9] = { realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand())  };
  realT d[9] = { realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand()), realT(rand())  };

  realT e[9] = {0.0};


  for( int i=0 ; i<9 ; ++i )
  {
    a[i] *= iRANDMAX;
    b[i] *= iRANDMAX;
    c[i] *= iRANDMAX;
    d[i] *= iRANDMAX;
  }

  R2TensorT<3> A,B,C,D,E;

  int count = 0;
  for( int i=0 ; i<3 ; ++i )
  {
    for( int j=0 ; j<3 ; ++j )
    {
      A(i,j) = a[count];
      B(i,j) = b[count];
      C(i,j) = c[count];
      D(i,j) = d[count++];

    }
  }


  realT t1 = getcputime();

  function1( a,b,c,d,e,num_iter );

  realT t2 = getcputime();

  function2( A,B,C,D,E,num_iter);

  realT t3 = getcputime();


  for( int i=0 ; i<9 ; ++i )
    GEOS_LOG(e[i]);
  GEOS_LOG("");


  for( int i=0 ; i<3 ; ++i )
    for( int j=0 ; j<3 ; ++j )
    {
      GEOS_LOG(E(i,j));
    }
  GEOS_LOG("");

  GEOS_LOG("baseline CPU time    = "<<t2-t1);
  GEOS_LOG("TensorClass CPU time = "<<t3-t2);


  return 0;


}


void function1( realT* __restrict__ a, realT* __restrict__ b, realT* __restrict__  c, realT* __restrict__  d, realT* __restrict__  e, const int n )
{

  for( int k = 0 ; k < n ; ++k )
  {

    for( int i=0 ; i<9 ; ++i )
    {
      e[i] += a[i] * b[i] + c[i] * d[i];
      a[i] *= a[i];
      d[i] *= c[i];
    }

  }



}

void function2( R2TensorT<3>& A, R2TensorT<3>& B, R2TensorT<3>& C, R2TensorT<3>& D, R2TensorT<3>& E, const int n )
{
  R2TensorT<3> TEMP, TEMP2;
  for( int k = 0 ; k < n ; ++k )
  {
/*
    TEMP = A;
    TEMP *= B;

    TEMP2 = C;
    TEMP2 *= D;

    E += TEMP;
    E += TEMP2;
 */
    E.AB_plus_CD(A,B,C,D);

    A *= A;
    D *= C;


  }



}


/**
 * @author Randy Settgast
 * @return cpu usage
 *
 * This function uses the rusage structure to query elapsed system time and user
 * time, and returns
 * the result.
 */
realT getcputime(void)
{
  struct timeval tim;
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);

  tim=ru.ru_utime;
  realT t=(realT)tim.tv_sec + (realT)tim.tv_usec / 1.0e6;

  tim=ru.ru_stime;
  t+=(realT)tim.tv_sec + (realT)tim.tv_usec / 1.0e6;
  return t;
}
