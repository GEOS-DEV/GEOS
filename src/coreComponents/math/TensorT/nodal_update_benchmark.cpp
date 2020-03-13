/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "TensorT.h"
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>

#include "common/ArrayT/ArrayT.h"

/// returns the amount of cpu time use for this process
realT getcputime( void );

void function1( realT * __restrict__ const xdisp,
                realT * __restrict__ const ydisp,
                realT * __restrict__ const zdisp,
                realT * __restrict__ const xvel,
                realT * __restrict__ const yvel,
                realT * __restrict__ const zvel,
                const realT * __restrict__ const xacc,
                const realT * __restrict__ const yacc,
                const realT * __restrict__ const zacc,
                realT * __restrict__ const xtot,
                realT * __restrict__ const ytot,
                realT * __restrict__ const ztot,
                realT * __restrict__ const dt,
                const int num_nodes,
                const int num_steps );

void function2( R1TensorT< 3 > * __restrict__ const Disp,
                R1TensorT< 3 > * __restrict__ const Vel,
                R1TensorT< 3 > * __restrict__ const Acc,
                R1TensorT< 3 > * __restrict__ Tot,
                realT * __restrict__ const dt,
                const int num_nodes,
                const int num_steps );

void function3( array1d< R1TensorT< 3 > > & Disp,
                array1d< R1TensorT< 3 > > & Vel,
                array1d< R1TensorT< 3 > > & Acc,
                R1TensorT< 3 > * __restrict__ Tot,
                realT * __restrict__ const dt,
                const int num_nodes,
                const int num_steps );


int main( int argc, char * argv[] )
{


  realT iRANDMAX = 1.0 / RAND_MAX;

  std::srand( 1234 );

  const int num_nodes = atoi( argv[1] );
  const int num_steps = atoi( argv[2] );

  realT dt[num_steps];


  R1TensorT< 3 > * const Disp = new R1TensorT< 3 >[num_nodes];
  R1TensorT< 3 > * const Vel = new R1TensorT< 3 >[num_nodes];
  R1TensorT< 3 > * const Acc = new R1TensorT< 3 >[num_nodes];
  R1TensorT< 3 > Tot;

  array1d< R1TensorT< 3 > > Disp2( num_nodes );
  array1d< R1TensorT< 3 > > Vel2( num_nodes );
  array1d< R1TensorT< 3 > > Acc2( num_nodes );
//  std::cout<<"sizeof( R1TensorT<3> ) = "<<sizeof( R1TensorT<3> )<<std::endl;
/*
   void* junk;
   posix_memalign( &junk, 16, num_nodes*sizeof(R1TensorT<3>));
   R1TensorT<3>* Disp = static_cast<R1TensorT<3>*>(junk);

   posix_memalign( &junk, 16, num_nodes*sizeof(R1TensorT<3>));
   R1TensorT<3>* Vel = static_cast<R1TensorT<3>*>(junk);

   posix_memalign( &junk, 16, num_nodes*sizeof(R1TensorT<3>));
   R1TensorT<3>* Acc = static_cast<R1TensorT<3>*>(junk);


   std::cout<<"sizeof(realT) = "<<sizeof(realT)<<std::endl;
   std::cout<<"__alignof__(realT) = "<<__alignof__(realT)<<std::endl;

   std::cout<<"sizeof(R1TensorT<3>) = "<<sizeof(R1TensorT<3>)<<std::endl;
   std::cout<<"__alignof__(R1TensorT<3>) =
      "<<__alignof__(R1TensorT<3>)<<std::endl;

   std::cout<<"sizeof(R2TensorT<3>) = "<<sizeof(R2TensorT<3>)<<std::endl;
   std::cout<<"__alignof__(R2TensorT<3>) =
      "<<__alignof__(R2TensorT<3>)<<std::endl;

   std::cout<<"sizeof(R2SymTensorT<3>) = "<<sizeof(R2SymTensorT<3>)<<std::endl;
   std::cout<<"__alignof__(R2SymTensorT<3>) =
      "<<__alignof__(R2SymTensorT<3>)<<std::endl;
 */


  realT * const xdisp = new realT[num_nodes];
  realT * const ydisp = new realT[num_nodes];
  realT * const zdisp = new realT[num_nodes];
  realT * const xvel  = new realT[num_nodes];
  realT * const yvel  = new realT[num_nodes];
  realT * const zvel  = new realT[num_nodes];
  realT * const xacc  = new realT[num_nodes];
  realT * const yacc  = new realT[num_nodes];
  realT * const zacc  = new realT[num_nodes];
  realT xtot = 0.0;
  realT ytot = 0.0;
  realT ztot = 0.0;


  for( int a=0; a<num_nodes; ++a )
  {
    xacc[a] = realT( rand()) * iRANDMAX;
    yacc[a] = realT( rand()) * iRANDMAX;
    zacc[a] = realT( rand()) * iRANDMAX;
    Acc[a]( 0 ) = xacc[a];
    Acc[a]( 1 ) = yacc[a];
    Acc[a]( 2 ) = zacc[a];
    Acc2[a]( 0 ) = xacc[a];
    Acc2[a]( 1 ) = yacc[a];
    Acc2[a]( 2 ) = zacc[a];

  }

  for( int i=0; i<num_steps; ++i )
  {
    dt[i] = realT( rand()); //* iRANDMAX;
  }

  int flag = atoi( argv[3] );
//  std::cout<<flag<<std::endl;

  // basic c-arrays
  realT t1 = getcputime();

  if( flag == 0 )
  {
    function1( xdisp, ydisp, zdisp,
               xvel, yvel, zvel,
               xacc, yacc, zacc,
               &xtot, &ytot, &ztot,
               dt, num_nodes, num_steps );
  }
  else if( flag == 1 )
  {
    function2( Disp, Vel, Acc, &Tot,
               dt, num_nodes, num_steps );
  }
  else
  {
    function3( Disp2, Vel2, Acc2, &Tot,
               dt, num_nodes, num_steps );
  }

  realT t2 = getcputime();


  for( int a=0; a<num_nodes; ++a )
  {
    xtot += xdisp[a];
    ytot += ydisp[a];
    ztot += zdisp[a];
  }
  for( int a=0; a<num_nodes; ++a )
  {
    Tot += Disp[a];
  }


  std::cout<<"\t\t\t\t"<<xtot<<' '<<ytot<<' '<<ztot<<std::endl;
  std::cout<<"\t\t\t\t"<<Tot( 0 )<<' '<<Tot( 1 )<<' '<<Tot( 2 )<<std::endl;
//  std::cout<<"baseline CPU time    = "<<t2-t1<<std::endl;
//  std::cout<<"Tensor CPU time      = "<<t3-t2<<std::endl<<std::endl;

  std::cout<<num_nodes<<' '<<t2-t1<<std::endl<<std::endl;

  delete [] xdisp;
  delete [] ydisp;
  delete [] zdisp;
  delete [] xvel;
  delete [] yvel;
  delete [] zvel;
  delete [] xacc;
  delete [] yacc;
  delete [] zacc;

  delete [] Acc;
  delete [] Vel;
  delete [] Disp;

  return 0;
}

inline
void function1( realT * __restrict__ const xdisp,
                realT * __restrict__ const ydisp,
                realT * __restrict__ const zdisp,
                realT * __restrict__ const xvel,
                realT * __restrict__ const yvel,
                realT * __restrict__ const zvel,
                const realT * __restrict__ const xacc,
                const realT * __restrict__ const yacc,
                const realT * __restrict__ const zacc,
                realT * __restrict__ const xtot,
                realT * __restrict__ const ytot,
                realT * __restrict__ const ztot,
                realT * __restrict__ const dt,
                const int num_nodes,
                const int num_steps )
{

  for( int i=0; i<num_steps; ++i )
  {
    const realT deltatime = dt[i];
    for( int a=0; a<num_nodes; ++a )
    {
      xvel[a] += xacc[a] * deltatime;
      yvel[a] += yacc[a] * deltatime;
      zvel[a] += zacc[a] * deltatime;
    }

    for( int a=0; a<num_nodes; ++a )
    {
      xdisp[a] += xvel[a] * deltatime;
      ydisp[a] += yvel[a] * deltatime;
      zdisp[a] += zvel[a] * deltatime;
    }
  }



}


inline
void function2( R1TensorT< 3 > * __restrict__ const Disp,
                R1TensorT< 3 > * __restrict__ const Vel,
                R1TensorT< 3 > * __restrict__ const Acc,
                R1TensorT< 3 > * __restrict__ Tot,
                realT * __restrict__ const dt,
                const int num_nodes,
                const int num_steps )
{
  // tensor class

  R1TensorT< 3 > Temp;
  for( int i=0; i<num_steps; ++i )
  {
    const realT deltatime = dt[i];

    for( int a=0; a<num_nodes; ++a )
    {
      Temp = Acc[a];
      Temp *= deltatime;
      Vel[a] += Temp;

//      Vel[a].plus_cA( deltatime, Acc[a] );

    }


    for( int a=0; a<num_nodes; ++a )
    {
      Temp = Vel[a];
      Temp *= deltatime;
      Disp[a] += Temp;

//      Disp[a].plus_cA( deltatime, Vel[a] );

    }

  }

}


inline
void function3( array1d< R1TensorT< 3 > > & Disp,
                array1d< R1TensorT< 3 > > & Vel,
                array1d< R1TensorT< 3 > > & Acc,
                R1TensorT< 3 > * __restrict__ Tot,
                realT * __restrict__ const dt,
                const int num_nodes,
                const int num_steps )
{
  // tensor class

  R1TensorT< 3 > Temp;
  for( int i=0; i<num_steps; ++i )
  {
    const realT deltatime = dt[i];

    for( int a=0; a<num_nodes; ++a )
    {
      Temp = Acc[a];
      Temp *= deltatime;
      Vel[a] += Temp;

//      Vel[a].plus_cA( deltatime, Acc[a] );

    }


    for( int a=0; a<num_nodes; ++a )
    {
      Temp = Vel[a];
      Temp *= deltatime;
      Disp[a] += Temp;

//      Disp[a].plus_cA( deltatime, Vel[a] );

    }

  }

}

/**
 * @return cpu usage
 *
 * This function uses the rusage structure to query elapsed system time and user
 * time, and returns
 * the result.
 */
realT getcputime( void )
{
  struct timeval tim;
  struct rusage ru;
  getrusage( RUSAGE_SELF, &ru );

  tim=ru.ru_utime;
  realT t=(realT)tim.tv_sec + (realT)tim.tv_usec / 1.0e6;

  tim=ru.ru_stime;
  t+=(realT)tim.tv_sec + (realT)tim.tv_usec / 1.0e6;
  return t;
}
