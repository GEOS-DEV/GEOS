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

/// returns the amount of cpu time use for this process
realT getcputime( void );


int main( int argc, char * argv[] )
{


  realT iRANDMAX = 1.0 / RAND_MAX;

  std::srand( 1234 );

  const int num_nodes = atoi( argv[1] );
  const int num_steps = atoi( argv[2] );


  R1TensorT< 3 > * const Disp = new R1TensorT< 3 >[num_nodes];
  R1TensorT< 3 > Tot;

  realT dt[num_steps];

  realT * const xdisp = new realT[num_nodes];
  realT * const ydisp = new realT[num_nodes];
  realT * const zdisp = new realT[num_nodes];
  realT xtot = 0.0;
  realT ytot = 0.0;
  realT ztot = 0.0;


  for( int a=0; a<num_nodes; ++a )
  {
    xdisp[a] = realT( rand()) * iRANDMAX;
    ydisp[a] = realT( rand()) * iRANDMAX;
    zdisp[a] = realT( rand()) * iRANDMAX;
    Disp[a]( 0 ) = xacc[a];
    Disp[a]( 1 ) = yacc[a];
    Disp[a]( 2 ) = zacc[a];
  }

  for( int i=0; i<num_steps; ++i )
  {
    dt[i] = realT( rand()) * iRANDMAX;
  }

  int flag = atoi( argv[3] );

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
  else
  {
    function2( Disp, Vel, Acc, &Tot,
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


  GEOSX_LOG( "\t\t\t\t"<<xtot<<' '<<ytot<<' '<<ztot );
  GEOSX_LOG( "\t\t\t\t"<<Tot( 0 )<<' '<<Tot( 1 )<<' '<<Tot( 2 ));
//  GEOSX_LOG("baseline CPU time    = "<<t2-t1);
//  GEOSX_LOG("Tensor CPU time      = "<<t3-t2<<std::endl);

  GEOSX_LOG( num_nodes<<' '<<t2-t1<<std::endl );

  delete [] xdisp;
  delete [] ydisp;
  delete [] zdisp;

  delete [] Disp;

  return 0;
}

inline
void CalculateGradient( realT Gradient[3][3],
                        const realT * const x,
                        const realT * const y,
                        const realT * const z,
                        const realT * const * const dNdX )

{

  for( int i=0; i<3; ++i )
    for( int j=0; j<3; ++j )
    {
      Gradient[i][j] = 0.0;
    }

  for( int a=0; a<8; ++a )
  {
    for( int i=0; i<3; ++i )
    {
      Gradient[0][i] += x[a] * dNdX[a][i];
      Gradient[1][i] += y[a] * dNdX[a][i];
      Gradient[2][i] += z[a] * dNdX[a][i];
    }
  }

}


void CalculateGradient( R2TensorT< nsdof > & Gradient,
                        const array1d< R1TensorT< nsdof > > & disp,
                        const array1d< R1TensorT< nsdof > > & dNdX )

{
  Gradient = 0.0;
  for( int a=0; a<8; ++a )
    Gradient.plus_dyadic_ab( disp( a ), dNdX( a ));

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
