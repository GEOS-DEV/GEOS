//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "TensorT.h"
#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>

#include "../../common/ArrayT/ArrayT.h"

/// returns the amount of cpu time use for this process
realT getcputime(void)  ;

void function1( realT* __restrict__ const xdisp,
                realT* __restrict__ const ydisp,
                realT* __restrict__ const zdisp,
                realT* __restrict__ const xvel ,
                realT* __restrict__ const yvel ,
                realT* __restrict__ const zvel ,
                const realT* __restrict__ const xacc ,
                const realT* __restrict__ const yacc ,
                const realT* __restrict__ const zacc ,
                realT* __restrict__ const xtot,
                realT* __restrict__ const ytot,
                realT* __restrict__ const ztot,
                realT* __restrict__ const dt ,
                const int num_nodes,
                const int num_steps );

void function2( R1TensorT<3>* __restrict__ const Disp,
                R1TensorT<3>* __restrict__ const Vel,
                R1TensorT<3>* __restrict__ const Acc,
                R1TensorT<3>* __restrict__ Tot,
                realT* __restrict__ const dt ,
                const int num_nodes,
                const int num_steps );

void function3( Array1dT<R1TensorT<3> >& Disp,
                Array1dT<R1TensorT<3> >& Vel,
                Array1dT<R1TensorT<3> >&  Acc,
                R1TensorT<3>* __restrict__ Tot,
                realT* __restrict__ const dt ,
                const int num_nodes,
                const int num_steps );


int main(int argc, char* argv[] )
{


  realT iRANDMAX = 1.0 / RAND_MAX;

  std::srand(1234);

  const int num_nodes = atoi(argv[1]);
  const int num_steps = atoi(argv[2]);

  realT dt[num_steps];


  R1TensorT<3>* const Disp = new R1TensorT<3>[num_nodes];
  R1TensorT<3>* const Vel = new R1TensorT<3>[num_nodes];
  R1TensorT<3>* const Acc = new R1TensorT<3>[num_nodes];
  R1TensorT<3> Tot;

  Array1dT< R1TensorT<3> > Disp2(num_nodes);
  Array1dT< R1TensorT<3> > Vel2(num_nodes);
  Array1dT< R1TensorT<3> > Acc2(num_nodes);
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
  std::cout<<"__alignof__(R1TensorT<3>) = "<<__alignof__(R1TensorT<3>)<<std::endl;

  std::cout<<"sizeof(R2TensorT<3>) = "<<sizeof(R2TensorT<3>)<<std::endl;
  std::cout<<"__alignof__(R2TensorT<3>) = "<<__alignof__(R2TensorT<3>)<<std::endl;

  std::cout<<"sizeof(R2SymTensorT<3>) = "<<sizeof(R2SymTensorT<3>)<<std::endl;
  std::cout<<"__alignof__(R2SymTensorT<3>) = "<<__alignof__(R2SymTensorT<3>)<<std::endl;
*/


  realT* const xdisp = new realT[num_nodes];
  realT* const ydisp = new realT[num_nodes];
  realT* const zdisp = new realT[num_nodes];
  realT* const xvel  = new realT[num_nodes];
  realT* const yvel  = new realT[num_nodes];
  realT* const zvel  = new realT[num_nodes];
  realT* const xacc  = new realT[num_nodes];
  realT* const yacc  = new realT[num_nodes];
  realT* const zacc  = new realT[num_nodes];
  realT xtot = 0.0 ;
  realT ytot = 0.0 ;
  realT ztot = 0.0 ;


  for( int a=0 ; a<num_nodes ; ++a )
  {
    xacc[a] = realT(rand()) * iRANDMAX;
    yacc[a] = realT(rand()) * iRANDMAX;
    zacc[a] = realT(rand()) * iRANDMAX;
    Acc[a](0) = xacc[a];
    Acc[a](1) = yacc[a];
    Acc[a](2) = zacc[a];
    Acc2[a](0) = xacc[a];
    Acc2[a](1) = yacc[a];
    Acc2[a](2) = zacc[a];

  }

  for( int i=0 ; i<num_steps ; ++i )
  {
    dt[i] = realT(rand()) ;//* iRANDMAX;
  }

  int flag = atoi(argv[3]) ;
//  std::cout<<flag<<std::endl;

  // basic c-arrays
  realT t1 = getcputime();

  if( flag == 0 )
  {
  function1(  xdisp, ydisp, zdisp,
              xvel , yvel , zvel ,
              xacc , yacc , zacc ,
              &xtot,  &ytot,  &ztot ,
              dt, num_nodes, num_steps );
  }
  else if( flag == 1)
  {
  function2(  Disp, Vel, Acc, &Tot,
              dt, num_nodes, num_steps );
  }
  else
  {
  function3(  Disp2, Vel2, Acc2, &Tot,
              dt, num_nodes, num_steps );
  }

  realT t2 = getcputime();


  for( int a=0 ; a<num_nodes ; ++a )
  {
    xtot += xdisp[a];
    ytot += ydisp[a];
    ztot += zdisp[a];
  }
  for( int a=0 ; a<num_nodes ; ++a )
  {
    Tot += Disp[a];
  }


  std::cout<<"\t\t\t\t"<<xtot<<' '<<ytot<<' '<<ztot<<std::endl;
  std::cout<<"\t\t\t\t"<<Tot(0)<<' '<<Tot(1)<<' '<<Tot(2)<<std::endl;
//  std::cout<<"baseline CPU time    = "<<t2-t1<<std::endl;
//  std::cout<<"Tensor CPU time      = "<<t3-t2<<std::endl<<std::endl;

  std::cout<<num_nodes<<' '<<t2-t1<<std::endl<<std::endl;

  delete [] xdisp;
  delete [] ydisp ;
  delete [] zdisp ;
  delete [] xvel  ;
  delete [] yvel  ;
  delete [] zvel  ;
  delete [] xacc  ;
  delete [] yacc  ;
  delete [] zacc  ;

  delete [] Acc;
  delete [] Vel;
  delete [] Disp;

  return 0;
}

inline
void function1( realT* __restrict__ const xdisp,
                realT* __restrict__ const ydisp,
                realT* __restrict__ const zdisp,
                realT* __restrict__ const xvel ,
                realT* __restrict__ const yvel ,
                realT* __restrict__ const zvel ,
                const realT* __restrict__ const xacc ,
                const realT* __restrict__ const yacc ,
                const realT* __restrict__ const zacc ,
                realT* __restrict__ const xtot,
                realT* __restrict__ const ytot,
                realT* __restrict__ const ztot,
                realT* __restrict__ const dt ,
                const int num_nodes,
                const int num_steps )
{

    for( int i=0 ; i<num_steps ; ++i )
    {
      const realT deltatime = dt[i];
      for( int a=0 ; a<num_nodes ; ++a )
      {
        xvel[a] += xacc[a] * deltatime;
        yvel[a] += yacc[a] * deltatime;
        zvel[a] += zacc[a] * deltatime;
      }

      for( int a=0 ; a<num_nodes ; ++a )
      {
        xdisp[a] += xvel[a] * deltatime;
        ydisp[a] += yvel[a] * deltatime;
        zdisp[a] += zvel[a] * deltatime;
      }
    }



}


inline
void function2( R1TensorT<3>* __restrict__ const Disp,
                R1TensorT<3>* __restrict__ const Vel,
                R1TensorT<3>* __restrict__ const Acc,
                R1TensorT<3>* __restrict__ Tot,
                realT* __restrict__ const dt,
                const int num_nodes,
                const int num_steps )
{
  // tensor class

  R1TensorT<3> Temp;
  for( int i=0 ; i<num_steps ; ++i )
  {
    const realT deltatime = dt[i];

    for( int a=0 ; a<num_nodes ; ++a )
    {
      Temp = Acc[a];
      Temp *= deltatime;
      Vel[a] += Temp;

//      Vel[a].plus_cA( deltatime, Acc[a] );

    }


    for( int a=0 ; a<num_nodes ; ++a )
    {
      Temp = Vel[a];
      Temp *= deltatime;
      Disp[a] += Temp;

//      Disp[a].plus_cA( deltatime, Vel[a] );

    }

  }

}


inline
void function3( Array1dT<R1TensorT<3> >& Disp,
                Array1dT<R1TensorT<3> >&  Vel,
                Array1dT<R1TensorT<3> >&  Acc,
                R1TensorT<3>* __restrict__ Tot,
                realT* __restrict__ const dt ,
                const int num_nodes,
                const int num_steps )
{
  // tensor class

  R1TensorT<3> Temp;
  for( int i=0 ; i<num_steps ; ++i )
  {
    const realT deltatime = dt[i];

    for( int a=0 ; a<num_nodes ; ++a )
    {
      Temp = Acc[a];
      Temp *= deltatime;
      Vel[a] += Temp;

//      Vel[a].plus_cA( deltatime, Acc[a] );

    }


    for( int a=0 ; a<num_nodes ; ++a )
    {
      Temp = Vel[a];
      Temp *= deltatime;
      Disp[a] += Temp;

//      Disp[a].plus_cA( deltatime, Vel[a] );

    }

  }

}

/**
 * @author Randy Settgast
 * @return cpu usage
 *
 * This function uses the rusage structure to query elapsed system time and user time, and returns
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

