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

/// returns the amount of cpu time use for this process
realT getcputime(void)  ;


int main(int argc, char* argv[] )
{


  realT iRANDMAX = 1.0 / RAND_MAX;

  std::srand(1234);

  const int num_nodes = atoi(argv[1]);
  const int num_steps = atoi(argv[2]);


  R1TensorT<3>* const Disp = new R1TensorT<3>[num_nodes];
  R1TensorT<3> Tot;

  realT dt[num_steps];

  realT* const xdisp = new realT[num_nodes];
  realT* const ydisp = new realT[num_nodes];
  realT* const zdisp = new realT[num_nodes];
  realT xtot = 0.0 ;
  realT ytot = 0.0 ;
  realT ztot = 0.0 ;


  for( int a=0 ; a<num_nodes ; ++a )
  {
    xdisp[a] = realT(rand()) * iRANDMAX;
    ydisp[a] = realT(rand()) * iRANDMAX;
    zdisp[a] = realT(rand()) * iRANDMAX;
    Disp[a](0) = xacc[a];
    Disp[a](1) = yacc[a];
    Disp[a](2) = zacc[a];
  }

  for( int i=0 ; i<num_steps ; ++i )
  {
    dt[i] = realT(rand()) * iRANDMAX;
  }

  int flag = atoi(argv[3]) ;

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
  else
  {
  function2(  Disp, Vel, Acc, &Tot,
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

  delete [] Disp;

  return 0;
}

inline
void CalculateGradient( realT Gradient[3][3] ,
                        const realT* const x,
                        const realT* const y,
                        const realT* const z,
                        const realT* const * const dNdX )

{

  for( int i=0 ; i<3 ; ++i )
    for( int j=0 ; j<3 ; ++j )
    {
      Gradient[i][j] = 0.0;
    }

  for( int a=0 ; a<8 ; ++a )
  {
    for( int i=0 ; i<3 ; ++i )
    {
      Gradient[0][i] += x[a] * dNdX[a][i];
      Gradient[1][i] += y[a] * dNdX[a][i];
      Gradient[2][i] += z[a] * dNdX[a][i];
    }
  }

}


void CalculateGradient( R2TensorT<nsdof>& Gradient ,
                        const Array1dT<R1TensorT<nsdof> >& disp,
                        const Array1dT<R1TensorT<nsdof> >& dNdX )

{
  Gradient = 0.0;
  for( int a=0 ; a<8 ; ++a )
    Gradient.plus_dyadic_ab( disp(a) , dNdX(a));

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

