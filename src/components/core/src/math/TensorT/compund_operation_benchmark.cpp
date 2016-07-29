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
#include <vector>

/// returns the amount of cpu time use for this process
realT getcputime(void)  ;



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
    a[i] *= iRANDMAX ;
    b[i] *= iRANDMAX ;
    c[i] *= iRANDMAX ;
    d[i] *= iRANDMAX ;
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
    std::cout<<e[i]<<std::endl;
  std::cout<<std::endl;


  for( int i=0 ; i<3 ; ++i )
    for( int j=0 ; j<3 ; ++j )
    {
      std::cout<<E(i,j)<<std::endl;
    }
  std::cout<<std::endl;

  std::cout<<"baseline CPU time    = "<<t2-t1<<std::endl;
  std::cout<<"TensorClass CPU time = "<<t3-t2<<std::endl;


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
