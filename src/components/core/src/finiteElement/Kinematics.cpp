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
/**
 * @file Kinematics.cpp
 * @author settgast1
 * @date Dec 20, 2010
 */


#include "Kinematics.h"
namespace geosx
{
void IncrementalKinematics( const R2TensorT<3>& A ,
                            R2SymTensorT<3>& Dadt ,
                            R2TensorT<3>& Rhat )
{
  R2SymTensorT<3> C_I ;
  R2SymTensorT<3> C_I2 ;
  //static R2SymTensorT<3> C_I3 ;
//  static R2TensorT<3> Amod ;
//  static R2TensorT<3> UhatA ;

  C_I.AijAkj_m_Aik_m_Aki(A);
  C_I2.AijAjk(C_I);
  //C_I3.AijBjk(C_I2,C_I);

  C_I *= -0.5;
  C_I2 *= 0.25;
  //C_I3 *= -1.0/6.0;

  Dadt = C_I;
  Dadt += C_I2;
  //Dadt += C_I3;


/*
  // ***** Precondition A for production of Rhat *****

  // Estimate Uhat-I store in C_I
  C_I2 *= 1.5;
  C_I3 *= 15.0/8.0;

  C_I += C_I2;
  C_I += C_I3;

  // Add I-Uhat to Amod
  Amod = C_I;
  Amod *= -1.0;

  // Make Uhat
  C_I.PlusIdentity(1.0);

  //store Uhat*A
  UhatA.AijBjk(C_I,A);

  // Construct Modified A
  Amod += UhatA;

  IncrementalRotation( Amod , Rhat );
*/
  IncrementalRotation( A , Rhat );


}

void IncrementalRotation( const R2TensorT<3>& A ,
                          R2TensorT<3>& Rot )
{
//  realT alpha[3];
  R1TensorT<3> alpha;

  alpha.eijkAjk(A);
//  alpha *=-1.0;

  realT Q = 0.25 * Dot(alpha,alpha);


  realT trA = A.Trace();
  realT trFhatinv_1 = ( 2.0 - trA );
  realT P = 0.25 * pow( trFhatinv_1, 2 );
//  static realT one_P = trA - 0.25 * pow(trA,2);

//  std::cout<<"P,Q = "<<P<<'\t'<<Q<<'\n';

  realT P2 = pow(P,2);
  realT P3 = P2*P;
//  static realT P4 = P3*P;

  realT Q2 = pow(Q,2);
//  static realT Q3 = Q2*Q;

  realT PpQ = P + Q;
  realT PpQ2 = pow(PpQ,2);
  realT PpQ3 = PpQ2*PpQ;

  realT term1 = sqrt( std::max( P + 3*P2*(1-PpQ)/PpQ2  - 2*P3*(1-PpQ)/PpQ3, 0.0 ) ) ;

  realT term2 = 0;

  if( fabs(Q) > 0.01 )
  {
    if( trFhatinv_1 > 0.0 )
      term2 = (1.0 - trFhatinv_1 / fabs(trFhatinv_1) * term1) / ( 4*Q );
    else
      term2 = (1.0) / ( 4*Q );

  }
  else
    term2 = 0.125 + Q * ( P2-12.0*(P-1.0) ) / (32 * P2)
                        + Q2*( (P-2.0)*(P2-10*P+32) ) / (64*P3);
//                        + Q3*( 1104 - 992*P + 376*P2 - 72*P3 +5*P4) / (512*P4) ;

  realT term3 = 0.5 * sqrt( ( P*Q*(3.0-Q) + P3 + Q2 ) / PpQ3 ) ;

//    cout<<term1<<' '<<term2<<' '<<term3<<endl;

    Rot.dyadic_aa(alpha);
    Rot *= term2;
    Rot.PlusIdentity(term1);

    alpha *= term3;

    Rot(0,1) += alpha(2);
    Rot(0,2) -= alpha(1);
    Rot(1,2) += alpha(0);
    Rot(1,0) -= alpha(2);
    Rot(2,0) += alpha(1);
    Rot(2,1) -= alpha(0);

/*
    Rot.t_data[1] += alpha.t_data[2];
    Rot.t_data[2] -= alpha.t_data[1];
    Rot.t_data[5] += alpha.t_data[0];
    Rot.t_data[3] -= alpha.t_data[2];
    Rot.t_data[6] += alpha.t_data[1];
    Rot.t_data[7] -= alpha.t_data[0];
*/

}






void CalculatePhantomGradient( R2TensorT<3>& Gradient ,
                               const int* bConnectivity,
                               const Array1dT<R1TensorT<3> >& disp,
                               const Array2dT<R1TensorT<3> >& dNdX )

{
  Gradient = 0.0;
  for( int side=1 ; side<=2 ; ++side )
    for( int a=1 ; a<=4 ; ++a )
      Gradient.plus_dyadic_ab( disp(bConnectivity[(side-1)*4+a-1]) , dNdX(side,a));

}

}
