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
 * @file IntegrationPointsT.cpp
 * @author settgast1
 * @date Dec 6, 2010
 */

#include "IntegrationRuleT.h"

IntegrationRuleT::IntegrationRuleT()
{
  // TODO Auto-generated constructor stub

}

IntegrationRuleT::~IntegrationRuleT()
{
  // TODO Auto-generated destructor stub
}


void IntegrationRuleT::CalculateShapeFunctionDerivatives( const array<R1Tensor>& X,
                                                        array<R1Tensor>& dNdX,
                                                        realT& detJ )
{
  const realT x[8] = { X(0)(0),
                        X(1)(0),
                        X(2)(0),
                        X(3)(0),
                        X(4)(0),
                        X(5)(0),
                        X(6)(0),
                        X(7)(0) };

  const realT y[8] = { X(0)(1),
                        X(1)(1),
                        X(2)(1),
                        X(3)(1),
                        X(4)(1),
                        X(5)(1),
                        X(6)(1),
                        X(7)(1) };

  const realT z[8] = { X(0)(2),
                        X(1)(2),
                        X(2)(2),
                        X(3)(2),
                        X(4)(2),
                        X(5)(2),
                        X(6)(2),
                        X(7)(2) };

  realT b[3][8];

  CalculateShapeFunctionDerivative(  y , z , b[0] ) ;
  CalculateShapeFunctionDerivative(  z , x , b[1] ) ;
  CalculateShapeFunctionDerivative(  x , y , b[2] ) ;

  detJ = 0.0;
  for( int a=0 ; a<8 ; ++a )
  {
    detJ += x[a]*b[0][a] + y[a]*b[1][a] + z[a]*b[2][a] ;
  }
  detJ /= 3.0;

  for( int a=0 ; a<8 ; ++a )
  {
    dNdX(a)(0) = b[0][a] / detJ;
    dNdX(a)(1) = b[1][a] / detJ;
    dNdX(a)(2) = b[2][a] / detJ;
  }

}


