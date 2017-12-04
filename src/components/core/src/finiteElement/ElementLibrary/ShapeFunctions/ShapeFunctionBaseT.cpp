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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * File: ShapeFunctionBaseT.h
 * Shape function Class
 *
 * created : RRS (09/14/2010)
 */

#include "ShapeFunctionBaseT.h"

ShapeFunctionBaseT::ShapeT() :
  rElement(0),
  pNodes(0)
{}

ShapeFunctionBaseT::~ShapeFunctionBaseT(void)
{}



void ShapeFunctionBaseT::CalculateJacobian(const int elem)
{
  J = 0.0;
  for( int i=1 ; i<=nsdof ; ++i )
    for( int j=1 ; j<=nsdof ; ++j )
      for( int a=1 ; a<=rElement->NumNodeElem() ; ++a )
        J(i,j) += pNodes->Xref( rElement->Connectivity(elem,a),i)
                  * dNdXi(a)(j);

}

void ContinuumShapeT::Calc_Shape_Deriv(const realT fac)
{

  realT factor = 1.0 / pow(static_cast<realT>(2),nsdof);

  // loop over "nodal indexed" shape functions and calculate the
  // shape function derivitives wrt the local coordinates
  for( int a=1 ; a<=rElement->NumNodeElem() ; ++a )
  {
    for( int i=1 ; i<=nsdof ; ++i )
    {
      dNdXi(a)(i) = factor * rElement->Xi_node(a,i);
      for( int j=1 ; j<=nsdof ; j++ )
        if( i!=j )
          dNdXi(a)(i) *= ( 1 + ip_coord_fac * rElement->Xi_node(ip,j) * rElement->Xi_node(a,j));

    }
    //std::cout<<a<<' '<<dNdXi(a)<<std::endl;
  }

  CalculateJacobian(elem);

  m_detJ(elem,ip) = J.Det();

  dNdX.AijBi(Jinv,dNdXi);


}


inline realT ShapeFunctionBaseT::Shape( const R1Tensor& Xi,
                                        const R1Tensor& Xi_node )
{
  realT N = 1.0 / pow(2.0,nsdof);

  for( int i=1 ; i<=nsdof ; ++i )
    N *= ( 1 +  Xi(i) * Xi_node(i) );

  return N;
}
