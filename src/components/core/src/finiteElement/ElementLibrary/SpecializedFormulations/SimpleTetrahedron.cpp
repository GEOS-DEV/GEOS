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
 * @file SimpleTetrahedron.cpp
 * @author Fu, Pengcheng
 * @date Apr 30, 2012
 */

#include "SimpleTetrahedron.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"


SimpleTetrahedron::SimpleTetrahedron():
  FiniteElement<3>(1,4,0)
{
  m_nodeOrdering.resize(4);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 2;
  m_nodeOrdering[3] = 3;

}

SimpleTetrahedron::~SimpleTetrahedron()
{
  // TODO Auto-generated destructor stub
}


/**
 * Reinitialize the finite element basis on a particular element.
 * We use the coordinates of the support points in real space to
 * construct the forward mapping from the parent coordinate system.  The
 * support points are assumed to follow a lexicographic ordering:
 * On the parent element, we loop over the x-coordinate fastest,
 * the y, then z (depending on the desired spatial dimension of the
 * element).
 */

void SimpleTetrahedron::reinit(const std::vector<R1TensorT<3> > &mapped_support_points)
{
  //See Chapter 16 of Adv. FEM of U Colorado by Carlos Felippa for detailed
  // formulation.
  //http://www.colorado.edu/engineering/cas/courses.d/AFEM.d/
  //Accessed in July 2012


  assert(mapped_support_points.size() == n_dofs);

  const unsigned int q = 0;

  const std::vector<R1TensorT<3> >& X = mapped_support_points;

  realT V, a[4], b[4], c[4];
  const realT sixth = 1.0 / 6.0;


  //a[1]=(y4-y2)*(z3-z2)-(y3-y2)*(z4-z2)
  a[0] = (X[3][1] - X[1][1]) * (X[2][2] - X[1][2]) - (X[2][1] - X[1][1]) * (X[3][2] - X[1][2]);
  //      b[1]=x32 z42 - x42 z32
  b[0] = (X[2][0] - X[1][0]) * (X[3][2] - X[1][2]) - (X[3][0] - X[1][0]) * (X[2][2] - X[1][2]);
  //      c[1]=x42 y32 - x32 y42
  c[0] = (X[3][0] - X[1][0]) * (X[2][1] - X[1][1]) - (X[2][0] - X[1][0]) * (X[3][1] - X[1][1]);
  //      a[2]=y31 z43 - y34 z13
  a[1] = (X[2][1] - X[0][1]) * (X[3][2] - X[2][2]) - (X[2][1] - X[3][1]) * (X[0][2] - X[2][2]);
  //      b[2]=x43 z31 - x13 z34
  b[1] = (X[3][0] - X[2][0]) * (X[2][2] - X[0][2]) - (X[0][0] - X[2][0]) * (X[2][2] - X[3][2]);
  //      c[2]=x31 y43 - x34 y13
  c[1] = (X[2][0] - X[0][0]) * (X[3][1] - X[2][1]) - (X[2][0] - X[3][0]) * (X[0][1] - X[2][1]);
  //      a[3]=y24 z14 - y14 z24
  a[2] = (X[1][1] - X[3][1]) * (X[0][2] - X[3][2]) - (X[0][1] - X[3][1]) * (X[1][2] - X[3][2]);
  //      b[3]=x14 z24 - x24 z14
  b[2] = (X[0][0] - X[3][0]) * (X[1][2] - X[3][2]) - (X[1][0] - X[3][0]) * (X[0][2] - X[3][2]);
  //      c[3]=x24 y14 - x14 y24
  c[2] = (X[1][0] - X[3][0]) * (X[0][1] - X[3][1]) - (X[0][0] - X[3][0]) * (X[1][1] - X[3][1]);
  //      a[4]=y13 z21 - y12 z31
  a[3] = (X[0][1] - X[2][1]) * (X[1][2] - X[0][2]) - (X[0][1] - X[1][1]) * (X[2][2] - X[0][2]);
  //      b[4]=x21 z13 - x31 z12
  b[3] = (X[1][0] - X[0][0]) * (X[0][2] - X[2][2]) - (X[2][0] - X[0][0]) * (X[0][2] - X[1][2]);
  //      c[4]=x13 y21 - x12 y31
  c[3] = (X[0][0] - X[2][0]) * (X[1][1] - X[0][1]) - (X[0][0] - X[1][0]) * (X[2][1] - X[0][1]);
  //
  //6V=x21 (y23 z34 - y34 z23) + x32 (y34 z12 - y12 z34) + x43 (y12 z23 - y23
  // z12),
  V = (X[1][0] - X[0][0]) * ((X[1][1] - X[2][1]) * (X[2][2] - X[3][2]) - (X[2][1] - X[3][1]) * (X[1][2] - X[2][2])) + (X[2][0] - X[1][0]) *
      ((X[2][1] - X[3][1]) * (X[0][2] - X[1][2]) - (X[0][1] - X[1][1]) * (X[2][2] - X[3][2])) + (X[3][0] - X[2][0]) *
      ((X[0][1] - X[1][1]) * (X[1][2] - X[2][2]) - (X[1][1] - X[2][1]) * (X[0][2] - X[1][2]));
  V *= sixth;

  data[q].jacobian_determinant = V;
  for( int iNd=0 ; iNd<4 ; ++iNd )
  {
    data[q].mapped_gradients[iNd](0) = a[iNd] * sixth / data[q].jacobian_determinant;
    data[q].mapped_gradients[iNd](1) = b[iNd] * sixth / data[q].jacobian_determinant;
    data[q].mapped_gradients[iNd](2) = c[iNd] * sixth / data[q].jacobian_determinant;

    //std::cout<<"data[q].mapped_gradients[a] =
    // "<<data[q].mapped_gradients[a]<<std::endl;
  }


}
