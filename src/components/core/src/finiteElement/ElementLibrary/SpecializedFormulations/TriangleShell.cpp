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
 * @file TriangleShell.cpp
 * @author Fu, Pengcheng
 * @date August 3, 2012
 */

#include "TriangleShell.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"


TriangleShell::TriangleShell():
FiniteElement<3>(1,3,0)
{
  m_nodeOrdering.resize(3);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 2;

}

TriangleShell::~TriangleShell()
{
  // TODO Auto-generated destructor stub
}

void TriangleShell::reinit(const std::vector<R1TensorT<3> > &mapped_support_points)
{

  assert(mapped_support_points.size() == n_dofs);
  //  printf("n=%d\n", n_dofs);

  //See Chapter 15 of Int. FEM of U Colorado by Carlos Felippa for detailed formulation.
  //http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/
  //Accessed in July 2012

  const std::vector<R1TensorT<3> >& Y = mapped_support_points;

  // transform

  realT l1 = sqrt(pow(Y[1][0] - Y[0][0], 2.0) + pow(Y[1][1] - Y[0][1], 2.0) + pow(Y[1][2] - Y[0][2], 2.0)); 
  realT l2 = sqrt(pow(Y[2][0] - Y[0][0], 2.0) + pow(Y[2][1] - Y[0][1], 2.0) + pow(Y[2][2] - Y[0][2], 2.0)); 

  realT theta = ((Y[1][0] - Y[0][0]) * (Y[2][0] - Y[0][0]) + (Y[1][1] - Y[0][1]) * (Y[2][1] - Y[0][1]) + (Y[1][2] - Y[0][2]) * (Y[2][2] - Y[0][2])) / l1 / l2;


  std::vector<R1TensorT<3> > X(n_dofs);

  X[0][0] = 0.0;
  X[0][1] = 0.0;

  X[1][0] = l1;
  X[1][1] = 0.0;

  X[2][0] = l2 * theta;
  X[2][1] = l2 * sqrt(1.0 - theta * theta);

  realT V;
  const realT half = 1.0 / 2.0;

  V = (X[1][0] * X[2][1] - X[2][0] * X[1][1]) + (X[2][0] * X[0][1] - X[0][0] * X[2][1]) + (X[0][0] * X[1][1] - X[1][0] * X[0][1]);
  V *= half;
  data[0].jacobian_determinant = V ;

  data[0].mapped_gradients[0](0) = (X[1][1] - X[2][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[0](1) = (X[2][0] - X[1][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[0](2) = 0.0;

  data[0].mapped_gradients[1](0) = (X[2][1] - X[0][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[1](1) = (X[0][0] - X[2][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[1](2) = 0.0;


  data[0].mapped_gradients[2](0) = (X[0][1] - X[1][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[2](1) = (X[1][0] - X[0][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[2](2) = 0.0;



}

