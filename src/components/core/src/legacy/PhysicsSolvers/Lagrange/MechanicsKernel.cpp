/*
 * MechanicsKernel.cpp
 *
 *  Created on: Oct 31, 2012
 *      Author: settgast1
 */

#include "MechanicsKernel.h"

#include "ElementLibrary/FiniteElement.h"

MechanicsKernel::MechanicsKernel()
{
  // TODO Auto-generated constructor stub

}

MechanicsKernel::~MechanicsKernel()
{
  // TODO Auto-generated destructor stub
}
/*
void MechanicsKernel::ElementStiffnessMatrix_SmallStrainLinearElastic( const MaterialBaseParameterDataT& matParams,
                                                                       const FiniteElement<3>& fe,
                                                                       const Array2dT<R1Tensor>& dNdX,
                                                                       const realT* const detJ,
                                                                       Epetra_SerialDenseMatrix& element_matrix )
{
  const int dim = 3;

  const realT G = matParams.ShearModulus;
  const realT lambda = matParams.Lame;

  element_matrix.Scale(0);

  R1Tensor dNdX_a;
  R1Tensor dNdX_b;

  for( unsigned int q=0 ; q<fe.n_quadrature_points() ; ++q )
  {
    for( unsigned int a=0 ; a<fe.dofs_per_element() ; ++a )
    {
      dNdX_a = dNdX(q,a);

      for( unsigned int b=0 ; b<fe.dofs_per_element() ; ++b )
      {
        dNdX_b = dNdX(q,b);

        element_matrix(a*dim+0,b*dim+0) = ( ( dNdX_a[1] * dNdX_b[1] + dNdX_a[2] * dNdX_b[2] ) * G + dNdX_a[0] * dNdX_b[0] * (2 * G + lambda) ) * detJ[q];
        element_matrix(a*dim+0,b*dim+1) = ( ( dNdX_a[1] * dNdX_b[0] ) * G + dNdX_a[0] * dNdX_b[1] * lambda ) * detJ[q];
        element_matrix(a*dim+0,b*dim+2) = ( ( dNdX_a[2] * dNdX_b[0] ) * G + dNdX_a[0] * dNdX_b[2] * lambda ) * detJ[q];

        element_matrix(a*dim+1,b*dim+0) = ( ( dNdX_a[0] * dNdX_b[1] ) * G + dNdX_a[1] * dNdX_b[0] * lambda ) * detJ[q];
        element_matrix(a*dim+1,b*dim+1) = ( ( dNdX_a[0] * dNdX_b[0] + dNdX_a[2] * dNdX_b[2] ) * G + dNdX_a[1] * dNdX_b[1] * (2 * G + lambda) ) * detJ[q];
        element_matrix(a*dim+1,b*dim+2) = ( ( dNdX_a[2] * dNdX_b[1] ) * G + dNdX_a[1] * dNdX_b[2] * lambda ) * detJ[q];

        element_matrix(a*dim+2,b*dim+0) = ( ( dNdX_a[0] * dNdX_b[2] ) * G + dNdX_a[2] * dNdX_b[0] * lambda ) * detJ[q];
        element_matrix(a*dim+2,b*dim+1) = ( ( dNdX_a[1] * dNdX_b[2] ) * G + dNdX_a[2] * dNdX_b[1] * lambda ) * detJ[q];
        element_matrix(a*dim+2,b*dim+2) = ( ( dNdX_a[0] * dNdX_b[0] + dNdX_a[1] * dNdX_b[1] ) * G + dNdX_a[2] * dNdX_b[2] * (2 * G + lambda) ) * detJ[q];


      }
    }
  }

}
*/
