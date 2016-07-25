/*
 * MechanicsKernel.h
 *
 *  Created on: Oct 31, 2012
 *      Author: settgast1
 */

#ifndef MECHANICSKERNEL_H_
#define MECHANICSKERNEL_H_

#include "Common/Common.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"



//class MaterialBaseParameterDataT;
template <int dim> class FiniteElement;


class MechanicsKernel
{
public:
  MechanicsKernel();
  virtual ~MechanicsKernel();

  /*
  void ElementStiffnessMatrix_SmallStrainLinearElastic( const MaterialBaseParameterDataT& matParams,
                                                        const FiniteElement<3>& fe,
                                                        const Array2dT<R1Tensor>& dNdX,
                                                        const realT* const detJ,
                                                        Epetra_SerialDenseMatrix& element_matrix );
*/
};

#endif /* MECHANICSKERNEL_H_ */
