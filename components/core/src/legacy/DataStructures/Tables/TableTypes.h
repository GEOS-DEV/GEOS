/*
 * TableTypes.h
 *
 *  Created on: Dec 8, 2013
 *      Author: johnson346
 */

#ifndef TABLETYPES_H_
#define TABLETYPES_H_

typedef Table<1, realT> Table1D;
typedef Table<2, realT> Table2D;
typedef Table<3, realT> Table3D;
typedef Table<4, realT> Table4D;

typedef Table<1, R1Tensor> VectorField1D;
typedef Table<2, R1Tensor> VectorField2D;
typedef Table<3, R1Tensor> VectorField3D;
typedef Table<4, R1Tensor> VectorField4D;

typedef Table<1, R2SymTensor> R2SymTensorField1D;
typedef Table<2, R2SymTensor> R2SymTensorField2D;
typedef Table<3, R2SymTensor> R2SymTensorField3D;
typedef Table<4, R2SymTensor> R2SymTensorField4D;

typedef Table<1, R2Tensor> R2TensorField1D;
typedef Table<2, R2Tensor> R2TensorField2D;
typedef Table<3, R2Tensor> R2TensorField3D;
typedef Table<4, R2Tensor> R2TensorField4D;

#endif /* TABLETYPES_H_ */
