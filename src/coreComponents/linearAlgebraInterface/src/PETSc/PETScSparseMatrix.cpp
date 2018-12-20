#include "PETScSparseMatrix.hpp"

PETScSparseMatrix::PETScSparseMatrix()
{
	// do nothing
}

/* Copy constructor */
PETScSparseMatrix::PETScSparseMatrix( PETScSparseMatrix const &in_matrix )
{
	MatDuplicate(in_matrix._mat, MAT_COPY_VALUES, &_mat);
}

/* Virtual destructor */
// virtual ~PETScSparseMatrix() = default;

// virtual ~PETScSparseMatrix()
// {
// 	if(_mat) MatDestroy(_mat);
// }

/* Create a square matrix from number of rows */
void PETScSparseMatrix::create( MPI_Comm const comm,
								int const m_nRowGlobal,
								int const nMaxEntriesPerRow )
{
	// set up matrix
	MatCreate(comm, &_mat);
	MatSetType(_mat, MATMPIAIJ);
	MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, m_nRowGlobal, m_nRowGlobal);
	MatMPIAIJSetPreallocation(_mat, nMaxEntriesPerRow, NULL, nMaxEntriesPerRow, NULL);
	MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	MatSetUp(_mat);

	// assemble
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	assembled = true;
}

 /* Create rectangular matrix */
 void PETScSparseMatrix::create( MPI_Comm const comm,
              					 int const m_nRowGlobal,
              					 int const m_nColGlobal,
              					 int const nMaxEntriesPerRow )
 {
	// set up matrix
	MatCreate(comm, &_mat);
	MatSetType(_mat, MATMPIAIJ);
	MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, m_nRowGlobal, m_nColGlobal);
	MatMPIAIJSetPreallocation(_mat, nMaxEntriesPerRow, NULL, nMaxEntriesPerRow, NULL);
	MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	MatSetUp(_mat);

	// assemble
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	assembled = true;
}

  /* Create square matrix */
  void PETScSparseMatrix::create( MPI_Comm const comm,
               					  int const m_nRowGlobal,
               					  std::vector<int> const nMaxEntriesPerRow )
  {
  	// set up matrix
	MatCreate(comm, &_mat);
	MatSetType(_mat, MATMPIAIJ);
	MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, m_nRowGlobal, m_nRowGlobal);
	MatMPIAIJSetPreallocation(_mat, 0, nMaxEntriesPerRow.data(), 0, nMaxEntriesPerRow.data());
	MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	MatSetUp(_mat);

	// assemble
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	assembled = true;
  }

  /* Create matrix */
  void PETScSparseMatrix::create( MPI_Comm const comm,
               int const m_nRowGlobal,
               int const m_nColGlobal,
               std::vector<int> const nMaxEntriesPerRow )
  {
  	// set up matrix
	MatCreate(comm, &_mat);
	MatSetType(_mat, MATMPIAIJ);
	MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, m_nRowGlobal, m_nColGlobal);
	MatMPIAIJSetPreallocation(_mat, 0, nMaxEntriesPerRow.data(), 0, nMaxEntriesPerRow.data());
	MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	MatSetUp(_mat);

	// assemble
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	assembled = true;
  }

 /* Create from a PETSc matrix */
 void PETScSparseMatrix::create( PETScSparseMatrix &matrix )
 {
 	MatDuplicate(matrix._mat, MAT_COPY_VALUES, &_mat);
 }

 /* Set all elements to zero, keep partition and sparsity pattern */
 void PETScSparseMatrix::zero()
 {
 	MatZeroEntries(_mat);
 }

 /* Empty function for Trilinos and PETSc implementation. Is required when the HYPRE library is used. */
 void PETScSparseMatrix::open()
 {
 	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
 }

 /* Assemble the matrix */
 void PETScSparseMatrix::close()
 {
 	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
 	assembled = true;
 }

//  /* Add values into iRow at columns cols. 
//   * iRow: global row index
//   * nCols: number of columns to modify
//   * values: values to add to prescribed locations
//   * cols: global column indices in which to add the values */
 void PETScSparseMatrix::add( int const iRow,
           					  int const nCols,
           					  double const *values,
           					  int const *cols )
 {
 	int rows[1] = {iRow};

 	MatSetValues(_mat, 1, rows, nCols, cols, values, ADD_VALUES);

 	// assemble
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	assembled = true;
 }

 /* Add value to position (iRow, iCol) */
 void PETScSparseMatrix::add( int const iRow,
           					  int const iCol,
           					  double const value )
 {
 	MatSetValue(_mat, iRow, iCol, value, ADD_VALUES);

 	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	assembled = true;
 }

//  /* Set values into iRow at columns cols. 
//   * iRow: global row index
//   * nCols: number of columns to modify
//   * values: values to add to prescribed locations
//   * cols: global column indices in which to add the values */
 void PETScSparseMatrix::set( int const iRow,
           					  int const nCols,
           					  double const *values,
           					  int const *cols )
  {
 	int rows[1];
 	rows[0] = iRow;

 	MatSetValues(_mat, 1, rows, nCols, cols, values, INSERT_VALUES);

 	// assemble
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	assembled = true;
 }

//  /* Set the value of (iRow, iCol) to value */
 void PETScSparseMatrix::set( int const iRow,
           					  int const iCol,
           					  double const value )
  {
 	MatSetValue(_mat, iRow, iCol, value, INSERT_VALUES);

 	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	assembled = true;
 }

//  /* Insert values into iRow at columns cols. 
//   * iRow: global row index
//   * nCols: number of columns to modify
//   * values: values to add to prescribed locations
//   * cols: global column indices in which to add the values */
 void PETScSparseMatrix::insert( int const iRow,
              					 int const nCols,
              					 double const *values,
              					 int const *cols )
  {
 	int rows[1];
 	rows[0] = iRow;

 	MatSetValues(_mat, 1, rows, nCols, cols, values, INSERT_VALUES);

 	// assemble
	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	assembled = true;
 }

/* Matrix-vector multiplication, stored in dst */
void PETScSparseMatrix::multiply( PETScVector &src,
            					  PETScVector &dst ) const
{
	// need to make sure dst is not NULL
	MatMult(_mat, src.getVec(), dst.getVec());
}

 /* Compute residual r = Ax - b */
void PETScSparseMatrix::residual( PETScVector &x,
            					  PETScVector &b,
            					  PETScVector &res ) const
{

	MatMult(_mat, x.getVec(), res.getVec());
	VecAXPY(res.getVec(), -1, b.getVec());

}

 /* Compute "gaxpy" alpha*A*x + beta*b */
 void PETScSparseMatrix::gaxpy( double alpha,
             					PETScVector  &x,
             					double beta,
             					PETScVector  &b,
             					PETScVector &res,
             					bool useTranspose)
 {
 	PETScVector x_(x);
 	PETScVector y_(x);
 	PETScVector b_(b);

 	x_.scale(alpha); // alpha*x_
 	b_.scale(beta); // beta*b_

 	MatMult(_mat, x_.getVec(), y_.getVec()); // alpha*A*x_ = y_
 	VecAXPY(b_.getVec(), 1, y_.getVec());

 	res = b_;
 }

 /* Multiply all elements by scalingFactor */
 void PETScSparseMatrix::scale( double scalingFactor )
 {
 	MatScale(_mat, scalingFactor);
 }

//  /* Left multiply with diagonal matrix with values in vec */
 void PETScSparseMatrix::leftScale( PETScVector &vec )
 {
 	MatDiagonalScale(_mat, vec.getVec(), NULL);
 }

 /* Right multiply with diagonal matrix with values in vec */
 void PETScSparseMatrix::rightScale( PETScVector &vec )
 {
 	MatDiagonalScale(_mat, NULL, vec.getVec());
 }

 /* Left multiply with diagonal matrix with values vecLeft and right multiply
    with diagonal matrix with values vecRight */
 void PETScSparseMatrix::leftRightScale( PETScVector &vecLeft,
                      PETScVector &vecRight )
 {
 	MatDiagonalScale(_mat, vecLeft.getVec(), vecRight.getVec());
 }

 /* Clear row and multiply diagonal by factor */
 void PETScSparseMatrix::clearRow( int const row,
                				   double const factor )
 {
 	int rows[1] = {row};
 	// get diagonal entry
 	MatZeroRows(_mat, 1, rows, factor, NULL, NULL);
 }

 /*
  * Get global row myRow
  * - numEntries: number of nonzeros 
  * - Values: array of values
  * - Indices: array of column indices */
 // void PETScSparseMatrix::getRow( int GlobalRow,
 //               					int &NumEntries,
 //               					double* Values,
 //               					int* Indices )
 // {
 // 	MatGetRow(_mat, GlobalRow, &NumEntries, &Indices, &Values);
 // }

 /*
  * Get global row myRow
  * - numEntries: number of nonzeros 
  * - vecValues: vector of values
  * - vecIndices: vector of column indices */
 // void PETScSparseMatrix::getRow( int GlobalRow,
 //              					int &NumEntries,
 //              					std::vector<double> &vecValues,
 //              					std::vector<int> &vecIndices );

 /*
  * Get local row myRow
  * - numEntries: number of nonzeros 
  * - Values: array of values
  * - Indices: array of column indices */
 // void PETScSparseMatrix::getLocalRow( int myRow,
 //                   					 int & NumEntries,
 //                   					 double * & Values,
 //                   					 int * & Indices );

//  /*
//   * Get local row myRow
//   * - numEntries: number of nonzeros 
//   * - vecValues: vector of values
//   * - vecIndices: vector of column indices */
//  // void PETScSparseMatrix::getLocalRow( int myRow,
//  //                   					int &NumEntries,
//  //                   					std::vector<double> &vecValues,
//  //                   					std::vector<int> &vecIndices );

//  /* Pointer to the underlying PETSc matrix */
//  // Mat* PETScSparseMatrix::getPointer() const;

 /* Global number of rows */
 int PETScSparseMatrix::globalRows() const
 {
	int num_rows;
	int num_cols;
	MatGetSize(_mat, &num_rows, &num_cols);
	return num_rows;
 }

 /* Global number of columns */
 int PETScSparseMatrix::globalCols() const
 {
 	int num_rows;
	int num_cols;
	MatGetSize(_mat, &num_rows, &num_cols);
	return num_cols;
 }

 /* Number of unique columns */
 // int PETScSparseMatrix::uniqueCols() const;

 /* Index of first global row owned by a process */
 int PETScSparseMatrix::ilower() const
 {
 	int firstrow;
 	int lastrow;
 	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
 	return firstrow; 
 }

 /* Index of last global row owned by a process */
int PETScSparseMatrix::iupper() const
{
 	int firstrow;
 	int lastrow;
 	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
 	return lastrow - 1; 
 } 

 /* Local number of rows */
int PETScSparseMatrix::myRows() const
{
	int firstrow;
	int lastrow;
	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
	return lastrow - firstrow; 
}

 /* Local number of columns */
 int PETScSparseMatrix::myCols() const
{
	int firstcol;
	int lastcol;
	MatGetOwnershipRangeColumn(_mat, &firstcol, &lastcol);
	return lastcol - firstcol; 
}

 /* Infinity norm */
 double PETScSparseMatrix::normInf() const
 {
	double normInf;
	MatNorm(_mat, NORM_INFINITY, &normInf);
	return normInf;
 }

 /* 1-norm */
 double PETScSparseMatrix::norm1() const
 {
	double norm1;
	MatNorm(_mat, NORM_1, &norm1);
	return norm1;
 }

 /* Frobenius norm */
 double PETScSparseMatrix::normFrobenius() const
 {
	double normFrob;
	MatNorm(_mat, NORM_FROBENIUS, &normFrob);
	return normFrob;
 }

 /* Check if the matrix is assembled */
 bool PETScSparseMatrix::isAssembled() const
 {
 	return assembled;
 }

  /* Print to terminal */
  void PETScSparseMatrix::print() const 
  {
  	MatView(_mat, PETSC_VIEWER_STDOUT_WORLD);
  }

  const Mat* PETScSparseMatrix::getPointer() const
  {
  	return &(_mat);
  }

  Mat PETScSparseMatrix::getMat()
  {
  	return _mat;
  }

