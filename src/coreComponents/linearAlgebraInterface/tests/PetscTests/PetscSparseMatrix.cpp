#include "PetscSparseMatrix.hpp"

PetscSparseMatrix::PetscSparseMatrix()
{
	// do nothing
}

/* Copy constructor */
PetscSparseMatrix::PetscSparseMatrix( PetscSparseMatrix const &in_matrix )
{
	MatDuplicate(in_matrix._mat, MAT_COPY_VALUES, &_mat);
}

/* Virtual destructor */
// virtual ~PetscSparseMatrix() = default;

// virtual ~PetscSparseMatrix()
// {
// 	if(_mat) MatDestroy(_mat);
// }

/* Create a square matrix from number of rows */
void PetscSparseMatrix::create( MPI_Comm const comm,
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

	// assembled = true;
}

 /* Create rectangular matrix */
 void PetscSparseMatrix::create( MPI_Comm const comm,
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
  void PetscSparseMatrix::create( MPI_Comm const comm,
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
  void PetscSparseMatrix::create( MPI_Comm const comm,
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
 void PetscSparseMatrix::create( PetscSparseMatrix &matrix )
 {
 	MatDuplicate(matrix._mat, MAT_COPY_VALUES, &_mat);
 }

 /* Set all elements to zero, keep partition and sparsity pattern */
 void PetscSparseMatrix::zero()
 {
 	MatZeroEntries(_mat);
 }

 /* Empty function for Trilinos and PETSc implementation. Is required when the HYPRE library is used. */
 void PetscSparseMatrix::open()
 {
 	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
 }

 /* Assemble the matrix */
 void PetscSparseMatrix::close()
 {
 	MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
 	MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);
 	assembled = true;
 }

//  /* Add values into iRow at columns cols. 
//   * iRow: global row index
//   * nCols: number of columns to modify
//   * values: values to add to prescribed locations
//   * cols: global column indices in which to add the values */
 void PetscSparseMatrix::add( int const iRow,
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

 // /* add values to sparse matrix located at indices */
 // void PetscSparseMatrix::add( array1d<int> const rowIndices,
 //           array1d<int> const colIndices,
 //           array2d<double> const values)
 // {

 // }

 /* Add value to position (iRow, iCol) */
 void PetscSparseMatrix::add( int const iRow,
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
 void PetscSparseMatrix::set( int const iRow,
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
 void PetscSparseMatrix::set( int const iRow,
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
 void PetscSparseMatrix::insert( int const iRow,
              					 int const nCols,
              					 double const *values,
              					 int const *cols )
  {
 	int rows[1];
 	rows[0] = iRow;

 	MatSetValues(_mat, 1, rows, nCols, cols, values, INSERT_VALUES);

 	// HANNAH: when do we assemble?
	// MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
	// MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

	// assembled = true;
 }

/* Matrix-vector multiplication, stored in dst */
void PetscSparseMatrix::multiply( PetscVector const &src,
            					  PetscVector &dst ) const
{
	// need to make sure dst is not NULL
	MatMult(_mat, src.getConstVec(), dst.getVec());
}

 /* Compute residual r = Ax - b */
void PetscSparseMatrix::residual( PetscVector const &x,
            					  PetscVector const &b,
            					  PetscVector &res ) const
{

	MatMult(_mat, x.getConstVec(), res.getConstVec());
	VecAXPY(res.getVec(), -1, b.getConstVec());

}

 /* Compute "gemv" y = alpha*A*x + beta*y */
 void PetscSparseMatrix::gemv( double const alpha,
             					PetscVector  const &x,
             					double const beta,
             					PetscVector  &y,
             					bool useTranspose)
 {
 	PetscVector x_(x);
 	PetscVector b_(x);

 	x_.scale(alpha); // alpha*x_
 	y.scale(beta); // beta*y

 	if (useTranspose){
 		MatMultTranspose(_mat, x_.getVec(), b_.getVec());
 	} else {
 		MatMult(_mat, x_.getVec(), b_.getVec()); // alpha*A*x_ = b_
 	}
 	VecAXPY(y.getVec(), 1, b_.getVec()); // alpha*A*x_ + beta*y = y
 }

 /* Multiply all elements by scalingFactor */
 void PetscSparseMatrix::scale( double const scalingFactor )
 {
 	 MatScale(_mat, scalingFactor);
 }

//  /* Left multiply with diagonal matrix with values in vec */
 void PetscSparseMatrix::leftScale( PetscVector const &vec )
 {
  	MatDiagonalScale(_mat, vec.getConstVec(), NULL);
 }

 /* Right multiply with diagonal matrix with values in vec */
 void PetscSparseMatrix::rightScale( PetscVector const &vec )
 {
  	MatDiagonalScale(_mat, NULL, vec.getConstVec());
 }

 /* Left multiply with diagonal matrix with values vecLeft and right multiply
    with diagonal matrix with values vecRight */
 void PetscSparseMatrix::leftRightScale( PetscVector const &vecLeft,
                      					 PetscVector const &vecRight )
 {
 	 MatDiagonalScale(_mat, vecLeft.getConstVec(), vecRight.getConstVec());
 }

 /* Clear row and multiply diagonal by factor */
 void PetscSparseMatrix::clearRow( int const row,
                				   double const factor )
 {
 	 int rows[1] = {row};
 	 // HANNAH: not exactly functional 
 	 double diag = 1;

   // get processor this belongs to
   int firstRow, lastRow;
   MatGetOwnershipRange(_mat, &firstRow, &lastRow);

   if(firstRow <= row && row < lastRow)
   {
    // get diagonal entry
     const double* vals;
     const int* inds;
     int numEntries;

    MatGetRow(_mat, row, &numEntries, &inds, &vals);

    for(int i = 0; i < numEntries; i++)
    {
      if(inds[i] == row) diag = vals[i];
    }
    MatRestoreRow(_mat, row, &numEntries, &inds, &vals);
   }
   // zero row and multiply diagonal by factor
    MatZeroRows(_mat, 1, rows, diag*factor, NULL, NULL);
 }

 /*
  * Get global row myRow
  * - numEntries: number of nonzeros 
  * - Values: array of values
  * - Indices: array of column indices */
 void PetscSparseMatrix::getRow(int GlobalRow,
               					int &NumEntries,
               					double* Values,
               					int* Indices ) const
 {
 	const double* vals;
 	const int* inds;
 	int numEntries;
  // 
	MatGetRow(_mat, GlobalRow, &numEntries, &inds, &vals);

    for(int i = 0; i < numEntries; i++)
    {
    	Values[i] = vals[i];
    	Indices[i] = inds[i];
    }
    NumEntries = numEntries;
    
	MatRestoreRow(_mat, GlobalRow, &numEntries, &inds, &vals);
 }

 /*
  * Get global row myRow
  * - numEntries: number of nonzeros 
  * - vecValues: vector of values
  * - vecIndices: vector of column indices */
 void PetscSparseMatrix::getRow( int GlobalRow,
              					int &NumEntries,
              					std::vector<double> &vecValues,
              					std::vector<int> &vecIndices ) const
 {
 	const double* vals;
 	const int* inds;
 	int numEntries;
	MatGetRow(_mat, GlobalRow, &numEntries, &inds, &vals);

 	vecIndices.assign(inds, inds + numEntries);
    vecValues.assign(vals, vals + numEntries);
 	NumEntries = numEntries;

 	MatRestoreRow(_mat, GlobalRow, &numEntries, &inds, &vals);
 }

 /*
  * Get local row myRow
  * - numEntries: number of nonzeros 
  * - Values: array of values
  * - Indices: array of column indices */
 void PetscSparseMatrix::getLocalRow( int myRow,
                   					 int &NumEntries,
                   					 double *Values,
                   					 int *Indices ) const
 {
 	// myRow -> globalRow
 	int firstRow, lastRow;
 	int GlobalRow;
 	MatGetOwnershipRange(_mat, &firstRow, &lastRow);
 	GlobalRow = firstRow + myRow;

 	const double* vals;
 	const int* inds;
 	int numEntries;
	MatGetRow(_mat, GlobalRow, &numEntries, &inds, &vals);

    for(int i = 0; i < numEntries; i++)
    {
    	Values[i] = vals[i];
    	Indices[i] = inds[i];
    }
    NumEntries = numEntries;
    
	MatRestoreRow(_mat, GlobalRow, &numEntries, &inds, &vals);
 }

 /*
  * Get local row myRow
  * - numEntries: number of nonzeros 
  * - vecValues: vector of values
  * - vecIndices: vector of column indices */
 void PetscSparseMatrix::getLocalRow( int myRow,
                   					int &NumEntries,
                   					std::vector<double> &vecValues,
                   					std::vector<int> &vecIndices ) const
 {
 	// myRow -> globalRow
 	int firstRow, lastRow;
 	int GlobalRow;
 	MatGetOwnershipRange(_mat, &firstRow, &lastRow);
 	GlobalRow = firstRow + myRow;

 	const double* vals;
 	const int* inds;
 	int numEntries;
	MatGetRow(_mat, GlobalRow, &numEntries, &inds, &vals);

	vecIndices.assign(inds, inds + numEntries);
    vecValues.assign(vals, vals + numEntries);
 	NumEntries = numEntries;

	MatRestoreRow(_mat, GlobalRow, &numEntries, &inds, &vals);
 }

 /* Global number of rows */
 int PetscSparseMatrix::globalRows() const
 {
	int num_rows;
	int num_cols;
	MatGetSize(_mat, &num_rows, &num_cols);
	return num_rows;
 }

 /* Global number of columns */
 int PetscSparseMatrix::globalCols() const
 {
 	int num_rows;
	int num_cols;
	MatGetSize(_mat, &num_rows, &num_cols);
	return num_cols;
 }

 /* Number of unique columns */
 // int PetscSparseMatrix::uniqueCols() const;

 /* Index of first global row owned by a process */
 int PetscSparseMatrix::ilower() const
 {
 	int firstrow;
 	int lastrow;
 	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
 	return firstrow; 
 }

 /* Index of last global row owned by a process */
int PetscSparseMatrix::iupper() const
{
 	int firstrow;
 	int lastrow;
 	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
 	return lastrow; 
 } 

/* get local row number from global row number if owned
     by process, return -1 otherwise */
int PetscSparseMatrix::rowMapLID( int const GID ) const
{
	int firstrow;
 	int lastrow;
 	MatGetOwnershipRange(_mat, &firstrow, &lastrow);

 	if (firstrow <= GID && GID < lastrow){
 		return GID - firstrow;
 	} else {
 		return -1;
 	}
}

 /* Local number of rows */
int PetscSparseMatrix::myRows() const
{
	int firstrow;
	int lastrow;
	MatGetOwnershipRange(_mat, &firstrow, &lastrow);
	return lastrow - firstrow; 
}

 /* Local number of columns */
 int PetscSparseMatrix::myCols() const
{
	int firstcol;
	int lastcol;
	MatGetOwnershipRangeColumn(_mat, &firstcol, &lastcol);
	return lastcol - firstcol; 
}

 /* Infinity norm */
 double PetscSparseMatrix::normInf() const
 {
	double normInf;
	MatNorm(_mat, NORM_INFINITY, &normInf);
	return normInf;
 }

 /* 1-norm */
 double PetscSparseMatrix::norm1() const
 {
	double norm1;
	MatNorm(_mat, NORM_1, &norm1);
	return norm1;
 }

 /* Frobenius norm */
 double PetscSparseMatrix::normFrobenius() const
 {
	double normFrob;
	MatNorm(_mat, NORM_FROBENIUS, &normFrob);
	return normFrob;
 }

 /* Check if the matrix is assembled */
 bool PetscSparseMatrix::isAssembled() const
 {
 	return assembled;
 }

  /* Print to terminal */
  void PetscSparseMatrix::print() const 
  {
  	MatView(_mat, PETSC_VIEWER_STDOUT_WORLD);
  }

  const Mat* PetscSparseMatrix::getPointer() const
  {
  	return &(_mat);
  }

	Mat* PetscSparseMatrix::getPointer_()
	{
		return &(_mat);
	}

  Mat PetscSparseMatrix::getMat()
  {
  	return _mat;
  }

	Mat PetscSparseMatrix::getConstMat() const
	{
		return _mat;
	}

  void PetscSparseMatrix::createWithLocalSize(int const localSize, int const maxEntriesPerRow, MPI_Comm const & comm)
	{

		// set up matrix
		MatCreate(comm, &_mat);
		MatSetType(_mat, MATMPIAIJ);
		MatSetSizes(_mat, localSize, localSize, PETSC_DETERMINE, PETSC_DETERMINE);
		MatMPIAIJSetPreallocation(_mat, maxEntriesPerRow, NULL, maxEntriesPerRow, NULL);
		MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		MatSetUp(_mat);

		// assemble
		MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

		assembled = true;

		printf("Hannah to do: createWithLocalSize()\n");
		// something wrong with number of columns

	}

  void PetscSparseMatrix::createWithGlobalSize(int const globalSize, int const maxEntriesPerRow, MPI_Comm const & comm)
	{
		// set up matrix
		MatCreate(comm, &_mat);
		MatSetType(_mat, MATMPIAIJ);
		MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, globalSize, globalSize);
		MatMPIAIJSetPreallocation(_mat, maxEntriesPerRow, NULL, maxEntriesPerRow, NULL);
		MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		MatSetUp(_mat);

		// assemble
		MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

		assembled = true;
	}

  void PetscSparseMatrix::createWithLocalSize(int const localRows, int const localCols, 
                           int const maxEntriesPerRow, MPI_Comm const & comm)
	{
		// set up matrix
		MatCreate(comm, &_mat);
		MatSetType(_mat, MATMPIAIJ);
		MatSetSizes(_mat, localRows, localCols, PETSC_DETERMINE, PETSC_DETERMINE);
		MatMPIAIJSetPreallocation(_mat, maxEntriesPerRow, NULL, maxEntriesPerRow, NULL);
		MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		MatSetUp(_mat);

		// assemble
		MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

		assembled = true;
	}

  void PetscSparseMatrix::createWithGlobalSize(int const globalRows, int const globalCols,
                            int const maxEntriesPerRow, MPI_Comm const & comm)
	{
		// set up matrix
		MatCreate(comm, &_mat);
		MatSetType(_mat, MATMPIAIJ);
		MatSetSizes(_mat, PETSC_DECIDE, PETSC_DECIDE, globalRows, globalCols);
		MatMPIAIJSetPreallocation(_mat, maxEntriesPerRow, NULL, maxEntriesPerRow, NULL);
		MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
		MatSetUp(_mat);

		// assemble
		MatAssemblyBegin(_mat, MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(_mat, MAT_FINAL_ASSEMBLY);

		assembled = true;
	}

	void PetscSparseMatrix::multiply(PetscSparseMatrix const & src, PetscSparseMatrix & dst) const
	{
		MatMatMult(_mat, src.getConstMat(), MAT_INITIAL_MATRIX, PETSC_DEFAULT, dst.getPointer_());
		// dst = C;
	}

	void PetscSparseMatrix::write(std::string const & filename) const
	{
		char filename_char[filename.length() + 1]; 
    strcpy(filename_char, filename.c_str()); 

		// set up PETSc viewer
		PetscViewer viewer;
		PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
		PetscViewerSetType(viewer, PETSCVIEWERBINARY);
		PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
		PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
		PetscViewerFileSetName(viewer, filename_char);
		MatView(_mat, viewer);

		printf("Hannah to do: write()\n");
	}

