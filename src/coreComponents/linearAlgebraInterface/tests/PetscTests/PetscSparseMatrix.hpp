#include "PetscVector.hpp"

#include <petscvec.h>
#include <petscmat.h>

class PetscSparseMatrix
{

public:

	/* Create an empty distributed matrix */
	PetscSparseMatrix();

	/* Copy constructor */
	PetscSparseMatrix(PetscSparseMatrix const &in_matrix);

	/* Virtual destructor */
	// virtual ~PetscSparseMatrix() = default;

	/* Create a square matrix from number of rows */
	void create( MPI_Comm const comm,
				       int const m_nRowGlobal,
				       int const nMaxEntriesPerRow );

  /* Create rectangular matrix */
  void create( MPI_Comm const comm,
               int const m_nRowGlobal,
               int const m_nColGlobal,
               int const nMaxEntriesPerRow);

   /* Create square matrix */
   void create( MPI_Comm const comm,
                int const m_nRowGlobal,
                std::vector<int> const nMaxEntriesPerRow );
  
   /* Create matrix */
   void create( MPI_Comm const comm,
                int const m_nRowGlobal,
                int const m_nColGlobal,
                std::vector<int> const nMaxEntriesPerRow );

  /* Create from a PETSc matrix */
  void create( PetscSparseMatrix &matrix );

  /* Set all elements to zero, keep partition and sparsity pattern */
  void zero();

  /* Empty function for Trilinos implementation. Is required when the HYPRE library is used. */
  void open();

  /* Assemble the matrix */
  void close();
  
  /* Add values into iRow at columns cols. 
   * iRow: global row index
   * nCols: number of columns to modify
   * values: values to add to prescribed locations
   * cols: global column indices in which to add the values */
  void add( int const iRow,
            int const nCols,
            double const *values,
            int const *cols );

  /* add values to sparse matrix located at indices */
  // void add( array1d<int> const rowIndices,
  //           array1d<int> const colIndices,
  //           array2d<double> const values);

  /* Add value to position (iRow, iCol) */
  void add( int const iRow,
            int const iCol,
            double const value );

  /* Set values into iRow at columns cols. 
   * iRow: global row index
   * nCols: number of columns to modify
   * values: values to add to prescribed locations
   * cols: global column indices in which to add the values */
  void set( int const iRow,
            int const nCols,
            double const *values,
            int const *cols );

  /* Set the value of (iRow, iCol) to value */
  void set( int const iRow,
            int const iCol,
            double const value );

  /* Insert values into iRow at columns cols. 
   * iRow: global row index
   * nCols: number of columns to modify
   * values: values to add to prescribed locations
   * cols: global column indices in which to add the values */
  void insert( int const iRow,
               int const nCols,
               double const *values,
               int const *cols );

  /* Matrix-vector multiplication, stored in dst */
  void multiply( PetscVector const &src,
                 PetscVector &dst ) const;

  /* Compute residual r = Ax - b */
  void residual( PetscVector  const &x,
                 PetscVector  const &b,
                 PetscVector &res ) const;

 /* Compute "gemv" y = alpha*A*x + beta*y */
 void gemv( double const alpha,
            PetscVector  const &x,
            double const beta,
            PetscVector  &y,
            bool useTranspose);

  /* Multiply all elements by scalingFactor */
  void scale( double const scalingFactor );

  /* Left multiply with diagonal matrix with values in vec */
  void leftScale( PetscVector const &vec );

  /* Right multiply with diagonal matrix with values in vec */
  void rightScale( PetscVector const &vec );

  /* Left multiply with diagonal matrix with values vecLeft and right multiply
     with diagonal matrix with values vecRight */
  void leftRightScale( PetscVector const &vecLeft,
                       PetscVector const &vecRight );

  /* Clear row and multiply diagonal by factor */
  void clearRow( int const row,
                 double const factor );

  /*
   * Get global row myRow
   * - numEntries: number of nonzeros 
   * - Values: array of values
   * - Indices: array of column indices */
  void getRow( int GlobalRow,
               int &NumEntries,
               double* Values,
               int* Indices ) const;

  /*
   * Get global row myRow
   * - numEntries: number of nonzeros 
   * - vecValues: vector of values
   * - vecIndices: vector of column indices */
  void getRow( int GlobalRow,
               int &NumEntries,
               std::vector<double> &vecValues,
               std::vector<int> &vecIndices ) const;

  /*
   * Get local row myRow
   * - numEntries: number of nonzeros 
   * - Values: array of values
   * - Indices: array of column indices */
  void getLocalRow( int myRow,
                    int &NumEntries,
                    double *Values,
                    int *Indices ) const;

  /*
   * Get local row myRow
   * - numEntries: number of nonzeros 
   * - vecValues: vector of values
   * - vecIndices: vector of column indices */
  void getLocalRow( int myRow,
                    int &NumEntries,
                    std::vector<double> &vecValues,
                    std::vector<int> &vecIndices ) const;

  /* Global number of rows */
  int globalRows() const;

  /* Global number of columns */
  int globalCols() const;

  /* Number of unique columns */
  int uniqueCols() const;

  /* Index of first global row owned by a process */
  int ilower() const;

  /* Index of last global row owned by a process */
  int iupper() const;

  /* Local number of rows */
  int myRows() const;

  /* Local number of columns */
  int myCols() const;

  /* get local row number from global row number if owned
     by process, return -1 otherwise */
  int rowMapLID( int const GID ) const;

  /* Infinity norm */
  double normInf() const;

  /* 1-norm */
  double norm1() const;

  /* Frobenius norm */
  double normFrobenius() const;

  /* Check if the matrix is assembled */
  bool isAssembled() const;

  /* Print to terminal */
  void print() const;

  /* Pointer to the underlying PETSc matrix */
  const Mat* getPointer() const;

  Mat* getPointer_();

  /* Get underlying PETSc object */
  Mat getMat();

  Mat getConstMat() const;

  /* create square matrix with local number of rows and columns */
  void createWithLocalSize(int const localSize, int const maxEntriesPerRow, MPI_Comm const & comm);

   /* create square matrix with global number of rows and columns */
   void createWithGlobalSize(int const globalSize, int const maxEntriesPerRow, MPI_Comm const & comm);

  /* create rectangular matrix with local number of rows and columns */
  void createWithLocalSize(int const localRows, int const localCols, 
                           int const maxEntriesPerRow, MPI_Comm const & comm);

   /* create rectangular matrix with global number of rows and columns */
   void createWithGlobalSize(int const globalRows, int const globalCols,
                             int const maxEntriesPerRow, MPI_Comm const & comm);

   /* matrix-matrix multiplicaiton */
   void multiply(PetscSparseMatrix const & src, PetscSparseMatrix & dst) const;

   /* write matrix to a file */
   void write(std::string const & filename) const;

protected:

	bool assembled = false;

  // underlying PETSc object
	Mat _mat;

};