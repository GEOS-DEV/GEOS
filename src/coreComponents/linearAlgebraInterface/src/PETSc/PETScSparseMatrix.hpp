#include "PETScVector.hpp"

#include <petscvec.h>
#include <petscmat.h>

class PETScSparseMatrix
{

public:

	/* Create an empty distributed matrix */
	PETScSparseMatrix();

	/* Copy constructor */
	PETScSparseMatrix(PETScSparseMatrix const &in_matrix);

	/* Virtual destructor */
	// virtual ~PETScSparseMatrix() = default;

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
   // void create( MPI_Comm const comm,
   //              int const m_nRowGlobal,
   //              int const m_nColGlobal,
   //              std::vector<int> const nMaxEntriesPerRow );

  /* Create from a PETSc matrix */
  void create( PETScSparseMatrix &matrix );

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
  void multiply( PETScVector &src,
                 PETScVector &dst ) const;

  /* Compute residual r = Ax - b */
  void residual( PETScVector  &x,
                 PETScVector  &b,
                 PETScVector &res ) const;

  /* Compute "gaxpy" alpha*A*x + beta*b */
  void gaxpy( double alpha,
              PETScVector  &x,
              double beta,
              PETScVector  &b,
              PETScVector &res,
              bool useTranspose=false );

  /* Multiply all elements by scalingFactor */
  void scale( double scalingFactor );

  /* Left multiply with diagonal matrix with values in vec */
  void leftScale( PETScVector  &vec );

  /* Right multiply with diagonal matrix with values in vec */
  void rightScale( PETScVector  &vec );

  /* Left multiply with diagonal matrix with values vecLeft and right multiply
     with diagonal matrix with values vecRight */
  void leftRightScale( PETScVector &vecLeft,
                       PETScVector &vecRight );

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
               int* Indices );

  /*
   * Get global row myRow
   * - numEntries: number of nonzeros 
   * - vecValues: vector of values
   * - vecIndices: vector of column indices */
  // void getRow( int GlobalRow,
  //              int &NumEntries,
  //              std::vector<double> &vecValues,
  //              std::vector<int> &vecIndices );

  /*
   * Get local row myRow
   * - numEntries: number of nonzeros 
   * - Values: array of values
   * - Indices: array of column indices */
  void getLocalRow( int myRow,
                    int & NumEntries,
                    double * & Values,
                    int * & Indices );

  /*
   * Get local row myRow
   * - numEntries: number of nonzeros 
   * - vecValues: vector of values
   * - vecIndices: vector of column indices */
  // void getLocalRow( int myRow,
  //                   int &NumEntries,
  //                   std::vector<double> &vecValues,
  //                   std::vector<int> &vecIndices );

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

  // Epetra_Map const & RowMap() const;
   
  // Epetra_Map const & ColMap() const;

  // Epetra_Map const & DomainMap() const;

  // int rowMapLID( int GID ) const;

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

  /* Get underlying PETSc object */
  Mat getMat();

protected:

	bool assembled = false;

  // underlying PETSc object
	Mat _mat;

};