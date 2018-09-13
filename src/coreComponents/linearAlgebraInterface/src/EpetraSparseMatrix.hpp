/**
 * @file EpetraSparseMatrix.hpp
 */

#ifndef EPETRASPARSEMATRIX_HPP_
#define EPETRASPARSEMATRIX_HPP_

#include "EpetraVector.hpp"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>
#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * \class EpetraSparseMatrix
 * \brief This class creates and provides basic support for the Epetra_CrsMatrix
 *        matrix object type used in Trilinos.
 */

class EpetraSparseMatrix
{
public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */
  EpetraSparseMatrix();

  /**
   * @brief Copy constructor.
   */
  EpetraSparseMatrix( EpetraSparseMatrix const &in_matrix );

  /**
   * @brief Virtual destructor.
   */
  virtual ~EpetraSparseMatrix() = default;
  //@}

  //! @name Create/Finalize/Reinitialize Methods
  //@{
  /**
   * @brief Create a square matrix from number of rows.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Maximum number of entries per row.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   *
   */
  void create( MPI_Comm const comm,
               globalIndex const m_nRowGlobal,
               integer const nMaxEntriesPerRow );

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Global number of rows.
   * \param m_nColGlobal Global number of columns.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   */
  void create( MPI_Comm const comm,
               globalIndex const m_nRowGlobal,
               globalIndex const m_nColGlobal,
               integer const nMaxEntriesPerRow = 0 );

  // TODO see if we need a vector of nnz
  //  /**
  //   * @brief Create a square matrix from number of unknowns.
  //   *
  //   * \param comm MPI communicator.
  //   * \param m_nRowGlobal Global number of unknowns.
  //   * \param nMaxEntriesPerRow Vector of maximum number of entries per row.
  //   */
  //  void create( MPI_Comm const comm,
  //               integer const m_nRowGlobal,
  //               std::vector<integer> const nMaxEntriesPerRow );
  //
  //  /**
  //   * @brief Create a square matrix from number of unknowns.
  //   *
  //   * \param comm MPI communicator.
  //   * \param m_nRowGlobal Global number of rows.
  //   * \param m_nColGlobal Global number of columns.
  //   * \param nMaxEntriesPerRow Vector of maximum number of entries per row.
  //   */
  //  void create( MPI_Comm const comm,
  //               globalIndex const m_nRowGlobal,
  //               globalIndex const m_nColGlobal,
  //               std::vector<integer> const nMaxEntriesPerRow );

  /**
   * @brief Create a square matrix from Epetra_Map.
   *
   * Prepare the matrix to be filled. Takes as inputs an existing Epetra_Map and
   * a hint on the maximum number of entries per row. This method is meant to be called
   * after the matrix has been declared from the empty constructor.
   *
   * \param input_map Epetra_Map.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   *
   */
  void create( Epetra_Map const &input_map,
               integer const nMaxEntriesPerRow );

  /**
   * @brief Create a rectangular matrix from two existing Epetra_Map, row and column maps.
   *
   * \param row_map Epetra_Map for rows.
   * \param col_map Epetra_Map for columns.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   */
  void create( Epetra_Map const &row_map,
               Epetra_Map const &col_map,
               integer const nMaxEntriesPerRow = 0 );

  /**
   * @brief Create a matrix from an existing Epetra_CrsGraph.
   *
   * TODO change that to whatever format the sparsity pattern will be.
   *
   * \param Epetra_CrsGraph existing graph.
   */
  void create( Epetra_CrsGraph &graph );

  /**
   * @brief Create a matrix from an existing Epetra_CrsMatrix.
   *
   * \param Epetra_CrsMatrix existing matrix.
   */
  void create( Epetra_CrsMatrix &matrix );

  /**
   * @brief Reinitialize the matrix.
   *
   * Keeps the parallel partitioning and the sparsity pattern but sets all elements to zero.
   *
   */
  void zero();

  /**
   * @brief Empty function for Trilinos implementation. Is required when the HYPRE library is used.
   *
   */
  void open();

  /**
   * @brief Assemble and compress the matrix.
   *
   * Compresses the matrix to CSR format with contiguous memory on each processor. Prevents from
   * adding new entries in the sparsity pattern but allows for modification of existing entries.
   *
   */
  void close();
  //@}

  //! @name Insertion/Replace/SumInto Methods
  //@{
  /**
   * @brief Add to row of elements.
   *
   * Adds the values <tt>values</tt> to row <tt>iRow</tt>, at locations specified
   * by <tt>cols</tt>. <tt>nCols</tt> is the number of entries (columns) in the row, and the
   * size of both <tt>values</tt> and <tt>cols</tt>.
   *
   * \param iRow Global row index.
   * \param nCols Number of columns to modify.
   * \param values Values to add to prescribed locations.
   * \param cols Global column indices in which to add the values.
   */
  void add( globalIndex const iRow,
            integer const nCols,
            real64 const *values,
            globalIndex const *cols );

  /**
   * @brief Add to one element.
   *
   * Adds the value <tt>value</tt> to location (<tt>iRow</tt>,<tt>iCol</tt>).
   *
   * \param iRow Global row index.
   * \param iCol Global column index.
   * \param value Value to add to prescribed locations.
   *
   */
  void add( globalIndex const iRow,
            globalIndex const iCol,
            real64 const value );

  /**
   * @brief Set row of elements.
   *
   * Sets the values <tt>values</tt> of row <tt>iRow</tt>, at locations specified
   * by <tt>cols</tt>. <tt>nCols</tt> is the number of entries (columns) in the row, and the
   * size of both <tt>values</tt> and <tt>cols</tt>.
   *
   * \param iRow Global row index.
   * \param nCols Number of columns to modify.
   * \param values Values to set in prescribed locations.
   * \param cols Global column indices in which to set the values.
   *
   */
  void set( globalIndex const iRow,
            integer const nCols,
            real64 const *values,
            globalIndex const *cols );

  /**
   * @brief Set one element.
   *
   * Sets the value of the location (<tt>iRow</tt>,<tt>iCol</tt>) to <tt>value</tt>.
   *
   * \param iRow Global row index.
   * \param iCol Global column index.
   * \param value Value to set at prescribed locations.
   *
   */
  void set( globalIndex const iRow,
            globalIndex const iCol,
            real64 const value );

  /**
   * @brief Insert to row of elements.
   *
   * Inserts the values <tt>values</tt> to row <tt>iRow</tt>, at locations specified
   * by <tt>cols</tt>. <tt>nCols</tt> is the number of entries (columns) in the row, and the
   * size of both <tt>values</tt> and <tt>cols</tt>.
   *
   * TODO remove the possibility to dynamically construct the sparsity pattern.
   *
   * \param iRow Global row index.
   * \param nCols Number of columns to modify.
   * \param values Values to add to prescribed locations.
   * \param cols Global column indices in which to add the values.
   */
  void insert( globalIndex const iRow,
               integer const nCols,
               real64 const *values,
               globalIndex const *cols );

  //@}

  //! @name Linear Algebra Methods
  //@{
  /**
   * @brief Matrix/Vector multiplication.
   */
  void multiply( EpetraVector const &src,
                 EpetraVector &dst ) const;

  /**
   * @brief Compute residual r = Ax - b.
   */
  void residual( EpetraVector const &x,
                 EpetraVector const &b,
                 EpetraVector &res ) const;

  /**
   * @brief Compute "gaxpy" res = alpha*A*x + beta*b.
   */
  void gaxpy( real64 alpha,
              EpetraVector const &x,
              real64 beta,
              EpetraVector const &b,
              EpetraVector &res,
              bool useTranspose=false );

  /**
   * @brief Multiply all elements by scalingFactor.
   */
  void scale( real64 scalingFactor );

  /**
   * @brief Pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   */
  void leftScale( EpetraVector const &vec );

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vec.
   */
  void rightScale( EpetraVector const &vec );

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vecRight
   * and pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   */
  void leftRightScale( EpetraVector const &vecLeft,
                       EpetraVector const &vecRight );

  /**
   * @brief Clear a row and multiplies the diagonal term by <tt>factor</tt>.
   */
  void clearRow( globalIndex const row,
                 real64 const factor );

  //@}

  //! @name Accessors Methods
  //@{

  /**
   * @brief Returns the row <tt>GlobalRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getRow( globalIndex GlobalRow,
               integer &NumEntries,
               real64* Values,
               globalIndex* Indices );

  /**
   * @brief Returns the row <tt>GlobalRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getRow( globalIndex GlobalRow,
               integer &NumEntries,
               std::vector<real64> &vecValues,
               std::vector<globalIndex> &vecIndices );

  /**
   * @brief Returns the row <tt>localRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getLocalRow( localIndex myRow,
                    integer &NumEntries,
                    real64* Values,
                    localIndex* Indices );

  /**
   * @brief Returns the row <tt>localRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getLocalRow( localIndex myRow,
                    integer &NumEntries,
                    std::vector<real64> &vecValues,
                    std::vector<localIndex> &vecIndices );

  /**
   * @brief Returns a pointer to the underlying matrix.
   */
  Epetra_CrsMatrix* getPointer() const;

  /**
   * @brief Returns the number of global rows.
   */
  globalIndex globalRows() const;

  /**
   * @brief Returns the number of global columns.
   */
  globalIndex globalCols() const;

  /**
   * @brief Returns the number of unique columns (can be used to check if matrix is square).
   */
  globalIndex uniqueCols() const;

  /**
   * @brief Returns the index of the first global row owned by that processor.
   */
  globalIndex ilower() const;

  /**
   * @brief Returns the index of the last global row owned by that processor.
   */
  globalIndex iupper() const;

  /**
   * @brief Returns the number of local rows.
   */
  int myRows() const;

  /**
   * @brief Returns the number of local columns.
   */
  int myCols() const;

  /**
   * @brief Returns the row map.
   */
  Epetra_Map const & RowMap() const;

  /**
   * @brief Returns the (usually overlapping) column map.
   */
  Epetra_Map const & ColMap() const;

  /**
   * @brief Returns the (1-to-1) domain column map.
   */
  Epetra_Map const & DomainMap() const;

  /**
   * @brief Wrapper for LID function. Returns the local map of the corresponding global index.
   * Returns -1 if the global row is not owned by the processor.
   */
  localIndex rowMapLID( globalIndex GID ) const;

  /**
   * @brief Returns the infinity norm of the matrix.
   */
  real64 normInf() const;

  /**
   * @brief Returns the one norm of the matrix.
   */
  real64 norm1() const;

  /**
   * @brief Returns the Frobenius norm of the matrix.
   */
  real64 normFrobenius() const;

  /**
   * @brief Returns true is the matrix has been assembled, false if not.
   */
  bool isAssembled() const;
  //@}

  //! @name I/O Methods
  //@{
  /**
   * @brief Print the matrix in Trilinos format to the terminal.
   */
  void print() const;
  //@}

private:

  /**
   * @brief Boolean value, true if the matrix had been finalized, false if not.
   */
  bool assembled = false;

  /**
   * @brief Pointer to the underlying Epetra_CrsMatrix.
   */
  std::unique_ptr<Epetra_CrsMatrix> m_matrix = nullptr;

};

} // namespace geosx

#endif /* EpetraSPARSEMATRIX_HPP_ */
