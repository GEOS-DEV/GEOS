/*
 * EpetraSparseMatrix.h
 *
 *  Created on: Jul 20, 2018
 *      Author: Matthias
 */

#ifndef EPETRASPARSEMATRIX_HPP_
#define EPETRASPARSEMATRIX_HPP_

#include <memory>
#include "EpetraVector.hpp"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
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
  void create(const MPI_Comm    comm,
              const globalIndex m_nRowGlobal,
              const integer     nMaxEntriesPerRow);

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Global number of rows.
   * \param m_nColGlobal Global number of columns.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   */
  void create(const MPI_Comm    comm,
              const globalIndex m_nRowGlobal,
              const globalIndex m_nColGlobal,
              const integer     nMaxEntriesPerRow = 0);

  /**
   * @brief Create a square matrix from number of unknowns.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Global number of unknowns.
   * \param nMaxEntriesPerRow Vector of maximum number of entries per row.
   */
  void create(const MPI_Comm             comm,
              const int                  m_nRowGlobal,
              const std::vector<integer> nMaxEntriesPerRow);

  /**
   * @brief Create a square matrix from number of unknowns.
   *
   * \param comm MPI communicator.
   * \param m_nRowGlobal Global number of rows.
   * \param m_nColGlobal Global number of columns.
   * \param nMaxEntriesPerRow Vector of maximum number of entries per row.
   */
  void create(const MPI_Comm             comm,
              const globalIndex          m_nRowGlobal,
              const globalIndex 		 m_nColGlobal,
              const std::vector<integer> nMaxEntriesPerRow);

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
  void create(const Epetra_Map &input_map,
              const integer     nMaxEntriesPerRow);

  /**
   * @brief Create a rectangular matrix from two existing Epetra_Map, row and column maps.
   *
   * \param row_map Epetra_Map for rows.
   * \param col_map Epetra_Map for columns.
   * \param nMaxEntriesPerRow Maximum number of entries per row.
   */
  void create(const Epetra_Map &row_map,
              const Epetra_Map &col_map,
              const integer     nMaxEntriesPerRow = 0);

  /**
   * @brief Create a matrix from an existing Epetra_FECrsMatrix.
   *
   * \param Epetra_CrsMatrix existing matrix.
   */
  void create(Epetra_CrsMatrix &matrix);

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
  void add(const globalIndex  iRow,
           const integer      nCols,
           const real64      *values,
           const globalIndex *cols);

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
  void add(const globalIndex    iRow,
           const globalIndex    iCol,
           const real64 value);

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
  void set(const globalIndex     iRow,
           const integer         nCols,
           const real64         *values,
           const globalIndex    *cols);

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
  void set(const globalIndex    iRow,
           const globalIndex    iCol,
           const real64 value);

  //@}

  //! @name Linear Algebra Methods
  //@{
  /**
   * @brief Matrix/Vector multiplication.
   */
  void apply(      EpetraVector &dst,
             const EpetraVector &src);
  //@}

  //! @name Accessors Methods
  //@{
  /**
   * @brief Returns the row <tt>GlobalRow</tt>. The number of non zeros in the row is <tt>NumEntries</tt>
   * , the values are sent to <tt>vecValues</tt> and the column indices in <tt>vecIndices</tt>.
   */
  void getRow(int                  GlobalRow,
              int                 &NumEntries,
              std::vector<real64> &vecValues,
              std::vector<int>    &vecIndices);


  /**
   * @brief Returns the number of global rows.
   */
  globalIndex globalRows();

  /**
   * @brief Returns the number of global columns.
   */
  globalIndex globalCols();

  /**
   * @brief Returns the number of unique columns (can be used to check if matrix is square).
   */
  globalIndex uniqueCols();

  /**
   * @brief Returns the number of local rows.
   */
  int myRows();

  /**
   * @brief Returns the number of local columns.
   */
  int myCols();

  /**
   * @brief Returns the row map.
   */
  const Epetra_Map RowMap();

  /**
   * @brief Returns the (usually overlapping) column map.
   */
  const Epetra_Map ColMap();

  /**
   * @brief Returns the (1-to-1) domain column map.
   */
  const Epetra_Map DomainMap();

  /**
   * @brief Returns true is the matrix has been assembled, false if not.
   */
  bool isAssembled();
  //@}

  //! @name I/O Methods
  //@{
  /**
   * @brief Print the matrix in Trilinos format to the terminal.
   */
  void print();
  //@}


private:

  /**
   * @brief Boolean value, true if the matrix had been finalized, false if not.
   */
  bool assembled=false;

  /**
   * @brief Pointer to the underlying Epetra_CrsMatrix.
   */
  std::unique_ptr<Epetra_CrsMatrix> matrix = nullptr;

};

}  // namespace geosx

#endif /* EpetraSPARSEMATRIX_HPP_ */
