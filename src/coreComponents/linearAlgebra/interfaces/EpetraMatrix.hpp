/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EpetraMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/EpetraVector.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

class Epetra_Map;
class Epetra_FECrsGraph;
class Epetra_FECrsMatrix;

namespace geosx
{

/**
 * \class EpetraMatrix
 * \brief This class creates and provides basic support for the Epetra_CrsMatrix
 *        matrix object type used in Trilinos.
 */
class EpetraMatrix : public MatrixBase<EpetraMatrix, EpetraVector>
{
public:

  using Base = MatrixBase<EpetraMatrix, EpetraVector>;

  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */
  EpetraMatrix();

  /**
   * @brief Copy constructor.
   *
   * Create new matrix from matrix <tt>src</tt>.
   */
  EpetraMatrix( EpetraMatrix const & src );

  /**
   * @brief Virtual destructor.
   */
  virtual ~EpetraMatrix() override;

  ///@}

  using Base::comm;
  using Base::isOpen;

  void createWithLocalSize( localIndex const localSize,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm ) final;

  void createWithGlobalSize( globalIndex const globalSize,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm ) final;

  void createWithLocalSize( localIndex const localRows,
                            localIndex const localCols,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm ) final;

  void createWithGlobalSize( globalIndex const globalRows,
                             globalIndex const globalCols,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm ) final;

  void open() final;

  void close() final;

  void set( real64 const value ) final;

  void zero() final;

  ///@}



  //! @name Linear Algebra Methods
  //@{
  /**
   * @brief Matrix/Vector multiplication.
   *
   * Compute <tt>Ax = b<tt>.
   *
   * @param src Input vector (x).
   * @param dst Output vector (b).
   *
   */
  virtual void multiply( EpetraVector const & src,
                         EpetraVector & dst ) const override;


  /**
   * @brief Matrix/Matrix multiplication.
   *
   * Compute <tt>this * B = C<tt>.
   *
   * @param src Input matrix (B).
   * @param dst Output matrix (C).
   * @param closeResult whether to close @p dst for additional entries.
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  void multiply( EpetraMatrix const & src,
                 EpetraMatrix & dst,
                 bool const closeResult = true ) const;

  /**
   * @brief Matrix/Matrix transpose multiplication.
   *
   * Compute <tt>this^T * B = C<tt>.
   *
   * @param src Input matrix (B).
   * @param dst Output matrix (C).
   * @param closeResult whether to close @p dst for additional entries.
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  void leftMultiplyTranspose( EpetraMatrix const & src,
                              EpetraMatrix & dst,
                              bool const closeResult = true ) const;

  /**
     * @brief Matrix/Matrix transpose multiplication.
     *
     * Compute <tt>B * this^T = C<tt>.
     *
     * @param src Input matrix (B).
     * @param dst Output matrix (C).
     * @param closeResult whether to close @p dst for additional entries.
     *
     * Note that the output matrix C should have the same
     * row-map as this.  If close() has already been called
     * on C, then C's sparsity pattern must already contain
     * the nonzero entries produced by the product this*B.
     */
    void rightMultiplyTranspose( EpetraMatrix const & src,
                                 EpetraMatrix & dst,
                                 bool const closeResult = true ) const;

  /**
   * @brief Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * @param alpha Scalar factor for added matvec product.
   * @param x Input vector.
   * @param beta Scalar factor for right hand side.
   * @param y Output vector.
   * @param useTranspose Boolean, set to true to use <tt>A^T</tt>.
   *
   */
  void gemv( real64 const alpha,
             EpetraVector const & x,
             real64 const beta,
             EpetraVector & y,
             bool useTranspose = false );

  /**
   * @brief Multiply all elements by scalingFactor.
   *
   * @param scalingFactor Scaling factor.
   *
   */
  void scale( real64 const scalingFactor );

  /**
   * @brief Pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   *
   * @param vec Vector to pre-multiply with.
   *
   */
  void leftScale( EpetraVector const & vec );

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vec.
   *
   * @param vec Vector to post-multiply with.
   *
   */
  void rightScale( EpetraVector const & vec );

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vecRight
   * and pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   *
   * @param vec vecLeft to pre-multiply with.
   * @param vec vecRight to post-multiply with.
   *
   */
  void leftRightScale( EpetraVector const & vecLeft,
                       EpetraVector const & vecRight );

  /**
   * @brief Clear a row, and optionally set diagonal element to <tt>diagValue</tt>.
   *
   * @param row globalIndex of the row to be cleared.
   * @param diagValue (Optional) set diagonal element to desired value.
   *
   */
  void clearRow( globalIndex const row,
                 real64 const diagValue = 0 );

  //@}

  //! @name Accessors Methods
  //@{

  /**
   * @brief Returns a copy of the data in row <tt>globalRow</tt>.
   * Note that the input arrays will be resized internally to fit the number of entries.
   */
  void getRowCopy( globalIndex globalRow,
                   array1d< globalIndex > & colIndices,
                   array1d< real64 > & values ) const;

  /**
   * @brief get diagonal element value on a given row
   * @param globalRow global row index
   * @return value of diagonal element on the row
   */
  real64 getDiagValue( globalIndex globalRow ) const;

  /**
   * @brief Returns a pointer to the underlying matrix.
   */
  Epetra_FECrsMatrix * unwrappedPointer() const;

  /**
   * @brief Returns the number of global rows.
   */
  globalIndex globalRows() const;

  /**
   * @brief Returns the number of global columns.
   */
  globalIndex globalCols() const;

  /**
 * @brief Return the local number of columns on each processor
 */
  localIndex localRows() const;

  /**
   * @brief Return the local number of columns on each processor
   */
  localIndex localCols() const;

  /**
   * @brief Returns the index of the first global row owned by that processor.
   */
  globalIndex ilower() const;

  /**
   * @brief Returns the next index after last global row owned by that processor.
   *
   * @note The intention is for [ilower; iupper) to be used as a half-open index range
   */
  globalIndex iupper() const;

  /**
   * @brief Returns the number of nonzeros in the local portion of the matrix
   */
  localIndex localNonzeros() const;

  /**
   * @brief Returns the total number of nonzeros in the matrix
   */
  globalIndex globalNonzeros() const;

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
  void print( std::ostream & os = std::cout ) const;

  /**
   * @brief Write the matrix to filename in a matlab-compatible format.
   *
   * Within octave / matlab:
   * >> load filename
   * >> M = spconvert(filename_root)
   */
  void write( string const & filename,
              bool const mtxFormat = true ) const;

  /**
   * @brief Map a global row index to local row index
   */
  localIndex getLocalRowID( globalIndex const index ) const;

  /**
   * @brief Map a local row index to global row index
   */
  globalIndex getGlobalRowID( localIndex const index ) const;

  //@}

private:
  /**
     * @brief Perform a matrix matrix product with Parallel Matrix
     */
    void multiply( bool const transA,
                   EpetraMatrix const & B,
                   bool const transB,
                   EpetraMatrix & C,
                   bool const closeResult ) const;

  /**
   * Boolean value, true if the matrix had been finalized, false if not.
   */
  bool m_assembled = false;

  /**
   * Pointer to the underlying Epetra_CrsMatrix.
   */
  std::unique_ptr< Epetra_FECrsMatrix > m_matrix;

  /*
   * Map representing the parallel partitioning of a source vector (x in y=Ax)
   */
  std::unique_ptr< Epetra_Map > m_src_map;

  /*
   * Map representing the parallel partitioning of a destination vector (y in y=Ax)
   */
  std::unique_ptr< Epetra_Map > m_dst_map;
};

/**
 * @brief Stream insertion operator for EpetraMatrix
 * @param os the output stream
 * @param matrix the matrix to be printed
 * @return reference to the output stream
 */
std::ostream & operator<<( std::ostream & os,
                           EpetraMatrix const & matrix );

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_*/
