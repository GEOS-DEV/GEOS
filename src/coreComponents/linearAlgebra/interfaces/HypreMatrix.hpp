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
 * @file HypreMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_HYPREMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_HYPREMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/HypreVector.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

// Just a placeholder to avoid to include two HYPRE header files
//#include "_hypre_IJ_mv.h"
//#include "_hypre_parcsr_mv.h"

// IJMatrix definition
struct hypre_IJMatrix_struct;
typedef struct hypre_IJMatrix_struct *HYPRE_IJMatrix;

// ParCSRMatrix definition
struct hypre_ParCSRMatrix_struct;
typedef struct hypre_ParCSRMatrix_struct *HYPRE_ParCSRMatrix;

namespace geosx
{

/**
 * \class HypreMatrix
 * \brief This class ...
 */
class HypreMatrix : public MatrixBase<HypreMatrix, HypreVector>
{
public:

  /// @name Constructor/Destructor Methods
  ///@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */

  HypreMatrix();

  /**
   * @brief Copy constructor.
   *
   * Create new matrix from matrix <tt>src</tt>.
   */
  HypreMatrix( HypreMatrix const & src );

  /**
   * @brief Virtual destructor.
   */
  ~HypreMatrix() final;

  ///@}

  void createWithLocalSize( localIndex const localSize,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm = MPI_COMM_WORLD ) final;

  void createWithGlobalSize( globalIndex const globalSize,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm = MPI_COMM_WORLD ) final;

  void createWithLocalSize( localIndex const localRows,
                            localIndex const localCols,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm = MPI_COMM_WORLD ) final;

  void createWithGlobalSize( globalIndex const globalRows,
                             globalIndex const globalCols,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm = MPI_COMM_WORLD ) final;

  void set( real64 const value ) final;

  void reset() final;

  void zero() final;

  void open() final;

  void close() final;

  /** @name Add/Set/Insert Methods
   *
   * TRILINOS logic:
   * The add and set methods assume entries already exist in the sparsity pattern.
   * Insert methods allow for dynamic allocation, but will temporarily use
   * extra memory if one attempts to insert multiple values to the same location.
   *
   * HYPRE logic:
   * The add and set methods can be used also if the sparsity pattern has not been
   * finalized. In Hypre the insert method is an alias for set
   *
   * Caution: In Trilinos these methods are not thread-safe.  //TODO: add thread safety
   */
  ///@{

  void add( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) final;

  void set( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) final;

  void insert( globalIndex const rowIndex,
               globalIndex const colIndex,
               real64 const value ) final;

  void add( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) final;

  void set( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) final;

  void insert( globalIndex const rowIndex,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const size ) final;

  void add( globalIndex const rowIndex,
            array1d< globalIndex > const & colIndices,
            array1d< real64 > const & values ) final;

  void set( globalIndex const rowIndex,
            array1d< globalIndex > const & colIndices,
            array1d< real64 > const & values ) final;

  void insert( globalIndex const rowIndex,
               array1d< globalIndex > const & colIndices,
               array1d< real64 > const & values ) final;

  void add( array1d< globalIndex > const & rowIndices,
            array1d< globalIndex > const & colIndices,
            array2d< real64 > const & values ) final;

  void set( array1d< globalIndex > const & rowIndices,
            array1d< globalIndex > const & colIndices,
            array2d< real64 > const & values ) final;

  void insert( array1d< globalIndex > const & rowIndices,
               array1d< globalIndex > const & colIndices,
               array2d< real64 > const & values ) final;

  void add( globalIndex const * rowIndices,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const numRows,
            localIndex const numCols ) final;

  void set( globalIndex const * rowIndices,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const numRows,
            localIndex const numCols ) final;

  void insert( globalIndex const * rowIndices,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const numRows,
               localIndex const numCols ) final;

  ///@}

  void multiply( HypreVector const & src,
                 HypreVector & dst ) const final;

  void multiply( HypreMatrix const & src,
                 HypreMatrix & dst,
                 bool const closeResult = true ) const final;

  void leftMultiplyTranspose( HypreMatrix const & src,
                              HypreMatrix & dst,
                              bool const closeResult = true ) const final;

  void rightMultiplyTranspose( HypreMatrix const & src,
                               HypreMatrix & dst,
                               bool const closeResult = true ) const final;

  void gemv( real64 const alpha,
             HypreVector const & x,
             real64 const beta,
             HypreVector & y,
             bool useTranspose=false ) const final;

  void scale( real64 const scalingFactor ) final;

  void leftScale( HypreVector const & vec ) final;

  void rightScale( HypreVector const &vec ) final;

  void leftRightScale( HypreVector const &vecLeft,
                       HypreVector const &vecRight ) final;

  void clearRow( globalIndex const row,
                 real64 const diagValue = 0 ) final;

  localIndex maxRowLength() const final;

  localIndex localRowLength( localIndex localRowIndex ) const final;

  localIndex globalRowLength( globalIndex globalRowIndex ) const final;

  void getRowCopy( globalIndex globalRow,
                   array1d< globalIndex > & colIndices,
                   array1d< real64 > & values ) const final;

  real64 getDiagValue( globalIndex globalRow ) const final;

  globalIndex globalRows() const final;

  globalIndex globalCols() const final;

  localIndex localRows() const final;

  localIndex localCols() const final;

  globalIndex ilower() const final;

  globalIndex iupper() const final;

  /**
   * @brief Returns the number of nonzeros in the local portion of the matrix
   */
  localIndex localNonzeros() const final;

  /**
   * @brief Returns the infinity norm of the matrix.
   */
  real64 normInf() const final;

  /**
   * @brief Returns the one norm of the matrix.
   */
  real64 norm1() const final;

  /**
   * @brief Returns the Frobenius norm of the matrix.
   */
  real64 normFrobenius() const final;

  /**
   * @brief Returns true is the matrix has been assembled, false if not.
   */
  bool isAssembled() const;
  //@}

  /**
   * @brief Returns true is the matrix has been closed, false if not.
   */
  bool isClosed() const;
  //@}

  void print( std::ostream & os = std::cout ) const final;

  void write( string const & filename,
              bool const mtxFormat = true ) const final;

  localIndex getLocalRowID( globalIndex const index ) const final;

  globalIndex getGlobalRowID( localIndex const index ) const final;

  ///@}

  /**
   * @brief Returns a pointer to the underlying HYPRE_IJMatrix object.
   */
  HYPRE_IJMatrix const * unwrappedPointer() const;

  HYPRE_IJMatrix * unwrappedPointer();

  operator HYPRE_IJMatrix()
  {
    return (HYPRE_IJMatrix) m_ij_mat;
  }

  operator HYPRE_ParCSRMatrix()
  {
    return (HYPRE_ParCSRMatrix) m_parcsr_mat;
  }

private:

  /**
   * @brief Returns the index of the first global col owned by that processor.
   */
  globalIndex jlower() const;

  /**
   * @brief Returns the next index after last global col owned by that processor.
   *
   * @note The intention is for [jlower; jupper) to be used as a half-open index range
   */
  globalIndex jupper() const;

  /**
   * @brief Perform a matrix matrix product with Parallel Matrix
   */
  void parCSRtoIJ( HYPRE_ParCSRMatrix const & parCSRMatrix );

  /**
   * Boolean value, true if the matrix sparsity pattern has been fixed.
   */
  bool m_is_pattern_fixed = false;

  /**
   * Boolean value, true if the matrix had been finalized, false if not.
   */
  bool m_is_ready_to_use = false;

  /**
   * Pointer to underlying HYPRE_IJMatrix type.
   */
  HYPRE_IJMatrix m_ij_mat = nullptr;

  /**
   * Pointer to underlying HYPRE_ParCSRMatrix type.
   */
  HYPRE_ParCSRMatrix m_parcsr_mat = nullptr;

};

/**
 * @brief Stream insertion operator for EpetraMatrix
 * @param os the output stream
 * @param matrix the matrix to be printed
 * @return reference to the output stream
 */
std::ostream & operator<<( std::ostream & os,
                           HypreMatrix const & matrix );

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_HYPREMATRIX_HPP_*/
