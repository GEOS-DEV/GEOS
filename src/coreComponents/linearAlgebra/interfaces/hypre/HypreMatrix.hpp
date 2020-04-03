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

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "HypreVector.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

// Just a placeholder to avoid to include two HYPRE header files
//#include "_hypre_IJ_mv.h"
//#include "_hypre_parcsr_mv.h"

// IJMatrix definition
struct hypre_IJMatrix_struct;
typedef struct hypre_IJMatrix_struct * HYPRE_IJMatrix;

// ParCSRMatrix definition
struct hypre_ParCSRMatrix_struct;
typedef struct hypre_ParCSRMatrix_struct * HYPRE_ParCSRMatrix;

namespace geosx
{

/**
 * \class HypreMatrix
 * \brief This class ...
 */
class HypreMatrix final : public virtual LinearOperator< HypreVector >,
  private MatrixBase< HypreMatrix, HypreVector >
{
public:

  /**
   * @name Constructor/Destructor Methods
   */
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
  ~HypreMatrix() override;

  ///@}

  /**
   * @name MatrixBase interface
   */
  ///@{

  using MatrixBase::createWithLocalSize;
  using MatrixBase::createWithGlobalSize;
  using MatrixBase::closed;
  using MatrixBase::assembled;
  using MatrixBase::insertable;
  using MatrixBase::modifiable;
  using MatrixBase::ready;

  virtual void createWithLocalSize( localIndex const localRows,
                                    localIndex const localCols,
                                    localIndex const maxEntriesPerRow,
                                    MPI_Comm const & comm = MPI_COMM_WORLD ) override;

  virtual void createWithGlobalSize( globalIndex const globalRows,
                                     globalIndex const globalCols,
                                     localIndex const maxEntriesPerRow,
                                     MPI_Comm const & comm = MPI_COMM_WORLD ) override;

  virtual void open() override;

  virtual void close() override;

  virtual bool created() const override;

  virtual void reset() override;

  virtual void set( real64 const value ) override;

  virtual void zero() override;

  virtual void add( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) override;

  virtual void set( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) override;

  virtual void insert( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value ) override;

  virtual void add( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) override;

  virtual void set( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) override;

  virtual void insert( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const size ) override;

  virtual void add( globalIndex const rowIndex,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void set( globalIndex const rowIndex,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void insert( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values ) override;

  virtual void add( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  virtual void set( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  virtual void insert( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  virtual void add( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  virtual void set( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  virtual void insert( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  virtual void add( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) override;

  virtual void set( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) override;

  virtual void insert( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols ) override;

  virtual void apply( HypreVector const & src,
                      HypreVector & dst ) const override;

  virtual void applyTranspose( Vector const & src,
                               Vector & dst ) const override;

  virtual void multiply( HypreMatrix const & src,
                         HypreMatrix & dst ) const override;

  virtual void leftMultiplyTranspose( HypreMatrix const & src,
                                      HypreMatrix & dst ) const override;

  virtual void rightMultiplyTranspose( HypreMatrix const & src,
                                       HypreMatrix & dst ) const override;

  virtual void multiplyRAP( HypreMatrix const & R,
                            HypreMatrix const & P,
                            HypreMatrix & dst ) const override;

  virtual void multiplyPtAP( HypreMatrix const & P,
                             HypreMatrix & dst ) const override;

  virtual void gemv( real64 const alpha,
                     HypreVector const & x,
                     real64 const beta,
                     HypreVector & y,
                     bool useTranspose = false ) const override;

  virtual void scale( real64 const scalingFactor ) override;

  virtual void leftScale( HypreVector const & vec ) override;

  virtual void rightScale( HypreVector const & vec ) override;

  virtual void leftRightScale( HypreVector const & vecLeft,
                               HypreVector const & vecRight ) override;

  virtual void transpose( HypreMatrix & dst ) const override;

  virtual real64 clearRow( globalIndex const row,
                           bool const keepDiag = false,
                           real64 const diagValue = 0.0 ) override;

  virtual localIndex maxRowLength() const override;

  virtual localIndex localRowLength( localIndex localRowIndex ) const override;

  virtual localIndex globalRowLength( globalIndex globalRowIndex ) const override;

  virtual void getRowCopy( globalIndex globalRow,
                           arraySlice1d< globalIndex > const & colIndices,
                           arraySlice1d< real64 > const & values ) const override;

  virtual real64 getDiagValue( globalIndex globalRow ) const override;

  virtual globalIndex numGlobalRows() const override;

  virtual globalIndex numGlobalCols() const override;

  virtual localIndex numLocalRows() const override;

  virtual localIndex numLocalCols() const override;

  virtual globalIndex ilower() const override;

  virtual globalIndex iupper() const override;

  virtual localIndex numLocalNonzeros() const override;

  virtual globalIndex numGlobalNonzeros() const override;

  virtual real64 normInf() const override;

  virtual real64 norm1() const override;

  virtual real64 normFrobenius() const override;

  virtual localIndex getLocalRowID( globalIndex const index ) const override;

  virtual globalIndex getGlobalRowID( localIndex const index ) const override;

  virtual MPI_Comm getComm() const override;

  virtual void print( std::ostream & os = std::cout ) const override;

  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

  /**
   * @brief Returns a pointer to the underlying HYPRE_IJMatrix object.
   */
  HYPRE_IJMatrix const & unwrapped() const;

  HYPRE_IJMatrix & unwrapped();

  HYPRE_ParCSRMatrix const & unwrappedParCSR() const;

  HYPRE_ParCSRMatrix & unwrappedParCSR();

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
   * Pointer to underlying HYPRE_IJMatrix type.
   */
  HYPRE_IJMatrix m_ij_mat = nullptr;

  /**
   * Pointer to underlying HYPRE_ParCSRMatrix type.
   */
  HYPRE_ParCSRMatrix m_parcsr_mat = nullptr;

};

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMATRIX_HPP_*/
