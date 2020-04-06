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
 * @file PetscMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "PetscVector.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

/**
 * See comment in PetscVector.hpp about the rationale for this declaration.
 */
struct _p_Mat;
typedef struct _p_Mat * Mat;

namespace geosx
{

/**
 * @class PetscMatrix
 * @brief This class creates and provides basic support for the Mat
 *        matrix object type used in PETSc.
 */
class PetscMatrix final : public virtual LinearOperator< PetscVector >,
  private MatrixBase< PetscMatrix, PetscVector >
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
  PetscMatrix();

  /**
   * @brief Copy constructor.
   *
   * Create new matrix from matrix <tt>src</tt>.
   */
  PetscMatrix( PetscMatrix const & src );

  /**
   * @brief Destructor.
   */
  ~PetscMatrix() override;

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
                                    MPI_Comm const & comm ) override;

  virtual void createWithGlobalSize( globalIndex const globalRows,
                                     globalIndex const globalCols,
                                     localIndex const maxEntriesPerRow,
                                     MPI_Comm const & comm ) override;

  virtual bool created() const override;

  virtual void reset() override;

  virtual void set( real64 const value ) override;

  virtual void zero() override;

  virtual void open() override;

  virtual void close() override;

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
                    arraySlice2d< real64 const, 1 > const & values ) override;

  virtual void set( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const, 1 > const & values ) override;

  virtual void insert( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, 1 > const & values ) override;

  virtual void add( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const, 0 > const & values ) override;

  virtual void set( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const, 0 > const & values ) override;

  virtual void insert( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const, 0 > const & values ) override;

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

  virtual void apply( PetscVector const & src,
                      PetscVector & dst ) const override;

  virtual void multiply( PetscMatrix const & src,
                         PetscMatrix & dst ) const override;

  virtual void applyTranspose( Vector const & src,
                               Vector & dst ) const override;

  virtual void leftMultiplyTranspose( PetscMatrix const & src,
                                      PetscMatrix & dst ) const override;

  virtual void rightMultiplyTranspose( PetscMatrix const & src,
                                       PetscMatrix & dst ) const override;

  virtual void multiplyRAP( PetscMatrix const & R,
                            PetscMatrix const & P,
                            PetscMatrix & dst ) const override;

  virtual void multiplyPtAP( PetscMatrix const & P,
                             PetscMatrix & dst ) const override;

  virtual void gemv( real64 const alpha,
                     PetscVector const & x,
                     real64 const beta,
                     PetscVector & y,
                     bool useTranspose = false ) const override;

  virtual void scale( real64 const scalingFactor ) override;

  virtual void leftScale( PetscVector const & vec ) override;

  virtual void rightScale( PetscVector const & vec ) override;

  virtual void leftRightScale( PetscVector const & vecLeft,
                               PetscVector const & vecRight ) override;

  virtual void transpose( PetscMatrix & dst ) const override;

  virtual real64 clearRow( globalIndex const row,
                           bool const keepDiag = false,
                           real64 const diagValue = 0.0 ) override;

  virtual localIndex maxRowLength() const override;

  virtual localIndex localRowLength( localIndex localRowIndex ) const override;

  virtual localIndex globalRowLength( globalIndex globalRowIndex ) const override;

  virtual real64 getDiagValue( globalIndex globalRow ) const override;

  virtual void getRowCopy( globalIndex globalRow,
                           arraySlice1d< globalIndex > const & colIndices,
                           arraySlice1d< real64 > const & values ) const override;

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
   * @brief Returns a pointer to the underlying matrix.
   */
  const Mat & unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying matrix.
   */
  Mat & unwrapped();

private:

  /// Underlying Petsc object.
  Mat m_mat;

  /// Indices of rows to be cleared on next close()
  array1d< globalIndex > m_rowsToClear;

  /// Diagonal values of rows to be set on next close()
  array1d< real64 > m_diagValues;

};

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_*/
