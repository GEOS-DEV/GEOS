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
 * @file PetscSparseMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/PetscVector.hpp"
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
 * @class PetscSparseMatrix
 * @brief This class creates and provides basic support for the Mat
 *        matrix object type used in PETSc.
 */
class PetscSparseMatrix final : public LinearOperator<PetscVector>,
                                private MatrixBase<PetscSparseMatrix, PetscVector>
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
  PetscSparseMatrix();

  /**
   * @brief Copy constructor.
   *
   * Create new matrix from matrix <tt>src</tt>.
   */
  PetscSparseMatrix( PetscSparseMatrix const & src );

  /**
   * @brief Destructor.
   */
  ~PetscSparseMatrix() override;

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

  void createWithLocalSize( localIndex const localRows,
                            localIndex const localCols,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm ) override;

  void createWithGlobalSize( globalIndex const globalRows,
                             globalIndex const globalCols,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm ) override;

  bool created() const override;

  void reset() override;

  void set( real64 const value ) override;

  void zero() override;

  void open() override;

  void close() override;

  void add( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) override;

  void set( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) override;

  void insert( globalIndex const rowIndex,
               globalIndex const colIndex,
               real64 const value ) override;

  void add( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) override;

  void set( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) override;

  void insert( globalIndex const rowIndex,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const size ) override;

  void add( globalIndex const rowIndex,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice1d<real64 const> const & values ) override;

  void set( globalIndex const rowIndex,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice1d<real64 const> const & values ) override;

  void insert( globalIndex const rowIndex,
               arraySlice1d<globalIndex const> const & colIndices,
               arraySlice1d<real64 const> const & values ) override;

  void add( arraySlice1d<globalIndex const> const & rowIndices,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice2d<real64 const, 1> const & values ) override;

  void set( arraySlice1d<globalIndex const> const & rowIndices,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice2d<real64 const, 1> const & values ) override;

  void insert( arraySlice1d<globalIndex const> const & rowIndices,
               arraySlice1d<globalIndex const> const & colIndices,
               arraySlice2d<real64 const, 1> const & values ) override;

  void add( arraySlice1d<globalIndex const> const & rowIndices,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice2d<real64 const, 0> const & values ) override;

  void set( arraySlice1d<globalIndex const> const & rowIndices,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice2d<real64 const, 0> const & values ) override;

  void insert( arraySlice1d<globalIndex const> const & rowIndices,
               arraySlice1d<globalIndex const> const & colIndices,
               arraySlice2d<real64 const, 0> const & values ) override;

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

  void insert( globalIndex const * rowIndices,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const numRows,
               localIndex const numCols ) override;

  void multiply( PetscVector const & src,
                 PetscVector & dst ) const override;

  void multiply( PetscSparseMatrix const & src,
                 PetscSparseMatrix & dst,
                 bool const closeResult = true ) const override;

  void leftMultiplyTranspose( PetscSparseMatrix const & src,
                              PetscSparseMatrix & dst,
                              bool const closeResult = true ) const override;

  void rightMultiplyTranspose( PetscSparseMatrix const & src,
                               PetscSparseMatrix & dst,
                               bool const closeResult = true ) const override;

  void gemv( real64 const alpha,
             PetscVector const & x,
             real64 const beta,
             PetscVector & y,
             bool useTranspose = false ) const override;

  void scale( real64 const scalingFactor ) override;

  void leftScale( PetscVector const & vec ) override;

  void rightScale( PetscVector const & vec ) override;

  void leftRightScale( PetscVector const & vecLeft,
                       PetscVector const & vecRight ) override;

  void clearRow( globalIndex const globalRow,
                 real64 const diagValue = 0 ) override;

  virtual localIndex maxRowLength() const override;

  virtual localIndex localRowLength( localIndex localRowIndex ) const override;

  virtual localIndex globalRowLength( globalIndex globalRowIndex ) const override;

  real64 getDiagValue( globalIndex globalRow ) const override;

  void getRowCopy( globalIndex globalRow,
                   array1d<globalIndex> & colIndices,
                   array1d<real64> & values ) const override;

  globalIndex numGlobalRows() const override;

  globalIndex numGlobalCols() const override;

  localIndex numLocalRows() const override;

  localIndex numLocalCols() const override;

  globalIndex ilower() const override;

  globalIndex iupper() const override;

  localIndex numLocalNonzeros() const override;

  globalIndex numGlobalNonzeros() const override;

  real64 normInf() const override;

  real64 norm1() const override;

  real64 normFrobenius() const override;

  localIndex getLocalRowID( globalIndex const index ) const override;

  globalIndex getGlobalRowID( localIndex const index ) const override;

  MPI_Comm getComm() const override;

  void print( std::ostream & os = std::cout ) const override;

  void write( string const & filename,
              MatrixOutputFormat const format ) const override;

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

};

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_*/
