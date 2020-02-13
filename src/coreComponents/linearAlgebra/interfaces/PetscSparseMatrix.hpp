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
class PetscSparseMatrix : public MatrixBase<PetscSparseMatrix, PetscVector>
{
public:

  using Base = MatrixBase<PetscSparseMatrix, PetscVector>;

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
  virtual ~PetscSparseMatrix() override;

  ///@}

  /**
   * @name MatrixBase interface
   */
  ///@{

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

  bool isCreated() const final;

  void reset() final;

  void set( real64 const value ) final;

  void zero() final;

  void open() final;

  void close() final;

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
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice1d<real64 const> const & values ) final;

  void set( globalIndex const rowIndex,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice1d<real64 const> const & values ) final;

  void insert( globalIndex const rowIndex,
               arraySlice1d<globalIndex const> const & colIndices,
               arraySlice1d<real64 const> const & values ) final;

  void add( arraySlice1d<globalIndex const> const & rowIndices,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice2d<real64 const, 1> const & values ) final;

  void set( arraySlice1d<globalIndex const> const & rowIndices,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice2d<real64 const, 1> const & values ) final;

  void insert( arraySlice1d<globalIndex const> const & rowIndices,
               arraySlice1d<globalIndex const> const & colIndices,
               arraySlice2d<real64 const, 1> const & values ) final;

  void add( arraySlice1d<globalIndex const> const & rowIndices,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice2d<real64 const, 0> const & values ) final;

  void set( arraySlice1d<globalIndex const> const & rowIndices,
            arraySlice1d<globalIndex const> const & colIndices,
            arraySlice2d<real64 const, 0> const & values ) final;

  void insert( arraySlice1d<globalIndex const> const & rowIndices,
               arraySlice1d<globalIndex const> const & colIndices,
               arraySlice2d<real64 const, 0> const & values ) final;

  virtual void add( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) final;

  virtual void set( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) final;

  void insert( globalIndex const * rowIndices,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const numRows,
               localIndex const numCols ) final;

  void multiply( PetscVector const & src,
                 PetscVector & dst ) const final;

  void multiply( PetscSparseMatrix const & src,
                 PetscSparseMatrix & dst,
                 bool const closeResult = true ) const final;

  void leftMultiplyTranspose( PetscSparseMatrix const & src,
                              PetscSparseMatrix & dst,
                              bool const closeResult = true ) const final;

  void rightMultiplyTranspose( PetscSparseMatrix const & src,
                               PetscSparseMatrix & dst,
                               bool const closeResult = true ) const final;

  void gemv( real64 const alpha,
             PetscVector const & x,
             real64 const beta,
             PetscVector & y,
             bool useTranspose = false ) const final;

  void scale( real64 const scalingFactor ) final;

  void leftScale( PetscVector const & vec ) final;

  void rightScale( PetscVector const & vec ) final;

  void leftRightScale( PetscVector const & vecLeft,
                       PetscVector const & vecRight ) final;

  void clearRow( globalIndex const globalRow,
                 real64 const diagValue = 0 ) final;

  real64 getDiagValue( globalIndex globalRow ) const final;

  void getRowCopy( globalIndex globalRow,
                   array1d<globalIndex> & colIndices,
                   array1d<real64> & values ) const final;

  globalIndex globalRows() const final;

  globalIndex globalCols() const final;

  globalIndex ilower() const final;

  globalIndex iupper() const final;

  localIndex localNonzeros() const final;

  globalIndex globalNonzeros() const final;

  real64 normInf() const final;

  real64 norm1() const final;

  real64 normFrobenius() const final;

  localIndex getLocalRowID( globalIndex const index ) const final;

  globalIndex getGlobalRowID( localIndex const index ) const final;

  localIndex localRows() const final;

  localIndex localCols() const final;

  MPI_Comm getComm() const final;

  void print( std::ostream & os = std::cout ) const final;

  void write( string const & filename,
              MatrixOutputFormat const format ) const final;

  /// @}

  /**
   * @brief Returns a pointer to the underlying matrix.
   */
  const Mat & unwrappedPointer() const;

  /**
   * @brief Returns a non-const pointer to the underlying matrix.
   */
  Mat & unwrappedPointer();

private:

  /// Underlying Petsc object.
  Mat m_mat;

};

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_*/
