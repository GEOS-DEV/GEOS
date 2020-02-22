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
 * @file EpetraVector.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/VectorBase.hpp"

class Epetra_FEVector;

namespace geosx
{

/**
 * @brief This class creates and provides basic support for the Epetra_FEVector
 *        vector object type used in Trilinos.  We use the FE version because
 *        Epetra_Vector support for long globalIDs is haphazard.
 */
class EpetraVector final : private VectorBase<EpetraVector>
{
public:

  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Empty vector constructor.
   * Create an empty (distributed) vector.
   */
  EpetraVector();

  /**
   * @brief Copy constructor.
   * @param src EpetraVector to be copied.
   */
  EpetraVector( EpetraVector const & src );

  /**
   * @brief Move constructor
   * @param src EpetraVector to move from
   */
  EpetraVector( EpetraVector && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src EpetraVector to be copied.
   */
  EpetraVector & operator=( EpetraVector const & src );

  /**
   * @brief Move assignment.
   * @param src EpetraVector to be moved from.
   */
  EpetraVector & operator=( EpetraVector && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~EpetraVector();

  ///@}

  /**
   * @name VectorBase interface
   */
  ///@{

  using VectorBase::closed;
  using VectorBase::ready;

  bool created() const override;

  void createWithLocalSize( localIndex const localSize, MPI_Comm const & comm ) override;

  void createWithGlobalSize( globalIndex const globalSize, MPI_Comm const & comm ) override;

  void create( arraySlice1d<real64 const> const & localValues, MPI_Comm const & comm ) override;

  void open() override;

  void close() override;

  virtual void reset() override;

  void set( globalIndex const globalRowIndex,
            real64 const value ) override;

  void add( globalIndex const globalRowIndex,
            real64 const value ) override;

  void set( globalIndex const * globalRowIndices,
            real64 const * values,
            localIndex size ) override;

  void add( globalIndex const * globalRowIndices,
            real64 const * values,
            localIndex const size ) override;

  void set( arraySlice1d<globalIndex const> const & globalRowIndices,
            arraySlice1d<real64 const> const & values ) override;

  void add( arraySlice1d<globalIndex const> const & globalRowIndices,
            arraySlice1d<real64 const> const & values ) override;

  void set( real64 const value ) override;

  void zero() override;

  void rand( unsigned const seed = 1984 ) override;

  void scale( real64 const scalingFactor ) override;

  real64 dot( EpetraVector const & vec ) const override;

  void copy( EpetraVector const & x ) override;

  void axpy( real64 const alpha,
             EpetraVector const & x ) override;

  void axpby( real64 const alpha,
              EpetraVector const & x,
              real64 const beta ) override;

  real64 norm1() const override;

  real64 norm2() const override;

  real64 normInf() const override;

  globalIndex globalSize() const override;

  localIndex localSize() const override;

  globalIndex ilower() const override;

  globalIndex iupper() const override;

  real64 get( globalIndex globalRow ) const override;

  void get( arraySlice1d<globalIndex const> const & globalRowIndices,
            array1d<real64> & values ) const override;

  localIndex getLocalRowID( globalIndex const globalRowIndex ) const override;

  globalIndex getGlobalRowID( localIndex const localRowIndex ) const override;

  real64 const * extractLocalVector() const override;

  real64 * extractLocalVector() override;

  MPI_Comm getComm() const override;

  void print( std::ostream & os = std::cout ) const override;

  void write( string const & filename,
              LAIOutputFormat const format ) const override;

  ///@}

  /**
   * @brief Returns a const pointer to the underlying Epetra object.
   */
  Epetra_FEVector const & unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying Epetra object.
   */
  Epetra_FEVector & unwrapped();

private:

  /// Unique pointer to underlying Epetra_FEVector object.
  std::unique_ptr<Epetra_FEVector> m_vector;
};

} // end geosx namespace

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_*/
