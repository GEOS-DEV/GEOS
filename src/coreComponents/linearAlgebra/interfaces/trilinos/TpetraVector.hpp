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
 * @file TpetraVector.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TPETRAVECTOR_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TPETRAVECTOR_HPP_

#include "linearAlgebra/interfaces/VectorBase.hpp"

#include <Tpetra_MultiVector_fwd.hpp>
#include <Tpetra_Vector_fwd.hpp>
#include <Tpetra_Map_fwd.hpp>

namespace geosx
{

/**
 * @brief Wrapper class for Trilinos/Tpetra's Vector class.
 */
class TpetraVector final : private VectorBase< TpetraVector >
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
  TpetraVector();

  /**
   * @brief Copy constructor.
   * @param src vector to be copied
   */
  TpetraVector( TpetraVector const & src );

  /**
   * @brief Move constructor
   * @param src vector to move from
   */
  TpetraVector( TpetraVector && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src vector to be copied
   * @return reference to this object
   */
  TpetraVector & operator=( TpetraVector const & src );

  /**
   * @brief Move assignment.
   * @param src vector to move from
   * @return reference to this object
   */
  TpetraVector & operator=( TpetraVector && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~TpetraVector();

  ///@}

  /**
   * @name VectorBase interface
   */
  ///@{

  using VectorBase::closed;
  using VectorBase::ready;

  virtual bool created() const override;

  virtual void createWithLocalSize( localIndex const localSize,
                                    MPI_Comm const & comm ) override;

  virtual void createWithGlobalSize( globalIndex const globalSize,
                                     MPI_Comm const & comm ) override;

  virtual void createWithLocalValues( arrayView1d< real64 > const & localValues,
                                      MPI_Comm const & comm ) override;

  virtual void open() override;

  virtual void close() override;

  virtual void reset() override;

  virtual void set( globalIndex const globalRowIndex,
                    real64 const value ) override;

  virtual void add( globalIndex const globalRowIndex,
                    real64 const value ) override;

  virtual void set( globalIndex const * globalRowIndices,
                    real64 const * values,
                    localIndex size ) override;

  virtual void add( globalIndex const * globalRowIndices,
                    real64 const * values,
                    localIndex const size ) override;

  virtual void set( arraySlice1d< globalIndex const > const & globalRowIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void add( arraySlice1d< globalIndex const > const & globalRowIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void set( real64 const value ) override;

  virtual void zero() override;

  virtual void rand( unsigned const seed = 1984 ) override;

  virtual void scale( real64 const scalingFactor ) override;

  virtual void reciprocal() override;

  virtual real64 dot( TpetraVector const & vec ) const override;

  virtual void copy( TpetraVector const & x ) override;

  virtual void axpy( real64 const alpha,
                     TpetraVector const & x ) override;

  virtual void axpby( real64 const alpha,
                      TpetraVector const & x,
                      real64 const beta ) override;

  virtual real64 norm1() const override;

  virtual real64 norm2() const override;

  virtual real64 normInf() const override;

  virtual globalIndex globalSize() const override;

  virtual localIndex localSize() const override;

  virtual globalIndex ilower() const override;

  virtual globalIndex iupper() const override;

  virtual real64 get( globalIndex globalRow ) const override;

  virtual void get( arraySlice1d< globalIndex const > const & globalIndices,
                    arraySlice1d< real64 > const & values ) const override;

  virtual localIndex getLocalRowID( globalIndex const globalRow ) const override;

  virtual globalIndex getGlobalRowID( localIndex const localRow ) const override;

  virtual real64 const * extractLocalVector() const override;

  virtual real64 * extractLocalVector() override;

  virtual void extract( arrayView1d< real64 > const & localVector ) const override;

  virtual MPI_Comm getComm() const override;

  virtual void print( std::ostream & os = std::cout ) const override;

  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

  /**
   * @brief Alias for Tpetra map template instantiation used by this class.
   */
  using Tpetra_Map = Tpetra::Map< int, globalIndex >;

  /**
   * @brief Alias for specific Tpetra vector template instantiation wrapped by this class.
   *
   * @note This uses Tpetra's default execution/memory space. When built with CUDA support,
   * this will be equal to Kokkos::Cuda, so we won't be able to create a host-only vector.
   * If we want both in the same executable, we'll have to make adjustments to our LAI approach.
   */
  using Tpetra_Vector = Tpetra::Vector< real64, int, globalIndex >;

  /**
   * @brief Alias for specific Tpetra::MultiVector vector template instantiation wrapped by this class.
   *
   * This is needed for correct specification of template arguments for solver classes (e.g. Belos)
   * which must be templated on MultiVector and not Vector to invoke the right explicit instantiations.
   */
  using Tpetra_MultiVector = Tpetra::MultiVector< real64, int, globalIndex >;

  /**
   * @brief Get the underlying Tpetra object.
   * @return reference to wrapped vector
   */
  Tpetra_Vector const & unwrapped() const;

  /**
   * @copydoc unwrapped() const
   */
  Tpetra_Vector & unwrapped();

private:

  /// Pointer to wrapped Tpetra object.
  std::unique_ptr< Tpetra_Vector > m_vector;

};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_INTERFACES_TPETRAVECTOR_HPP_
