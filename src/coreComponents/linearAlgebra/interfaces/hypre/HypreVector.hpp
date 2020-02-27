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
 * @file HypreVector.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREVECTOR_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREVECTOR_HPP_

#include "linearAlgebra/interfaces/VectorBase.hpp"

// Just a placeholder to avoid to include two HYPRE header files
// #include "HYPRE_IJ_mv.h"
// #include "HYPRE_parcsr_mv.h"

// IJVector definition
struct hypre_IJVector_struct;
typedef struct hypre_IJVector_struct *HYPRE_IJVector;

// ParVector definition
struct hypre_ParVector_struct;
typedef struct hypre_ParVector_struct *HYPRE_ParVector;

namespace geosx
{

/**
 * @brief This class creates and provides basic support for the HYPRE_ParVector
 *        vector object type used in Hypre using the linear-algebraic system
 *        interface (IJ interface).
 */
class HypreVector final : private VectorBase<HypreVector>
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
  HypreVector();

  /**
   * @brief Copy constructor.
   * @param src vector to be copied
   *
   */
  HypreVector( HypreVector const & src );

  /**
   * @brief Move constructor.
   * @param src vector to be moved
   */
  HypreVector( HypreVector && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src HypreVector to be copied.
   */
  HypreVector & operator=( HypreVector const & src );

  /**
   * @brief Move assignment.
   * @param src HypreVector to be moved.
   */
  HypreVector & operator=( HypreVector && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~HypreVector();

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

  void reset() override;

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

  real64 dot( HypreVector const & vec ) const override;

  void copy( HypreVector const & x ) override;

  void axpy( real64 const alpha,
             HypreVector const & x ) override;

  void axpby( real64 const alpha,
              HypreVector const & x,
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
            arraySlice1d<real64> const & values ) const override;

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
   * @brief Returns a const pointer to the underlying HYPRE_IJVector object.
   */
  HYPRE_IJVector const & unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying HYPRE_IJVector object.
   */
  HYPRE_IJVector & unwrapped();

  /**
   * @brief Returns a const pointer to the underlying HYPRE_IJVector object.
   */
  HYPRE_ParVector const & unwrappedParVector() const;

  /**
   * @brief Returns a non-const pointer to the underlying HYPRE_IJVector object.
   */
  HYPRE_ParVector & unwrappedParVector();

private:

  /**
   * Pointer to underlying HYPRE_IJVector type.
   */
  HYPRE_IJVector m_ij_vector;

  /**
   * Pointer to underlying HYPRE_ParVector type.
   */
  HYPRE_ParVector m_par_vector;

};

}// end namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREVECTOR_HPP_*/
