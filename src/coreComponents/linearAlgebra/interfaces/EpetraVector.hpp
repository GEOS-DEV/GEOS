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
#include "mpiCommunications/MpiWrapper.hpp"

class Epetra_FEVector;
class Epetra_Map;

namespace geosx
{

/**
 * \class EpetraVector
 * \brief This class creates and provides basic support for the Epetra_FEVector
 *        vector object type used in Trilinos.  We use the FE version because
 *        Epetra_Vector support for long globalIDs is haphazard.
 */
class EpetraVector
{
public:
  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Empty vector constructor.
   *
   * Create an empty (distributed) vector.
   */
  EpetraVector();

  /**
   * @brief Copy constructor.
   *
   * \param src EpetraVector to be copied.
   *
   */
  EpetraVector( EpetraVector const & src );

  /**
   * @brief Destructor.
   */
  ~EpetraVector();
  //@}


  //! @name Create Methods
  //@{

  /**
   * @brief Create a vector based on a previous vector.
   *
   * \param vector an already formed EpetraVector.
   *
   */
  void create( EpetraVector const & src );

  /**
   * @brief Create a vector based on local number of elements.
   *
   * Create a vector based on local number of elements.  Global size is
   * the sum across processors.  For specifying a global size and having
   * automatic partitioning, see createGlobal().
   *
   * \param localSize local number of elements.
   *
   */
  void createWithLocalSize( localIndex const localSize, MPI_Comm const & comm = MPI_COMM_WORLD );

  /**
   * @brief Create a vector based on global number of elements.
   *
   * Create a vector based on global number of elements. Every processors
   * gets the same number of local elements except proc 0, which gets any
   * remainder elements as well if the split can't be done evenly.
   *
   * \param globalSize Global number of elements.
   *
   */
  void createWithGlobalSize( globalIndex const globalSize, MPI_Comm const & comm = MPI_COMM_WORLD );

  /**
   * @brief Construct parallel vector from a local array.
   *
   * Create a vector from local data
   *
   * \param localValues local data to put into vector
   *
   */
  void create( array1d<real64> const & localValues, MPI_Comm const & comm = MPI_COMM_WORLD );

  //@}
  //! @name Open / close
  //@{

  /**
   * @brief Empty function for Trilinos implementation. May be required by other libraries.
   *
   */
  void open();

  /**
   * @brief Assemble vector
   *
   * Performs parallel communication to scatter assembled entries to appropriate locations
   */
  void close();

  //@}
  //! @name Add/Set Methods
  //@{

  /**
   * @brief Set vector value.
   *
   * Set vector value at given element.
   *
   * \param globalRow global row index
   * \param value Value to add at given row.
   *
   */
  void set( globalIndex const globalRow,
            real64 const value );

  /**
   * @brief Add into vector value.
   *
   * Add into vector value at given row.
   *
   * \param globalRow global row.
   * \param value Values to add in given row.
   *
   */
  void add( globalIndex const globalRow,
            real64 const value );

  /**
   * @brief Set vector values.
   *
   * Set vector values at given elements.
   *
   * \param globalIndices global row indices.
   * \param values Values to add in given rows.
   * \param size Number of elements
   *
   */
  void set( globalIndex const * globalIndices,
            real64 const * values,
            localIndex size );

  /**
   * @brief Add vector values.
   *
   * Add vector values at given elements.
   *
   * \param globalIndices global row indices.
   * \param values Values to add in given rows.
   * \param size Number of elements
   *
   */
  void add( globalIndex const * globalIndices,
            real64 const * values,
            localIndex size );

  /**
   * @brief Set vector values using array1d
   *
   * Set vector values at given elements.
   *
   * \param globalIndices global row indices.
   * \param values Values to add in given rows.
   *
   */
  void set( array1d<globalIndex> const & globalIndices,
            array1d<real64> const & values );


  /**
   * @brief Add into vector values using array1d
   *
   * Add into vector values at given rows.
   *
   * \param globalIndices global rows indices
   * \param values Values to add in given rows.
   *
   */
  void add( array1d<globalIndex> const & globalIndices,
            array1d<real64> const & values );

  /**
   * @brief Set all elements to a constant value.
   *
   * \param value Values to set vector elements to.
   *
   */
  void set( real64 const value );

  /**
   * @brief Set vector elements to zero.
   *
   */
  void zero();

  /**
   * @brief Set vector elements to random entries.
   *
   */
  void rand();

  //@}

  //! @name Algebraic Operations
  //@{

  /**
   * @brief Multiply all elements by scalingFactor.
   *
   * \param scalingFactor Scaling Factor.
   *
   */
  void scale( real64 const scalingFactor );

  /**
   * @brief Dot product with the vector vec.
   *
   * \param vec EpetraVector to dot-product with.
   *
   */
  real64 dot( EpetraVector const &vec ) const;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>x</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * \param x EpetraVector to copy.
   *
   */
  void copy( EpetraVector const &x );

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>alpha*x + y</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param alpha Scaling factor for added vector.
   * \param x EpetraVector to add.
   *
   */
  void axpy( real64 const alpha,
             EpetraVector const &x );

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>alpha*x + beta*y</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param alpha Scaling factor for added vector.
   * \param x EpetraVector to add.
   * \param beta Scaling factor for self vector.
   *
   */
  void axpby( real64 const alpha,
              EpetraVector const &x,
              real64 const beta );

  /**
   * @brief 1-norm of the vector.
   *
   */
  real64 norm1() const;

  /**
   * @brief 2-norm of the vector.
   *
   */
  real64 norm2() const;

  /**
   * @brief Infinity-norm of the vector.
   *
   */
  real64 normInf() const;

  //@}

  //! @name Accessor Methods
  //@{

  /**
   * @brief Returns the global of the vector.
   */
  globalIndex globalSize() const;

  /**
   * @brief Returns the local size of the vector.
   */
  localIndex localSize() const;

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
   * @brief Returns value globalRow of the vector. TODO: Not yet implemented, since not built-in
   */
  real64 get( globalIndex globalRow ) const;

  /**
   * @brief Returns array of values at globalIndices of the vector. TODO: Not yet implemented, since not built-in
   */
  void get( array1d<globalIndex> const & globalIndices,
            array1d<real64> & values ) const;

  /**
   * @brief Returns a const pointer to the underlying Epetra object.
   */
  Epetra_FEVector const * unwrappedPointer() const;

  /**
   * @brief Returns a non-const pointer to the underlying Epetra object.
   */
  Epetra_FEVector * unwrappedPointer();

  //@}

  //! @name I/O Methods
  //@{

  /**
   * @brief Print the vector in Trilinos format to the terminal.
   */
  void print( std::ostream & os = std::cout ) const;

  /**
   * @brief Write the vector to a matlab-compatible file
   */
  void write( string const & filename,
              bool const mtxFormat = true ) const;

  /**
   * Map a global row index to local row index
   */
  localIndex getLocalRowID( globalIndex const index ) const;

  /**
   * Map a local row index to global row index
   */
  globalIndex getGlobalRowID( localIndex const index ) const;

  /**
   * Extract a view of the local portion of the array
   */
  real64 const * extractLocalVector() const;

  /**
   * Extract a view of the local portion of the array
   */
  real64 * extractLocalVector();

  //@}

private:

  /**
   * Unique pointer to underlying Epetra_FEVector type.
   */
  std::unique_ptr<Epetra_FEVector> m_vector;

  /**
   * @brief Create a vector from an Epetra_Map.
   *
   * Create a vector from an Epetra_Map.  Allows for maximum flexibility
   * for advanced users.
   *
   * \param map Input Epetra Map.
   */
  void create( Epetra_Map const & map );
};

/**
 * @brief Stream insertion operator for EpetraVector
 * @param os the output stream
 * @param vec the vector to be printed
 * @return reference to the output stream
 */
std::ostream & operator<<( std::ostream & os, EpetraVector const & vec );

} // end geosx namespace

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_*/
