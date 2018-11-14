/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file EpetraVector.hpp
 *  Created on: Jul 24, 2018
 *  Author: Matthias Cremon
 *
 */

#ifndef EPETRAVECTOR_HPP_
#define EPETRAVECTOR_HPP_

#include "InterfaceTypes.hpp"

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * \class EpetraVector
 * \brief This class creates and provides basic support for the Epetra_Vector
 *        vector object type used in Trilinos.
 */
class EpetraVector
{
public:

  using lid = trilinosTypes::lid;
  using gid = trilinosTypes::gid;

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
   * \param vector Vector to be copied.
   *
   */
  EpetraVector( EpetraVector const & vector );

  /**
   * @brief Virtual destructor.
   */
  virtual ~EpetraVector() = default;
  //@}

  //! @name Create Methods
  //@{we

  /**
   * @brief Construct vector from array.
   *
   * Create a vector from an Epetra_Map and an array of values.
   *
   * \param map Input Epetra Map.
   * \param values Array of values to populate vector.
   *
   */
  void create( Epetra_Map const &map,
               real64  *values );

  /**
   * @brief Construct vector from array.
   *
   * Create a vector from a size and an array of values.
   *
   * \param size Global number of elements.
   * \param values Array of values to populate vector.
   *
   */
  void create( gid const size,
               real64 *values );

  /**
   * @brief Set vector value.
   *
   * Set vector value at given element.
   *
   * \param element elements index.
   * \param value Values to add in given element.
   *
   */
  void set( gid element,
            real64 value );

  /**
   * @brief Set vector values.
   *
   * Set vector values at given elements.
   *
   * \param elements elements indices.
   * \param values Values to add in given rows.
   *
   */
  void set( array1d<gid> elements,
            array1d<real64> values );

  /**
   * @brief Add into vector value.
   *
   * Add into vector value at given element.
   *
   *
   * \param element elements index.
   * \param value Values to add in given element.
   *
   */
  void add( gid element,
            real64 value );

  /**
   * @brief Add into vector values.
   *
   * Add into vector values at given elements.
   *
   *
   * \param elements elements indices.
   * \param values Values to add in given rows.
   *
   */
  void add( array1d<gid> const elements,
            array1d<real64> const values );


  void add( lid const numIndices,
            gid const * const elements,
            real64 const * const values );

  /**
   * @brief Construct vector from std::vector.
   *
   * Create a vector from an std vector.
   *
   * \param vec Std vector to cast as EpetraVector
   *
   */
  void create( std::vector<real64> &vec );
  //@}

  //! @name Linear Algebra Methods
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
   * \param vec Vector to dot-product with.
   * \param dst Result.
   *
   */
  void dot( EpetraVector const &vec,
            real64 &dst );

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>x</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * \param x Vector to add.
   *
   */
  void copy( EpetraVector const &x );

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>alpha*x + y</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param alpha Scaling factor for added vector.
   * \param x Vector to add.
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
   * \param x Vector to add.
   * \param beta Scaling factor for self vector.
   *
   */
  void axpby( real64 const alpha,
              EpetraVector const &x,
              real64 const beta );

  /**
   * @brief 1-norm of the vector.
   *
   * \param 1-norm of the vector.
   *
   */
  void norm1( real64 &result ) const;

  /**
   * @brief 2-norm of the vector.
   *
   * \param 2-norm of the vector.
   *
   */
  void norm2( real64 &result ) const;

  /**
   * @brief Infinity-norm of the vector.
   *
   * \param Inf-norm of the vector.
   *
   */
  void normInf( real64 &result ) const;

  //@}

  //! @name Accessor Methods
  //@{

  /**
   * @brief Returns the global of the vector.
   */
  gid globalSize() const;

  /**
   * @brief Returns the local size of the vector.
   */
  trilinosTypes::lid localSize() const;


  /**
   *
   * @return Returns a pointer to the const data in the vector
   */
  real64 const * getValues() const;

  /**
   *
   * @return Returns a pointer to the data in the vector
   */
  real64 * getValues();

  /**
   * @brief Returns element i of the vector.
   */
  real64 getElement(gid i) const;

  /**
   * @brief Returns a const pointer to the underlying Epetra_Vector.
   */
  const Epetra_Vector* getPointer() const;

  /**
   * @brief Returns a non-const pointer to the underlying Epetra_Vector.
   */
  Epetra_Vector* getPointer();

  //@}

  //! @name I/O Methods
  //@{

  /**
   * @brief Print the vector in Trilinos format to the terminal.
   */
  void print( std::ostream & outputStream ) const;

  //@}

protected:
  // Unique pointer to underlying Epetra_Vector type.
  std::unique_ptr<Epetra_Vector> m_vector = nullptr;
};

}

#endif /* EPETRAVECTOR_HPP_ */
