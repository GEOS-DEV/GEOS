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
 * @file TrilinosVector.hpp
 *  Created on: Jul 24, 2018
 *  Author: Matthias Cremon
 *
 */

#ifndef TRILINOSVECTOR_HPP_
#define TRILINOSVECTOR_HPP_

#include "InterfaceTypes.hpp"

#include <Epetra_Vector.h>        // TODO: come back and remove nonessential includes
#include <Epetra_FEVector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

#include "common/DataTypes.hpp"

namespace geosx
{

//namespace Tpetra  // placeholder for future Tpetra implementation
//{}

namespace Epetra  // for now use Epetra (specifically Epetra_FEVector)
{

/**
 * \class Vector
 * \brief This class creates and provides basic support for the Epetra_FEVector
 *        vector object type used in Trilinos.  We use the FE version because
 *        Epetra_Vector support for long globalIDs is haphazard.
 */
class Vector
{
public:
  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Empty vector constructor.
   *
   * Create an empty (distributed) vector.
   */
  Vector();

  /**
   * @brief Copy constructor.
   *
   * \param vector Vector to be copied.
   *
   */
  Vector( Vector const & vector );

  /**
   * @brief Virtual destructor.
   */
  virtual ~Vector() = default;
  //@}


  //! @name Create Methods
  //@{

  /**
   * @brief Create a vector from an Epetra_Map.
   *
   * Create a vector from an Epetra_Map.  Allows for maximum flexibility
   * for advanced Trilinos users.
   *
   * \param map Input Epetra Map.
   */
  void create( Epetra_Map const &map);

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
  void createWithLocalSize( trilinosTypes::lid const localSize );

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
  void createWithGlobalSize( trilinosTypes::gid const globalSize );

  /**
   * @brief Construct parallel vector from a local array.
   *
   * Create a vector from local data
   *
   * \param localValues local data to put into vector
   *
   */
  void create( array1d<real64> const localValues );


  //! @name Modify Methods
  //@{

  /**
   * @brief Set vector value.
   *
   * Set vector value at given element.
   *
   * \param globalIndex global row index
   * \param value Value to add at given row.
   *
   */
  void set( trilinosTypes::gid const globalIndex,
            real64 const value );

  /**
   * @brief Set vector values.
   *
   * Set vector values at given elements.
   *
   * \param globalIndices global row indices.
   * \param values Values to add in given rows.
   *
   */
  void set( array1d<trilinosTypes::gid> const globalIndices,
            array1d<real64> const values );

  /**
   * @brief Add into vector value.
   *
   * Add into vector value at given row.
   *
   * \param globalIndex global row.
   * \param value Values to add in given row.
   *
   */
  void add( trilinosTypes::gid const globalIndex,
            real64 const value );

  /**
   * @brief Add into vector values.
   *
   * Add into vector values at given rows.
   *
   * \param globalIndices global rows indices
   * \param values Values to add in given rows.
   *
   */
  void add( array1d<trilinosTypes::gid> const globalIndices,
            array1d<real64> const values );
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
   * \param vec Vector to dot-product with.
   * \param dst Result.
   *
   */
  void dot( Vector const &vec,
            real64 &dst );

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>x</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * \param x Vector to copy.
   *
   */
  void copy( Vector const &x );

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
             Vector const &x );

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
              Vector const &x,
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
  trilinosTypes::gid globalSize() const;

  /**
   * @brief Returns the local size of the vector.
   */
  trilinosTypes::lid localSize() const;

  /**
   * @brief Returns element i of the vector.
   */
  real64 getElement(trilinosTypes::gid i) const;

  /**
   * @brief Returns a const pointer to the underlying Epetra_Vector.
   */
  const Epetra_FEVector* getPointer() const;

  /**
   * @brief Returns a non-const pointer to the underlying Epetra_Vector.
   */
  Epetra_FEVector* getPointer();

  //@}

  //! @name I/O Methods
  //@{

  /**
   * @brief Print the vector in Trilinos format to the terminal.
   */
  void print() const;

  //@}

protected:
  // Unique pointer to underlying Epetra_FEVector type.
  std::unique_ptr<Epetra_FEVector> m_vector = nullptr;
};

} // end Epetra
} // end geosx

#endif /* TRILINOSVECTOR_HPP_ */
