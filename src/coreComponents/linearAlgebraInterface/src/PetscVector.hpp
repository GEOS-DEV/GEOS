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
 * @file PetscVector.hpp
 */

#ifndef LAI_PETSCVECTOR_HPP_
#define LAI_PETSCVECTOR_HPP_

#include "InterfaceTypes.hpp"
#include <petscvec.h>

#include "common/DataTypes.hpp"
#include <vector>
#include <string> // Hannah: need these?
#include <time.h>

namespace geosx
{

/**
 * \class PetscVector
 * \brief This class creates and provides basic support for Vec 
 *        vector object type used in PETSc.
 */
class PetscVector
{

 public:
  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Empty vector constructor.
   *
   * Create an empty (distributed) vector.
   */
  PetscVector();

    /**
   * @brief Copy constructor.
   *
   * \param vec Vector to be copied.
   *
   */
  PetscVector(PetscVector const & vec);

  /* Construct from Petsc vector */
  PetscVector(Vec vec); // Hannah: do I need this?

  /**
   * @brief Virtual destructor.
   */
  virtual ~PetscVector() = default;
  //@}

  //! @name Create Methods
  //@{

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
  // Hannah: to do

  //@}
  //! @name Open / close
  //@{

  /**
   * @brief Empty function for PETSc implementation. May be required by other libraries.
   *
   */
  void open();

  /**
   * @brief Assemble vector
   *
   * // Hannah: blurb here
   */
  void close();

  //@}
  //! @name Add/Set Methods
  //@

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
   * Add into vector value at given element.
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
  void scale(real64 const scalingFactor);

  /**
   * @brief Dot product with the vector vec.
   *
   * \param vec EpetraVector to dot-product with.
   *
   */
  real64 dot( PetscVector const &vec );
  // Hannah: fix spacing later

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>x</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * \param x PetscVector to add.
   *
   */
  void copy(PetscVector const &x);

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>alpha*x + y</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param alpha Scaling factor for added vector.
   * \param x Vector to add.
   *
   */
  void axpy(real64 const alpha, PetscVector const &x);

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
  void axpby(PetscScalar const alpha, PetscVector &x, real64 const beta);

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
  real64 normInf() const

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
   * @brief Returns a single element.
   */
  real64 get(globalIndex globalRow) const;

  /**
   * @brief Returns array of values at globalIndices of the vector. TODO: Not yet implemented, since not built-in
   */
  void get( array1d<globalIndex> const & globalIndices,
            array1d<real64> & values ) const;
  // Hannah: to do

  /**
   * @brief Returns a const pointer to the underlying Vec.
   */
  const Vec* unwrappedPointer() const;
  // Hannah: Vec const * ?

  /**
   * @brief Returns a non-const pointer to the underlying Vec.
   */
  Vec* unwrappedPointer();

  /* Returns vector */
  Vec getConstVec() const;

  /* Returns vector */
  Vec getVec();

  //@}

  //! @name I/O Methods
  //@{

  /**
   * @brief Print the vector in PETSc format to the terminal.
   */
  void print() const;

  /**
   * @brief Write the vector to a matlab-compatible file
   */
  void write( string const & filename ) const;
  // Hannah: to do 

  //@}

 protected:
  
  // Underlying Petsc Vec type
  Vec _vec;
};

} // end geosx namespace

#endif /* PETSCVECTOR_HPP_ */
