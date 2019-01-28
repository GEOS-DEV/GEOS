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
 *  Created on: Jan 28, 2019
 *  Author: Hannah Morgan
 *
 */

#ifndef PETSCVECTOR_HPP_
#define PETSCVECTOR_HPP_

#include "InterfaceTypes.hpp"
#include <petscvec.h>
#include "common/DataTypes.hpp"

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
   * \param vector Vector to be copied.
   *
   */
  PetscVector(PetscVector const & vec);

  /* Construct from Petsc vector */
  PetscVector(Vec vec);

  /**
   * @brief Virtual destructor.
   */
  virtual ~PetscVector() = default;
  //@}

  //! @name Create Methods
  //@{we

  /**
   * @brief Construct vector from array.
   *
   * Create a vector from a size and an array of values.
   *
   * \param size Global number of elements.
   * \param values Array of values to populate vector.
   *
   */
  void create(int const size, real64 *values);

  /**
   * @brief Construct vector from std::vector.
   *
   * Create a vector from an std vector.
   *
   * \param vec Std vector to cast as PetscVector
   *
   */
  void create(std::vector<real64> &vec);

  /**
   * @brief Set vector value.
   *
   * Set vector value at given element.
   *
   * \param element elements index.
   * \param value Values to add in given element.
   *
   */
  void set(int element, real64 value);

  // /* set values of vector elements */
  // void set(array1d<int> elements, array1d<real64> values);

  /**
   * @brief Add into vector value.
   *
   * Add into vector value at given element.
   *
   * (TODO This needs to use integers for some reason! No longlong).
   *
   * \param element elements index.
   * \param value Values to add in given element.
   *
   */
  void add(int element, real64 value);

  // /* add values to vector elements */
  // void add(array1d<int> elements, array1d<real64> values);

  //@}

  //! @name Linear Algebra Methods
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
   * \param vec Vector to dot-product with.
   * \param dst Result.
   *
   */
  void dot(PetscVector const vec, real64 *dst);

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>x</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * \param x Vector to add.
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
  void axpby(real64 const alpha, PetscVector &x, real64 const beta);

  /**
   * @brief 1-norm of the vector.
   *
   * \param 1-norm of the vector.
   *
   */
  void norm1(real64 &result) const;

  /**
   * @brief 2-norm of the vector.
   *
   * \param 2-norm of the vector.
   *
   */
  void norm2(real64 &result) const;

  /**
   * @brief Infinity-norm of the vector.
   *
   * \param Inf-norm of the vector.
   *
   */
  void normInf(real64 &result) const;

  /**
   * @brief Returns the global of the vector.
   */
  int globalSize() const;

  /**
   * @brief Returns the local size of the vector.
   */
  int localSize() const;

  /**
   * @brief Returns a const pointer to the underlying Vec.
   */
  const Vec* getPointer() const;

  /**
   * @brief Returns a non-const pointer to the underlying Vec.
   */
  const Vec* getPointer();

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

 protected:
  
  // Underlying Petsc Vec type
  Vec _vec;
};

}

#endif /* PETSCVECTOR_HPP_ */
