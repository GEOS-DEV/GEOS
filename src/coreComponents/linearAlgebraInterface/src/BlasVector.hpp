/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file BlasVector.hpp
 */

#ifndef CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASVECTOR_HPP_
#define CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASVECTOR_HPP_

#include "BlasMatrix.hpp"

namespace geosx
{

/**
 * \class BlasVector
 * \brief This class creates and provides basic support for for manipulating
 *        a BLAS-style vectors
 */

class BlasVector : public BlasMatrix
{
public:
   //----------------------------------------------------------------------------
   //! @name Constructor/Destructor Methods
   //@{

  /**
   * @brief Default constructor.
   */
  BlasVector();

   /**
    * @brief Shaped constructor; defines a variable-sized vector.
    *
    * \param IN
    * length - vector size.
    *
    * @note All values are initialized to 0.
    */
    BlasVector(localIndex length);

   //@}

   //----------------------------------------------------------------------------
   //! @name Mathematical methods
   //@{

   /**
    * @brief Vector-Vector sum;
    * \a this = scalarX*<tt>x</tt> + \a this.
    *
    * Computes (scalarX*<tt>x</tt> + \a this) and overwrites the result on the
    * current \a this, with optional scaling.
    *
    * \param IN
    * <tt>x</tt> - Dense vector.
    * \param [IN]
    * scalarx - Optional scalar to multiply with <tt>x</tt>.
    *
    * @warning
    * Assumes that <tt>x</tt> and \a this have the same size.
    *
    * @note This function calls BLAS function DAXPY using Trilinos/Epetra BLAS
    *       Wrapper Class if built with Trilinos.
    */
   void vectorAdd(BlasVector& x, real64 const scalarX=1.);

   //@}

   //----------------------------------------------------------------------------
   //! @name Data Accessor methods
   //@{

   using BlasMatrix::operator(); //overload base-class function operator()

   /**
    * @brief Element access function.
    *
    * The parentheses operator returns the element in the ith vector entry.
    */
   inline real64 & operator ()( localIndex Index);

   /**
    * @brief Element access function.
    *
    * The parentheses operator returns the element in the ith vector entry.
    */
   inline real64 const & operator ()( localIndex Index) const;

   /**
    * @brief Returns vector length.
    */
   localIndex length() const;

   //@}

};

// Inline methods
inline real64 & BlasVector::operator ()( localIndex index )
{
  GEOS_ASSERT_MSG( 0 <= index &&
                   index <= m_numRows,
                   "Requested value out of bounds" );
  return m_values[index];
}

inline real64 const & BlasVector::operator ()( localIndex index) const
                                               {
  GEOS_ASSERT_MSG( 0 <= index &&
                   index <= m_numRows,
                   "Requested value out of bounds" );
  return m_values[index];
}

} // namespace geosx


#endif /* CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLASVECTOR_HPP_ */
