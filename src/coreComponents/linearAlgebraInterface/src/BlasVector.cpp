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
 * @file BlasVector.cpp
 */

// Include the corresponding header file.
#include <BlasVector.hpp>

// Put everything under the geosx namespace.
namespace geosx
{

//----------------------------------------------Constructor/destructor methods---
BlasVector::BlasVector()
{

}

BlasVector::BlasVector( localIndex length)
{
  this->resize( length, 1 );
}

// vector-vector sum (optional scaling)
void BlasVector::vectorAdd(BlasVector& x, const double scalarX)
{

  GEOS_ASSERT_MSG(x.getNumRows() == m_numRows && x.getNumCols() == m_numCols,
                  "Vector dimensions not compatible for sum");
  this->matrixAdd(x, scalarX);

  return;
}

//------------------------------------------------------Data Accessor methods---

// Returns number of matrix rows.
localIndex BlasVector::length() const
{
  return m_numRows;
}

} // end geosx namespace

