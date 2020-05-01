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
 * @file Polynomial.cpp
 */

#include "Polynomial.hpp"

/*
 * Constructor.  Takes a vector of coefficients
 */

Polynomial :: Polynomial( const std::vector< double > _coefficients )
  :
  m_coefficients( _coefficients )
{}


/*
 * Destructor.
 */

Polynomial :: ~Polynomial()
{}


/*
 * Return polynomial degree
 */

int Polynomial :: Degree ()
{
  return static_cast< int >(m_coefficients.size());
}


/*
 * Use Horner's method to evaluate the polynomial
 * function value p(x).
 */

double Polynomial :: Value ( const double x ) const
{
  std::vector< double >::const_reverse_iterator
    it     = m_coefficients.rbegin(),
    end_it = m_coefficients.rend();

  double val = 0;
  for(; it != end_it; ++it )
    val = *it + val*x;

  return val;
}


/*
 * Use Horner's method to evaluate the polynomial
 * function derivative p'(x).
 */

double Polynomial :: Deriv ( const double x ) const
{
  std::vector< double >::const_reverse_iterator
    it     = m_coefficients.rbegin(),
    end_it = m_coefficients.rend();

  double value = *it;
  double deriv = 0;
  ++it;

  for(; it != end_it; ++it )
  {
    deriv = value + deriv*x;
    value = *it + value*x;
  }
  return deriv;
}


/*
 * Use Horner's method to recursively evaluate
 * the polynomial value p(x) and derivative p'(x).
 * Here, both are computed simultaneously to
 * maximize efficiency.
 */

void Polynomial :: Evaluate ( const double x,
                              double & value,
                              double & deriv )
{
  std::vector< double >::reverse_iterator
    it     = m_coefficients.rbegin(),
    end_it = m_coefficients.rend();

  value = *it;
  deriv = 0;
  ++it;

  for(; it != end_it; ++it )
  {
    deriv = value + deriv*x;
    value = *it + value*x;
  }
}
