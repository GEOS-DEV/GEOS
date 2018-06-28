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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file Polynomial.cpp
 * @author white203
 * @date Apr 11, 2010
 */

#include "Polynomial.hpp"

/*
 * Constructor.  Takes a vector of coefficients
 */

Polynomial :: Polynomial(const std::vector<double> _coefficients)
  :
  m_coefficients(_coefficients)
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
  return static_cast<int>(m_coefficients.size());
}


/*
 * Use Horner's method to evaluate the polynomial
 * function value p(x).
 */

double Polynomial :: Value (const double x) const
{
  std::vector<double>::const_reverse_iterator
    it     = m_coefficients.rbegin(),
    end_it = m_coefficients.rend();

  double val = 0;
  for( ; it != end_it ; ++it)
    val = *it + val*x;

  return val;
}


/*
 * Use Horner's method to evaluate the polynomial
 * function derivative p'(x).
 */

double Polynomial :: Deriv (const double x) const
{
  std::vector<double>::const_reverse_iterator
    it     = m_coefficients.rbegin(),
    end_it = m_coefficients.rend();

  double value = *it;
  double deriv = 0;
  ++it;

  for( ; it != end_it ; ++it)
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

void Polynomial :: Evaluate (const double x,
                             double &value,
                             double &deriv)
{
  std::vector<double>::reverse_iterator
    it     = m_coefficients.rbegin(),
    end_it = m_coefficients.rend();

  value = *it;
  deriv = 0;
  ++it;

  for( ; it != end_it ; ++it)
  {
    deriv = value + deriv*x;
    value = *it + value*x;
  }
}
