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
 * @file Polynomial.hpp
 */

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

//#include "legacy/Common/Common.h"
#include <vector>

/*! A class to represent polynomial objects of arbitrary order. */

class Polynomial
{
public:

  Polynomial( const std::vector< double > _coefficients );

  Polynomial( Polynomial const & ) = default;

  ~Polynomial();

  int  Degree();

  double Value( const double x ) const;
  double Deriv( const double x ) const;

  void Evaluate( const double x,
                 double & value,
                 double & deriv );

private:

  std::vector< double > m_coefficients;
};

#endif
