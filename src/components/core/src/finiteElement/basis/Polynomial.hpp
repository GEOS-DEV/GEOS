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
 * @file Polynomial.h
 * @author white230
 */

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

//#include "legacy/Common/Common.h"
#include <vector>

/*! A class to represent polynomial objects of arbitrary order. */

class Polynomial
{
public:

  Polynomial(const std::vector<double> _coefficients);

  Polynomial( Polynomial const & ) = default;

  ~Polynomial();

  int  Degree();

  double Value(const double x) const;
  double Deriv(const double x) const;

  void Evaluate(const double x,
                double &value,
                double &deriv);

private:

  std::vector<double> m_coefficients;
};

#endif
