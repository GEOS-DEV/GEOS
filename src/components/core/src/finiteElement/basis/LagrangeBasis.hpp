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

#ifndef LAGRANGE_BASIS_H
#define LAGRANGE_BASIS_H

//#include "legacy/Common/Common.h"
#include <cassert>

#include "BasisBase.hpp"
#include "Polynomial.hpp"
//#include "legacy/Utilities/StructuredGridUtilities.h"

/** A class to define a parent finite element space consisting
 *  of Lagrangian polynomial basis function. The class is
 *  templated on the spatial dimension. */

namespace geosx
{

template<int dim>
class LagrangeBasis : public BasisBase
{
public:
  static string CatalogName()
  {
    string name = "LagrangeBasis";
    name.append(std::to_string(dim));
    return name;
  }

  LagrangeBasis(void){}

  LagrangeBasis( const int degree );

  int size() const override final;

  double value( const int index,
                const R1Tensor &point ) const override final;

  R1Tensor gradient( const int index,
                     const R1Tensor &point ) const override final;

  R1Tensor support_point( const int index ) override final;

  void ReadXML( xmlWrapper::xmlNode const & targetNode ) override final;

private:

  int m_degree;
  int n_shape_functions;

  std::vector<Polynomial> m_polynomials;
};
}
#endif
