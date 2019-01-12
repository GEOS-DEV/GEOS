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

#ifndef GAUSS_QUADRATURE_H
#define GAUSS_QUADRATURE_H

//#include "legacy/Common/Common.h"
#include "finiteElement/quadrature/QuadratureBase.hpp"
//#include "legacy/Utilities/StructuredGridUtilities.h"

#include <cassert>
namespace geosx
{

template<int dim>
class GaussQuadrature : public QuadratureBase
{
public:

  static string CatalogName()
  {
    string name = "GaussQuadrature";
    name.append( std::to_string( dim ) );
    return name;
  }

  GaussQuadrature() = default;
  ~GaussQuadrature() override;

  int size() const override final;
  R1Tensor integration_point( const int index ) const override final;
  double integration_weight( const int index ) const override final;

  void ReadXML( xmlWrapper::xmlNode const & xmlNode ) override final;

private:

  int m_degree;
  int m_n_gauss_points;

  std::vector<double> m_points_1d;
  std::vector<double> m_weights_1d;

};
}

#endif
