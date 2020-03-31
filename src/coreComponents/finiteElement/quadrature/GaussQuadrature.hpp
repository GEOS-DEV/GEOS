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

#ifndef GAUSS_QUADRATURE_H
#define GAUSS_QUADRATURE_H

//#include "legacy/Common/Common.h"
#include "finiteElement/quadrature/QuadratureBase.hpp"
//#include "legacy/Utilities/StructuredGridUtilities.h"

#include <cassert>
namespace geosx
{

template< int dim >
class GaussQuadrature : public QuadratureBase
{
public:

  static string CatalogName()
  {
    string name = "GaussQuadrature";
    name.append( std::to_string( dim ) );
    return name;
  }

  GaussQuadrature( std::string const & name, Group * const parent );

  virtual ~GaussQuadrature() override;

  virtual void PostProcessInput() override;

  int size() const override final;
  R1Tensor integration_point( const int index ) const override final;
  double integration_weight( const int index ) const override final;

  struct viewKeyStruct
  {
    static constexpr auto degreeString = "degree";
  } viewKeys;

private:

  int m_degree;
  int m_n_gauss_points;

  std::vector< double > m_points_1d;
  std::vector< double > m_weights_1d;

};
}

#endif
