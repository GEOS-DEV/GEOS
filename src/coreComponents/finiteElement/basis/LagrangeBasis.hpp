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

template< int dim >
class LagrangeBasis : public BasisBase
{
public:
  static string CatalogName()
  {
    string name = "LagrangeBasis";
    name.append( std::to_string( dim ));
    return name;
  }

  LagrangeBasis( std::string const & name, Group * const parent );

  virtual ~LagrangeBasis() override;

  virtual void PostProcessInput() override;

  int size() const override final;

  double value( const int index,
                const R1Tensor & point ) const override final;

  R1Tensor gradient( const int index,
                     const R1Tensor & point ) const override final;

  R1Tensor support_point( const int index ) override final;

  struct viewKeyStruct
  {
    static constexpr auto degreeString = "degree";
  } viewKeys;

private:

  int m_degree;
  int n_shape_functions;

  std::vector< Polynomial > m_polynomials;
};
}
#endif
