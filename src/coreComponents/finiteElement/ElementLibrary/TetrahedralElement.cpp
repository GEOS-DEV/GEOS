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
 * @file Tetrahedron.cpp
 */

#include "TetrahedralElement.hpp"

#include "finiteElement/basis/LagrangeBasis.hpp"
#include "finiteElement/quadrature/GaussQuadrature.hpp"

namespace geosx
{

TetrahedralElement::TetrahedralElement( BasisBase const & GEOSX_UNUSED_PARAM( basis ),
                                        QuadratureBase const & GEOSX_UNUSED_PARAM( quadrature ),
                                        const int GEOSX_UNUSED_PARAM( num_zero_energy_modes ) ):
//  FiniteElementBase( dim, quadrature.size(), basis.size(), num_zero_energy_modes)
  FiniteElementBase( 3, 1, 4, 0)
{
  m_nodeOrdering.resize(4);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 2;
  m_nodeOrdering[3] = 3;

}

TetrahedralElement::~TetrahedralElement()
{}


REGISTER_CATALOG_ENTRY( FiniteElementBase, TetrahedralElement, BasisBase const &, QuadratureBase const &, const int )

}
