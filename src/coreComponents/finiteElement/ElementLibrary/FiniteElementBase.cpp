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
 * @file FiniteElementBase.cpp
 */

#include "FiniteElementBase.h"

namespace geosx
{
FiniteElementBase::FiniteElementBase( const int dim,
                                      const int num_q_points,
                                      const int num_dofs,
                                      const int num_zero_energy_modes ):
  m_nodeOrdering(),
  n_q_points(num_q_points),
  n_dofs(num_dofs),
  m_zero_energy_modes(num_zero_energy_modes),
  m_dim(dim)
{
  data.resize(n_q_points);
  for(auto q=0 ; q<n_q_points ; ++q)
  {
//    data[q].parent_q_point = quadrature.integration_point(q);
    data[q].parent_q_weight = 1;

    data[q].parent_values.resize(n_dofs);
    data[q].parent_gradients.resize(n_dofs, 3);
    data[q].mapped_gradients.resize(n_dofs, 3);

  }
}

FiniteElementBase::~FiniteElementBase()
{
  // TODO Auto-generated destructor stub
}

}
