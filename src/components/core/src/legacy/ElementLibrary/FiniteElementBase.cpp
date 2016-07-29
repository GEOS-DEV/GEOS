/*
 * FiniteElementBase.cpp
 *
 *  Created on: Nov 20, 2012
 *      Author: settgast
 */

#include "FiniteElementBase.h"

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
  for(unsigned q=0; q<n_q_points; ++q)
  {
//    data[q].parent_q_point = quadrature.integration_point(q);
    data[q].parent_q_weight = 1;

    data[q].parent_values.resize(n_dofs);
    data[q].parent_gradients.resize(n_dofs);
    data[q].mapped_gradients.resize(n_dofs);

  }
}

FiniteElementBase::~FiniteElementBase()
{
  // TODO Auto-generated destructor stub
}

