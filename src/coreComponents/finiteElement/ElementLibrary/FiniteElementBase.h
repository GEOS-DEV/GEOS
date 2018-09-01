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

/*
 * FiniteElementBase.h
 *
 *  Created on: Nov 20, 2012
 *      Author: settgast
 */

#ifndef FINITEELEMENTBASE_H_
#define FINITEELEMENTBASE_H_

#include <assert.h>
#include "common/DataTypes.hpp"

namespace geosx
{
class FiniteElementBase
{
public:
  FiniteElementBase( const int dim,
                     const int num_q_points,
                     const int num_dofs,
                     const int num_zero_energy_modes );

  virtual ~FiniteElementBase();



  virtual void reinit( const array1d<R1TensorT<3> > &mapped_support_points) = 0;


//  virtual void zero_energy_mode_control( const array1d<R1Tensor>& dNdx,
//                                         const realT& volume,
//                                         const array1d<R1Tensor>& x,
//                                         const array1d<R1Tensor>& vel,
//                                         const realT& dampcoef,
//                                         const realT& stiffcoef,
//                                         const realT& rho,
//                                         const realT& modulus,
//                                         const realT& dt,
//                                         array1d<R1Tensor>& Qstiffness,
//                                         array1d<R1Tensor>& force ) {}

  virtual void zero_energy_mode_control( const array1d<R1Tensor>&,
                                         const realT&,
                                         const array1d<R1Tensor>&,
                                         const array1d<R1Tensor>&,
                                         const realT&,
                                         const realT&,
                                         const realT&,
                                         const realT&,
                                         const realT&,
                                         array1d<R1Tensor>&,
                                         array1d<R1Tensor>&  ) {}


  double value(const int shape_index,
               const int q_index) const
  {
    assert(q_index < n_q_points);
    assert(shape_index < n_dofs);
    return data[q_index].parent_values[shape_index];
  }

  std::vector<double> const & values( const int q_index ) const
  {
    assert(q_index < n_q_points);
    return data[q_index].parent_values;
  }

  R1Tensor gradient( const localIndex shape_index,
                     const localIndex q_index ) const
  {
    assert(q_index < n_q_points);
    assert(shape_index < n_dofs);
    return data[q_index].mapped_gradients[shape_index];
  }

  double JxW(const localIndex q_index) const
  {
    assert(q_index < n_q_points);
    return data[q_index].jacobian_determinant *
           data[q_index].parent_q_weight;
  }


  int Dim()                             { return m_dim; }
  int n_quadrature_points() const  { return n_q_points;  }
  int dofs_per_element() const     { return n_dofs;  }
  inline int zero_energy_modes() const  { return m_zero_energy_modes; }

  std::string m_type;

protected:
  array1d<integer> m_nodeOrdering;
  int n_q_points;
  int n_dofs;
  int m_zero_energy_modes;

  struct QuadraturePointData
  {
    R1TensorT<3>                 parent_q_point;
    double                       parent_q_weight;
    std::vector<double>          parent_values;
    std::vector<R1TensorT<3> >   parent_gradients;

    std::vector<R1TensorT<3> >   mapped_gradients;
    double                       jacobian_determinant;
  };

  std::vector<QuadraturePointData> data;

private:
  int m_dim;


  FiniteElementBase();
  FiniteElementBase( const FiniteElementBase& );


};
}
#endif /* FINITEELEMENTBASE_H_ */
