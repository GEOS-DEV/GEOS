/*
 * FiniteElementBase.h
 *
 *  Created on: Nov 20, 2012
 *      Author: settgast
 */

#ifndef FINITEELEMENTBASE_H_
#define FINITEELEMENTBASE_H_

#include "Common/Common.h"
#include <assert.h>

class FiniteElementBase
{
public:
  FiniteElementBase( const int dim,
                     const int num_q_points,
                     const int num_dofs,
                     const int num_zero_energy_modes );

  virtual ~FiniteElementBase();



  virtual void reinit( const std::vector<R1TensorT<3> > &mapped_support_points) = 0;


//  virtual void zero_energy_mode_control( const Array1dT<R1Tensor>& dNdx,
//                                         const realT& volume,
//                                         const Array1dT<R1Tensor>& x,
//                                         const Array1dT<R1Tensor>& vel,
//                                         const realT& dampcoef,
//                                         const realT& stiffcoef,
//                                         const realT& rho,
//                                         const realT& modulus,
//                                         const realT& dt,
//                                         Array1dT<R1Tensor>& Qstiffness,
//                                         Array1dT<R1Tensor>& force ) {}

  virtual void zero_energy_mode_control( const Array1dT<R1Tensor>& ,
                                         const realT& ,
                                         const Array1dT<R1Tensor>& ,
                                         const Array1dT<R1Tensor>& ,
                                         const realT& ,
                                         const realT& ,
                                         const realT& ,
                                         const realT& ,
                                         const realT& ,
                                         Array1dT<R1Tensor>& ,
                                         Array1dT<R1Tensor>&  ) {}


  double value(const unsigned shape_index,
               const unsigned q_index) const
  {
    assert(q_index < n_q_points);
    assert(shape_index < n_dofs);
    return data[q_index].parent_values[shape_index];
  }

  std::vector<double> const & values( const unsigned int q_index ) const
  {
    assert(q_index < n_q_points);
    return data[q_index].parent_values;
  }

  R1Tensor gradient( const unsigned shape_index,
                 const unsigned q_index ) const
  {
    assert(q_index < n_q_points);
    assert(shape_index < n_dofs);
    return data[q_index].mapped_gradients[shape_index];
  }

  double JxW(const unsigned q_index) const
  {
    assert(q_index < n_q_points);
    return data[q_index].jacobian_determinant *
           data[q_index].parent_q_weight;
  }


  int Dim()                             { return m_dim; }
  unsigned n_quadrature_points() const  { return n_q_points;  }
  unsigned dofs_per_element() const     { return n_dofs;  }
  inline int zero_energy_modes() const  { return m_zero_energy_modes; }

  std::string m_type;

protected:
  iArray1d m_nodeOrdering;
  unsigned n_q_points;
  unsigned n_dofs;
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

#endif /* FINITEELEMENTBASE_H_ */
