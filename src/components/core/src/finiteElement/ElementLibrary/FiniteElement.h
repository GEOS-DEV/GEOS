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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef FINITE_ELEMENT_H
#define FINITE_ELEMENT_H

/**
 * @file FiniteElement.h
 * @author white230
 */

#include "finiteElement/basis/BasisBase.hpp"
#include "finiteElement/quadrature/QuadratureBase.hpp"
#include "FiniteElementBase.h"

/**
 * Class representing a generic finite element.  Its constructor
 * takes a specific interpolation basis and quadrature rule in
 * order to define a complete element.
 *
 * The class assume that the
 * mapping from parent coordinates to real coordinates is
 * iso-parametric, and therefore the same basis is used for both
 * interpolation and mapping.
 *
 * The class also defines a generic interface for accessing finite
 * element data.  In the future, more sophisticated element
 * definitions that do not fit within the current class can be
 * defined through derived classes.
 */
namespace geosx
{
template <int dim>
class FiniteElement : public FiniteElementBase
{
public:

  FiniteElement( const int num_q_points,
                 const int num_dofs,
                 const int num_zero_energy_modes = 0 );

  FiniteElement(BasisBase const & basis,
                QuadratureBase const & quadrature,
                const int num_zero_energy_modes = 0 );

  virtual ~FiniteElement(){}


  virtual void reinit(const array<R1TensorT<3> > &mapped_support_points);



};



template <int dim>
FiniteElement<dim> :: FiniteElement(const int num_q_points,
                                    const int num_dofs,
                                    const int num_zero_energy_modes):
  FiniteElementBase( dim, num_q_points, num_dofs, num_zero_energy_modes)
{}


/**
 * Constructor.  Takes an interpolation basis and quadrature rule,
 * and pre-computes all static finite element data.  Any data
 * that depends on the mapped configuration of the element, however,
 * is left uninitialized until reinit(...) is called.
 */

template <int dim>
FiniteElement<dim> :: FiniteElement(BasisBase const &basis,
                                    QuadratureBase const &quadrature,
                                    const int num_zero_energy_modes ):
  FiniteElementBase( dim, quadrature.size(), basis.size(), num_zero_energy_modes)
{

  data.resize(n_q_points);
  for(auto q=0 ; q<n_q_points ; ++q)
  {
    data[q].parent_q_point = quadrature.integration_point(q);
    data[q].parent_q_weight = quadrature.integration_weight(q);

    data[q].parent_values.resize(n_dofs);
    data[q].parent_gradients.resize(n_dofs);
    data[q].mapped_gradients.resize(n_dofs);

    for(auto i=0 ; i<n_dofs ; ++i)
    {
      data[q].parent_values[i]    = basis.value(i,data[q].parent_q_point);
      data[q].parent_gradients[i] = basis.gradient(i,data[q].parent_q_point);
    }
  }
}



/**
 * Reinitialize the finite element basis on a particular element.
 * We use the coordinates of the support points in real space to
 * construct the forward mapping from the parent coordinate system.  The
 * support points are assumed to follow a lexicographic ordering:
 * On the parent element, we loop over the x-coordinate fastest,
 * the y, then z (depending on the desired spatial dimension of the
 * element).
 */

template <int dim>
void FiniteElement<dim> :: reinit(const array<R1TensorT<3> > &mapped_support_points)
{
  assert(mapped_support_points.size() == n_dofs);

  R2TensorT<3> jacobian;
  R2TensorT<3> inv_jacobian;

  for(auto q=0 ; q<n_q_points ; ++q)
  {

    jacobian = 0;
    for(auto a=0 ; a<n_dofs ; ++a)
    {
      jacobian.plus_dyadic_ab( mapped_support_points[a], data[q].parent_gradients[a] );
    }

    if( dim==2 )
    {
      jacobian(2,2) = 1;
    }

    data[q].jacobian_determinant = inv_jacobian.Inverse(jacobian);

    for(auto i=0 ; i<n_dofs ; ++i)
    {
      data[q].mapped_gradients[i].AijBi( inv_jacobian, data[q].parent_gradients[i] );
    }
  }
}


}



#endif
