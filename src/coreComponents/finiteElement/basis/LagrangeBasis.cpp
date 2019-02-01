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

#include "LagrangeBasis.hpp"

#include "meshUtilities/StructuredGridUtilities.hpp"
/*
 * Constructor.
 */
namespace geosx
{

using namespace dataRepository;

template <int dim>
LagrangeBasis<dim>::LagrangeBasis(std::string const & name, ManagedGroup * const parent)
  :
  BasisBase(name, parent),
  m_degree(0),
  n_shape_functions(0)
{
  RegisterViewWrapper( viewKeyStruct::degreeString, &m_degree, 0 )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Basis degree");
}

template <int dim>
LagrangeBasis<dim>::~LagrangeBasis()
{}

/*
 * Number of degrees of freedom in basis
 */

template <int dim>
int LagrangeBasis<dim>::size() const
{
  return n_shape_functions;
}


/*
 * Evaluate the basis at a particular point in parent coordinates
 */
template <int dim>
double LagrangeBasis<dim>::value(const int       index,
                                   const R1Tensor &point) const
{
  std::vector<int> indices(dim);
  StructuredGrid::map_index<dim>(index,m_degree+1,indices);

  double val = 1;
  for(int d=0 ; d<dim ; ++d)
    val *= m_polynomials[indices[d]].Value(point[d]);

  return val;
}


/*
 * Evaluate the gradient at a particular point in parent coordinates.
 */
template <int dim>
R1Tensor LagrangeBasis<dim>::gradient(const int index,
                                        const R1Tensor &point) const
{
  std::vector<int> indices(dim);
  StructuredGrid::map_index<dim>(index,m_degree+1,indices);

  R1Tensor grad;

  for(int r=0 ; r<dim ; ++r)
    grad[r] = 1.0;

  double pvalue,pderiv;

  for(int c=0 ; c<dim ; ++c)
  {
    pvalue = m_polynomials[indices[c]].Value(point[c]);
    pderiv = m_polynomials[indices[c]].Deriv(point[c]);
    for(int r=0 ; r<dim ; ++r)
    {
      grad[r] *= ( r==c ? pderiv : pvalue );
    }
  }

  return grad;
}


/*
 * Get the support point for a particular basis function.  For the
 * Lagrange basis, these points are equidistant within the parent
 * element.
 */

template <int dim>
R1Tensor LagrangeBasis<dim>::support_point(const int index)
{
  R1Tensor pt;

  if(m_degree == 0)
  {
    for(int i=0 ; i<dim ; ++i)
      pt(i) = 0.5;
  }
  else
  {
    std::vector<int> indices(dim);
    StructuredGrid::map_index<dim>(index,m_degree+1,indices);

    const double h = 1.0/m_degree;

    for(int i=0 ; i<dim ; ++i)
      pt(i) = h*indices[i];
  }

  return pt;
}

template <int dim>
void LagrangeBasis<dim>::PostProcessInput()
{
  n_shape_functions = StructuredGrid::dimpower<dim>(m_degree+1);

  std::vector<std::vector<double> > coeff( m_degree+1, std::vector<double>(m_degree+1) );

  switch(m_degree)
  {
  case 0:
    coeff[0][0] =  1.0;
    break;

  case 1:
    coeff[0][0] =  1.0;
    coeff[0][1] = -1.0;

    coeff[1][0] =  0.0;
    coeff[1][1] =  1.0;
    break;

  case 2:
    coeff[0][0] =  1.0;
    coeff[0][1] = -3.0;
    coeff[0][2] =  2.0;

    coeff[1][0] =  0.0;
    coeff[1][1] =  4.0;
    coeff[1][2] = -4.0;

    coeff[2][0] =  0.0;
    coeff[2][1] = -1.0;
    coeff[2][2] =  2.0;
    break;

  case 3:
    coeff[0][0] =  1.0;
    coeff[0][1] = -11.0/2.0;
    coeff[0][2] =  9.0;
    coeff[0][3] = -9.0/2.0;

    coeff[1][0] =  0.0;
    coeff[1][1] =  9.0;
    coeff[1][2] = -45.0/2.0;
    coeff[1][3] =  27.0/2.0;

    coeff[2][0] =  0.0;
    coeff[2][1] = -9.0/2.0;
    coeff[2][2] =  18.0;
    coeff[2][3] = -27.0/2.0;

    coeff[3][0] =  0.0;
    coeff[3][1] =  1.0;
    coeff[3][2] = -9.0/2.0;
    coeff[3][3] =  9.0/2.0;
    break;

  case 4:
    coeff[0][0] =  1.0;
    coeff[0][1] = -25.0/3.0;
    coeff[0][2] =  70.0/3.0;
    coeff[0][3] = -80.0/3.0;
    coeff[0][4] =  32.0/3.0;

    coeff[1][0] =  0.0;
    coeff[1][1] =  16.0;
    coeff[1][2] = -208.0/3.0;
    coeff[1][3] =  96.0;
    coeff[1][4] = -128.0/3.0;

    coeff[2][0] =  0.0;
    coeff[2][1] = -12.0;
    coeff[2][2] =  76.0;
    coeff[2][3] = -128.0;
    coeff[2][4] =  64.0;

    coeff[3][0] =  0.0;
    coeff[3][1] =  16.0/3.0;
    coeff[3][2] = -112.0/3.0;
    coeff[3][3] =  224.0/3.0;
    coeff[3][4] = -128.0/3.0;

    coeff[4][0] =  0.0;
    coeff[4][1] = -1.0;
    coeff[4][2] =  22.0/3.0;
    coeff[4][3] = -16.0;
    coeff[4][4] =  32.0/3.0;
    break;

  default:
    assert(m_degree<5);
  }

  for(int n=0 ; n<m_degree+1 ; ++n)
    m_polynomials.push_back(Polynomial(coeff[n]));
}

/*
 * Explicit instantiations.
 */

//template class LagrangeBasis<1>;
//template class LagrangeBasis<2>;
//template class LagrangeBasis<3>;

//REGISTER_CATALOG_ENTRY( BasisBase, LagrangeBasis<1>,void )
namespace { cxx_utilities::CatalogEntryConstructor<BasisBase, LagrangeBasis<1>, std::string const &, ManagedGroup * const > catEntry_LagrangeBasis1; }
namespace { cxx_utilities::CatalogEntryConstructor<BasisBase, LagrangeBasis<2>, std::string const &, ManagedGroup * const > catEntry_LagrangeBasis2; }
namespace { cxx_utilities::CatalogEntryConstructor<BasisBase, LagrangeBasis<3>, std::string const &, ManagedGroup * const > catEntry_LagrangeBasis3; }
}
