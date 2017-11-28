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
#include "LagrangeBasis.hpp"

#include "MeshUtilities/StructuredGridUtilities.hpp"
/*
 * Constructor.
 */
namespace geosx
{

template <int dim>
LagrangeBasis<dim> :: LagrangeBasis(const int degree)
  :
  m_degree(degree),
  n_shape_functions(StructuredGrid::dimpower<dim>(degree+1))
{
  std::vector<std::vector<double> > coeff( degree+1, std::vector<double>(degree+1) );

  switch(degree)
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
    assert(degree<5);
  }

  for(int n=0 ; n<degree+1 ; ++n)
    m_polynomials.push_back(Polynomial(coeff[n]));
}


/*
 * Number of degrees of freedom in basis
 */

template <int dim>
int LagrangeBasis<dim> :: size() const
{
  return n_shape_functions;
}


/*
 * Evaluate the basis at a particular point in parent coordinates
 */


template <int dim>
double LagrangeBasis<dim> :: value(const int       index,
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
R1Tensor LagrangeBasis<dim> :: gradient(const int index,
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
R1Tensor LagrangeBasis<dim> :: support_point(const int index)
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
void LagrangeBasis<dim>::ReadXML( xmlWrapper::xmlNode const & targetNode )
{
  m_degree = targetNode.attribute("degree").as_int(1);
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
namespace { cxx_utilities::CatalogEntryConstructor<BasisBase,LagrangeBasis<1> > catEntry_LagrangeBasis1; }
namespace { cxx_utilities::CatalogEntryConstructor<BasisBase,LagrangeBasis<2> > catEntry_LagrangeBasis2; }
namespace { cxx_utilities::CatalogEntryConstructor<BasisBase,LagrangeBasis<3> > catEntry_LagrangeBasis3; }
}
