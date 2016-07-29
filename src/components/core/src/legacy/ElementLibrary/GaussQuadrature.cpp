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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "ElementLibrary/GaussQuadrature.h"

/*
 * Constructor. Compute weights and unit integration points
 * for a Gauss-Legendre quadrature rule of degree p > 0, where
 * p corresponds to the number of integration points in each coordinate
 * direction. 
 */

template <int dim>
GaussQuadrature<dim> :: GaussQuadrature(const unsigned degree):
  m_degree(degree),
  m_n_gauss_points(StructuredGrid::dimpower<dim>(degree))
{
  assert(degree > 0);

  m_points_1d.resize(degree);
  m_weights_1d.resize(degree);

  const double tolerance = 5e-16;

  const unsigned m = (degree+1)/2;

  for(unsigned i=1; i<=m; ++i)
  {
    double z = std::cos(M_PI*(i-0.25)/(degree+0.5));

    double pp,p1,p2,p3;

    do
    {
      p1 = 1.0;
      p2 = 0.0;
      for (unsigned int j=0;j<degree;++j)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*j+1.0)*z*p2-j*p3)/(j+1);
      }
      pp = degree*(z*p1-p2)/(z*z-1);
      z = z-p1/pp;
    }
    while(std::fabs(p1/pp) > tolerance);

    double x = 0.5*z;

    m_points_1d[i-1]      = 0.5-x;
    m_points_1d[degree-i] = 0.5+x;

    double w = 1.0/((1.0-z*z)*pp*pp);

    m_weights_1d[i-1]      = w;
    m_weights_1d[degree-i] = w;
  }
}


/*
 * Destructor.
 */

template <int dim>
GaussQuadrature<dim> :: ~GaussQuadrature()
{}


/*
 * Get number of integration points.
 */

template <int dim>
unsigned GaussQuadrature<dim> :: size()
{
  return m_n_gauss_points;
}


/*
 * Get integration point on unit cell.
 */ 

template <int dim>
R1Tensor GaussQuadrature<dim> :: integration_point(const unsigned index)
{
  std::vector<unsigned> indices(dim);
  StructuredGrid::map_index<dim>(index,m_degree,indices);

  R1Tensor point;
  for(unsigned d=0; d<dim; ++d)
    point[d] = m_points_1d[indices[d]];

  return point;
}


/* 
 * Get integration weight value on unit cell
 */

template <int dim>
double GaussQuadrature<dim> :: integration_weight(const unsigned index)
{
  std::vector<unsigned> indices(dim);
  StructuredGrid::map_index<dim>(index,m_degree,indices);

  double weight = 1.0;
  for(unsigned d=0; d<dim; ++d)
    weight *= m_weights_1d[indices[d]];

  return weight;
}


/*
 * Explicit instantiations.
 */

template class GaussQuadrature<1>;
template class GaussQuadrature<2>;
template class GaussQuadrature<3>;

