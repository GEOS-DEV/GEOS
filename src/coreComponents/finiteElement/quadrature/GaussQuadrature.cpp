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

#include "finiteElement/quadrature/GaussQuadrature.hpp"
#include "meshUtilities/StructuredGridUtilities.hpp"
/*
 * Constructor. Compute weights and unit integration points
 * for a Gauss-Legendre quadrature rule of degree p > 0, where
 * p corresponds to the number of integration points in each coordinate
 * direction.
 */

/*
 * Destructor.
 */
namespace geosx
{

using namespace dataRepository;

template< int dim >
GaussQuadrature< dim >::GaussQuadrature( std::string const & name, Group * const parent )
  :
  QuadratureBase( name, parent ),
  m_degree( 0 ),
  m_n_gauss_points( 0 )
{
  registerWrapper( viewKeyStruct::degreeString, &m_degree )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Quadrature degree" );
}


template< int dim >
GaussQuadrature< dim >::~GaussQuadrature()
{}


/*
 * Get number of integration points.
 */

template< int dim >
int GaussQuadrature< dim >::size() const
{
  return m_n_gauss_points;
}

/*
 * Get integration point on unit cell.
 */

template< int dim >
R1Tensor GaussQuadrature< dim >::integration_point( const int index ) const
{
  std::vector< int > indices( dim );
  StructuredGrid::map_index< dim >( index, m_degree, indices );

  R1Tensor point;
  for( int d = 0; d < dim; ++d )
    point[d] = m_points_1d[indices[d]];

  return point;
}

/*
 * Get integration weight value on unit cell
 */

template< int dim >
double GaussQuadrature< dim >::integration_weight( const int index ) const
{
  std::vector< int > indices( dim );
  StructuredGrid::map_index< dim >( index, m_degree, indices );

  double weight = 1.0;
  for( int d = 0; d < dim; ++d )
    weight *= m_weights_1d[indices[d]];

  return weight;
}

template< int dim >
void GaussQuadrature< dim >::PostProcessInput()
{
  m_n_gauss_points = StructuredGrid::dimpower< dim >( m_degree );

  GEOSX_ASSERT_GT( m_degree, 0 );

  m_points_1d.resize( m_degree );
  m_weights_1d.resize( m_degree );

  const double tolerance = 5e-16;

  const int m = ( m_degree + 1 ) / 2;

  for( int i = 1; i <= m; ++i )
  {
    double z = std::cos( M_PI * ( i - 0.25 ) / ( m_degree + 0.5 ) );

    double pp, p1, p2, p3;

    do
    {
      p1 = 1.0;
      p2 = 0.0;
      for( int j = 0; j < m_degree; ++j )
      {
        p3 = p2;
        p2 = p1;
        p1 = ( ( 2.0 * j + 1.0 ) * z * p2 - j * p3 ) / ( j + 1 );
      }
      pp = m_degree * ( z * p1 - p2 ) / ( z * z - 1 );
      z = z - p1 / pp;
    }
    while( std::fabs( p1 / pp ) > tolerance );

    double x = 0.5 * z;

    m_points_1d[i - 1] = 0.5 - x;
    m_points_1d[m_degree - i] = 0.5 + x;

    double w = 1.0 / ( ( 1.0 - z * z ) * pp * pp );

    m_weights_1d[i - 1] = w;
    m_weights_1d[m_degree - i] = w;
  }
}

namespace
{
dataRepository::CatalogEntryConstructor< QuadratureBase, GaussQuadrature< 1 >, std::string const &, Group * const > catEntry_GaussQuadrature1;
dataRepository::CatalogEntryConstructor< QuadratureBase, GaussQuadrature< 2 >, std::string const &, Group * const > catEntry_GaussQuadrature2;
dataRepository::CatalogEntryConstructor< QuadratureBase, GaussQuadrature< 3 >, std::string const &, Group * const > catEntry_GaussQuadrature3;
}

}
