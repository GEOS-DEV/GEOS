/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file RTheta.cpp
 */

#include "RTheta.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
using namespace dataRepository;

RTheta::RTheta( const string & name, Group * const parent ):
  BoundedPlanarObject( name, parent ),
  m_center{ 0.0, 0.0, 0.0 },
  m_tolerance()
{
  registerWrapper( viewKeyStruct::centerString(), &m_center ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "(x,y,z) coordinates of the center of the RTheta" );

  // registerWrapper( viewKeyStruct::radiusString(), &m_radius ).
  //   setInputFlag( InputFlags::OPTIONAL ).
  //   setDescription( "Coefficients of the RTheta function." );

  registerWrapper( viewKeyStruct::toleranceString(), &m_tolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 1e-5 ).
    setDescription( "Tolerance to determine if a point sits on the RTheta or not. "
                    "It is relative to the maximum dimension of the RTheta." );

}

RTheta::~RTheta()
{}

void RTheta::postProcessInput()
{
  // Make sure that you have an orthonormal basis.
  LvArray::tensorOps::normalize< 3 >( m_normal );

}

bool RTheta::isCoordInObject( real64 const ( &coord ) [3] ) const
{
  bool isInside = true;

  //gets the vector from coord to the center
  real64 dummy[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( coord );
  LvArray::tensorOps::subtract< 3 >( dummy, m_center );

  // 1. Check if point is on the plane of the RTheta
  if( std::abs( LvArray::tensorOps::AiBi< 3 >( dummy, m_normal ) ) > m_tolerance )
  {
    isInside = false;
  }

  // 2. Get angle with the plane's length vector
  real64 dummyNorm = sqrt(LvArray::tensorOps::l2NormSquared< 3 >( dummy ));
  real64 cosTheta = LvArray::tensorOps::AiBi< 3 >( dummy, getLengthVector() )/(dummyNorm + 1e-15); //assume lengthVector is unitary
  real64 theta = acos(cosTheta);

  // 2. Check if it is inside the RTheta
  if( LvArray::tensorOps::l2NormSquared< 3 >( dummy ) > getRadius(theta)*getRadius(theta) )
  {
    isInside = false;
  }

  return isInside;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, RTheta, string const &, Group * const )

} /* namespace geosx */
