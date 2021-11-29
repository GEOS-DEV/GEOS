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
 * @file Cylinder.cpp
 * @brief Generate a cylindrical geometry.
 * @param Center point of one (upper or lower) face of the cylinder
 * @param Center point of the other face of the cylinder
 * @param Radius of the cylinder
 */

#include "Cylinder.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
using namespace dataRepository;

Cylinder::Cylinder( const string & name, Group * const parent ):
  SimpleGeometricObjectBase( name, parent ),
  m_point1{ 0.0, 0.0, 0.0 },
  m_point2{ 0.0, 0.0, 0.0 },
  m_radius{ 0.0 }
{
  registerWrapper( viewKeyStruct::point1String(), &m_point1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Center point of one (upper or lower) face of the cylinder" );

  registerWrapper( viewKeyStruct::point2String(), &m_point2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Center point of the other face of the cylinder" );

  registerWrapper( viewKeyStruct::radiusString(), &m_radius ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Radius of the cylinder" );

  registerWrapper( viewKeyStruct::innerRadiusString(), &m_innerRadius ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Inner radius of the anulus" );

}

Cylinder::~Cylinder()
{}


bool Cylinder::isCoordInObject( real64 const ( &coord ) [3] ) const
{
  bool rval = false;

  real64 axisVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_point2 );
  LvArray::tensorOps::subtract< 3 >( axisVector, m_point1 );
  real64 const height = LvArray::tensorOps::normalize< 3 >( axisVector );

  real64 coord_minus_point1[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( coord );
  LvArray::tensorOps::subtract< 3 >( coord_minus_point1, m_point1 );

  real64 projection[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( axisVector );
  LvArray::tensorOps::scale< 3 >( projection, LvArray::tensorOps::AiBi< 3 >( axisVector, coord_minus_point1 ) );

  real64 distance[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( coord_minus_point1 );
  LvArray::tensorOps::subtract< 3 >( distance, projection );

  if( LvArray::tensorOps::l2Norm< 3 >( distance ) < m_radius &&
      LvArray::tensorOps::l2Norm< 3 >( distance ) >= m_innerRadius &&
      LvArray::tensorOps::l2Norm< 3 >( projection ) < height )
  {
    rval = true;
  }

  return rval;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, Cylinder, string const &, Group * const )

} /* namespace geosx */
