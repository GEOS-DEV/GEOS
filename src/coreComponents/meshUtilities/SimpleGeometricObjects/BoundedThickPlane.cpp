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

/**
 * @file BoundedThickPlane.cpp
 */

#include "BoundedThickPlane.hpp"

namespace geosx
{
using namespace dataRepository;

BoundedThickPlane::BoundedThickPlane( const std::string& name, Group * const parent ):
    ThickPlane( name, parent ),
  m_lengthVector{0.0,0.0,0.0},
  m_widthVector{0.0,0.0,0.0}
{
  registerWrapper( viewKeyStruct::mLengthVectorString, &m_lengthVector, false )->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Tangent vector defining the orthonormal basis along with the normal.");

  registerWrapper( viewKeyStruct::mWidthVectorString, &m_widthVector, false )->
       setInputFlag(InputFlags::REQUIRED)->
       setDescription("Tangent vector defining the orthonormal basis along with the normal.");

  registerWrapper( viewKeyStruct::dimensionsString, &m_dimensions, false )->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Length and width of the bounded plane");
}

BoundedThickPlane::~BoundedThickPlane()
{}

void BoundedThickPlane::PostProcessInput()
{
  m_thickness *= 0.5; // actually store the half-thickness
  GEOS_ERROR_IF(m_thickness <= 0, "Error: the plane appears to have zero or negative thickness");

  // Make sure that you have an orthonormal basis.
  m_normal.Normalize();
  GEOS_ERROR_IF(std::fabs(m_normal.L2_Norm()-1.0) > 1e-15, "Error: could not properly normalize input normal.");
  m_lengthVector.Normalize();
  GEOS_ERROR_IF(std::fabs(m_lengthVector.L2_Norm()-1.0) > 1e-15, "Error: could not properly normalize input lengthVector.");
  m_widthVector.Normalize();
  GEOS_ERROR_IF(std::fabs(m_widthVector.L2_Norm()-1.0) > 1e-15, "Error: could not properly normalize input widthVector.");

  //check if they are all orthogonal
  R1Tensor vector = Cross(m_lengthVector, m_widthVector);
  real64 scalarProd = Dot(m_normal, vector);
  GEOS_ERROR_IF(std::fabs(scalarProd) - 1 > 1e-15, "Error: the 3 vectors provided do not form an orthonormal basis!");
  GEOS_ERROR_IF(m_dimensions.size() != 2, "Error: Need to provide both length and width!");
}


bool BoundedThickPlane::IsCoordInObject( const R1Tensor& coord ) const
{
  real64 normalDistance = 0.0;
  for(int i=0; i<3; ++i)
  {
    normalDistance += m_normal[i]*(coord[i]-m_origin[i]);
  }

  return std::fabs(normalDistance) <= m_thickness;


}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, BoundedThickPlane, std::string const &, Group * const )

} /* namespace geosx */

