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

/*
 * ThickPlane.cpp
 */

#include "ThickPlane.hpp"

namespace geosx
{
using namespace dataRepository;

ThickPlane::ThickPlane( const std::string& name, ManagedGroup * const parent ):
  SimpleGeometricObjectBase( name, parent ),
  m_origin{0.0,0.0,0.0},
  m_normal{0.0,0.0,1.0},
  m_thickness{0.0}
{
  RegisterViewWrapper( viewKeyStruct::originString, &m_origin, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Origin point (x,y,z) of the plane (basically, any point on the plane)");

  RegisterViewWrapper( viewKeyStruct::normalString, &m_normal, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)");

  RegisterViewWrapper( viewKeyStruct::thicknessString, &m_thickness, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("The total thickness of the plane (with half to each side)");
}

ThickPlane::~ThickPlane()
{}


void ThickPlane::PostProcessInput()
{
  m_thickness *= 0.5; // actually store the half-thickness
  GEOS_ERROR_IF(m_thickness <= 0, "Error: the plane appears to have zero or negative thickness");

  m_normal.Normalize();
  GEOS_ERROR_IF(std::fabs(m_normal.L2_Norm()-1.0) > 1e-15, "Error: could not properly normalize input normal.");
}


bool ThickPlane::IsCoordInObject( const R1Tensor& coord ) const
{
  real64 normalDistance = 0.0;
  for(int i=0; i<3; ++i)
  {
    normalDistance += m_normal[i]*(coord[i]-m_origin[i]);
  }

  return std::fabs(normalDistance) <= m_thickness;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, ThickPlane, std::string const &, ManagedGroup * const )

} /* namespace geosx */

