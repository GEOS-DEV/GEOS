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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * Box.cpp
 *
 *  Created on: Aug 4, 2017
 *      Author: settgast
 */

#include "Box.hpp"

namespace geosx
{
using namespace dataRepository;

Box::Box( const std::string& name, ManagedGroup * const parent ):
  SimpleGeometricObjectBase( name, parent ),
  m_min{0.0,0.0,0.0},
  m_max{0.0,0.0,0.0},
  m_strikeAngle{0.0},
  m_boxCenter{0.0,0.0,0.0},
  m_cosStrike{0.0},
  m_sinStrike{0.0}
{
  RegisterViewWrapper( viewKeyStruct::xMinString, &m_min, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Minimum (x,y,z) coordinates of the box");

  RegisterViewWrapper( viewKeyStruct::xMaxString, &m_max, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Maximum (x,y,z) coordinates of the box");

  RegisterViewWrapper( viewKeyStruct::strikeAngleString, &m_strikeAngle, false )->
    setApplyDefaultValue(-90.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("The strike angle of the box");

  RegisterViewWrapper( viewKeyStruct::boxCenterString, &m_boxCenter, false );
  RegisterViewWrapper( viewKeyStruct::cosStrikeString, &m_cosStrike, false );
  RegisterViewWrapper( viewKeyStruct::sinStrikeString, &m_sinStrike, false );
}

Box::~Box()
{}




void Box::PostProcessInput()
{
  m_boxCenter = m_min;
  m_boxCenter += m_max;
  m_boxCenter *= 0.5;

  m_strikeAngle += 90; // Counterclockwise from x-axis
  if (std::fabs(m_strikeAngle) > 1e-20)
  {
    GEOS_ERROR_IF( (m_max[0]-m_min[0]) < (m_max[1]-m_min[1]),
                   "Error: When a strike angle is specified, the box is supposed to represent a plane normal to the "
                   "y direction. This box seems to be too thick.");

    m_cosStrike = std::cos(m_strikeAngle / 180 *M_PI);
    m_sinStrike = std::sin(m_strikeAngle / 180 *M_PI);
  }

}

//void Box::ReadXML( xmlWrapper::xmlNode const & xmlNode )
//{
//  pugi::xml_attribute xmlatt =
// targetNode.attribute(subDocNode.getStringKey().c_str());
//  m_degree = targetNode.attribute("degree").as_int(1);
//  as_type( xmlVal, xmlatt.value(), defVal );
//
//  m_min = xmlWrapper::as_type( xmlNode, "xMin", {-1e99,-1e99,-1e99});
//  m_max = xmlWrapper::as_type( xmlNode, "xMax", { 1e99, 1e99, 1e99});
//  m_strikeAngle = hdn.GetAttributeOrDefault<realT>("strikeAngle", -90.0); //
// from North
//  m_strikeAngle += 90; // Counterclockwise from x-axis
//  if (std::fabs(m_strikeAngle) > 1e-20)
//  {
//    if ((m_max[0]-m_min[0]) < (m_max[1]-m_min[1]))
//      throw GPException("Error: When a strike angle is specified, the box is
// supposed to represent a plane normal to the y direction. This box seems to be
// too thick.");
//
//    m_cosStrike = std::cos(m_strikeAngle / 180 *3.1415926535);
//    m_sinStrike = std::sin(m_strikeAngle / 180 *3.1415926535);
//    m_boxCenter = m_min;
//    m_boxCenter += m_max;
//    m_boxCenter *= 0.5;
//  }
//
// }

bool Box::IsCoordInObject( const R1Tensor& coord ) const
{
  bool rval = false;
  if (std::fabs(m_strikeAngle) < 1e-20)
  {
    if( coord <= m_max && coord >= m_min )
    {
      rval = true;
    }
  }
  else
  {
    R1Tensor coordR, coord0(coord);
    coord0 -= m_boxCenter;
    coordR[0] = coord0[0] * m_cosStrike + coord0[1] * m_sinStrike;
    coordR[1] = -coord0[0] * m_sinStrike + coord0[1] * m_cosStrike;
    coordR[2] = coord0[2];
//    if( coordR <= (m_max-m_boxCenter) && coordR >= (m_min-m_boxCenter) )
//    {
//      rval = true;
//    }
    coordR += m_boxCenter;
    if( coordR <= (m_max) && coordR >= (m_min) )
    {
      rval = true;
    }

  }

  return rval;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, Box, std::string const &, ManagedGroup * const )

} /* namespace geosx */
