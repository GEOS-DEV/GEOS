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
  SimpleGeometricObjectBase( name, parent )
{}

Box::~Box()
{
  // TODO Auto-generated destructor stub
}

void Box::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName(this->CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("A simple box object");

  docNode->AllocateChildNode( viewKeys.xmin.Key(),
                              viewKeys.xmin.Key(),
                              -1,
                              "R1Tensor",
                              "R1Tensor",
                              "Lower corner of box",
                              "Lower corner of box",
                              "-1e99,-1e99,-1e99",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.xmax.Key(),
                              viewKeys.xmax.Key(),
                              -1,
                              "R1Tensor",
                              "R1Tensor",
                              "Upper corner of box",
                              "Upper corner of box",
                              "1e99,1e99,1e99",
                              "",
                              0,
                              1,
                              0 );
}

void Box::ReadXML_PostProcess()
{
  m_min = *(this->getData<R1Tensor>(viewKeys.xmin));
  m_max = *(this->getData<R1Tensor>(viewKeys.xmax));
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
