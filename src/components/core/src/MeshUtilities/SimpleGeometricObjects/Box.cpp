/*
 * Box.cpp
 *
 *  Created on: Aug 4, 2017
 *      Author: settgast
 */

#include "Box.hpp"

namespace geosx
{

Box::Box():
  SimpleGeometricObjectBase()
{
  // TODO Auto-generated constructor stub

}

Box::~Box()
{
  // TODO Auto-generated destructor stub
}

void Box::ReadXML( xmlWrapper::xmlNode const & xmlNode )
{
//  pugi::xml_attribute xmlatt =
// targetNode.attribute(subDocNode.getStringKey().c_str());
//  m_degree = targetNode.attribute("degree").as_int(1);
//  as_type( xmlVal, xmlatt.value(), defVal );

  m_min = xmlWrapper::as_type( xmlNode, "xMin", {-1e99,-1e99,-1e99});
  m_max = xmlWrapper::as_type( xmlNode, "xMax", { 1e99, 1e99, 1e99});
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

}

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

REGISTER_CATALOG_ENTRY0( SimpleGeometricObjectBase, Box )

} /* namespace geosx */
