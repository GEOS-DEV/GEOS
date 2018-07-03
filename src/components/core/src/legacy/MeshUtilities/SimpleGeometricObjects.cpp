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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * SimpleGeometricObjects.cpp
 *
 *  Created on: Dec 4, 2012
 *      Author: settgast1
 */

#include "SimpleGeometricObjects.h"

#include "../../codingUtilities/Functions.hpp"
#include "ObjectManagers/FunctionManager.h"


// type strings
const std::string SimpleGeometricObjectBase::BoxStr = "Box";
const std::string SimpleGeometricObjectBase::CylinderStr = "Cylinder";
const std::string SimpleGeometricObjectBase::SphereStr = "Sphere";
const std::string SimpleGeometricObjectBase::EllipsoidStr = "Ellipsoid";
const std::string SimpleGeometricObjectBase::CylinderBy2EndsStr = "CylinderBy2Ends";
const std::string SimpleGeometricObjectBase::NotStr = "Not";
const std::string SimpleGeometricObjectBase::IntersectionStr = "Intersection";
const std::string SimpleGeometricObjectBase::UnionStr = "Union";
const std::string SimpleGeometricObjectBase::AndStr = "And"; // another way to
                                                             // specify
                                                             // intersection
const std::string SimpleGeometricObjectBase::OrStr = "Or"; // another way to
                                                           // specify union
const std::string SimpleGeometricObjectBase::InvertStr = "Invert"; // another
                                                                   // way to
                                                                   // specify
                                                                   // not

const std::string SimpleGeometricObjectBase::TransformStr = "Transform";
const std::string SimpleGeometricObjectBase::GeometryFunctionStr = "GeometryFunction";



SimpleGeometricObjectBase* SimpleGeometricObjectBase::Allocate( const Types type )
{
  // need to replace this with a proper factory
  SimpleGeometricObjectBase* object = NULL;
  if( type==box )
    object = new Box;
  else if (type == cylinder)
    object = new Cylinder;
  else if (type == sphere)
    object = new Sphere;
  else if (type == cylinderby2ends)
    object = new CylinderBy2Ends;
  else if (type == ellipsoid)
    object = new Ellipsoid;
  else if (type == intersectionGeometry)
    object = new IntersectionGeometry;
  else if (type == unionGeometry)
    object = new UnionGeometry;
  else if (type == notGeometry)
    object = new NotBooleanGeometry;
  else if (type == transformGeometry)
    object = new TransformGeometry;
  else if (type == geometryFunction)
    object = new GeometryFunction;
  else
    throw GPException("SimpleGeometricObjectBase::Allocate: Unrecognized type");
  return object;
}


SimpleGeometricObjectBase* SimpleGeometricObjectBase::Allocate( TICPP::HierarchicalDataNode* hdn){
  std::string geomTypeStr = hdn->Heading();
  SimpleGeometricObjectBase::Types type = fromString<Types>(geomTypeStr);
  SimpleGeometricObjectBase* object = SimpleGeometricObjectBase::Allocate( type );
  object->ReadXML(*hdn);
  return object;
}

void SimpleGeometricObjectBase::Deallocate( SimpleGeometricObjectBase* object )
{
  delete object;
}



void Box::ReadXML( TICPP::HierarchicalDataNode& hdn )
{
  m_min = hdn.GetAttributeTensor("xmin");
  m_max = hdn.GetAttributeTensor("xmax");
  m_strikeAngle = hdn.GetAttributeOrDefault<realT>("strikeAngle", -90.0); // from
                                                                          // North
  m_strikeAngle += 90; // Counterclockwise from x-axis
  if (std::fabs(m_strikeAngle) > 1e-20)
  {
    if ((m_max[0]-m_min[0]) < (m_max[1]-m_min[1]))
      throw GPException(
              "Error: When a strike angle is specified, the box is supposed to represent a plane normal to the y direction. This box seems to be too thick.");

    m_cosStrike = std::cos(m_strikeAngle / 180 *3.1415926535);
    m_sinStrike = std::sin(m_strikeAngle / 180 *3.1415926535);
    m_boxCenter = m_min;
    m_boxCenter += m_max;
    m_boxCenter *= 0.5;
  }

}

bool Box::IsCoordInObject( const R1Tensor& coord )
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
    if( coordR <= (m_max-m_boxCenter) && coordR >= (m_min-m_boxCenter) )
    {
      rval = true;
    }
  }

  return rval;
}



void Cylinder::ReadXML( TICPP::HierarchicalDataNode& hdn )
{


  m_refPoint = hdn.GetAttributeTensor("point1");
  m_axis = hdn.GetAttributeTensor("point2");
  m_axis -= m_refPoint;
  m_length = m_axis.Normalize();
  m_radius = hdn.GetAttributeValue<realT>("radius");
}

bool Cylinder::IsCoordInObject( const R1Tensor& coord )
{
  bool rval = false;

  R1Tensor coord_minus_ref;
  coord_minus_ref = coord;
  coord_minus_ref -= m_refPoint;

  R1Tensor projection;
  projection  = m_axis;
  projection *= Dot( m_axis, coord_minus_ref );  //Fu:  Used to be *= -Dot...  I
                                                 // didn't understand and I
                                                 // thought it's wrong.

  R1Tensor distance;
  distance  = coord_minus_ref;
  distance -= projection;


  if( distance.L2_Norm()<m_radius && projection.L2_Norm()<m_length )
  {
    rval = true;
  }

  return rval;
}

void CylinderBy2Ends::ReadXML( TICPP::HierarchicalDataNode& hdn )
{
  m_refPoint = hdn.GetAttributeTensor("point1");
  m_axis = hdn.GetAttributeTensor("point2");
  m_axis -= m_refPoint;
  m_length = m_axis.Normalize();
  m_radius = hdn.GetAttributeValue<realT>("radius");
}

bool CylinderBy2Ends::IsCoordInObject( const R1Tensor& coord )
{
  bool rval = false;
  realT projectionLength;

  R1Tensor coord_minus_ref;
  coord_minus_ref = coord;
  coord_minus_ref -= m_refPoint;

  R1Tensor projection;
  projection  = m_axis;
  projectionLength = Dot(m_axis, coord_minus_ref);
  projection *= projectionLength;

  R1Tensor distance;
  distance  = coord_minus_ref;
  distance -= projection;



  if( distance.L2_Norm()<m_radius && projectionLength >= 0.0 && projectionLength < m_length )
  {
    rval = true;
  }

  return rval;
}


void Sphere::ReadXML( TICPP::HierarchicalDataNode& hdn )
{

  m_refPoint = hdn.GetAttributeTensor("center");
  m_radius = hdn.GetAttributeValue<realT>("radius");
  m_radiusSqrd = m_radius*m_radius;
}

bool Sphere::IsCoordInObject( const R1Tensor& coord )
{
  bool rval = false;

  R1Tensor coord_minus_ref;
  coord_minus_ref = coord;
  coord_minus_ref -= m_refPoint;

  realT distanceSqrd = Dot( coord_minus_ref, coord_minus_ref );

  if( distanceSqrd<m_radiusSqrd )
  {
    rval = true;
  }

  return rval;
}


void Ellipsoid::ReadXML( TICPP::HierarchicalDataNode& hdn )
{

  m_refPoint = hdn.GetAttributeTensor("center");
  m_rx = hdn.GetAttributeValue<realT>("rx");
  m_ry = hdn.GetAttributeValue<realT>("ry");
  m_rz = hdn.GetAttributeValue<realT>("rz");
}

bool Ellipsoid::IsCoordInObject( const R1Tensor& coord )
{
  bool rval = false;

  R1Tensor coord_minus_ref;
  coord_minus_ref = coord;
  coord_minus_ref -= m_refPoint;

  coord_minus_ref[0] /= m_rx;
  coord_minus_ref[1] /= m_ry;
  coord_minus_ref[2] /= m_rz;

  realT distanceSqrd = Dot( coord_minus_ref, coord_minus_ref );

  if( distanceSqrd<1 )
  {
    rval = true;
  }

  return rval;
}


/////////////////////////

// Transformations

void TransformGeometry::ReadXML( TICPP::HierarchicalDataNode& hdn )
{
  UnionGeometry::ReadXML(hdn);
  m_U = hdn.GetAttributeR2Tensor("U");
  m_R = hdn.GetAttributeTensor("r");
  m_Uinv = m_U.Inverse();
}

bool TransformGeometry::IsCoordInObject( const R1Tensor& coord )
{

  R1Tensor transformedCoord;
  transformedCoord.AijBi( m_Uinv, coord - m_R);
  bool rval = UnionGeometry::IsCoordInObject(transformedCoord);

  return rval;
}


/////////////////////////

// Geometry function
/* <GeometryFunction function="myLevelSet"/>
 *
 */
void GeometryFunction::ReadXML( TICPP::HierarchicalDataNode& hdn )
{
  m_functionName      = hdn.GetAttributeString("function");

  // get function
  if(!(m_functionName.empty()))
  {
    FunctionManager& fm = FunctionManager::Instance();
    m_function = &(fm.GetFunction(m_functionName));
  }

}

bool GeometryFunction::IsCoordInObject( const R1Tensor& coord )
{
  realT value =  (*m_function)(coord.Data()[0]);
  return value > 0;
}
