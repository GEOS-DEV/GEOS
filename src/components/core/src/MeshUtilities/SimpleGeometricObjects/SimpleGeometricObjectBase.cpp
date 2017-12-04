/*
 * SimpleGeometricObjects.cpp
 *
 *  Created on: Dec 4, 2012
 *      Author: settgast1
 */

//#include "codingUtilities/Functions.hpp"
//#include "managers/FunctionManager.h"
#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{


SimpleGeometricObjectBase::SimpleGeometricObjectBase( std::string const & name,
                                                      ManagedGroup * const parent ):
  ManagedGroup( name, parent )
{}


SimpleGeometricObjectBase::~SimpleGeometricObjectBase()
{}


SimpleGeometricObjectBase::CatalogInterface::CatalogType& SimpleGeometricObjectBase::GetCatalog()
{
  static SimpleGeometricObjectBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


//
//// type strings
//const std::string SimpleGeometricObjectBase::BoxStr = "Box";
//const std::string SimpleGeometricObjectBase::CylinderStr = "Cylinder";
//const std::string SimpleGeometricObjectBase::SphereStr = "Sphere";
//const std::string SimpleGeometricObjectBase::EllipsoidStr = "Ellipsoid";
//const std::string SimpleGeometricObjectBase::CylinderBy2EndsStr =
// "CylinderBy2Ends";
//const std::string SimpleGeometricObjectBase::NotStr = "Not";
//const std::string SimpleGeometricObjectBase::IntersectionStr = "Intersection";
//const std::string SimpleGeometricObjectBase::UnionStr = "Union";
//const std::string SimpleGeometricObjectBase::AndStr = "And"; // another way to
// specify intersection
//const std::string SimpleGeometricObjectBase::OrStr = "Or"; // another way to
// specify union
//const std::string SimpleGeometricObjectBase::InvertStr = "Invert"; // another
// way to specify not
//
//const std::string SimpleGeometricObjectBase::TransformStr = "Transform";
//const std::string SimpleGeometricObjectBase::GeometryFunctionStr =
// "GeometryFunction";
//
//
//
//SimpleGeometricObjectBase* SimpleGeometricObjectBase::Allocate( const Types
// type )
//{
//	// need to replace this with a proper factory
//  SimpleGeometricObjectBase* object = NULL;
//  if( type==box )
//    object = new Box;
//  else if (type == cylinder)
//    object = new Cylinder;
//  else if (type == sphere)
//    object = new Sphere;
//  else if (type == cylinderby2ends)
//    object = new CylinderBy2Ends;
//  else if (type == ellipsoid)
//    object = new Ellipsoid;
//  else if (type == intersectionGeometry)
//    object = new IntersectionGeometry;
//  else if (type == unionGeometry)
//    object = new UnionGeometry;
//  else if (type == notGeometry)
//    object = new NotBooleanGeometry;
//  else if (type == transformGeometry)
//    object = new TransformGeometry;
//  else if (type == geometryFunction)
//    object = new GeometryFunction;
//  else
//    throw GPException("SimpleGeometricObjectBase::Allocate: Unrecognized
// type");
//  return object;
//}
//
//
//SimpleGeometricObjectBase* SimpleGeometricObjectBase::Allocate(
// TICPP::HierarchicalDataNode* hdn){
//	std::string geomTypeStr = hdn->Heading();
//	SimpleGeometricObjectBase::Types type = fromString<Types>(geomTypeStr);
//	SimpleGeometricObjectBase* object = SimpleGeometricObjectBase::Allocate(
// type );
//	object->ReadXML(*hdn);
//	return object;
//}
//
//void SimpleGeometricObjectBase::Deallocate( SimpleGeometricObjectBase* object
// )
//{
//  delete object;
//}

}
