/*
 * SimpleGeometricObjects.h
 *
 *  Created on: Dec 4, 2012
 *      Author: settgast1
 */

#ifndef SIMPLEGEOMETRICOBJECTS_H_
#define SIMPLEGEOMETRICOBJECTS_H_

#include "../legacy/IO/ticpp/HierarchicalDataNode.h.old"
#include "Common/Common.h"
#include "Utilities/StringUtilities.h"

class Function;

class SimpleGeometricObjectBase
{

public:
  enum Types
  {
    box = 0,
    cylinder = 1,
    sphere = 2,
    cylinderby2ends = 3,
    ellipsoid = 4,
    intersectionGeometry = 5,
    unionGeometry = 6,
    notGeometry = 7,
    transformGeometry = 8,
    geometryFunction = 9,
    numGeometricObjectTypes = 10
  };
  // The default cylinder object use point1 as the center and point2 as a point on one face.  The new type uses the centers of the two faces.

  // The default cylinder object use point1 as the center and point2 as a point on one face.  The new type uses the centers of the two faces.

  // Type strings
  static const std::string BoxStr;
  static const std::string CylinderStr;
  static const std::string SphereStr;
  static const std::string EllipsoidStr;
  static const std::string CylinderBy2EndsStr;
  static const std::string IntersectionStr;
  static const std::string UnionStr;
  static const std::string InvertStr;
  static const std::string AndStr;
  static const std::string OrStr;
  static const std::string NotStr;
  static const std::string TransformStr;
  static const std::string GeometryFunctionStr;

  virtual ~SimpleGeometricObjectBase()
  {
  }

  virtual void ReadXML( TICPP::HierarchicalDataNode& hdn ) = 0;

  virtual bool IsCoordInObject( const R1Tensor& coord ) = 0;

  static SimpleGeometricObjectBase* Allocate( const Types type );

  static SimpleGeometricObjectBase* Allocate( TICPP::HierarchicalDataNode* hdn );

  static void Deallocate( SimpleGeometricObjectBase* object );

  static Types IntToType( const int input )
  {
    Types rval( numGeometricObjectTypes );
    if( input < (int) numGeometricObjectTypes && input >= 0 )
    {
      rval = (Types) input;
    }
    /*
     if( input == box )
     rval = box;
     if( input == cylinder )
     rval = cylinder;
     */

    return rval;
  }
};

// override the fromString template for SimpleGeometricObjectBase::Types
template<>
inline SimpleGeometricObjectBase::Types fromString<SimpleGeometricObjectBase::Types>( std::string theString )
{
  if( theString == SimpleGeometricObjectBase::BoxStr )
  {
    return SimpleGeometricObjectBase::box;
  }
  else if( theString == SimpleGeometricObjectBase::CylinderStr )
  {
    return SimpleGeometricObjectBase::cylinder;
  }
  else if( theString == SimpleGeometricObjectBase::CylinderBy2EndsStr )
  {
    return SimpleGeometricObjectBase::cylinderby2ends;
  }
  else if( theString == SimpleGeometricObjectBase::SphereStr )
  {
    return SimpleGeometricObjectBase::sphere;
  }
  else if( theString == SimpleGeometricObjectBase::EllipsoidStr )
  {
    return SimpleGeometricObjectBase::ellipsoid;
  }
  else if( theString == SimpleGeometricObjectBase::IntersectionStr || theString == SimpleGeometricObjectBase::AndStr )
  {
    return SimpleGeometricObjectBase::intersectionGeometry;
  }
  else if( theString == SimpleGeometricObjectBase::UnionStr || theString == SimpleGeometricObjectBase::OrStr )
  {
    return SimpleGeometricObjectBase::unionGeometry;
  }
  else if( theString == SimpleGeometricObjectBase::NotStr || theString == SimpleGeometricObjectBase::InvertStr )
  {
    return SimpleGeometricObjectBase::notGeometry;
  }
  else if( theString == SimpleGeometricObjectBase::TransformStr )
  {
    return SimpleGeometricObjectBase::transformGeometry;
  }
  else if( strIsInt( theString ) )
  { // backwards compatability
    int ii = fromString<int>( theString );
    return SimpleGeometricObjectBase::IntToType( ii );
  }
  else
  {
    std::cout << SimpleGeometricObjectBase::CylinderBy2EndsStr << std::endl;
    throw GPException( "Error unrecognized geometric object type: " + theString + "." );
  }

  // should never get here
  throw GPException( "Error unrecognized geometric object type: " + theString + "." );
  return SimpleGeometricObjectBase::numGeometricObjectTypes;
}

class Box : public SimpleGeometricObjectBase
{
public:
  void ReadXML( TICPP::HierarchicalDataNode& hdn );

  bool IsCoordInObject( const R1Tensor& coord );
  private:
  R1Tensor m_min;
  R1Tensor m_max;
  realT m_strikeAngle;
  R1Tensor m_boxCenter;
  realT m_cosStrike, m_sinStrike;

};

class Cylinder : public SimpleGeometricObjectBase
{

public:
  void ReadXML( TICPP::HierarchicalDataNode& hdn );

  bool IsCoordInObject( const R1Tensor& coord );
  private:
  R1Tensor m_refPoint;
  R1Tensor m_axis;
  realT m_radius;
  realT m_length;

};

class CylinderBy2Ends : public SimpleGeometricObjectBase
{
public:
  void ReadXML( TICPP::HierarchicalDataNode& hdn );

  bool IsCoordInObject( const R1Tensor& coord );
  private:
  R1Tensor m_refPoint;
  R1Tensor m_axis;
  realT m_radius;
  realT m_length;

};

class Sphere : public SimpleGeometricObjectBase
{

public:
  void ReadXML( TICPP::HierarchicalDataNode& hdn );

  bool IsCoordInObject( const R1Tensor& coord );
  private:
  R1Tensor m_refPoint;
  realT m_radius;
  realT m_radiusSqrd;

};

class Ellipsoid : public SimpleGeometricObjectBase
{

public:
  void ReadXML( TICPP::HierarchicalDataNode& hdn );

  bool IsCoordInObject( const R1Tensor& coord );
  private:
  R1Tensor m_refPoint;
  realT m_rx;
  realT m_ry;
  realT m_rz;

};

class BooleanGeometry : public SimpleGeometricObjectBase
{
public:
  ~BooleanGeometry()
  {
    for( unsigned i = 0 ; i < objectPointers.size() ; ++i )
    {
      delete objectPointers[i];
    }
  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& hdn )
  {
    objectPointers.clear();
    for( TICPP::HierarchicalDataNode* childNode = hdn.Next( true ) ; childNode ;
        childNode = hdn.Next() )
    {
      SimpleGeometricObjectBase* objectPtr = Allocate( childNode );
      objectPointers.push_back( objectPtr );
    }
  }
  virtual bool IsCoordInObject( const R1Tensor& coord )=0;
  protected:
  std::vector<SimpleGeometricObjectBase*> objectPointers;

};

class IntersectionGeometry : public BooleanGeometry
{

public:
  virtual bool IsCoordInObject( const R1Tensor& coord )
  {
    bool rv = objectPointers.size() > 0; // returns false if empty - prob should just throw an error

    for( unsigned i = 0 ; i < objectPointers.size() ; ++i )
    {
      rv = objectPointers[i]->IsCoordInObject( coord );
      if( !rv )
        break;
    }
    return rv;
  }
private:

};

class UnionGeometry : public BooleanGeometry
{

public:
  virtual bool IsCoordInObject( const R1Tensor& coord )
  {
    bool rv = false; // returns false if empty

    for( unsigned i = 0 ; i < objectPointers.size() ; ++i )
    {
      rv = objectPointers[i]->IsCoordInObject( coord );
      if( rv )
        break;
    }
    return rv;
  }
private:

};

class NotBooleanGeometry : public BooleanGeometry
{

public:
  void ReadXML( TICPP::HierarchicalDataNode& hdn )
  {
    BooleanGeometry::ReadXML( hdn );
    if( objectPointers.size() != 1 )
    {
      throw GPException( "Error: NotObject must include one and only one geometric object." );
    }
  }

  bool IsCoordInObject( const R1Tensor& coord )
  {
    bool rv = objectPointers[0]->IsCoordInObject( coord );
    return !rv;
  }
private:

};

//////////////

// Transformations

// transforms the geometric object(s) from x to x' = Ux+r

class TransformGeometry : public UnionGeometry
{

public:
  virtual void ReadXML( TICPP::HierarchicalDataNode& hdn );
  virtual bool IsCoordInObject( const R1Tensor& coord );
  protected:
  R2Tensor m_U;
  R2Tensor m_Uinv;
  R1Tensor m_R;
};

//////////////

// Geometry function

// returns true (inside) if function is positive, false (outside) if negative

class GeometryFunction : public SimpleGeometricObjectBase
{

public:
  virtual void ReadXML( TICPP::HierarchicalDataNode& hdn );
  virtual bool IsCoordInObject( const R1Tensor& coord );
  protected:

  Function* m_function;
  std::string m_functionName;
};

#endif /* SIMPLEGEOMETRICOBJECTS_H_ */
