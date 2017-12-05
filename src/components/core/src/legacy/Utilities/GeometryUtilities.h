/*
 * GeometryUtilities.h
 *
 *  Created on: Nov 19, 2011
 *      Author: settgast1
 */

#ifndef GEOMETRYUTILITIES_H_
#define GEOMETRYUTILITIES_H_

//#include "legacy/Common/Common.h"
#include "Utilities.h"

namespace GeometryUtilities
{

realT Area_2DPolygon(const array<R1TensorT<2> >& pointsLocalProjection);

realT CentroidAndVolume_3DTetrahedron(const R1Tensor& x0, const R1Tensor& x1, const R1Tensor& x2,
                                      const R1Tensor& x3, R1Tensor& x);

realT Centroid_2DPolygon(const array<R1TensorT<2> >& points, R1TensorT<2>& center);

realT Centroid_3DPolygon(const lArray1d& pointsEntries,
                         const array<R1Tensor>& points,
                         R1Tensor& center,
                         R1Tensor& normal);

realT Centroid_3DPolygon(const lArray1d& pointsEntries,
                         const array<R1Tensor>& pointsCommonPlanes,
                         R1Tensor& center);

realT Centroid_3DPolygon(const lArray1d& pointsEntries,
                         const array<R1Tensor>& pointReferences,
                         const array<R1Tensor>& pointDisplacements,
                         R1Tensor& center,
                         R1Tensor& normal);

realT Centroid_3DPolygon(const lArray1d& pointsEntries,
                         const array<R1Tensor>& pointReferences,
                         const array<R1Tensor>& pointDisplacements,
                         R1Tensor& center);

void FaceNormal3DPolygonReferenceConfig(const lArray1d& pointsEntries,
                                        const array<R1Tensor>& pointReferences,
                                        R1Tensor& normal );

void FindProjectionInParentSpace(const R1Tensor& xPoint, const R1Tensor& normal,
                                 const array<R1Tensor>& xNodes, R1Tensor& soln, realT N[4],
                                 const bool initialGuess = false,
                                 const realT tol = 1e-6);

realT ProjectPointToPlane( const R1Tensor& point,
                           const R1Tensor& xplane,
                           const R1Tensor& nplane,
                           R1Tensor& pointproj);

realT OrthogonalVectorComponent( const R1Tensor& vtarget,
                                 const R1Tensor& vorthogonal,
                                 R1Tensor& vresidual);

void ProjectPointToPlaneAlongUnitVector(const R1Tensor& point,
                                        const R1Tensor& xplane,
                                        const R1Tensor& nplane,
                                        const R1Tensor& vproj,
                                        R1Tensor& pointproj);

void ProjectPointToLineSegment( const R1Tensor& x0,
                                const R1Tensor& x1,
                                const R1Tensor& point,
                                realT& ndist,
                                realT& udist,
                                realT& segmentLength,
                                R1Tensor& pointproj);

localIndex LineIntersection( const R1TensorT<2>& uv1,
                             const R1TensorT<2>& uv1p,
                             const R1TensorT<2>& uv2,
                             const R1TensorT<2>& uv2p,
                             R1TensorT<2>& intersection,
                             R1TensorT<2>& intersection2,
                             const realT small = 0);

bool PointInTetrahedron(const R1Tensor& x,
                        const R1Tensor& x0,
                        const R1Tensor& x1,
                        const R1Tensor& x2,
                        const R1Tensor& x3);

bool PointInFace(const array<R1TensorT<2> >& uv,
                 const R1TensorT<2>& uv0,
                 const realT tol = 1e-10);

bool PointInPolyhedron(const R1Tensor& point,
                       const localIndex elementIndex,
                       const lArray2d& toFaces,
                       const lArray2d& toNodes,
                       const array<lArray1d>& faceToNodes,
                       const array<R1Tensor>& referencePosition,
                       const array<R1Tensor>& displacement,
                       const bool planarFaces = true);

bool PointInTriangle(
  const R1TensorT<2>& uva, const R1TensorT<2>& uvb,
  const R1TensorT<2>& uvc, const R1TensorT<2>& uv0,
  const realT tol = 1e-10);

array<R1TensorT<2> >::size_type CoincidentIntersection_2DPolygons(
  const array<R1TensorT<2> >& a1,
  const array<R1TensorT<2> >& a2,
  array<R1TensorT<2> >& aa
  );

void Intersections_2DPolygons(const array<R1TensorT<2> >& a1,
                              const array<R1TensorT<2> >& a2,
                              array<R1TensorT<2> >& intersections);

void OrderPoints_2DPolygon(array<R1TensorT<2> >& polygonPoints, const realT tol );

template< int T_size >
void RemoveCoincidentOrderedPoints(array<R1TensorT<T_size> >& polygonPoints, const realT tol = 1e-6);

realT Intersection_2DPolygons(const array<R1TensorT<2> >& a1,
                              const array<R1TensorT<2> >& a2);

realT Intersection_2DPolygons(const array<R1TensorT<2> >& a1,
                              const array<R1TensorT<2> >& a2,
                              const realT& positionTolerance,
                              array<R1TensorT<2> >& polygonPoints);

bool VectorInPlane(const R1Tensor& normal,
                   const R1Tensor& point,
                   const realT penetration,
                   R1Tensor& e1);

bool VectorsInPlane(const R1Tensor& normal,
                    R1Tensor& e1,
                    R1Tensor& e2);

void CartesianPointsProjectedToPlanarPoints(const lArray1d& pointsEntries,
                                            const array<R1Tensor>& pointReferences,
                                            const array<R1Tensor>& pointDisplacements,
                                            const R1Tensor& e1,
                                            const R1Tensor& e2,
                                            array<R1TensorT<2> >& localPoints,
                                            R1TensorT<2>& min,
                                            R1TensorT<2>& max);

void CartesianPointsProjectedToPlanarPoints(  const array<R1Tensor>& points,
                                              const R1Tensor& e1,
                                              const R1Tensor& e2,
                                              array<R1TensorT<2> >& localPoints,
                                              R1TensorT<2>& min,
                                              R1TensorT<2>& max);

inline void PlanarPointProjectedToCartesianCoordinates(const R1TensorT<2>& localPoint,
                                                       const R1Tensor& e1,
                                                       const R1Tensor& e2,
                                                       R1Tensor& point)
{
  point = e1;
  point *= localPoint(0);
  R1Tensor tmp(e2);
  tmp *= localPoint(1);
  point += tmp;
}

realT DistanceToPlane(const R1Tensor& normal,
                      const R1Tensor& pointOnPlane,
                      const R1Tensor& pointToQuery);

void ClipByPlane_3DPolygon(const R1Tensor& normal,
                           const R1Tensor& pointOnPlane,
                           const array<R1Tensor>& points,
                           const R1Tensor& e1,
                           const R1Tensor& e2,
                           const array<real64>& penetration,
                           array<R1TensorT<2> >& inContact,
                           array<R1TensorT<2> >& outContact,
                           array<real64>& inDistance,
                           array<real64>& outDistance,
                           const bool reverse = false);

void ClipByPlane_3DPolygon(const R1Tensor& normal,
                           const R1Tensor& pointOnPlane,
                           const array<R1Tensor>& points,
                           const R1Tensor& e1,
                           const R1Tensor& e2,
                           const array<real64>& penetration,
                           array<array<R1TensorT<2> > >& inContact,
                           array<array<R1TensorT<2> > >& outContact,
                           const bool reverse = false);

void ClipByPlane_3DPolygon(const R1Tensor& normal,
                           const R1Tensor& pointOnPlane,
                           const array<R1Tensor>& points,
                           array<array<R1TensorT<2> > >& inContact,
                           array<array<R1TensorT<2> > >& outContact,
                           const bool reverse = false);

realT Overlap2DPolygons( const array<R1TensorT<2> >& a1,
                         const array<R1TensorT<2> >& a2,
                         const realT& positionTolerance,
                         array<R1TensorT<2> >& polygonPoints);

void Intersections2DPolygon( const array<R1TensorT<2> >& a1,
                             const array<R1TensorT<2> >& a2,
                             array<R1TensorT<2> >& intersections);

localIndex LineIntersection( const R1TensorT<2>& uv1,
                             const R1TensorT<2>& uv1p,
                             const R1TensorT<2>& uv2,
                             const R1TensorT<2>& uv2p,
                             R1TensorT<2>& intersection,
                             const realT small);


//****** GEOLOGY OPERATIONS **************************************************
void Strike( const R1Tensor& nt,
             const R1Tensor& upt,
             R1Tensor& strikeVec);

void Strike( const R1Tensor& east,
             const R1Tensor& north,
             const realT strikeAngle,
             R1Tensor& strike);

realT StrikeAngle( const R1Tensor& nt,
                   const R1Tensor& upt,
                   const R1Tensor& strikeVec);

///Get the up-dip vector for the given up and normal vectors
/**
   @param[in] nt Unit normal to the fault surface (not checked to be unit)
   @param[in] upt Unit vector in the direction of vertical
 */
inline realT Dip( const R1Tensor& nt, const R1Tensor& upt) { return acos(Dot(nt, upt));}

void Dip( const R1Tensor& up,
          const R1Tensor& strike,
          const realT dipAngle,
          R1Tensor& dip);

///Get the transform matrix from plane frame to space frame
/**
   @param[in] normalToPlane Vector normal to the fault surface
   @param[in] up Vector in the direction of vertical
   @param[out] T Transform tensor
 */
void TransformPlaneFrame(const R1Tensor& normalToPlane,
                         const R1Tensor& up,
                         R2Tensor& T);

///Get the transform matrix from strike frame to space frame
/**
   @param[in] normalToPlane Vector normal to the fault surface
   @param[in] up Vector in the direction of vertical
   @param[out] T Transform tensor
 */
void TransformStrikeFrame(const R1Tensor& normalToPlane,
                          const R1Tensor& up,
                          R2Tensor& T);

///Get the transform matrix from strike frame to space frame
/**
   @param[in] normalToPlane Vector normal to the fault surface
   @param[in] strike Vector in the direction of strike
   @param[out] T Transform tensor
 */
void TransformStrikeFrame2(const R1Tensor& strike,
                           const R1Tensor& up,
                           R2Tensor& T);

///Get the transform matrix from WPP frame to space frame
/**
   @param[in] north Vector in the direction of north
   @param[in] up Vector in the direction of vertical
   @param[out] T Transform tensor
 */
void TransformWPPFrame(const R1Tensor& north,
                       const R1Tensor& up,
                       R2Tensor& T);

///Transform the R2 tensor from a frame to another
/**
   @param[in] transform Transform matrix to prime frame
   @param[in] tensor Current matrix
   @param[out] tensorPrime Matrix transformed to the prime frame
 */
void TransformTensorFrame(const R2Tensor& transform,
                          const R2Tensor& tensor,
                          R2Tensor& tensorPrime);

///Transform the R2 tensor from a frame to another
/**
   @param[in] transform Transform matrix FROM prime frame
   @param[in] tensor Current matrix
   @param[out] tensorPrime Matrix transformed to the prime frame
 */
void TransformTensorFrameTranspose(const R2Tensor& transform,
                                   const R2Tensor& tensor,
                                   R2Tensor& tensorPrime);

void StrikeDipRotationMatrix(const realT strike,
                             const realT dip,
                             bool invert,
                             R2Tensor& rotationMatrix);

template<typename T>
void MapFromRegion(const T& p0, const T& p1, const realT fct0,
                   const realT fct1, T& p,
                   const bool tangential = false)
{
  //let's inform this weighted by their current contributions to common planes
  const realT coef0 = (isZero(fct0) && isZero(fct1)) ? 0.5 :
                      (tangential ? (-fct0 / (fct0 + fct1)) : (fct0 / (fct0 + fct1)));

  const realT coef1 = 1 - coef0;

  p = p0;
  p *= coef0;

  T tmp = p1;
  tmp *= coef1;
  p += tmp;

}

template<typename T>
void MapToRegion(const realT fctNormal, const realT fct0,
                 const realT fct1, const T& p, T& p0, T& p1,
                 const bool tangential = false)
{
  const realT alpha0 = tangential ? (-fctNormal / fct0) : (fctNormal / fct0);
  const realT alpha1 = fctNormal / fct1;

  T t;

  t = p;
  t *= alpha0;
  p0 += t;

  t = p;
  t *= alpha1;
  p1 += t;
}


/**
 * @author Scott Johnson
 * @brief Structure to hold local point data for sorting
 */
class PointStruct
{
public:
  /**
   * @brief Constructor
   * @author Scott Johnson
   * @param[in] point Point
   */
  PointStruct(const R1TensorT<2>& point): m_angle(0), m_value(0)
  {
    this->m_point = point;
  }

  /**
   * @brief Comparison operator
   * @author Scott Johnson
   * @param a Object 1
   * @param b Object 2
   * @return Comparison flag
   */
  static bool Compare( const PointStruct& a, const PointStruct& b )
  {
    return isEqual(a.m_angle, b.m_angle) ? a.m_value > b.m_value : a.m_angle > b.m_angle;
  }

  /**
   * @brief Lowest operator
   * @author Scott Johnson
   * @param a Object 1
   * @param b Object 2
   * @return Comparison flag
   */
  static bool Lowest( const PointStruct& a, const PointStruct& b )
  {
    return isEqual(a.y(), b.y()) ? a.x() < b.x() : a.y() < b.y();
  }

  static int CounterClockwise(const PointStruct& a, const PointStruct& b, const PointStruct& c)
  {
    return (b.x() - a.x())*(c.y() - a.y()) - (b.y() - a.y())*(c.x() - a.x());
  }

  inline static realT Cross(const PointStruct& a, const PointStruct& b)
  {
    return (a.x()*b.y()) - (a.y()*b.x());
  }

  inline static realT Cross(const R1TensorT<2> a, const R1TensorT<2> b)
  {
    return (a.Data()[0]*b.Data()[1]) - (a.Data()[1]*b.Data()[0]);
  }

  inline realT x() const
  {
    return m_point.Data()[0];
  }

  inline realT y() const
  {
    return m_point.Data()[1];
  }

  ///Point
  R1TensorT<2> m_point;

  ///Dot product value
  realT m_angle;

  ///Total value
  realT m_value;
};

} /* namespace GeometryUtilities */
#endif /* GEOMETRYUTILITIES_H_ */
