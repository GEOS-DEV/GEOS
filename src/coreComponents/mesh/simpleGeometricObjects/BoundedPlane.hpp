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
 * @file BoundedPlane.hpp
 */

#ifndef GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_BOUNDEDPLANE_HPP_
#define GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_BOUNDEDPLANE_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

/**
 * @class BoundedPlane
 * @brief Class to represent a geometric box in GEOSX.
 */
class BoundedPlane : public SimpleGeometricObjectBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  BoundedPlane( const string & name,
                Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~BoundedPlane() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "BoundedPlane"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;

  /**
   * @brief Find the bounds of the plane.
   */
  void findRectangleLimits();

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the normal to the plane.
   * @return the normal vector
   */
  R1Tensor & getNormal() {return m_normal;}

  /**
   * @copydoc getNormal()
   */
  R1Tensor const & getNormal() const {return m_normal;}

  /**
   * @brief Get the origin of the plane.
   * @return the origin of the plane
   */
  R1Tensor & getCenter() {return m_origin;}

  /**
   * @copydoc getCenter()
   */
  R1Tensor const & getCenter() const {return m_origin;}

  /**
   * @brief Get one of the tangent vectors defining the orthonormal basis along with the normal.
   * @return the tangent vector
   */
  R1Tensor & getWidthVector() {return m_widthVector;}

  /**
   * @copydoc getWidthVector()
   */
  R1Tensor const & getWidthVector() const {return m_widthVector;}

  /**
   * @brief Get one of the tangent vectors defining the orthonormal basis along with the normal.
   * @return the length vector
   */
  R1Tensor & getLengthVector() {return m_lengthVector;}

  /**
   * @copydoc getLengthVector()
   */
  R1Tensor const & getLengthVector() const {return m_lengthVector;}


protected:

  /**
   * @brief This function provides the capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void postProcessInput() override final;

private:

  /// Origin point (x,y,z) of the plane (basically, any point on the plane)
  R1Tensor m_origin;
  /// Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)
  R1Tensor m_normal;
  /// Length vector in the orthonormal basis along with the normal
  R1Tensor m_lengthVector;
  /// Width vector in the orthonormal basis along with the normal
  R1Tensor m_widthVector;
  /// Dimensions of the bounded plane
  array1d< real64 > m_dimensions;
  /// Length and width of the bounded plane
  array2d< real64 > m_points;
  /// tolerance to determine if a point sits on the plane or not
  real64 m_tolerance;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * originString() { return "origin"; }
    static constexpr char const * normalString() { return "normal"; }
    static constexpr char const * dimensionsString() { return "dimensions"; }
    static constexpr char const * mLengthVectorString() { return "lengthVector"; }
    static constexpr char const * mWidthVectorString() { return "widthVector"; }
    static constexpr char const * toleranceString() { return "tolerance"; }
  };

  /// @endcond

};
} /* namespace geosx */

#endif /* GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_BOUNDEDPLANE_HPP_*/
