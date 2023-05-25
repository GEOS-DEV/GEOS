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
 * @file Circle.hpp
 */

#ifndef GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_CIRCLE_HPP_
#define GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_CIRCLE_HPP_

#include "SimpleGeometricObjectBase.hpp"
#include "BoundedPlanarObject.hpp"

namespace geos
{

/**
 * @class Circle
 * @brief Class to represent a geometric disc in GEOSX.
 */
class Circle : public BoundedPlanarObject
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
  Circle( const string & name,
          Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~Circle() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "Circle"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the center of the circle.
   * @return the center of the circle
   */
  virtual R1Tensor & getCenter() override final {return m_center;}

  /**
   * @copydoc getCenter()
   */
  virtual R1Tensor const & getCenter() const override final {return m_center;}

protected:

  /**
   * @brief This function provides the capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void postProcessInput() override final;

private:

  /// center of the circle in (x,y,z) coordinates
  R1Tensor m_center;
  /// Dimensions of the circle's radius
  real64 m_radius;
  /// tolerance to determine if a point sits on the circle or not
  real64 m_tolerance;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * centerString() { return "center"; }
    static constexpr char const * radiusString() { return "radius"; }
    static constexpr char const * toleranceString() { return "tolerance"; }
  };

  /// @endcond

};
} /* namespace geosx */

#endif /* GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_CIRCLE_HPP_*/
