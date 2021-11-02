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
 * @file Cylinder.hpp
 *
 */

#ifndef GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
#define GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

/**
 * @class Cylinder
 * @brief Class to represent a geometric cylinder in GEOSX.
 */
class Cylinder : public SimpleGeometricObjectBase
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
  Cylinder( const string & name,
            Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~Cylinder() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "Cylinder"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;


private:

  /// Center point of one (upper or lower) face of the cylinder
  R1Tensor m_point1;
  /// Center point of the other face of the cylinder
  R1Tensor m_point2;
  /// Radius of the cylinder
  real64 m_radius = 0.0;

  real64 m_innerRadius = 0.0;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr char const * point1String() { return "point1"; }
    static constexpr char const * point2String() { return "point2"; }
    static constexpr char const * radiusString() { return "radius"; }
    static constexpr char const * innerRadiusString() { return "innerRadius"; }
  };

  /// @endcond

};
} /* namespace geosx */

#endif /* GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
        */
