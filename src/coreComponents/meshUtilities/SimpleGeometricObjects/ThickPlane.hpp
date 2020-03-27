/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ThickPlane.hpp
 */

#ifndef GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_THICKPLANE_HPP_
#define GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_THICKPLANE_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class ThickPlane : public SimpleGeometricObjectBase
{
public:
  ThickPlane( const std::string & name,
              Group * const parent );

  virtual ~ThickPlane() override;

  static string CatalogName() { return "ThickPlane"; }

  bool IsCoordInObject( const R1Tensor & coord ) const override final;

  /*
   * Accessors
   */
  // normal vector
  R1Tensor & getNormal() {return m_normal;}

  R1Tensor const & getNormal() const {return m_normal;}

  // origin of the plane
  R1Tensor & getCenter() {return m_origin;}

  R1Tensor const & getCenter() const {return m_origin;}

protected:
  virtual void PostProcessInput() override final;

private:

  R1Tensor m_origin;
  R1Tensor m_normal;
  real64 m_thickness;

  struct viewKeyStruct
  {
    static constexpr auto originString = "origin";
    static constexpr auto normalString = "normal";
    static constexpr auto thicknessString = "thickness";
  };

};
} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_THICKPLANE_HPP_
        */
