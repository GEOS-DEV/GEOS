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
 * @file BoundedPlane.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOUNDEDPLANE_HPP_
#define SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOUNDEDPLANE_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class BoundedPlane : public SimpleGeometricObjectBase
{
public:
  BoundedPlane( const std::string & name,
                Group * const parent );

  virtual ~BoundedPlane() override;

  static string CatalogName() { return "BoundedPlane"; }

  bool IsCoordInObject( const R1Tensor & coord ) const override final;

  void findRectangleLimits();

  /*
   * Accessors
   */
  // normal vector
  R1Tensor & getNormal() {return m_normal;}

  R1Tensor const & getNormal() const {return m_normal;}

  // origin of the plane
  R1Tensor & getCenter() {return m_origin;}

  R1Tensor const & getCenter() const {return m_origin;}

  // width vector
  R1Tensor & getWidthVector() {return m_widthVector;}

  R1Tensor const & getWidthVector() const {return m_widthVector;}

  // length vector
  R1Tensor & getLengthVector() {return m_lengthVector;}

  R1Tensor const & getLengthVector() const {return m_lengthVector;}


protected:
  virtual void PostProcessInput() override final;

private:
  R1Tensor m_origin;
  R1Tensor m_normal;
  R1Tensor m_lengthVector;
  R1Tensor m_widthVector;
  array1d< real64 >   m_dimensions;
  array1d< R1Tensor > m_points;

  struct viewKeyStruct
  {
    static constexpr auto originString = "origin";
    static constexpr auto normalString = "normal";
    static constexpr auto dimensionsString    = "dimensions";
    static constexpr auto mLengthVectorString = "lengthVector";
    static constexpr auto mWidthVectorString  = "widthVector";
  };


};
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOUNDEDPLANE_HPP_
        */
