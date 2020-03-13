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
 * @file Cylinder.hpp
 *
 */

#ifndef GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
#define GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class Cylinder : public SimpleGeometricObjectBase
{
public:
  Cylinder( const std::string & name,
            Group * const parent );

  virtual ~Cylinder() override;

  static string CatalogName() { return "Cylinder"; }

  bool IsCoordInObject( const R1Tensor & coord ) const override final;


private:

  R1Tensor m_point1;
  R1Tensor m_point2;
  realT m_radius = 0.0;

  struct viewKeyStruct
  {
    static constexpr auto point1String = "point1";
    static constexpr auto point2String = "point2";
    static constexpr auto radiusString = "radius";

  };


};
} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
        */
