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
 * @file BoundedThickPlane.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOUNDEDTHICKPLANE_HPP_
#define SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOUNDEDTHICKPLANE_HPP_

#include "ThickPlane.hpp"

namespace geosx
{

class BoundedThickPlane : public ThickPlane
{
public:
  BoundedThickPlane( const std::string& name,
              Group * const parent );

  virtual ~BoundedThickPlane() override;

  static string CatalogName() { return "BoundedThickPlane"; }

  bool IsCoordInObject( const R1Tensor& coord ) const override final;

protected:
  virtual void PostProcessInput() override final;

private:
  R1Tensor            m_lengthVector;
  R1Tensor            m_widthVector;
  array1d < real64 >  m_dimensions;

  struct viewKeyStruct
  {
    static constexpr auto dimensionsString    = "dimensions";
    static constexpr auto mLengthVectorString = "lengthVector";
    static constexpr auto mWidthVectorString  = "widthVector";
  };


};
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOUNDEDTHICKPLANE_HPP_
        */

