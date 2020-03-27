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

/*
 * @file Box.hpp
 *
 */

#ifndef GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_
#define GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

class Box : public SimpleGeometricObjectBase
{
public:
  Box( const std::string & name,
       Group * const parent );

  virtual ~Box() override;

  static string CatalogName() { return "Box"; }

  bool IsCoordInObject( const R1Tensor & coord ) const override final;

protected:
  virtual void PostProcessInput() override final;

private:

  R1Tensor m_min;
  R1Tensor m_max;
  realT m_strikeAngle=0.0;
  R1Tensor m_boxCenter={0.0, 0.0, 0.0};
  realT m_cosStrike=0.0;
  real64 m_sinStrike=0.0;

  struct viewKeyStruct
  {
    static constexpr auto xMinString = "xMin";
    static constexpr auto xMaxString = "xMax";
    static constexpr auto strikeAngleString = "strike";
    static constexpr auto boxCenterString = "center";
    static constexpr auto cosStrikeString = "cosStrike";
    static constexpr auto sinStrikeString = "sinStrike";

    dataRepository::ViewKey xmin = { "xMin" };
    dataRepository::ViewKey xmax = { "xMax" };
  } viewKeys;


};
} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_
        */
