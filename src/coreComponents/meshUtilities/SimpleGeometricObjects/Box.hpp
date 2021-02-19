/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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

/**
 * @class Box
 * @brief Class to represent a geometric box in GEOSX.
 */
class Box : public SimpleGeometricObjectBase
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
  Box( const string & name,
       Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~Box() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string catalogName() { return "Box"; }

  ///@}

  bool isCoordInObject( real64 const ( &coord ) [3] ) const override final;

protected:

  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void postProcessInput() override final;

private:

  /// Mininum (x,y,z) coordinates of the box
  R1Tensor m_min;
  /// Maximum (x,y,z) coordinates of the box
  R1Tensor m_max;
  /// Strike angle of the box
  real64 m_strikeAngle=0.0;
  /// Coordinates of the center of the box
  R1Tensor m_boxCenter={0.0, 0.0, 0.0};
  /// Cosine of the strike angle of the box
  real64 m_cosStrike=0.0;
  /// Sine of the strike angle of the box
  real64 m_sinStrike=0.0;

  /// @cond DO_NOT_DOCUMENT

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

  /// @endcond

};
} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_BOX_HPP_
        */
