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

/**
 * @file InternalWellboreGenerator.hpp
 */

#ifndef GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP
#define GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP

#include "common/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "InternalMeshGenerator.hpp"

namespace geosx
{

/**
 * @class InternalWellboreGenerator
 * @brief The InternalWellboreGenerator class is a class generating internal wellbore mesh.
 */
class InternalWellboreGenerator : public InternalMeshGenerator
{
public:

  /**
   * @brief Main constructor for InternalWellboreGenerator.
   * @param[in] name of the InternalWellboreGenerator
   * @param[in] parent point to the parent Group of the InternalWellboreGenerator
   */
  InternalWellboreGenerator( const string & name, Group * const parent );

  ~InternalWellboreGenerator() override = default;

  /**
   * @brief Return the name of the InternalWellboreGenerator in object Catalog.
   * @return string that contains the key name to InternalWellboreGenerator in the Catalog
   */
  static string catalogName() { return "InternalWellbore"; }

  void generateMesh( DomainPartition & domain ) override final;

protected:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * radiusString() { return "radius"; }
    constexpr static char const * thetaString() { return "theta"; }
    constexpr static char const * rOutString() { return "rOut"; }
    constexpr static char const * rElemsString() { return "nr"; }
    constexpr static char const * tElemsString() { return "nt"; }
    constexpr static char const * rBiasString() { return "rBias"; }
    constexpr static char const * trajectoryString() { return "trajectory"; }
  };
  /// @endcond

  void postProcessInput() override final;

private:

  /// Wellbore radius
  real64 m_radius;

  /// Farfield distance
  real64 m_rOut;

  /// Maximum tangent angle in degrees
  real64 m_theta;

  /// Number of elements in radial direction
  real64 m_rElems;

  /// Number of elements in tangent direction
  real64 m_tElems;

  /// Bias in the radial direction
  real64 m_rBias;

  /// Trajectory defined by coordinates of the wellbore centers
  array2d< real64 > m_trajectory;

};

} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP */
