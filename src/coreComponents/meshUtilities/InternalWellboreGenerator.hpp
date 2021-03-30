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


  virtual void reduceNumNodesForPeriodicBoundary( integer ( &numNodes )[3] ) override final;

  virtual void setNodeGlobalIndicesOnPeriodicBoundary( int (& index)[3],
                                                       real64 ( &minExtent )[3],
                                                       real64 ( &maxExtent )[3],
                                                       arraySlice1d< real64 const > const & X,
                                                       real64 const tol ) override final;

  virtual void setConnectivityForPeriodicBoundaries( integer const i,
                                                     integer const j,
                                                     integer const k,
                                                     integer const iBlock,
                                                     integer const jBlock,
                                                     integer const kBlock,
                                                     int (& globalIJK)[3],
                                                     int const (&numElemsInDirForBlock)[3],
                                                     integer const (&numNodesInDir)[3],
                                                     int const (&firstElemIndexInPartition)[3],
                                                     localIndex ( &nodeOfBox )[8] ) override final;

  virtual void coordinateTransformation( NodeManager & nodeManager ) override final;

  virtual inline bool isCartesian() const override final
  {
    return false;
  }


protected:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct : public InternalMeshGenerator::viewKeyStruct
  {
    constexpr static char const * radiusString() { return "radius"; }
    constexpr static char const * thetaString() { return "theta"; }
    constexpr static char const * rOutString() { return "rOut"; }
    constexpr static char const * rElemsString() { return "nr"; }
    constexpr static char const * tElemsString() { return "nt"; }
    constexpr static char const * rBiasString() { return "rBias"; }
    constexpr static char const * trajectoryString() { return "trajectory"; }
    constexpr static char const * cartesianOuterBoundaryString() { return "useCartesianOuterBoundary"; }
    constexpr static char const * autoSpaceRadialElemsString() { return "autoSpaceRadialElems"; }
  };
  /// @endcond

  void postProcessInput() override final;

private:

  /// Trajectory defined by coordinates of the wellbore centers
  array2d< real64 > m_trajectory;

  int m_cartesianOuterBoundary;

  bool m_isFullAnnulus;

  array1d< int > m_autoSpaceRadialElems;

  array1d< real64 > & m_radialCoords;

};

} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP */
