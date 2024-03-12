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
 * @file RegionMultiphaseStatistics.hpp
 */

#ifndef GEOS_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_REGIONMULTIPHASESTATISTICS_HPP_
#define GEOS_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_REGIONMULTIPHASESTATISTICS_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"

#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

class CompositionalMultiphaseBase;

/**
 * @class RegionMultiphaseStatistics
 *
 * Task class allowing for the computation of aggregate statistics in compositional multiphase simulations
 */
class RegionMultiphaseStatistics : public FieldStatisticsBase< CompositionalMultiphaseBase >
{
public:
  enum PropertyNameType : integer
  {
    Pressure,
    Temperature,
    PoreVolume,
    VolumeFraction,
    Mass
  };

public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  RegionMultiphaseStatistics( const string & name,
                              Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "RegionMultiphaseStatistics"; }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**@}*/

private:

  using Base = FieldStatisticsBase< CompositionalMultiphaseBase >;

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for the list of region names
    constexpr static char const * regionNamesString() { return "regionNames"; }
    /// String for the list of region indices
    constexpr static char const * regionIdentifiersString() { return "regionIdentifiers"; }
    /// String for the list of properties to report
    constexpr static char const * propertyNamesString() { return "propertyNames"; }
    /// String for the suffix of the field name
    constexpr static char const * fieldNameString() { return "_region"; }
    /// String for the file name
    constexpr static char const * fileNameString() { return "values"; }
  };

  void postProcessInput() override;

  void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @brief Compute some statistics on the regions
   * @param[in] time current time
   * @param[in] mesh the mesh level object
   * @param[in] regionNames the array of target region names
   */
  void computeRegionStatistics( real64 const time,
                                MeshLevel & mesh,
                                arrayView1d< string const > const & regionNames ) const;

  void initializeFile() const;

  template< typename ARRAY >
  std::ostream & writeArray( std::ostream & os, ARRAY const & array ) const;

  struct RegionStatistics;
  struct RegionStatisticsKernel;

  // The names of regions
  string_array m_regionNames;

  // The ids of regions
  array1d< real64 > m_regionIdentifiers;

  // The names of properties to output
  string_array m_propertyNames;

  // The names converted into integers
  array1d< PropertyNameType > m_propertyNameTypes;
};

ENUM_STRINGS( RegionMultiphaseStatistics::PropertyNameType,
              "pressure",
              "temperature",
              "poreVolume",
              "volumeFraction",
              "mass" );

} /* namespace geos */

#endif /* GEOS_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_REGIONMULTIPHASESTATISTICS_HPP_ */
