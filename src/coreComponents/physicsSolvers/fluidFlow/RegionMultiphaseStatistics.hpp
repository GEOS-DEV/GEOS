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
    PoreVolume,
    Saturation
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
  };

   struct RegionStatistics
   {
    explicit RegionStatistics(integer regionCount);
    /// average region pressure
    real64 averagePressure;
    /// minimum region pressure
    real64 minPressure;
    /// maximum region pressure
    real64 maxPressure;

    /// minimum region delta pressure
    real64 minDeltaPressure;
    /// maximum region delta pressure
    real64 maxDeltaPressure;

    /// average region temperature
    real64 averageTemperature;
    /// minimum region temperature
    real64 minTemperature;
    /// maximum region temperature
    real64 maxTemperature;

    /// total region pore volume
    real64 totalPoreVolume;
    /// total region uncompacted pore volume
    real64 totalUncompactedPoreVolume;
    /// phase region phase pore volume
    array1d< real64 > phasePoreVolume;

    /// region phase mass (trapped and non-trapped, immobile and mobile)
    array1d< real64 > phaseMass;
    /// trapped region phase mass
    array1d< real64 > trappedPhaseMass;
    /// immobile region phase mass
    array1d< real64 > immobilePhaseMass;
    /// region component mass
    array2d< real64 > componentMass;
   };

  void postProcessInput() override;

  void registerDataOnMesh( Group & meshBodies ) override;

  // The names of regions
  string_array m_regionNames;

  // The ids of regions
  array1d<localIndex> m_regionIdentifiers;

  // The names of properties to output
  string_array m_propertyNames;

  // The names converted into integers
  array1d< PropertyNameType > m_propertyNameTypes;
};

ENUM_STRINGS( RegionMultiphaseStatistics::PropertyNameType,
              "pressure",
              "poreVolume",
              "saturation" );

} /* namespace geos */

#endif /* GEOS_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_REGIONMULTIPHASESTATISTICS_HPP_ */
