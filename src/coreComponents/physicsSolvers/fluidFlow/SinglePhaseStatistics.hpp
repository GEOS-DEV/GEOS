/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseStatistics.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"

namespace geos
{

class SinglePhaseBase;

/**
 * @class SinglePhaseStatistics
 *
 * Task class allowing for the computation of aggregate statistics in single-phase simulations
 */
class SinglePhaseStatistics : public FieldStatisticsBase< SinglePhaseBase >
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  SinglePhaseStatistics( const string & name,
                         Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "SinglePhaseStatistics"; }

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
  class RegionStatistics : public dataRepository::Group
  {
public:
    RegionStatistics( string const & name,
                      Group * const parent )
      : Group( name, parent )
    {
      registerWrapper( viewKeyStruct::averagePressureString(), &m_averagePressure ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "average region pressure" );

      registerWrapper( viewKeyStruct::minPressureString(), &m_minPressure ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "minimum region pressure" );

      registerWrapper( viewKeyStruct::maxPressureString(), &m_maxPressure ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "maximum region pressure" );

      registerWrapper( viewKeyStruct::minDeltaPressureString(), &m_minDeltaPressure ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "minimum region delta pressure" );

      registerWrapper( viewKeyStruct::maxDeltaPressureString(), &m_maxDeltaPressure ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "maximum region delta pressure" );

      registerWrapper( viewKeyStruct::maxDeltaPressureString(), &m_totalMass ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "fluid mass" );


      registerWrapper( viewKeyStruct::averageTemperatureString(), &m_averageTemperature ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "average region temperature" );

      registerWrapper( viewKeyStruct::minTemperatureString(), &m_minTemperature ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "minimum region temperature" );

      registerWrapper( viewKeyStruct::maxTemperatureString(), &m_maxTemperature ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "maximum region temperature" );


      registerWrapper( viewKeyStruct::totalPoreVolumeString(), &m_totalPoreVolume ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "total region pore volume" );

      registerWrapper( viewKeyStruct::totalUncompactedPoreVolumeString(), &m_totalUncompactedPoreVolume ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "total region uncompacted pore volume" );

    }

    struct viewKeyStruct
    {
      constexpr static char const * averagePressureString() { return "averagePressure"; }
      constexpr static char const * minPressureString() { return "minPressure"; }
      constexpr static char const * maxPressureString() { return "maxPressure"; }

      constexpr static char const * minDeltaPressureString() { return "minDeltaPressure"; }
      constexpr static char const * maxDeltaPressureString() { return "maxDeltaPressure"; }

      constexpr static char const * totalMassString() { return "totalMass"; }

      constexpr static char const * averageTemperatureString() { return "averageTemperature"; }
      constexpr static char const * minTemperatureString() { return "minTemperature"; }
      constexpr static char const * maxTemperatureString() { return "maxTemperature"; }

      constexpr static char const * totalPoreVolumeString() { return "totalPoreVolume"; }
      constexpr static char const * totalUncompactedPoreVolumeString() { return "totalUncompactedPoreVolume"; }
    };


private:
    /// average region pressure
    real64 m_averagePressure;
    /// minimum region pressure
    real64 m_minPressure;
    /// maximum region pressure
    real64 m_maxPressure;

    /// minimum region delta pressure
    real64 m_minDeltaPressure;
    /// maximum region delta pressure
    real64 m_maxDeltaPressure;

    // fluid mass
    real64 m_totalMass;

    /// average region temperature
    real64 m_averageTemperature;
    /// minimum region temperature
    real64 m_minTemperature;
    /// maximum region temperature
    real64 m_maxTemperature;

    /// total region pore volume
    real64 m_totalPoreVolume;
    /// total region uncompacted pore volume
    real64 m_totalUncompactedPoreVolume;
    /// phase region phase pore volume
    array1d< real64 > m_phasePoreVolume;

  };

  /**@}*/

private:

  using Base = FieldStatisticsBase< SinglePhaseBase >;

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for the region statistics
    constexpr static char const * regionStatisticsString() { return "regionStatistics"; }
  };


  /**
   * @brief Compute some statistics on the reservoir (average field pressure, etc)
   * @param[in] mesh the mesh level object
   * @param[in] regionNames the array of target region names
   */
  void computeRegionStatistics( real64 const time,
                                MeshLevel & mesh,
                                arrayView1d< string const > const & regionNames );


  void registerDataOnMesh( Group & meshBodies ) override;

};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_ */
