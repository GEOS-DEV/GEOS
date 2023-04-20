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
 * @file CompositionalMultiphaseStatistics.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICS_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"

namespace geos
{

class CompositionalMultiphaseBase;

/**
 * @class CompositionalMultiphaseStatistics
 *
 * Task class allowing for the computation of aggregate statistics in compositional multiphase simulations
 */
class CompositionalMultiphaseStatistics : public FieldStatisticsBase< CompositionalMultiphaseBase >
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  CompositionalMultiphaseStatistics( const string & name,
                                     Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "CompositionalMultiphaseStatistics"; }

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
    /// String for the flag deciding the computation of the CFL numbers
    constexpr static char const * computeCFLNumbersString() { return "computeCFLNumbers"; }
    /// String for the flag deciding the computation of the region statistics
    constexpr static char const * computeRegionStatisticsString() { return "computeRegionStatistics"; }
    /// String for the region statistics
    constexpr static char const * regionStatisticsString() { return "regionStatistics"; }
    /// String for the relperm threshold
    constexpr static char const * relpermThresholdString() { return "relpermThreshold"; }
  };

  class RegionStatistics : public dataRepository::Group
  {
public:
    RegionStatistics( string const & name,
                      Group * const parent )
      : Group( name, parent )
    {
      //setInputFlags( dataRepository::InputFlags::REQUIRED );

      registerWrapper( viewKeyStruct::phasePoreVolumeString(), &m_phasePoreVolume ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "Phase region phase pore volume" );


      registerWrapper( viewKeyStruct::phaseMassString(), &m_phaseMass ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "Region phase mass (trapped and non-trapped, immobile and mobile)" );

      registerWrapper( viewKeyStruct::trappedPhaseMassString(), &m_trappedPhaseMass ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "Trapped region phase mass" );

      registerWrapper( viewKeyStruct::immobilePhaseMassString(), &m_immobilePhaseMass ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "Immobile region phase mass" );

      registerWrapper( viewKeyStruct::dissolvedComponentMassString(), &m_dissolvedComponentMass ).
        setApplyDefaultValue( 0 ).
        //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
        setDescription( "Dissolved region component mass" );

    }

    struct viewKeyStruct
    {
      constexpr static char const * phasePoreVolumeString() { return "phasePoreVolume"; }
      constexpr static char const * phaseMassString() { return "phaseMass"; }
      constexpr static char const * trappedPhaseMassString() { return "trappedPhaseMass"; }
      constexpr static char const * immobilePhaseMassString() { return "immobilePhaseMass"; }
      constexpr static char const * dissolvedComponentMassString() { return "dissolvedComponentMass"; }
    };

    void init( integer const numPhases, integer const numComps )
    {
      m_phasePoreVolume.resizeDimension< 0 >( numPhases );
      m_phaseMass.resizeDimension< 0 >( numPhases );
      m_trappedPhaseMass.resizeDimension< 0 >( numPhases );
      m_immobilePhaseMass.resizeDimension< 0 >( numPhases );
      m_dissolvedComponentMass.resizeDimension< 0, 1 >( numPhases, numComps );
    }

private:
    RegionStatistics() = delete;



    // /// average region pressure
    // real64 averagePressure;
    // /// minimum region pressure
    // real64 minPressure;
    // /// maximum region pressure
    // real64 maxPressure;

    // /// minimum region delta pressure
    // real64 minDeltaPressure;
    // /// maximum region delta pressure
    // real64 maxDeltaPressure;

    // /// average region temperature
    // real64 averageTemperature;
    // /// minimum region temperature
    // real64 minTemperature;
    // /// maximum region temperature
    // real64 maxTemperature;

    // /// total region pore volume
    // real64 totalPoreVolume;
    // /// total region uncompacted pore volume
    // real64 totalUncompactedPoreVolume;

    /// phase region phase pore volume
    array1d< real64 > m_phasePoreVolume;

    /// region phase mass (trapped and non-trapped, immobile and mobile)
    array1d< real64 > m_phaseMass;
    /// trapped region phase mass
    array1d< real64 > m_trappedPhaseMass;
    /// immobile region phase mass
    array1d< real64 > m_immobilePhaseMass;
    /// dissolved region component mass
    array2d< real64 > m_dissolvedComponentMass;
  };

  /**
   * @brief Compute some statistics on the reservoir (average field pressure, etc)
   * @param[in] mesh the mesh level object
   * @param[in] regionNames the array of target region names
   */
  void computeRegionStatistics( MeshLevel & mesh,
                                arrayView1d< string const > const & regionNames );

  /**
   * @brief Compute CFL numbers
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   */
  void computeCFLNumbers( real64 const & dt,
                          DomainPartition & domain );

  void postProcessInput() override;

  void registerDataOnMesh( Group & meshBodies ) override;

  /// Flag to decide whether CFL numbers are computed or not
  integer m_computeCFLNumbers;

  /// Flag to decide whether region statistics are computed or not
  integer m_computeRegionStatistics;

  /// Threshold to decide whether a phase is considered "mobile" or not
  real64 m_relpermThreshold;

};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICS_HPP_ */
