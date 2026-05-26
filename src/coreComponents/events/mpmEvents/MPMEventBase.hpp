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
 * @file MPMEventBase.hpp
 */

#ifndef GEOSX_MPMEVENTBASE_HPP_
#define GEOSX_MPMEVENTBASE_HPP_

#include "dataRepository/Group.hpp"

namespace geos
{

/**
 * @class MPMEventBase
 *
 * This class implements the abstract base mpm event for the solid mechanics material point method solver
 */
class MPMEventBase : public dataRepository::Group
{
public:
  /**
   * Constructor
   * @param time at which event starts
   * @param interval time interval between during which the event is performed
   */
  MPMEventBase( string const & name,
                Group * const parent );

  /// Destructor
  virtual ~MPMEventBase() override;


  /**
   * @copydoc dataRepository::Group::createChild()
   *
   * An event may have an arbitrary number of sub-events defined as children in the input xml.
   * e.g.:
   * @code{.unparsed}
   * <Events>
   *   <PeriodicEvent name="base_event" ...>
   *     <PeriodicEvent name="sub_event" .../>
   *     ...
   *   </PeriodicEvent>
   * </Events>
   * @endcode
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;


  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void expandObjectCatalogs() override;


  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "MPMEventBase"; }

  virtual string getCatalogName() const = 0;

  virtual void postInputInitialization() override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * startTimeString() { return "startTime"; }
    static constexpr char const * endTimeString() { return "endTime"; }
    static constexpr char const * delayString() { return "delay"; }
    static constexpr char const * durationString() { return "duration"; }
    static constexpr char const * dependenciesString() { return "dependencies"; }
    static constexpr char const * hasStartedString() { return "hasStarted"; }
    static constexpr char const * isCompleteString() { return "isComplete"; }

    dataRepository::ViewKey startTime = { startTimeString() };
    dataRepository::ViewKey endTime = { endTimeString() };
    dataRepository::ViewKey isComplete = { isCompleteString() };
  } viewKeys;
  /// @endcond


  /// Catalog interface
  using CatalogInterface = dataRepository::CatalogInterface< MPMEventBase, string const &, Group * const >;
  /// @copydoc dataRepository::Group::getCatalog()
  static CatalogInterface::CatalogType & getCatalog();

 void setStartTime( real64 startTime ) { m_startTime = startTime; }
  void setEndTime( real64 endTime ) { m_endTime = endTime; }
 
  real64 getStartTime() const { return m_startTime; }
  real64 getEndTime() const { return m_endTime; }
  real64 getTimeInterval() const { return m_endTime - m_startTime; }
  real64 getDelay() const { return m_delay; }
  real64 getDuration() const { return m_duration; }

  int hasStarted() const { return m_hasStarted; }
  void setHasStarted( int hasStarted ) { m_hasStarted = hasStarted; }
  int isComplete() const { return m_isComplete; }
  void setIsComplete( int isComplete ) { m_isComplete = isComplete; }

  bool hasDependencies() const { return m_dependencies.size() > 0; }
  string_array getDependencies() { return m_dependencies; }

protected:
  // Event variables
  real64 m_startTime;
  real64 m_endTime;

  real64 m_delay;
  real64 m_duration;
  string_array m_dependencies;

  int m_hasStarted;
  int m_isComplete;
};

} /* namespace geos */

#endif /* GEOSX_MPMEVENTBASE_HPP_ */
