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
 * @file TemperatureProfileMPMEvent.hpp
 */

#ifndef GEOSX_TEMPERATUREPROFILE_MPMEVENT_HPP_
#define GEOSX_TEMPERATUREPROFILE_MPMEVENT_HPP_

#include "MPMEventBase.hpp"
#include "physicsSolvers/solidMechanics/MPMSolverEnums.hpp"

namespace geos
{

/**
 * @class TemperatureProfileMPMEvent
 *
 * This class implements a uniform temperature-table event for the solid mechanics material point method solver
 */
class TemperatureProfileMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  TemperatureProfileMPMEvent( const string & name,
                              Group * const parent );

  /// Destructor
  virtual ~TemperatureProfileMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "TemperatureProfile"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * temperatureTableString() { return "temperatureTable"; }
    static constexpr char const * interpolationTypeString() { return "interpolationType"; }
  } TemperatureProfileMPMEventViewKeys;
  /// @endcond

  arrayView2d< real64 const > getTemperatureTable() const
  {
    return m_temperatureTable;
  }

  mpm::InterpolationOption getInterpolationType() const
  {
    return m_interpolationType;
  }

protected:
  virtual void postInputInitialization() override final;

  array2d< real64 > m_temperatureTable;
  mpm::InterpolationOption m_interpolationType;
};

} /* namespace geos */

#endif /* GEOSX_TEMPERATUREPROFILE_MPMEVENT_HPP_ */