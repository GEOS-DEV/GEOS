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
 * @file DeformationUpdateMPMEvent.hpp
 */

#ifndef GEOSX_DEFORMATIONUPDATE_MPMEVENT_HPP_
#define GEOSX_DEFORMATIONUPDATE_MPMEVENT_HPP_

#include "MPMEventBase.hpp"
#include "physicsSolvers/solidMechanics/MPMSolverEnums.hpp"

namespace geos
{

/**
 * @class DeformationUpdateMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class DeformationUpdateMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  DeformationUpdateMPMEvent( const string & name,
                             Group * const parent );

  /// Destructor
  virtual ~DeformationUpdateMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "DeformationUpdate"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * relativeDeformationString() { return "relativeDeformation"; }

    static constexpr char const * prescribedFTableString() { return "prescribedFTable"; }
    static constexpr char const * prescribedBoundaryFTableString() { return "prescribedBoundaryFTable"; }
    static constexpr char const * stressControlString() { return "stressControl"; }

    static constexpr char const * fTableInterpTypeString() { return "fTableInterpType"; }
    static constexpr char const * stressTableInterpTypeString() { return "stressTableInterpType"; }
    
    static constexpr char const * fTableString() { return "fTable"; }
    static constexpr char const * stressTableString() { return "stressTable"; }

  } DeformationUpdateMPMEventViewKeys;
  /// @endcond

  int getRelativeDeformation() const { return m_relativeDeformation; }

  int getPrescribedBoundaryFTable() const { return m_prescribedBoundaryFTable; }
  int getPrescribedFTable() const { return m_prescribedFTable; }
  array1d< int > getStressControl() const { return m_stressControl; }

  mpm::InterpolationOption getFTableInterpolation() const { return m_fTableInterpType; }
  mpm::InterpolationOption getStressTableInterpolation() const { return m_stressTableInterpType; }

  arrayView2d< real64 const > getFTable() const { return m_fTable; }
  arrayView2d< real64 const > getStressTable() const { return m_stressTable; }

protected:
  virtual void postInputInitialization() override final;

  // Event variables
  int m_relativeDeformation;

  int m_prescribedFTable;
  int m_prescribedBoundaryFTable;
  array1d< int > m_stressControl;

  mpm::InterpolationOption  m_fTableInterpType;
  mpm::InterpolationOption  m_stressTableInterpType;

  array2d< real64 > m_fTable;
  array2d< real64 > m_stressTable;
};

} /* namespace geos */

#endif /* GEOSX_DEFORMATIONUPDATE_MPMEVENT_HPP_ */