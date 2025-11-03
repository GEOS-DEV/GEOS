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
 * @file CohesiveZoneReferenceMPMEvent.hpp
 */

#ifndef GEOSX_COHESIVEZONEREFERENCE_MPMEVENT_HPP_
#define GEOSX_COHESIVEZONEREFERENCE_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

#include "physicsSolvers/solidMechanics/SolidMechanicsMPM.hpp"

namespace geos
{

/**
 * @class CohesiveZoneReferenceMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class CohesiveZoneReferenceMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  CohesiveZoneReferenceMPMEvent( const string & name,
                  Group * const parent );

  /// Destructor
  virtual ~CohesiveZoneReferenceMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "CohesiveZoneReference"; }

  virtual string getCatalogName() const override { return catalogName(); }

 /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {

  } CohesiveZoneReferenceMPMEventViewKeys;
  /// @endcond

  int getCZVolumeNormalization() const { return m_czVolumeNormalization; } 
  int getComputeNormalsAndPositions() const { return m_computeNormalsAndPositions; } 
  SolidMechanicsMPM::NormalsAndPositionsMethodOption getNormalsAndPositionsMethod() const { return m_normalsAndPositionsMethod; }

protected:
  virtual void postInputInitialization() override final;

  int m_czVolumeNormalization;
  int m_computeNormalsAndPositions;
  SolidMechanicsMPM::NormalsAndPositionsMethodOption m_normalsAndPositionsMethod;
};

} /* namespace geos */

#endif /* GEOSX_COHESIVEZONEREFERENCE_MPMEVENT_HPP_ */
