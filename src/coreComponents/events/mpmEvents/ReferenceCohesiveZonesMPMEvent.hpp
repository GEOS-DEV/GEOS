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
 * @file ReferenceCohesiveZoneMPMEvent.hpp
 */

#ifndef GEOSX_REFERENCECOHESIVEZONES_MPMEVENT_HPP_
#define GEOSX_REFERENCECOHESIVEZONES_MPMEVENT_HPP_

#include "MPMEventBase.hpp"

#include "physicsSolvers/solidMechanics/MPMSolverEnums.hpp"

namespace geos
{

/**
 * @class ReferenceCohesiveZonesMPMEvent
 *
 * This class implements the material swap mpm event for the solid mechanics material point method solver
 */
class ReferenceCohesiveZonesMPMEvent : public MPMEventBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  ReferenceCohesiveZonesMPMEvent( const string & name,
                                 Group * const parent );

  /// Destructor
  virtual ~ReferenceCohesiveZonesMPMEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "ReferenceCohesiveZones"; }

  virtual string getCatalogName() const override { return catalogName(); }

  string_array const & getRegionNames() { return m_regionNames; }

  //  /**
  //   * @brief Create a new CohesiveZoneBase object as a child of this group.
  //   * @param childKey catalog key of the new CohesiveZoneBase derived type to create
  //   * @param childName name of the new CohesiveZoneBase object
  //   * @return pointer to the created CohesiveZoneBase object
  //   */
  // virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  { 
    static constexpr char const * regionNamesString() { return "regionNames"; }
    static constexpr char const * constitutiveModelsString() { return "constitutiveModels"; }
    static constexpr char const * czTagsString() { return "czTags"; }
    static constexpr char const * czVolumeNormalizationString() { return "czVolumeNormalization"; }
    static constexpr char const * computeNormalsAndPositionsString() { return "computeNormalsAndPositions"; }
    static constexpr char const * normalsAndPositionsMethodString() { return "normalsAndPositionsMethod"; }
    static constexpr char const * czSurfaceDisplacementUpdateString() { return "czSurfaceDisplacementUpdate"; }
  } CohesiveZoneMPMEventViewKeys;
  /// @endcond

  int getCZVolumeNormalization() const { return m_czVolumeNormalization; }
  int getComputeNormalsAndPositions() const { return m_computeNormalsAndPositions; }
  mpm::NormalsAndPositionsMethodOption getNormalsAndPositionsMethod() const { return m_normalsAndPositionsMethod; }
  mpm::CohesiveSurfaceDisplacementUpdateOption getCZSurfaceDisplacementUpdate() const { return m_czSurfaceDisplacementUpdate; }

private:

protected:
  virtual void postInputInitialization() override final;

  string_array m_regionNames;

  int m_czVolumeNormalization;
  int m_computeNormalsAndPositions;
  mpm::NormalsAndPositionsMethodOption m_normalsAndPositionsMethod;
  mpm::CohesiveSurfaceDisplacementUpdateOption m_czSurfaceDisplacementUpdate;
};

} /* namespace geos */

#endif /* GEOSX_REFERENCECOHESIVEZONES_MPMEVENT_HPP_ */
