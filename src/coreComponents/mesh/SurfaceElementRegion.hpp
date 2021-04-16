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
 * @file SurfaceElementRegion.hpp
 *
 */

#ifndef GEOSX_MESH_SURFACEELEMENTREGION_HPP_
#define GEOSX_MESH_SURFACEELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"
#include "codingUtilities/EnumStrings.hpp"

namespace geosx
{


class EdgeManager;

/**
 * @class SurfaceElementRegion
 *
 * The SurfaceElementRegion class contains the functionality to support the concept of a SurfaceElementRegion in the
 * element hierarchy. SurfaceElementRegion derives from ElementRegion and has an entry in the ObjectManagerBase
 * catalog.
 */
class SurfaceElementRegion : public ElementRegionBase
{
public:

  /**
   * @enum SurfaceSubRegionType
   *
   * The options for the surface subregion type
   */
  enum class SurfaceSubRegionType : integer
  {
    faceElement,
    embeddedElement
  };

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name the name of the object in the data hierarchy.
   * @param parent a pointer to the parent group in the data hierarchy.
   */
  SurfaceElementRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  SurfaceElementRegion() = delete;

  /**
   * @brief Default destructor.
   */
  virtual ~SurfaceElementRegion() override;

  ///@}

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Get the key name for the SurfaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string catalogName()
  { return "SurfaceElementRegion"; }

  virtual const string getCatalogName() const override final
  { return SurfaceElementRegion::catalogName(); }

  ///@}


  /**
   * @name Generation of the face element mesh region
   */
  ///@{

  virtual void generateMesh( Group & ) override;

  /**
   * @brief This function generates and adds entries to the face/fracture mesh.
   * @param time_np1 rupture time
   * @param edgeManager pointer to the EdgeManager object.
   * @param faceManager pointer to the FaceManager object.
   * @param originalFaceToEdges face-to-edge map before the rupture.
   * @param subRegionName the name of the FaceElementSubRegion to insert the new entries.
   * @param faceIndices the local indices of the new faces that define the face element.
   * @return the local index of the new FaceElement entry.
   */
  localIndex addToFractureMesh( real64 const time_np1,
                                EdgeManager * const edgeManager,
                                FaceManager const * const faceManager,
                                ArrayOfArraysView< localIndex const > const & originalFaceToEdges,
                                string const & subRegionName,
                                localIndex const faceIndices[2] );

  ///@}


  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get default aperture value.
   * @return default aperture value
   */
  real64 getDefaultAperture() const { return m_defaultAperture; }

  /**
   * @brief Get subRegion type.
   * @return subRegion type
   */
  SurfaceSubRegionType subRegionType() const { return m_subRegionType; }


  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// @return subRegion type string
    static constexpr char const * subRegionTypeString() { return "subRegionType"; }

    /// @return Fracture set string
    static constexpr char const * fractureSetString() { return "fractureSet"; }

    /// @return Default fracture aperture
    static constexpr char const * defaultApertureString() { return "defaultAperture"; }

    /// @return Rupture time string
    static constexpr char const * ruptureTimeString() { return "ruptureTime"; }
  };

protected:
  virtual void initializePreSubGroups() override;

private:

  SurfaceSubRegionType m_subRegionType;

  real64 m_defaultAperture;

};

ENUM_STRINGS( SurfaceElementRegion::SurfaceSubRegionType, "faceElement", "embeddedElement" )

} /* namespace geosx */

#endif /* CORECOMPONENTS_MESH_SurfaceElementRegion_HPP_ */
