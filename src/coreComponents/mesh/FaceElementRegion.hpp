/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FaceElementRegion.hpp
 *
 */

#ifndef GEOSX_MESH_FACEELEMENTREGION_HPP_
#define GEOSX_MESH_FACEELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{

class EdgeManager;

/**
 * @class FaceElementRegion
 *
 * The FaceElementRegion class contains the functionality to support the concept of a FaceElementRegion in the element
 * hierarchy. FaceElementRegion derives from ElementRegion and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class FaceElementRegion : public ElementRegionBase
{
public:
  /**
   * @brief constructor
   * @param name The name of the object in the data hierarchy.
   * @param parent Pointer to the parent group in the data hierarchy.
   */
  FaceElementRegion( string const & name, Group * const parent );

  FaceElementRegion() = delete;
  virtual ~FaceElementRegion() override;

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string CatalogName()
  { return "FaceElementRegion"; }

  virtual const string getCatalogName() const override final
  { return FaceElementRegion::CatalogName(); }


  virtual void GenerateMesh( Group const * ) override {}

  /**
   * @brief This function generates and adds entries to the face/fracture mesh
   * @param faceManager A pointer to the FaceManager object.
   * @param subRegionName The name of the FaceElementSubRegion to insert the new entries.
   * @param faceIndices The local indices of the new faces that define the face element.
   * @return The local index of the new FaceElement entry.
   */
  localIndex AddToFractureMesh( real64 const time_np1,
                                EdgeManager * const edgeManager,
                                FaceManager const * const faceManager,
                                ArrayOfArraysView< localIndex const > const & originalFaceToEdges,
                                string const & subRegionName,
                                localIndex const faceIndices[2] );


  real64 getDefaultAperture() const { return m_defaultAperture; }


  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    static constexpr auto fractureSetString = "fractureSet";
    static constexpr auto defaultApertureString = "defaultAperture";
    constexpr static auto ruptureTimeString = "ruptureTime";

  };

protected:
  virtual void InitializePreSubGroups( Group * const ) override;


private:
  /// The
  real64 m_defaultAperture;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_FACEELEMENTREGION_HPP_ */
