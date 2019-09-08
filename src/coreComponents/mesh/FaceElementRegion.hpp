/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file FaceElementRegion.hpp
 *
 */

#ifndef CORECOMPONENTS_MESH_FACEELEMENTREGION_HPP_
#define CORECOMPONENTS_MESH_FACEELEMENTREGION_HPP_

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
  localIndex AddToFractureMesh( EdgeManager * const edgeManager,
                                FaceManager const * const faceManager,
                                array1d< array1d<localIndex> > const & originalFaceToEdges,
                                string const & subRegionName,
                                localIndex const faceIndices[2] );


  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    static constexpr auto fractureSetString = "fractureSet";
  };


private:

};

} /* namespace geosx */

#endif /* CORECOMPONENTS_MESH_FACEELEMENTREGION_HPP_ */
