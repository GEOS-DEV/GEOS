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
 * @file ParticleRegion.hpp
 *
 */

#ifndef GEOSX_MESH_PARTICLEREGION_HPP_
#define GEOSX_MESH_PARTICLEREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{
class EdgeManager;
class EmbeddedSurfaceGenerator;

/**
 * @class ParticleRegion
 *
 * The ParticleRegion class contains the functionality to support the concept of a ParticleRegion in the element
 * hierarchy. ParticleRegion derives from ElementRegionBase and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class ParticleRegion : public ElementRegionBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{


  /**
   * @brief Constructor.
   * @param name the name of the object in the data hierarchy.
   * @param parent a pointer to the parent group in the data hierarchy.
   */
  ParticleRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  ParticleRegion() = delete;

  /**
   * @brief Destructor.
   */
  virtual ~ParticleRegion() override;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string catalogName()
  { return "ParticleRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual const string getCatalogName() const override final
  { return ParticleRegion::catalogName(); }

  ///@}

  /**
   * @name Generation of the particle region
   */
  ///@{

  /**
   * @brief Add a cellBlockRegion name to the list.
   * @param cellBlockName string containing the cell block region name.
   */
  void addCellBlockName( string const & cellBlockName )
  {
    m_cellBlockNames.emplace_back( cellBlockName );
  }

  virtual void generateMesh( Group & cellBlocks ) override;


  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// @return String key for the cell block names
    static constexpr char const * sourceCellBlockNamesString() {return "cellBlocks"; }
  };


private:

  // Cell block names
  string_array m_cellBlockNames;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_PARTICLEREGION_HPP_ */
