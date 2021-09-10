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
 * @file GeometricObjectManager.hpp
 */

#ifndef GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_GEOMETRICOBJECTMANAGER_HPP_
#define GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_GEOMETRICOBJECTMANAGER_HPP_

#include "dataRepository/Group.hpp"


namespace geosx
{

/**
 * GeometricObjectManager
 * @brief Manager of the simple geometric objects
 */
class GeometricObjectManager : public dataRepository::Group
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  GeometricObjectManager( string const & name,
                          Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~GeometricObjectManager() override;

  ///@}

  /**
   * @brief @return a pointer to this GeometricObjectManager
   */
  static GeometricObjectManager & getInstance();

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  virtual void expandObjectCatalogs() override;

private:
  GeometricObjectManager() = delete;
  static GeometricObjectManager * m_instance;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_SIMPLEGEOMETRICOBJECTS_GEOMETRICOBJECTMANAGER_HPP_ */
