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
 * @file GeometricObjectManager.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_GEOMETRICOBJECTMANAGER_HPP_
#define GEOSX_PHYSICSSOLVERS_GEOMETRICOBJECTMANAGER_HPP_

#include "dataRepository/Group.hpp"


namespace geosx
{

/**
 * @class GeometricObjectManager
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
  GeometricObjectManager( std::string const & name,
                          Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~GeometricObjectManager() override;

  ///@}

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /**
   * @brief This function is used to expand any catalogs in the data structure.
   */
  virtual void ExpandObjectCatalogs() override;

private:
  GeometricObjectManager() = delete;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_GEOMETRICOBJECTMANAGER_HPP_ */
