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

/*
 * GeometricObjectManager.hpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOMETRICOBJECTMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOMETRICOBJECTMANAGER_HPP_

#include "dataRepository/Group.hpp"


namespace geosx
{

class GeometricObjectManager : public dataRepository::Group
{
public:
  GeometricObjectManager( std::string const & name,
                          Group * const parent );

  virtual ~GeometricObjectManager() override;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

private:
  GeometricObjectManager() = delete;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOMETRICOBJECTMANAGER_HPP_ */
