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


#ifndef GEOSX_PHYSICSSOLVERS_MESHMANAGER_HPP_
#define GEOSX_PHYSICSSOLVERS_MESHMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{
class SolverBase;

class MeshManager : public dataRepository::Group
{
public:
  MeshManager( std::string const & name,
               Group * const parent );

  virtual ~MeshManager() override;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

  void GenerateMeshes( DomainPartition * const domain );
  void GenerateMeshLevels( DomainPartition * const domain );

private:
  MeshManager() = delete;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MESHMANAGER_HPP_ */
