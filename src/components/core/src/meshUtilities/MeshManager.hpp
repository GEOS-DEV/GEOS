// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * MeshManager.hpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_MESHMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_MESHMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{
class SolverBase;

class MeshManager : public dataRepository::ManagedGroup
{
public:
  MeshManager( std::string const & name,
               ManagedGroup * const parent );

  virtual ~MeshManager() override;

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  void GenerateMeshes( DomainPartition * const domain );
  void GenerateMeshLevels( DomainPartition * const domain );

private:
  MeshManager() = delete;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_MESHMANAGER_HPP_ */
