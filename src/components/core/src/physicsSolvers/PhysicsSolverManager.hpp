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
 * PhysicsSolverManager.hpp
 *
 *  Created on: Sep 7, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"

namespace pugi
{
class xml_node;
}

namespace geosx
{
class SolverBase;

class PhysicsSolverManager : public dataRepository::ManagedGroup
{
public:
  PhysicsSolverManager( std::string const & name,
                        ManagedGroup * const parent );

  virtual ~PhysicsSolverManager() override;

  virtual void FillDocumentationNode() override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

private:
  PhysicsSolverManager() = delete;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_ */
