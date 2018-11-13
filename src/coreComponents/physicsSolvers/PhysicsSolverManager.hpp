/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_

#include "../systemSolverInterface/LinearSystemRepository.hpp"
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

  struct viewKeyStruct
  {
      constexpr static auto gravityVectorString = "gravityVector";
      constexpr static auto blockSystemRepositoryString = "blockSystemRepository";
  } viewKeys;

  struct groupKeyStruct
  {
  } groupKeys;

  R1Tensor const & gravityVector() const { return m_gravityVector; }
  R1Tensor       & gravityVector()       { return m_gravityVector; }



private:
  PhysicsSolverManager() = delete;

  R1Tensor m_gravityVector;

  /// this is a block structured linear system object used to hold the system
  systemSolverInterface::LinearSystemRepository m_blockSystemRepository;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_PHYSICSSOLVERMANAGER_HPP_ */
