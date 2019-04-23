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

/*
 * ComponentFluxBoundaryCondition.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_COMPONENTFLUXBOUNDARYCONDITION_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_COMPONENTFLUXBOUNDARYCONDITION_HPP_

#include "FieldSpecificationBase.hpp"

namespace geosx
{

class ComponentFluxBoundaryCondition : public FieldSpecificationBase
{
public:
  ComponentFluxBoundaryCondition( string const & name, dataRepository::ManagedGroup *const parent );
  ComponentFluxBoundaryCondition() = delete;
  virtual ~ComponentFluxBoundaryCondition();

  static string CatalogName() { return "ComponentFlux"; }

  virtual const string getCatalogName() const 
  {
    return ComponentFluxBoundaryCondition::CatalogName();
  }

};



} /* namespace geosx */

#endif /*
          SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_COMPONENTFLUXBOUNDARYCONDITION_HPP_
        */
