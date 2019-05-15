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
 * SourceFluxBoundaryCondition.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_SOURCEFLUXBOUNDARYCONDITION_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_SOURCEFLUXBOUNDARYCONDITION_HPP_

#include "FieldSpecificationBase.hpp"

namespace geosx
{

class SourceFluxBoundaryCondition : public FieldSpecificationBase
{
public:
  SourceFluxBoundaryCondition( string const & name, dataRepository::ManagedGroup *const parent );
  SourceFluxBoundaryCondition() = delete;
  virtual ~SourceFluxBoundaryCondition() override;

  virtual void InitializePreSubGroups( ManagedGroup * const ) override;

  static string CatalogName() { return "SourceFlux"; }

  virtual const string getCatalogName() const override
  {
    return SourceFluxBoundaryCondition::CatalogName();
  }

};



} /* namespace geosx */

#endif /*
          SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_SOURCEFLUXBOUNDARYCONDITION_HPP_
        */
