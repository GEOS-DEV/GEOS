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
 * DirichletBoundaryCondition.hpp
 *
 *  Created on: Jun 2, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_DIRICHLETBOUNDARYCONDITION_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_DIRICHLETBOUNDARYCONDITION_HPP_

#include "BoundaryConditionBase.hpp"

namespace geosx
{

class DirichletBoundaryCondition : public BoundaryConditionBase
{
public:
  DirichletBoundaryCondition( string const & name, dataRepository::ManagedGroup *const parent );
  DirichletBoundaryCondition() = delete;
  virtual ~DirichletBoundaryCondition();

  static string CatalogName() { return "Dirichlet"; }



};



} /* namespace geosx */

#endif /*
          SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_DIRICHLETBOUNDARYCONDITION_HPP_
        */
