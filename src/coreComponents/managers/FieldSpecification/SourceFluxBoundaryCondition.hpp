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
  SourceFluxBoundaryCondition( string const & name, dataRepository::Group *const parent );
  SourceFluxBoundaryCondition() = delete;
  virtual ~SourceFluxBoundaryCondition() override;

  virtual void InitializePreSubGroups( Group * const ) override;

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
