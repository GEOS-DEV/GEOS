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
 * @file DirichletBoundaryCondition.hpp
 */

#ifndef GEOSX_MANAGERS_FIELDSPECIFICATION_DIRICHLETBOUNDARYCONDITION_HPP_
#define GEOSX_MANAGERS_FIELDSPECIFICATION_DIRICHLETBOUNDARYCONDITION_HPP_

#include "FieldSpecificationBase.hpp"

namespace geosx
{

class DirichletBoundaryCondition : public FieldSpecificationBase
{
public:
  DirichletBoundaryCondition( string const & name, dataRepository::Group * const parent );
  DirichletBoundaryCondition() = delete;
  virtual ~DirichletBoundaryCondition();

  static string CatalogName() { return "Dirichlet"; }

  virtual const string getCatalogName() const
  {
    return DirichletBoundaryCondition::CatalogName();
  }

};



} /* namespace geosx */

#endif /* GEOSX_MANAGERS_FIELDSPECIFICATION_DIRICHLETBOUNDARYCONDITION_HPP_ */
