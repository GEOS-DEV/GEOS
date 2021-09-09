/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DirichletBoundaryCondition.hpp
 */

#ifndef GEOSX_FIELDSPECIFICATION_DIRICHLETBOUNDARYCONDITION_HPP_
#define GEOSX_FIELDSPECIFICATION_DIRICHLETBOUNDARYCONDITION_HPP_

#include "FieldSpecificationBase.hpp"

namespace geosx
{

/**
 * @class DirichletBoundaryCondition
 * A class to manage Dirichlet boundary conditions
 */
class DirichletBoundaryCondition : public FieldSpecificationBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group( string const & name, Group * const parent )
  DirichletBoundaryCondition( string const & name, dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  DirichletBoundaryCondition() = delete;

  /**
   * @brief destructor
   */
  virtual ~DirichletBoundaryCondition();

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "Dirichlet"; }

  virtual const string getCatalogName() const
  {
    return DirichletBoundaryCondition::catalogName();
  }

};

} /* namespace geosx */

#endif /* GEOSX_FIELDSPECIFICATION_DIRICHLETBOUNDARYCONDITION_HPP_ */
