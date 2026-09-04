/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MixedMimeticDiscretizationManager.hpp
 */

#ifndef GEOS_MIXEDMIMETIC_MIXEDMIMETICDISCRETIZATIONMANAGER_HPP_
#define GEOS_MIXEDMIMETIC_MIXEDMIMETICDISCRETIZATIONMANAGER_HPP_

#include "dataRepository/Group.hpp"

namespace geos
{

class MixedMimeticDiscretization;

/**
 * @class MixedMimeticDiscretizationManager
 *
 * Contains the mixed mimetic discretizations, registered under the NumericalMethods block.
 */
class MixedMimeticDiscretizationManager : public dataRepository::Group
{
public:

  MixedMimeticDiscretizationManager() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the MixedMimeticDiscretizationManager in the data repository
   * @param parent the parent group of this group.
   */
  MixedMimeticDiscretizationManager( string const & name, Group * const parent );

  virtual ~MixedMimeticDiscretizationManager() override;

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  virtual void expandObjectCatalogs() override;

  /**
   * @brief Get a mixed mimetic discretization by name, const version.
   * @param name the name of the discretization
   * @return the discretization object
   */
  MixedMimeticDiscretization const & getMixedMimeticDiscretization( string const & name ) const;

  /**
   * @brief Get a mixed mimetic discretization by name.
   * @param name the name of the discretization
   * @return the discretization object
   */
  MixedMimeticDiscretization & getMixedMimeticDiscretization( string const & name );

};

} // namespace geos

#endif //GEOS_MIXEDMIMETIC_MIXEDMIMETICDISCRETIZATIONMANAGER_HPP_
