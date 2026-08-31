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
 * @file MixedVEMManager.hpp
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMMANAGER_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMMANAGER_HPP_

#include "dataRepository/Group.hpp"

namespace geos
{

class MixedVEMDiscretization;

/**
 * @class MixedVEMManager
 *
 * Class managing the mixed virtual element discretizations declared in the input file.
 */
class MixedVEMManager : public dataRepository::Group
{
public:

  MixedVEMManager() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the MixedVEMManager in the data repository
   * @param parent the parent group of this group.
   */
  MixedVEMManager( string const & name, Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~MixedVEMManager() override;

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  virtual void expandObjectCatalogs() override;

  /**
   * @brief Return the discretization associated with the provided name.
   * @param[in] name the provided name
   * @return the discretization associated with the provided name
   */
  MixedVEMDiscretization const & getDiscretization( string const & name ) const;

  /**
   * @copydoc getDiscretization( string const & ) const
   */
  MixedVEMDiscretization & getDiscretization( string const & name );

};

} // namespace geos

#endif // GEOS_MIXEDVEM_MIXEDVEMMANAGER_HPP_
