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
 * @file OutputManager.hpp
 */

#ifndef GEOSX_MANAGERS_OUTPUTS_OUTPUTMANAGER_HPP_
#define GEOSX_MANAGERS_OUTPUTS_OUTPUTMANAGER_HPP_

#include "dataRepository/Group.hpp"


namespace geosx
{
namespace dataRepository
{
namespace keys
{}
}

/**
 * @class OutputManager
 *
 * An class for managing output types
 */
class OutputManager : public dataRepository::Group
{
public:
  /**
   * @brief Main constructor.
   * @param name The name of the object in the data repository.
   * @param parent The parent of this object in the data repository.
   **/
  OutputManager( std::string const & name,
                 Group * const parent );

  /// Destructor
  virtual ~OutputManager() override;

  /**
   * @brief Construct and retrieve a pointer to a registered group by catalog key.
   * @param childKey A catalog key denoting the group to create.
   * @param childName The name of the newly created group.
   * @return A pointer to the newly created group.
   */
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    dataRepository::ViewKey time = { "time" };
  } viewKeys;
  /// @endcond
};


} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_OUTPUTMANAGER_HPP_ */
