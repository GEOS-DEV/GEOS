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
 * @file OutputManager.hpp
 */

#ifndef GEOSX_FILEIO_OUTPUTS_OUTPUTMANAGER_HPP_
#define GEOSX_FILEIO_OUTPUTS_OUTPUTMANAGER_HPP_

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
  /// @copydoc geosx::dataRepository::Group::Group( string const & name, Group * const parent )
  OutputManager( string const & name,
                 Group * const parent );

  /// Destructor
  virtual ~OutputManager() override;

  /// @copydoc geosx::dataRepository::Group::createChild( string const & childKey, string const & childName )
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    dataRepository::ViewKey time = { "time" };
  } viewKeys;
  /// @endcond
};


} /* namespace geosx */

#endif /* GEOSX_FILEIO_OUTPUTS_OUTPUTMANAGER_HPP_ */
