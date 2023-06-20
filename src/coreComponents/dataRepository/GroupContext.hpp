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
 * @file GroupContext.hpp
 */

#ifndef GEOS_DATAREPOSITORY_GROUPCONTEXT_HPP_
#define GEOS_DATAREPOSITORY_GROUPCONTEXT_HPP_

#include "DataContext.hpp"
#include "WrapperBase.hpp"
#include "Group.hpp"

namespace geos
{
namespace dataRepository
{


/// Helps to know where a Group is in the hierarchy.
/// See DataContext class for more info.
class GroupContext : public DataContext
{
public:

  /**
   * @brief Construct a new GroupContext object
   * @param group The reference to the Group related to this GroupContext.
   */
  GroupContext( Group & group );

  /**
   * @brief Destroy the GroupContext object
   */
  virtual ~GroupContext() {}

  /**
   * @return the reference to the Group related to this GroupContext.
   */
  Group & getGroup() const;

  /**
   * @return the group path with the file & line of the first parent for which this information exists.
   */
  virtual string toString() const;

protected:

  /**
   * @brief Construct a new GroupContext object
   * @param group The reference to the Group related to this GroupContext.
   * @param objectName Target object name.
   */
  GroupContext( Group & group, string const & objectName );

  /// The reference to the Group related to this GroupContext.
  Group & m_group;

};

/// Dedicated implementation of GroupContext for Wrapper.
/// See DataContext class for more info.
class WrapperContext final : public GroupContext
{
public:

  /**
   * @brief Construct a new WrapperContext object
   * @param wrapper the target Wrapper object
   */
  WrapperContext( WrapperBase & wrapper );

  /**
   * @brief Destroy the WrapperContext object
   */
  virtual ~WrapperContext() {}

  /**
   * @return the parent group DataContext followed by the wrapper name.
   */
  virtual string toString() const;

};


} /* namespace dataRepository */
} /* namespace geos */

#endif /* GEOS_DATAREPOSITORY_GROUPCONTEXT_HPP_ */
