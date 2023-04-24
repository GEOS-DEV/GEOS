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
 * @file DataContext.hpp
 */

#ifndef GEOS_DATAREPOSITORY_SOURCECONTEXT_HPP_
#define GEOS_DATAREPOSITORY_SOURCECONTEXT_HPP_

#include "common/DataTypes.hpp"
#include "common/Logger.hpp"
#include "xmlWrapper.hpp"

namespace geos
{

namespace dataRepository
{

class Group;
class WrapperBase;

/// This abstract class stores data that helps to retrieve an object:
/// - from which position in a file (if applicable), see DataFileContext,
/// - where it is located in the data hierarchy, see GroupContext and WrapperContext.
/// Typically, the target object contain an unique_ptr< DataContext > instance of this class.
class DataContext
{
public:

  /**
   * @brief Construct a new DataContext object.
   * @param objectName the target object name
   * @param isDataFileContext true if this Context is a DataFileContext (see isDataFileContext for more infos)
   */
  DataContext( string const & objectName, bool isDataFileContext );

  /**
   * @return A string that mention all the known informations to retrieve from where the target
   * object comes from.
   */
  virtual string toString() const = 0;

  string getObjectName() const
  { return m_objectName; }

  /**
   * @brief In some cases, we need to know if a DataContext is from a file. It means that it
   * is a more user-friendly information compared to other DataContext classes.
   * @return true if the context is from a file.
   */
  bool isDataFileContext() const
  { return m_isDataFileContext; }

  /**
   * @brief Insert toString() result in a stream.
   */
  friend std::ostream & operator<<( std::ostream & os, const DataContext & dt );

protected:

  /// see getObjectName()
  string const m_objectName;

  /// see isDataFileContext()
  bool const m_isDataFileContext;

};

/// Helps to know where a Group is in the hierarchy.
/// See DataContext class for more infos.
class GroupContext : public DataContext
{
public:

  /**
   * @brief Construct a new GroupContext object
   * @param group see getGroup()
   */
  GroupContext( Group & group );

  /**
   * @brief Get the target Group object (wrapper's parent in case of a WrapperContext).
   */
  Group & getGroup() const;

  /**
   * @copydoc DataContext::toString()
   */
  virtual string toString() const;

protected:

  /**
   * @brief Construct a new GroupContext object
   * @param group see getGroup()
   * @param objectName Target object name.
   */
  GroupContext( Group & group, string const & objectName );

  /// see getGroup()
  Group & m_group;

};

/// Helps to know the source context of a Wrapper in the hierarchy, or in the source file, if possible.
/// See DataContext class for more infos.
class WrapperContext final : public GroupContext
{
public:

  /**
   * @brief Construct a new WrapperContext object
   */
  WrapperContext( WrapperBase & wrapper );

  /**
   * @copydoc DataContext::toString()
   */
  virtual string toString() const;

};

/// Helps to know from where a Group or a Wrapper has been declared in its source file (a xml typically).
/// See DataContext class for more infos.
class DataFileContext final : public DataContext
{
public:

  /**
   * @brief Construct the file context of a Group from an xml node.
   */
  DataFileContext( Group & group, xmlWrapper::xmlNodePos const & nodePos, string const & nodeTagName );
  /**
   * @brief Construct the file context of a Group from an xml attribute.
   */
  DataFileContext( WrapperBase & wrapper, xmlWrapper::xmlAttributePos const & attPos );

  virtual string toString() const;

  /**
   * @brief Get the type name in the source file (XML node tag name / attribute name).
   */
  string getTypeName() const
  { return m_typeName; }

  /**
   * @brief Get the source file path where the target object has been declared.
   */
  string getFilePath() const
  { return m_filePath; }

  /**
   * @brief Get the line (starting from 1) where the target object has been declared in the source
   * file.
   */
  size_t getLine() const
  { return m_line; }

  /**
   * @brief Get the character offset in the line (starting from 1) where the target object has been declared in the source
   * file.
   */
  size_t getOffsetInLine() const
  { return m_offsetInLine; }

  /**
   * @brief Get the character offset (starting from 0) in the source file path where the target object has been
   * declared.
   */
  size_t getOffset() const
  { return m_offset; }

protected:

  /// see getTypeName()
  string const m_typeName;
  /// see getFilePath()
  string const m_filePath;
  /// see getLine()
  size_t const m_line;
  /// see getLineOffset()
  size_t const m_offsetInLine;
  /// see getOffset()
  size_t const m_offset;

};

} /* namespace dataRepository */

} /* namespace geos */

#endif /* GEOS_DATAREPOSITORY_INPUTFLAGS_HPP_ */
