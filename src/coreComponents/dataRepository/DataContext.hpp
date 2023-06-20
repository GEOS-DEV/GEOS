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

#ifndef GEOS_DATAREPOSITORY_DATACONTEXT_HPP_
#define GEOS_DATAREPOSITORY_DATACONTEXT_HPP_

#include "common/DataTypes.hpp"
#include "common/Logger.hpp"
#include "xmlWrapper.hpp"
#include "common/Format.hpp"
#include "Group.hpp"
#include "WrapperBase.hpp"

namespace geos
{

namespace dataRepository
{

/// DataContext is an abstract class storing contextual information on an object:
/// - Either its line position in a file (if applicable, implementation in DataFileContext),
/// - or its location in the data hierarchy (implementation in GroupContext and WrapperContext).
/// Typically, the target object contains an unique_ptr< DataContext > instance of this class.
class DataContext
{
public:

  /**
   * @brief Construct a new DataContext object.
   * @param targetName the target object name
   * @param isDataFileContext true if this Context is a DataFileContext
   */
  DataContext( string const & targetName, bool isDataFileContext );

  /**
   * @brief Destroy the DataContext object
   */
  virtual ~DataContext() {}

  /**
   * @return A string that mention all the known informations to retrieve from where the target
   * object comes from.
   */
  virtual string toString() const = 0;

  /**
   * @return Get the target object name
   */
  string getTargetName() const
  { return m_targetName; }

  /**
   * @brief Flag on availability of file information. Used to provide more user-friendly information.
   * @return true if the context is from a file.
   */
  bool isDataFileContext() const
  { return m_isDataFileContext; }

  /**
   * @brief Insert contextual information in the provided stream.
   */
  friend std::ostream & operator<<( std::ostream & os, const DataContext & dt );

protected:

  /// @see getObjectName()
  string const m_targetName;

  /// @see isDataFileContext()
  bool const m_isDataFileContext;

};

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

/// Stores information to retrieve where a Group or Wrapper has been declared in the input source file (e.g. XML)
class DataFileContext final : public DataContext
{
public:

  /**
   * @brief Construct the file context of a Group from an xml node.
   * @param group the target Group object
   * @param nodePos the target object xml node position
   * @param nodeTagName the xml node tag name.
   */
  DataFileContext( Group & group, xmlWrapper::xmlNodePos const & nodePos, string const & nodeTagName );
  /**
   * @brief Construct the file context of a Group from an xml attribute.
   * @param wrapper the target Wrapper object
   * @param attPos the target object xml attribute position
   */
  DataFileContext( WrapperBase & wrapper, xmlWrapper::xmlAttributePos const & attPos );

  /**
   * @brief Destroy the DataFileContext object
   */
  virtual ~DataFileContext() {}

  /**
   * @return the target object name followed by the the file and line declaring it.
   */
  virtual string toString() const;

  /**
   * @return the type name in the source file (XML node tag name / attribute name).
   */
  string getTypeName() const
  { return m_typeName; }

  /**
   * @return the source file path where the target object has been declared.
   */
  string getFilePath() const
  { return m_filePath; }

  /**
   * @return the line (starting from 1) where the target object has been declared in the source file.
   */
  size_t getLine() const
  { return m_line; }

  /**
   * @return the character offset in the line (starting from 1) where the target object has been
   * declared in the source file.
   */
  size_t getOffsetInLine() const
  { return m_offsetInLine; }

  /**
   * @return the character offset (starting from 0) in the source file path where the target object has been
   * declared.
   */
  size_t getOffset() const
  { return m_offset; }

protected:

  /// @see getTypeName()
  string const m_typeName;
  /// @see getFilePath()
  string const m_filePath;
  /// @see getLine()
  size_t const m_line;
  /// @see getLineOffset()
  size_t const m_offsetInLine;
  /// @see getOffset()
  size_t const m_offset;

};


} /* namespace dataRepository */

} /* namespace geos */



/// Formatter to be able to directly use a DataContext as a GEOS_FMT() argument.
template<>
struct GEOS_FMT_NS::formatter< geos::dataRepository::DataContext >
{
  /**
   * @brief Format the specified DataContext to a string.
   * @param dataContext the DataContext object to format
   * @param ctx formatting state consisting of the formatting arguments and the output iterator
   * @return iterator to the output buffer
   */
  auto format( geos::dataRepository::DataContext const & dataContext, format_context & ctx )
  {
    return format_to( ctx.out(), dataContext.toString() );
  }

  /**
   * @brief Method to parse a dataContext from a string. Not implemented!
   * @param ctx formatting state consisting of the formatting arguments and the output iterator
   * @return iterator to the output buffer (leaved unchanged)
   */
  constexpr auto parse( format_parse_context & ctx )
  { return ctx.begin(); }
};

#endif /* GEOS_DATAREPOSITORY_DATACONTEXT_HPP_ */
