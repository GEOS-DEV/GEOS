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
   */
  DataContext( string const & targetName );

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
   * @return The result of toString() properly suffixed by the name of a contained object.
   */
  virtual string toString( string const & innerObjectName ) const
  { return toString() + '/' + innerObjectName; }

  /**
   * @return Get the target object name
   */
  string getTargetName() const
  { return m_targetName; }

  virtual string getTargetNameInPath( bool & fileLineFound ) const
  { fileLineFound = false; return m_targetName; }

  /**
   * @brief Insert contextual information in the provided stream.
   */
  friend std::ostream & operator<<( std::ostream & os, const DataContext & dt );

protected:

  /// @see getObjectName()
  string const m_targetName;

};

/// Stores information to retrieve where a target object has been declared in the input source
/// file (e.g. XML).
class DataFileContext final : public DataContext
{
public:

  /**
   * @brief Construct the file context of a Group from an xml node.
   * @param targetNode the target object xml node
   * @param nodePos the target object xml node position
   */
  DataFileContext( xmlWrapper::xmlNode const & targetNode, xmlWrapper::xmlNodePos const & nodePos );
  /**
   * @brief Construct the file context of a Group from an xml node.
   * @param targetNode the xml node containing the xml attribute
   * @param att the target object xml attribute
   * @param attPos the target object xml attribute position
   */
  DataFileContext( xmlWrapper::xmlNode const & targetNode, xmlWrapper::xmlAttribute const & att,
                   xmlWrapper::xmlAttributePos const & attPos );

  /**
   * @brief Destroy the DataFileContext object
   */
  virtual ~DataFileContext() {}

  /**
   * @return the target object name followed by the the file and line declaring it.
   */
  virtual string toString() const;
  /**
   * @copydoc DataContext::toString()
   */
  string toString( string const & innerObjectName ) const override
  { return toString() + ", attribute " + innerObjectName; }

  virtual string getTargetNameInPath( bool & foundNearestLine ) const override;

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
