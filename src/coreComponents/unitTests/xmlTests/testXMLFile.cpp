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

// Source includes
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/FieldSpecificationBase.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/SourceContext.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>
#include <algorithm>

using namespace geos;
using namespace geos::dataRepository;
using namespace xmlWrapper;

// Tests if the xml file parsing works with one file, one file with nested includes, and multiple files
TEST( testXML, testXMLFile )
{
  geos::ProblemManager & problemManager = geos::getGlobalState().getProblemManager();
  problemManager.parseCommandLineInput();
  problemManager.parseInputFile();

  // Check that we've read the full XML with all nested includes by inspecting boundary conditions
  EXPECT_TRUE( problemManager.getFieldSpecificationManager().hasGroup< FieldSpecificationBase >( "v0" ) );
}

/**
 * @brief Returns the xml nodes and attribute in hierarchy that GEOS is supposed to read in
 * the Group hierarchy.
 * @param document root xml document (which potentially contains includes).
 * @param targetNode the function will be called on this node and its children.
 * @param elementNames the list of xml elements names with the form "GroupName" or
 * "GroupName/WrapperName". The node name is supposed to be the name attribute value, it it exists,
 * or the tag name.
 */
void getElementsRecursive( xmlDocument const & document, xmlNode const & targetNode,
                           std::set< string > & elementNames )
{
  // Store here every node name that only exist in the xml (and not in the Group hierarchy).
  static const std::set< string > xmlOnlyNodes {
    "ElementRegions"
  };

  // The Group name will be the name attribute value, or the node tag name if the name attribute
  // doesn't exist.
  string groupName = [&]() {
    xmlAttribute nameAtt = targetNode.attribute( "name" );
    return nameAtt ? string( nameAtt.value() ) : string( targetNode.name() );
  }();

  if( xmlOnlyNodes.find( groupName ) == xmlOnlyNodes.end() )
  {
    elementNames.emplace( groupName );

    for( xmlAttribute att : targetNode.attributes() )
    {
      if( !isFileMetadataAttribute( att.name() ) )
      {
        elementNames.emplace( groupName + '/' + att.name() );
      }
    }
  }

  for( xmlNode subNode : targetNode.children() )
  {
    getElementsRecursive( document, subNode, elementNames );
  }
}
/**
 * @brief Verify if the FileContext data correctly locates from where the object comes from (in the
 * correct document/include, at the correct line and at the correct character offset).
 * @param document root xml document (which potentially contains includes).
 * @param fileContext the FileContext to verify.
 */
void verifyFileContext( FileContext const & fileContext,
                        xmlDocument const & document )
{
  string const & strToVerify = fileContext.getTypeName();
  string const & errInfos = "Verifying " + strToVerify + " in " + fileContext.toString();

  // verifying if all FileContext data have been found
  EXPECT_FALSE( fileContext.getFilePath().empty() ) << errInfos;
  EXPECT_FALSE( fileContext.getObjectName().empty() ) << errInfos;
  EXPECT_FALSE( fileContext.getTypeName().empty() ) << errInfos;
  EXPECT_NE( fileContext.getOffset(), xmlDocument::npos ) << errInfos;
  EXPECT_NE( fileContext.getLine(), xmlDocument::npos ) << errInfos;
  EXPECT_NE( fileContext.getOffsetInLine(), xmlDocument::npos ) << errInfos;

  // will crash if the source buffer specified by fileContext.getFilePath() isn't available
  string const * buffer = document.getOriginalBuffer( fileContext.getFilePath() );
  EXPECT_NE( buffer, nullptr ) << errInfos;

  size_t curLine = 1;
  bool lineFound = false;

  // Does fileContext.getOffset() locates the object?
  EXPECT_LT( fileContext.getOffset() + strToVerify.size(), buffer->size() ) << errInfos;
  EXPECT_EQ( strToVerify, buffer->substr( fileContext.getOffset(), strToVerify.size() ) ) << errInfos;

  // Were trying to reach the line return by FileContext::getLine()
  for( size_t offset=0;
       offset < buffer->size() && curLine < fileContext.getLine();
       ++offset )
  {
    if( (*buffer)[offset] == '\n' )
    {
      curLine++;
    }

    // the theorical line has been reach!
    if( curLine == fileContext.getLine() )
    {
      // Does fileContext.getLine() and fileContext.getOffsetInLine() locates the object?
      EXPECT_LT( offset + fileContext.getOffsetInLine() + strToVerify.size(),
                 buffer->size() ) << errInfos;
      EXPECT_EQ( strToVerify,
                 buffer->substr( offset + fileContext.getOffsetInLine(), strToVerify.size() ) ) << errInfos;
      lineFound = true;
    }
  }
  // does the fileContext line has been reached?
  EXPECT_TRUE( lineFound );
}
/**
 * @brief Verifies if the specified group, its children, and the associated wrappers have a
 * FileContext object that correctly locate from where the objects where declared in the
 * source file.
 * @param document root xml document (which potentially contains includes).
 * @param group The group to (recursively) verify.
 * @param verifCount a set that will be filled with the names, with the form
 * "GroupName" or "GroupName/WrapperName".
 */
void verifyGroupFileContextRecursive( xmlDocument const & document, Group const & group,
                                      std::set< string > & verifications )
{
  // GEOS_LOG( "Verifying "<< group.getName());
  if( group.getSourceContext().isFileContext() )
  {
    verifyFileContext( dynamic_cast< FileContext const & >( group.getSourceContext() ),
                       document );
    verifications.emplace( group.getName() );
  }

  for( auto const & wrapperIterator : group.wrappers() )
  {
    WrapperBase const * wrapper = wrapperIterator.second;
    if( wrapper->getSourceContext().isFileContext() )
    {
      verifyFileContext( dynamic_cast< FileContext const & >( wrapper->getSourceContext() ),
                         document );
      verifications.emplace( group.getName() + '/' + wrapper->getName() );
    }
  }

  for( auto subGroup : group.getSubGroups() )
  {
    verifyGroupFileContextRecursive( document, *subGroup.second, verifications );
  }
}

/**
 * @brief Returns the element that exists in setB but not in setA.
 */
std::set< string > getDifference( std::set< string > & setA,
                                  std::set< string > & setB )
{
  std::vector< string > result( setA.size() + setB.size() );

  auto it = std::set_difference( setA.begin(), setA.end(),
                                 setB.begin(), setB.end(),
                                 result.begin() );
  result.resize( it - result.begin() );

  return std::set< string >( result.begin(), result.end() );
}

// Tests
// - if the line information of each nodes and attributes can be retrieved,
// - if the resulting Group & Wrapper hierarchy matches with the input xml documents and includes hierarchy.
TEST( testXML, testXMLFileLines )
{
  xmlDocument xmlDocument;
  ProblemManager & problemManager = getGlobalState().getProblemManager();

  {
    problemManager.parseCommandLineInput();
    Group & commandLine = problemManager.getGroup( problemManager.groupKeys.commandLine );
    string const & inputFileName = commandLine.getReference< string >( problemManager.viewKeys.inputFileName );
    xmlDocument.load_file( inputFileName.c_str(), true );
    problemManager.parseXMLDocument( xmlDocument );
  }

  GEOS_LOG( "Loaded files : " );
  for( auto const & buffer: xmlDocument.getOriginalBuffers() )
  {
    GEOS_LOG( "    " << buffer.first << " (" << buffer.second.size() << " chars)" );
  }

  std::set< string > expectedElements;
  getElementsRecursive( xmlDocument, xmlDocument.root().child( "Problem" ), expectedElements );

  std::set< string > verifiedElements;
  verifyGroupFileContextRecursive( xmlDocument, problemManager, verifiedElements );

  std::set< string > notFound = getDifference( expectedElements, verifiedElements );
  EXPECT_TRUE( notFound.empty() ) << "Infos : There should not exists xml element that were not in "
                                     "the Group hierarchy.\nNot in Group hierarchy : {"
                                  << stringutilities::join( notFound, "," ) << "}";

  std::set< string > notExpected = getDifference( verifiedElements, expectedElements );
  EXPECT_TRUE( notExpected.empty() ) << "Infos : There should not exists an object in the Group "
                                        "hierarchy that contains a FileContext but which were not "
                                        "declared in the Xml.\nNot in XML hierarchy : {"
                                     << stringutilities::join( notExpected, "," ) << "}";
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  GeosxState state( basicSetup( argc, argv, true ) );

  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
