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
 * @file InputExtensionGroup.hpp
 */

#ifndef GEOS_DATAREPOSITORY_INPUTEXTENSIONGROUP_HPP_
#define GEOS_DATAREPOSITORY_INPUTEXTENSIONGROUP_HPP_

#include "codingUtilities/Patterns.hpp"
#include "ExecutableGroup.hpp"
#include "StubGroup.hpp"
#include "KeyNames.hpp"
#include "InputExtension.hpp"
#include "InputParsing.hpp"
#include "xmlInputExtension.hpp"

namespace geos
{

namespace dataRepository
{

template < typename Base, typename DocNode = inputExtension::default_document_node >
class InputExtensionGroup : public Base
{
public:
  using super = Base;
  using document_node = DocNode;
  using callback_coordinator = typename inputExtension::InputExtender< document_node >::super;

  static_assert((std::is_base_of<Group, Base>::value), "Base must be a subclass of Group.");
  using super::super;

  using Group::getInputExtensionRules; // ensure visibility of the virtual function

  virtual std::vector< inputExtension::Rule< > > getInputExtensionRules( callback_coordinator & callbackCoordinator ) const override final
  {
    return inputExtensionRules( callbackCoordinator );
  }

  virtual std::vector< inputExtension::Rule<> > inputExtensionRules( callback_coordinator & ) const = 0;
};
namespace internal
{
template < typename Base, typename DocNode = inputExtension::default_document_node >
class InputRemovalGroup : public Base
{
public:
  using super = Base;
  using document_node = DocNode;
  using remove_callback_coordinator = typename inputExtension::internal::InputRemover< document_node >::super;

  static_assert((std::is_base_of<Group, Base>::value), "Base must be a subclass of Group.");
  using super::super;

  std::vector< inputExtension::Rule< > > getInputRemovalRules( remove_callback_coordinator & callbackCoordinator )
  {
    return inputRemovalRules( callbackCoordinator );
  }

  virtual std::vector< inputExtension::Rule<> > inputRemovalRules( remove_callback_coordinator & ) const = 0;
};

} // namespace internal

template < typename Base, typename... StubWrapperInfos >
class DeprecatedGroup : public StubGroup< InputExtensionGroup< Base >, StubWrapperInfos... >
{
public:
  using super = StubGroup< InputExtensionGroup< Base >, StubWrapperInfos... >;
  using typename super::document_node;
  using typename super::super::callback_coordinator;
  using parsing_context = inputExtension::StaticParsingContext< document_node,  inputExtension::SingleNodePolicy >;

  using super::super;

  enum Flags
  {
    RuleApplied,
    NodeAdded,
    AttributeAdded
  };

  // Explicit constructor which appends "[DEPRECATED]" to the input wrapper descriptions
  explicit DeprecatedGroup( std::string const & name, Group * const parent, StubWrapperInfos... ws )
    : super(name, parent, modifyWrapperInfo(ws)...)
  {}

  virtual std::vector< inputExtension::Rule< > > const & getInputExtensionRules( callback_coordinator & callbackCoordinator ) const override
  {
    std::vector< inputExtension::Rule< > > & extensionRules = this->inputExtensionRules();
    for (auto & rule : extensionRules)
    {
      rule.flags.insert( Flags::RuleApplied ); // Flag for rules
    }
    std::stringstream notification;
    callbackCoordinator.registerCallback( Flags::RuleApplied, [&notification]( inputExtension::Rule< > & rule, parsing_context const & parsingContext, document_node & extensionNode )
    {
      std::unordered_map<std::string, std::string> formatArgs;
      decltype(auto) originalNode = parsingContext.getNode( 0 );
      // external_node_pos_type & nodePos = parsingContext.getPosition(0);
      formatArgs["old_element"] = originalNode.tag();
      int elementCounter = 0;
      std::vector< string > newElements;
      StaticTreeIteration<>::processTree( originalNode, [&]( document_node & node )
      {
        formatArgs[GEOS_FMT("new_element_{}",elementCounter++)] = node.tag();
        newElements.push_back( node.tag( ) );
      } );
      formatArgs["new_elements"] = stringutilities::join( newElements.begin(), newElements.end(), ", " );
      notification << "DEPRECATION WARNING at node '" << originalNode.name()
                  // << "' from file=" << nodePos.filePath
                  << ", path=" << originalNode.path();
                  // << ", line=" << nodePos.line << ":\n";
      notification << GEOS_VFMT( "The '{old_element}' element(s) is/are deprecated and will be removed in a future version. Please transition to the '{new_elements}' element(s) instead.\nThis input configuration is being injected into the simulation as a replacement:\n", GEOS_FMT_ARG_MAP( formatArgs ) );
      extensionNode.print(notification, "  ", pugi::format_default, pugi::encoding_utf8);
      notification << "\n";
      for( auto & note : rule.devNotes )
      {
        notification << "NOTE: " << GEOS_FMT( note, GEOS_FMT_ARG_MAP( formatArgs ) );
      }
    });
    return extensionRules;
  }

  virtual void postInputInitialization() override
  {
    super::deregisterAllRecursive( );
  }

private:
  // Helper function to modify a single WrapperInfo
  template < typename WrapperInfo >
  static WrapperInfo modifyWrapperInfo(WrapperInfo const & info)
  {
    WrapperInfo modifiedInfo = info;
    modifiedInfo.description += " [DEPRECATED]";
    return modifiedInfo;
  }

  // Helper function to modify each WrapperInfo in the pack
  template < typename... Ws >
  static auto modifyWrapperInfo( Ws const &... ws )
  {
    return std::make_tuple( modifyWrapperInfo(ws)... );
  }
};

template < typename Base, typename... StubWrapperInfos >
class DeprecatedExecutableGroup : public DeprecatedGroup< Base, StubWrapperInfos... >
{
  static_assert( (std::is_base_of< ExecutableGroup, Base >::value), "The Base group must inherit from ExecutableGroup" );
public:
  using super = DeprecatedGroup< Base, StubWrapperInfos... >;
  using typename super::super::super::callback_coordinator;
  using typename super::super::super::document_node;
  using parsing_context = inputExtension::StaticParsingContext< document_node,  inputExtension::SingleNodePolicy >;


  enum Flags
  {
    Executable // Flag for executable nodes
  };

  virtual std::vector< inputExtension::Rule< > > const & getInputExtensionRules( callback_coordinator & callbackCoordinator ) const override
  {
    // Register a callback for handling executable nodes
    callbackCoordinator.registerCallback( Flags::Executable, [this]( inputExtension::Node< > & extensionNode, parsing_context const & parsingContext, document_node & node )
    {
      std::string dataRepoPath = xmlPathToDataRepoPath(node);
      m_pathsForExecution.push_back(dataRepoPath);
    });

    return super::getInputExtensionRules( callbackCoordinator );
  }

  virtual void postInputInitialization() override
  {
    super::postInputInitialization(); // Call base post-initialization
    initializeExecutableReplacements(); // Initialize executable replacements
  }

  virtual bool execute(real64 const time_n, real64 const dt, integer const cycleNumber, integer const eventCounter, real64 const eventProgress, DomainPartition& domain) override final
  {
    bool executed = false;
    for (auto & replacement : m_executableReplacements)
    {
      executed |= replacement->execute(time_n, dt, cycleNumber, eventCounter, eventProgress, domain);
    }
    return executed;
  }

  virtual void cleanup(real64 const time_n, integer const cycleNumber, integer const eventCounter, real64 const eventProgress, DomainPartition& domain) override final
  {
    for (auto & replacement : m_executableReplacements)
    {
      replacement->cleanup(time_n, cycleNumber, eventCounter, eventProgress, domain);
    }
  }

protected:
  void initializeExecutableReplacements()
  {
    for (const auto & path : m_pathsForExecution)
    {
      auto execGroup = this->template getGroupByPath<ExecutableGroup>(path);
      if ( execGroup )
      {
        m_executableReplacements.push_back(execGroup);
      }
    }
  }

private:
  mutable std::vector<std::string> m_pathsForExecution; // Paths to the executable groups created
  std::vector<std::shared_ptr<ExecutableGroup>> m_executableReplacements; // References to executable groups for direct invocation
};

template < typename Base, typename NodeType = inputExtension::default_document_node, typename... StubWrapperInfos >
class ObsoleteGroup : public StubGroup< Base, StubWrapperInfos... >
{
public:
  using super = StubGroup< Base, StubWrapperInfos... >;
  using callback_coordinator = typename inputExtension::InputExtender< NodeType >::super;

  ObsoleteGroup( string const & name, Group * const parent, StubWrapperInfos... ws ):
    super( name, parent, ws... )
  { }

  virtual void getInputExtensionRules( callback_coordinator & callbackCoordinator ) override final
  {
    std::ostringstream notification;
    // Generate error notification for the presence of an obsolete group in the input XML.
    notification << "ERROR: Obsolete Group Encountered!\n"
                  << "The group '" << this->getCatalogName() << "' is obsolete and should not be present in the input file.\n"
                  // << "Found in file: " << nodePos.filePath << "\n"
                  // << "Line: " << nodePos.line << "\n"
                  << "Please remove or replace this group according to the current schema guidelines.";
    GEOS_THROW( notification.str(), InputError );
  }
};

template < typename Document >
void reconstructIncludesAndSave( Document & document, const std::string & outputDirectory )
{
  using doc_node = typename Document::node_type;

  std::map< std::string, std::vector< doc_node > > documentMap;

  // Traverse the document and distribute nodes into separate documentPaths
  doc_node root = document.getFirstChild();
  StaticTreeIteration< >::processTree( root, [&]( doc_node & node )
  {
    std::string filepath = node.attribute( inputParsing::filePathString ).value();
    documentMap[ filepath ].push_back( node );
  } );

  // Function to strip attributes
  auto stripAttributes = [] ( doc_node & node )
  {
    node.remove_attribute( inputParsing::filePathString );
    node.remove_attribute( inputParsing::charOffsetString );
  };

  // Helper function to create or get the corresponding node in the new document
  auto getOrCreateNode = []( doc_node & parent, const std::string & name ) -> doc_node
  {
    for ( doc_node & child : parent.children() )
    {
      if ( child.name() == name )
      {
        return child;
      }
    }
    return parent.append_child( name.c_str() );
  };

  // Now reconstruct the files and reintroduce Include elements
  for ( const auto & entry : documentMap )
  {
    const std::string & filepath = entry.first;
    const std::vector< doc_node > & nodes = entry.second;

    // Create a new document for each original file
    Document doc;
    doc_node newRoot = doc.appendChild( geos::dataRepository::keys::ProblemManager );

    for (const doc_node & node : nodes)
    {
      // potentially recreate the full tree structure for each node in the new document
      std::vector< std::string > nodePath;
      for (doc_node n = node; n != root; n = n.parent())
      {
        nodePath.push_back(n.name());
      }
      std::reverse(nodePath.begin(), nodePath.end());

      doc_node current = newRoot;
      for (const std::string & name : nodePath)
      {
        current = getOrCreateNode(current, name);
      }

      doc_node newChild = current.append_copy(node);
      stripAttributes(newChild);
    }

    // TODO (wrt) : add -preprocessed suffix

    // Save the file
    std::string fullpath = joinPath( outputDirectory, filepath );
    doc.saveFile(fullpath.c_str());
  }
}

}

}

#endif