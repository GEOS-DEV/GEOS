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
 * @file InputExtensions.hpp
 */

#ifndef GEOS_DATAREPOSITORY_INPUTEXTENSIONS_HPP_
#define GEOS_DATAREPOSITORY_INPUTEXTENSIONS_HPP_

#include "common/DataTypes.hpp"
#include "codingUtilities/Patterns.hpp"
#include "xmlWrapper.hpp"
#include <functional>
#include <optional>

namespace geos
{

namespace dataRepository
{

namespace inputExtension
{

using default_document = xmlWrapper::xmlDocument;
using default_document_node = typename default_document::node_type;

// Define a single monolithic trait to check for all required 'using' statements
template <typename, typename = void >
struct is_valid_document : std::false_type {};

template <typename T>
struct is_valid_document<T, std::void_t<
  typename T::read_result_type,
  typename T::node_type,
  typename T::node_pos_type,
  typename T::attribute_type,
  typename T::attribute_pos_type
>> : std::true_type {};


// Policy for handling single document node.
template< typename DocNode = default_document_node >
struct SingleNodePolicy
{
  using container_type = std::reference_wrapper< DocNode >;
  SingleNodePolicy( std::initializer_list< std::reference_wrapper< DocNode > > nodes )
    : m_container(*nodes.begin())
  {
    if (nodes.size() != 1)
    {
      GEOS_ERROR("SingleNodePolicy expects exactly one node.");
    }
  }
  DocNode & getNode( int GEOS_UNUSED_PARAM( index ) )
  {
    return m_container;
  }
  size_t size( ) const { return 1; }
  container_type m_container;
};

// Policy for handling multiple document nodes.
template< typename DocNode = default_document_node >
struct NodeSetPolicy
{
  using container_type = std::vector< std::reference_wrapper< DocNode > >;
  NodeSetPolicy( std::initializer_list< std::reference_wrapper< DocNode > > nodes )
    : m_container( nodes.begin(), nodes.end() )
  { }
  DocNode & getNode( int index )
  {
    return m_container.at(index);
  }
  size_t size() const { return m_container.size(); }
  container_type m_container;
};

// Refactored StaticParsingContext with context policy.
template <
  typename DocNode = default_document_node,
  template <typename> class ContextPolicy = SingleNodePolicy
>
class StaticParsingContext
{
public:
  using policy = ContextPolicy< DocNode >;

  StaticParsingContext( std::initializer_list< std::reference_wrapper< DocNode > > nodes )
    : m_context(nodes)
  {}

  size_t scope() const { return m_context.size( ); }

  DocNode & getNode( int index ) { return m_context.getNode( index ); }

  template <typename LAMBDA>
  void forInContext(LAMBDA&& lambda)
  {
    for (int ii = 0; ii < scope(); ++ii)
    {
      lambda( getNode(ii) );
    }
  }

private:
  policy m_context;
};

namespace singletonParsingContext
{
  template <typename DocNode = default_document_node, template <typename> class ContextPolicy = SingleNodePolicy>
  auto isValid()
  {
    return [] ( StaticParsingContext< DocNode, ContextPolicy > & context ) { return context.scope() == 1; };
  };

  template <typename DocNode = default_document_node, template <typename> class ContextPolicy = SingleNodePolicy>
  auto isNamed(const std::string & name)
  {
    return [name] ( StaticParsingContext< DocNode, ContextPolicy > & context )
    {
      return (context.scope() == 1) && context.getNode(0).name() == name;
    };
  };

  template <typename DocNode = default_document_node, template <typename> class ContextPolicy = SingleNodePolicy, typename... Names>
  auto hasAttributes( Names... names )
  {
    return [names...] ( StaticParsingContext< DocNode, ContextPolicy > & context )
    {
      return (context.scope() == 1) && (context.getNode(0).has_attribute(names) && ...);
    };
  };

  template < typename DocNode = default_document_node, template <typename> class ContextPolicy = SingleNodePolicy >
  auto root()
  {
    return [] ( StaticParsingContext< DocNode, ContextPolicy > & context ) -> decltype(auto)
    {
      if (context.scope() != 1)
      {
        GEOS_ERROR("This root resolver relies on a unitary parsing context!");
      }
      return context.getNode(0).root();
    };
  };

  template < typename DocNode = default_document_node, template <typename> class ContextPolicy = SingleNodePolicy >
  auto thisNode()
  {
    return [] ( StaticParsingContext< DocNode, ContextPolicy > & context ) -> decltype(auto)
    {
      if (context.scope() != 1)
      {
        GEOS_ERROR("This resolver relies on a unitary parsing context!");
      }
      return context.getNode(0);
    };
  };

  template < typename DocNode = default_document_node, template <typename> class ContextPolicy = SingleNodePolicy >
  auto thisNodesParent()
  {
    return [] ( StaticParsingContext< DocNode, ContextPolicy > & context ) -> decltype(auto)
    {
      if (context.scope() != 1)
      {
        GEOS_ERROR("This parent resolver relies on a unitary parsing context!");
      }
      return context.getNode(0).parent();
    };
  };

  template < typename DocNode = default_document_node, template <typename> class ContextPolicy = SingleNodePolicy >
  auto thisNodesChild( std::string const & name )
  {
    return [name] ( StaticParsingContext< DocNode, ContextPolicy > & context ) -> decltype(auto)
    {
      if (context.scope() != 1)
      {
        GEOS_ERROR("This parent resolver relies on a unitary parsing context!");
      }
      for( auto & child : context.getNode(0).children() )
      {
        if ( child.name() == name )
        {
          return child;
        }
      }
      GEOS_ERROR( GEOS_FMT( "Child node '{}' for extension not found!", name ) );
      return context.getNode(0); // never reached, prevents compilation failure
    };
  };
}

template < typename DocNode = default_document_node, typename FlagType = int >
struct Attribute
{
  std::optional< string > source{};
  string destination{};
  std::function< string( const string & ) > transform = []( const string & ){ return string{}; };
  std::function< DocNode ( StaticParsingContext< DocNode > & ) > sourceNode = singletonParsingContext::thisNode< DocNode >();
  std::set< FlagType > flags{};
};

template < typename DocNode = default_document_node, typename FlagType = int >
struct Node
{
  string identifier = {};
  std::vector< Attribute< DocNode, FlagType > > attributes{};
  std::vector< Node< DocNode > > subNodes{};
  std::function< std::vector< Attribute< DocNode, FlagType > >( StaticParsingContext< DocNode > & ) > dynamicAttributes = [](const auto&){ return std::vector< Attribute< DocNode, FlagType > >(); };
  std::function< std::vector< Node< DocNode, FlagType > >( StaticParsingContext< DocNode > & ) > dynamicSubNodes = [](const auto&){ return std::vector< Node< DocNode, FlagType > >(); };
  std::set< FlagType > flags{};
};

template < typename DocNode = default_document_node, typename FlagType = int >
struct Rule
{
  std::vector< Node< DocNode, FlagType > > subTrees {};
  std::function< bool( StaticParsingContext< DocNode > & ) > appliesWhen = singletonParsingContext::isValid< DocNode >();
  std::function< DocNode ( StaticParsingContext< DocNode > & ) > determineRoot = singletonParsingContext::thisNodesParent< DocNode >();
  std::vector< string > devNotes{};
  std::set< FlagType > flags{};
};

namespace mapAttribute
{
  template < typename DocNode = default_document_node, typename FlagType = int >
  Attribute< DocNode, FlagType > newAttribute( string const & newName, string const & newValue )
  {
    Attribute< DocNode, FlagType > attribute;
    attribute.destination = newName;
    attribute.transform = [newValue](string const & GEOS_UNUSED_PARAM( oldValue ) ) { return newValue; };
    return attribute;
  }

  // Helper function for identity transformation between source and destination names
  template < typename DocNode = default_document_node, typename FlagType = int >
  Attribute< DocNode, FlagType > identity(const string & oldName, const string & newName)
  {
    return { oldName, newName, [](const string & value) { return value; } };
  }

  template < typename DocNode = default_document_node, typename FlagType = int >
  Attribute< DocNode, FlagType > append( const string & oldName, const string & newName, const string & suffix )
  {
    return { oldName, newName, [suffix](const string & oldValue) { return oldValue + suffix; } };
  }

  // Helper function for appending node ID to the attribute value, with different source and destination names
  template < typename T,  typename DocNode = default_document_node, typename FlagType = int >
  Attribute< DocNode, FlagType > appendCatalogName(const string & oldName, const string & newName)
  {
    return append( oldName, newName, T::getCatalogName() );
  }

  template < typename DocNode = default_document_node, typename FlagType = int >
  Attribute< DocNode, FlagType > atomAsList( const string & oldName, const string & newName )
  {
    return { oldName, newName, [](const string & value) { return "{ " + value + " }"; } };
  }
}

// template < typename RuleFor, typename LAMBDA, typename DocNode = default_document_node >
// std::vector< RuleFor > forEachInAttributeList( const string & listAttribute, LAMBDA && lambda )
// {
//   using static_parsing_context = StaticParsingContext< DocNode >;

//   return [&]( static_parsing_context const & context )
//   {
//     GEOS_ERROR_IF( ! singletonParsingContext::isValid< DocNode >()( context ), "This dynamic attribute generator generator is intended for use in unitary parsing contexts!" );
//     DocNode contextNode = context.getNode( 0 );
//     std::vector< RuleFor > rulesFor;
//     auto attribute = contextNode.attribute( listAttribute.c_str() );
//     if ( ! attribute )
//     {
//       return rulesFor;
//     }
//     // TODO (wrt) : use our actual parsing mechanisms for arrays to deal with e.g. nested lists
//     std::vector< string > items{};
//     for ( const auto & item : items )
//     {
//       stringutilities::trimSpaces( item );
//       rulesFor.push_back( lambda( item ) );
//     }
//     return rulesFor;
//   };
// }

template < typename DocNode = default_document_node, typename FlagType = int >
class InputExtender : public FlagCallbackCoordinator< FlagType,
                                                      std::tuple< Rule< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >,
                                                      std::tuple< Node< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >,
                                                      std::tuple< Attribute< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & > >
{
public:
  using super = FlagCallbackCoordinator< FlagType,
                                         std::tuple< Rule< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >,
                                         std::tuple< Node< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >,
                                         std::tuple< Attribute< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & > >;

  InputExtender( std::set< std::string > const & mergable )
    : super( )
    , m_mergable( mergable )
  { }

  void applyExtensions( StaticParsingContext< DocNode > & staticParsingContext, std::vector< Rule< DocNode, FlagType > > & extensions )
  {
    for ( auto & extension : extensions )
    {
      if ( extension.appliesWhen( staticParsingContext ) )
      {
        decltype(auto) extensionRoot = extension.determineRoot( staticParsingContext );
        createNodes( extensionRoot, staticParsingContext, extension.subTrees );
        this->template executeCallbacks< Rule< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >( extension.flags, extension, staticParsingContext, extensionRoot );
      }
    }
  }

private:
  void createNodes( DocNode & parent, StaticParsingContext< DocNode > & staticParsingContext, std::vector< Node< DocNode > > & subnodes )
  {
    for ( auto & node : subnodes )
    {
      // Lambda to handle operations on newNode
      auto handleNewNode = [&]( auto & newNode )
      {
        // Determine static and dynamic attributes to add to the node
        auto dynamicAttributes = node.dynamicAttributes( staticParsingContext );
        std::vector< Attribute< DocNode, FlagType > > allAttributes = node.attributes;
        allAttributes.insert( allAttributes.end(), dynamicAttributes.begin(), dynamicAttributes.end() );
        createAttributes( newNode, staticParsingContext, allAttributes );

        // Determine static and dynamic subnodes to instantiate under the node
        auto dynamicNodes = node.dynamicSubNodes( staticParsingContext );
        std::vector< Node< DocNode > > allSubNodes = node.subNodes;
        allSubNodes.insert( allSubNodes.end(), dynamicNodes.begin(), dynamicNodes.end() );

        // recur
        createNodes( newNode, staticParsingContext, allSubNodes );
        this->template executeCallbacks< Node< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >( node.flags, node, staticParsingContext, newNode );
      };

      bool create = m_mergable.find(node.identifier) == m_mergable.end();

      // Conditionally create newNode and handle it
      if ( create )
      {
        // need decltype auto to allow return-by-reference
        decltype(auto) newNode = parent.append_child( node.identifier.c_str() );
        handleNewNode( newNode );
      }
      else
      {
        handleNewNode( parent );
      }
    }
  }

  void createAttributes( DocNode & extensionNode, StaticParsingContext< DocNode > & staticParsingContext, std::vector< Attribute< DocNode, FlagType > > & attributes )
  {
    DocNode & newNode = extensionNode;
    for ( auto & attribute : attributes )
    {
      decltype(auto) newAttr = newNode.append_attribute( attribute.destination.c_str() );
      string sourceValue = "";
      if ( attribute.source.has_value() )
      {
        // Otherwise, fetch the source attribute value from the node and apply the transformation
        sourceValue = attribute.sourceNode( staticParsingContext ).attribute( attribute.source.value().c_str() ).value();
      }
      newAttr.set_value( attribute.transform( sourceValue ).c_str() );
      this->template executeCallbacks< Attribute< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >( attribute.flags, attribute, staticParsingContext, newNode );
    }
  }

  std::set< string > const m_mergable;
};

// namespace internal
// {

// template < typename DocNode = default_document_node, typename FlagType = int >
// class InputRemover : public FlagCallbackCoordinator< FlagType,
//                                                      std::tuple< Rule< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >,
//                                                      std::tuple< Node< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >,
//                                                      std::tuple< Attribute< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & > >
// {
// public:
//   using super = FlagCallbackCoordinator< FlagType,
//                                          std::tuple< Rule< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >,
//                                          std::tuple< Node< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & >,
//                                          std::tuple< Attribute< DocNode, FlagType > &, StaticParsingContext< DocNode > &, DocNode & > >;
//   void applyRemovals(StaticParsingContext<DocNode> & staticParsingContext, const std::vector< Rule< DocNode, FlagType > > & removalRules )
//   {
//     for (const auto & rule : removalRules)
//     {
//       if ( rule.appliesWhen( staticParsingContext ) )
//       {
//         DocNode extensionRoot = rule.determineRoot( staticParsingContext );
//         removeNodes( extensionRoot, staticParsingContext, rule.subTrees );
//         this->executeCallbacks( rule.flags, rule, staticParsingContext, extensionRoot );
//       }
//     }
//   }

// protected:
//   void removeNodes(DocNode & extensionNode, StaticParsingContext< DocNode > & staticParsingContext, const std::vector< Node< DocNode > > & nodes)
//   {
//     DocNode & parent = extensionNode;
//     for (const auto & node : nodes)
//     {
//       removeNodes( node, staticParsingContext, node.subNodes );
//       DocNode childNode = parent.find_child(node.identifier.c_str());
//       if ( childNode )
//       {
//         parent.remove_child(childNode);
//       }
//     }
//   }
// };

// } // namespace internal

} // namespace inputExtension

} // namespace dataRepository

} // namespace geos

#endif