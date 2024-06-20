#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../InputProcessing.hpp"

using namespace geos::dataRepository;
using namespace geos::dataRepository::inputProcessing;
using ::testing::Return;

struct DocumentAttributePos
{
};

struct DocumentNodePos
{
  std::string filePath;
  int offset;
  int line;
  int offsetInLine;
};

class DocumentAttribute
{
private:
  std::string m_key;
  std::optional< std::string > m_value;
public:
  DocumentAttribute(const std::string & k) : m_key(k), m_value() {}
  DocumentAttribute(const std::string & k, const std::string & v) : m_key(k), m_value(v) {}
  std::string name() const { return m_key; }
  std::string value() const { return m_value.value(); }
  void set_value( std::string const & value ) { m_value = value; }
  bool empty() const { return !m_value.has_value(); }
  operator bool() const { return m_value.has_value(); }
};

class DocumentNode
{
private:
  std::string m_name;
  DocumentNode * m_parent;
  std::vector< DocumentNode > m_children;
  std::map< std::string, DocumentAttribute > m_attributes;

public:
  DocumentNode() : m_name(""), m_parent(nullptr) { }
  DocumentNode(const std::string & name) : m_name(name), m_parent(nullptr) { }
  DocumentNode(const std::string & name, DocumentNode * parent ) : m_name(name), m_parent(parent) {  }

  DocumentNode & append_child( const std::string & name )
  {
    m_children.emplace_back( name, this );
    return m_children.back();
  }

  DocumentAttribute & append_attribute( const std::string & name )
  {
    m_attributes.emplace( std::make_pair( name, DocumentAttribute( name, "" ) ) );
    return m_attributes.at( name );
  }

  DocumentAttribute & attribute( const std::string & name )
  {
    auto it = m_attributes.find(name);
    if (it != m_attributes.end())
    {
      return it->second;
    }
    static DocumentAttribute empty("", "");
    return empty;  // Return an empty attribute if not found
  }

  DocumentAttribute const & attribute( const std::string & name ) const
  {
    auto it = m_attributes.find(name);
    if (it != m_attributes.end())
    {
      return it->second;
    }
    static DocumentAttribute empty("", "");
    return empty;  // Return an empty attribute if not found
  }

  auto attributes( ) const
  {
    std::vector< DocumentAttribute > atts;
    for( auto att : m_attributes )
    {
      atts.push_back(att.second);
    }
    return atts;
  }

  std::string name() const { return m_name; }
  DocumentNode & parent() const { return *m_parent; }
  std::string path() const { return m_name; }

  const std::vector< DocumentNode > & children() const { return m_children; }
  std::vector< DocumentNode > & children() { return m_children; }

  friend bool operator==(const DocumentNode& a, const DocumentNode& b)
  {
    return a.name() == b.name();
  }
  friend bool operator!=(const DocumentNode& a, const DocumentNode& b)
  {
    return ! (a == b);
  }
  friend bool operator<(const DocumentNode& a, const DocumentNode& b)
  {
    return a.name() < b.name();
  }
};

class Document
{
private:
  DocumentNode & m_root;

public:
  using read_result_type = bool;
  using node_type = DocumentNode;
  using node_pos_type = DocumentNodePos;
  using attribute_type = DocumentAttribute;
  using attribute_pos_type = DocumentAttributePos;

  Document( DocumentNode & root ) : m_root( root ) {}

  DocumentNode & appendChild(const std::string & name)
  {
    return m_root.append_child(name);
  }

  DocumentNode & getRoot() { return m_root; }

  std::string getFilePath() { return "__internal__"; }

  DocumentNodePos getNodePosition( DocumentNode & docNode ) const
  {
    return DocumentNodePos{"__internal__",0,0,0};
  }
};

class DocumentContext
{
public:
  DocumentContext(const DocumentNode&, const DocumentNodePos&) {}
  DocumentContext(const DocumentNode&, const DocumentAttribute&, const DocumentNodePos&) {}
};

class TestWrapperBase
{
protected:
  std::string m_name;
public:
  TestWrapperBase( std::string const & name ) : m_name( name ) {}
  virtual bool initFromInput( DocumentNode const &, DocumentNodePos const & ) = 0;
  std::string const & getName() const
  {
    return m_name;
  }
  virtual void makePolymorphic() = 0;
 };

template < typename T >
class TestWrapper : public TestWrapperBase
{
public:
  TestWrapper( std::string const & name ) : TestWrapperBase(name), m_data() {}
  TestWrapper( std::string const & name, T data ) : TestWrapperBase(name), m_data( data ) {}
  TestWrapper( std::string const & name, T & data ) : TestWrapperBase(name), m_data( data ) {}
  virtual bool initFromInput( DocumentNode const & docNode, DocumentNodePos const & docNodePos )
  {
    m_data = docNode.attribute(m_name).value();
    return true;
  }
  T const & reference() const { return m_data; }
  T & reference() { return m_data; }
  virtual void makePolymorphic() override { }
private:
  T m_data;
};

class TestGroup
{
private:
  std::string m_name;
  TestGroup * m_parent;
  std::vector< std::shared_ptr< TestGroup > > m_children;
  std::vector< std::unique_ptr< TestWrapperBase > > m_wrappers;
  std::unique_ptr< DataContext > m_fileContext;

  // setting up mock wrappers for specific tests
  void addTestWrappers( )
  {
    if( m_name == "ConfigNode" )
    {
      addWrapper( std::make_unique< TestWrapper< std::string > >("name") );
      addWrapper( std::make_unique< TestWrapper< std::string > >("status") );
      addWrapper( std::make_unique< TestWrapper< std::string > >("mode") );
    }
    else if ( m_name == "DetailNode" )
    {
      addWrapper( std::make_unique< TestWrapper< std::string > >("name") );
      addWrapper( std::make_unique< TestWrapper< std::string > >("information") );
    }
  }

public:
  std::vector< inputExtension::Rule< DocumentNode > > m_inputExtensions;

  explicit TestGroup( const std::string & name ) : m_name(name), m_parent(nullptr) { addTestWrappers(); }
  explicit TestGroup( const std::string & name, TestGroup * parent ) : m_name(name), m_parent(parent) { addTestWrappers(); }

  const std::string & getName() const
  {
    return m_name;
  }

  std::shared_ptr< TestGroup > createChild( const std::string & docNodeName, const std::string & dataNodeName )
  {
    m_children.emplace_back(std::make_unique<TestGroup>(docNodeName, this));
    return m_children.back();
  }

  std::vector< std::shared_ptr< TestGroup > > & children()
  {
    return m_children;
  }

  std::shared_ptr< TestGroup > getGroupPointer( const std::string & childName )
  {
    for (auto & child : m_children)
    {
      if (child->getName() == childName)
      {
        return child;
      }
    }
    return std::unique_ptr< TestGroup >( nullptr );
  }

  void addWrapper( std::unique_ptr<TestWrapperBase> wrapper )
  {
    m_wrappers.push_back( std::move(wrapper) );
  }

  bool hasWrapper( std::string const & wrapperName )
  {
    for( auto & wrapper : m_wrappers )
    {
      if( wrapper->getName() == wrapperName )
      {
        return true;
      }
    }
    return false;
  }

  TestWrapperBase & getWrapperBase( std::string const & wrapperName )
  {
    for( auto & wrapper : m_wrappers )
    {
      if( wrapper->getName() == wrapperName )
      {
        return *wrapper;
      }
    }
    static TestWrapper< std::string > empty("");
    return empty;
  }

  template < typename T >
  TestWrapper< T > & getWrapper( std::string const & wrapperName )
  {
    for( auto & wrapper : m_wrappers )
    {
      if( wrapper->getName() == wrapperName )
      {
        return * dynamic_cast< TestWrapper< T > * >( wrapper.get() );
      }
    }
    static TestWrapper< T > empty("empty");
    return empty;
  }

  void registerDataContext( std::unique_ptr< DataContext > context )
  {
    m_fileContext = std::move( context );
  }

  DataContext const & getDataContext() const
  {
    return *m_fileContext;
  }

  std::string dumpInputOptions() const
  {
    return "";
  }

  void forWrappers(std::function<void(TestWrapperBase&)> func)
  {
    for (auto & wrapper : m_wrappers)
    {
      func(*wrapper);
    }
  }

  std::vector< inputExtension::Rule< DocumentNode > > getInputExtensionRules( inputExtension::InputExtender< DocumentNode > & )
  {
    return m_inputExtensions;
  }
};


// Test Fixture
template < typename PhaseType >
class InputProcessingPhaseTest : public ::testing::Test
{
protected:
  DocumentNode m_docroot;
  Document m_document;
  TestGroup m_group;
  std::set< std::string > m_mergableNodes;

  InputProcessingPhaseTest()
    : m_docroot("Problem")
    , m_document(m_docroot)
    , m_group("Problem")
    , m_mergableNodes()
  {
    m_mergableNodes.insert("Problem");
  }

  void ExecutePhase()
  {
    PhaseType phase( m_document, m_mergableNodes );
    phase.execute( m_group, m_docroot );
  }

  virtual void testChecks() = 0;
};

class DeclarationTest : public InputProcessingPhaseTest< Declaration< Document, TestGroup > >
{
  void SetUp() override
  {
    // mostly just tests things actually expand the document tree, and that multiple extensions on the same group are applied
    inputExtension::Rule< DocumentNode > extensionRule;
    extensionRule.determineRoot = inputExtension::singletonParsingContext::thisNode< DocumentNode >();
    extensionRule.appliesWhen = inputExtension::singletonParsingContext::isNamed< DocumentNode >( "Problem" );

    inputExtension::Node< DocumentNode > subNode;
    subNode.identifier = "SubNode";

    extensionRule.subTrees.push_back( subNode );
    m_group.m_inputExtensions.push_back( extensionRule );

    inputExtension::Rule< DocumentNode > subnodeExtensionRule;
    subnodeExtensionRule.determineRoot = inputExtension::singletonParsingContext::thisNodesChild< DocumentNode >( "SubNode" );
    subnodeExtensionRule.appliesWhen = inputExtension::singletonParsingContext::isNamed< DocumentNode >( "Problem" );

    inputExtension::Node< DocumentNode > subsubNode;
    subsubNode.identifier = "SubSubNode";

    subnodeExtensionRule.subTrees.push_back( subsubNode );
    m_group.m_inputExtensions.push_back( subnodeExtensionRule );
  }
protected:
  void testChecks() override
  {
    // Add specific checks for Declaration phase
    geos::StaticTreeIteration< >::processTree( m_docroot, []( DocumentNode & node )
    {
      std::cout << node.name() << std::endl;
    } );
    EXPECT_EQ( m_docroot.children().size(), 1 );
    EXPECT_EQ( m_docroot.children()[0].children().size(), 1 );
  }
};

TEST_F(DeclarationTest, TestDeclaration)
{
  ExecutePhase();
  testChecks();
}

// Define a rule configuration function for a specific rule
void configureConfigNodeRule( inputExtension::Rule< DocumentNode > & rule )
{
  // Define a terse rule for node expansion with dynamic attributes and subnodes
  rule.determineRoot = inputExtension::singletonParsingContext::thisNode< DocumentNode >();
  rule.appliesWhen = inputExtension::singletonParsingContext::isNamed< DocumentNode >( "Problem" );
  inputExtension::Node< DocumentNode > configNode;
  configNode.identifier = "ConfigNode";
  configNode.attributes.push_back( inputExtension::mapAttribute::newAttribute< DocumentNode >( "mode", "active" ) );
  configNode.dynamicAttributes = [](  inputExtension::StaticParsingContext< DocumentNode > & context )
  {
    std::vector< inputExtension::Attribute< DocumentNode > > attrs;
    attrs.push_back( inputExtension::mapAttribute::newAttribute<DocumentNode>("status", "enabled"));
    return attrs;
  };
  configNode.dynamicSubNodes = []( inputExtension::StaticParsingContext< DocumentNode > & context)
  {
    std::vector< inputExtension::Node< DocumentNode > > subNodes;
    inputExtension::Node< DocumentNode > detailNode;
    detailNode.identifier = "DetailNode";
    detailNode.attributes.push_back( inputExtension::mapAttribute::newAttribute< DocumentNode >( "information", "detail" ) );
    subNodes.push_back( detailNode );
    return subNodes;
  };
  rule.subTrees.push_back( configNode );
}

// Register the rule using the provided macro
REGISTER_SUGAR_EXTENSION_RULE( "Problem", DocumentNode, ConfigNodeRule, configureConfigNodeRule );

class TerseSyntaxTest : public InputProcessingPhaseTest< TerseSyntax< Document, TestGroup > >
{
protected:
  void SetUp() override
  { }
protected:
  void testChecks() override
  {
    geos::StaticTreeIteration< >::processTree( m_docroot, []( DocumentNode & node )
    {
      std::cout << node.name() << std::endl;
    } );
    // Add specific checks for Declaration phase
    EXPECT_EQ(m_docroot.children().size(), 1);
    auto & expandedNode = m_docroot.children()[0];
    EXPECT_EQ(expandedNode.name(), "ConfigNode");
    ASSERT_STREQ(expandedNode.attribute("mode").value().c_str(), "active");
    ASSERT_STREQ(expandedNode.attribute("status").value().c_str(), "enabled");
    EXPECT_EQ(expandedNode.children().size(), 1);
    ASSERT_STREQ(expandedNode.children()[0].name().c_str(), "DetailNode");
    ASSERT_STREQ(expandedNode.children()[0].attribute("information").value().c_str(), "detail");
  }
};

TEST_F( TerseSyntaxTest, TestTerseSyntax )
{
  ExecutePhase();
  testChecks();
}

class DefinitionTest : public ::testing::Test
{
protected:
  DocumentNode m_docroot;
  Document m_document;
  TestGroup m_group;
  std::set< std::string > m_mergableNodes;

  DefinitionTest()
    : m_docroot("Problem")
    , m_document(m_docroot)
    , m_group("Problem")
    , m_mergableNodes()
  {
    m_mergableNodes.insert("Problem");
    prepareInitialTree();
  }

  void prepareInitialTree()
  {
    // Simulate the expansions and modifications made in DeclarationTest and TerseSyntaxTest
    m_docroot.append_child("ConfigNode").append_attribute("mode").set_value("active");
    auto & configNode = m_docroot.children().back();
    configNode.append_attribute("name").set_value("configNode");
    configNode.append_attribute("status").set_value("enabled");
    configNode.append_child("DetailNode").append_attribute("information").set_value("detail");
    auto & detailNode = configNode.children().back();
    detailNode.append_attribute("name").set_value("detailNode");
  }

  void ExecutePhase()
  {
    Definition< Document, TestGroup, TestWrapperBase > phase( m_document, m_mergableNodes );
    std::cout << "preexisting document nodes: " << std::endl;
    geos::StaticTreeIteration< >::processTree( m_docroot, [] ( DocumentNode & node )
    {
      std::cout << node.name() << std::endl;
    } );
    phase.execute( m_group, m_docroot );
  }

  void testChecks()
  {
    std::cout << "document nodes:" << std::endl;
    geos::StaticTreeIteration< >::processTree( m_docroot, []( DocumentNode & node )
    {
      std::cout << node.name() << std::endl;
    } );
    std::cout << "groups:" << std::endl;
    geos::StaticTreeIteration< >::processTree( m_group, []( TestGroup & node )
    {
      std::cout << node.getName() << std::endl;
    } );
    EXPECT_TRUE( validateInternalDataStructures( m_group, m_docroot ) );
  }

  bool validateInternalDataStructures( TestGroup & group, const DocumentNode & docRoot )
  {
    // Start by validating the root node
    if ( group.getName() != docRoot.name() )
    {
      std::cout << "group name mismatch" << std::endl;
      return false;  // Root node names must match
    }
    return validateChildren( group, docRoot );
  }

  bool validateChildren( TestGroup & group, const DocumentNode & docNode )
  {
    // Check each child node in the DocumentNode against the TestGroup's children
    auto & groupChildren = group.children();
    const auto & docChildren = docNode.children();

    if (groupChildren.size() != docChildren.size())
    {
      std::cout << "child size mismatch" << std::endl;
      return false;  // The number of children must match
    }

    for ( size_t ii = 0; ii < docChildren.size(); ii++ )
    {
      const DocumentNode & docChild = docChildren[ii];
      const std::shared_ptr< TestGroup > & groupChild = groupChildren[ii];

      if ( groupChild->getName() != docChild.name() )
      {
        std::cout << "group name mismatch" << std::endl;
        return false;  // Child node names must match
      }

      // Recursively validate the children of each node
      if ( ! validateChildren( *groupChild, docChild ) )
      {
        return false;
      }

      // verify attributes are set
      for( auto & attr : docChild.attributes() )
      {
        if( ! groupChild->hasWrapper( attr.name() ) )
        {
          std::cout << "wrapper " << attr.name() << " doesn't exist" << std::endl;
          return false;
        }
        if( groupChild->getWrapper< std::string >( attr.name() ).reference() != attr.value() )
        {
          std::cout << "wrapper " << attr.name() << " value mismatch " << groupChild->getWrapper< std::string >( attr.name() ).reference() << " != " << attr.value() << std::endl;
          return false;
        }
      }
    }

    return true;
  }

};

TEST_F(DefinitionTest, TestDefinition)
{
  ExecutePhase();
  testChecks();
}
