#include <gtest/gtest.h>
#include "../InputExtension.hpp"
#include "../xmlWrapper.hpp"

using namespace geos;
using namespace geos::dataRepository::inputExtension;

class InputExtensionsTest : public ::testing::Test
{
protected:
  xmlWrapper::xmlDocument doc;
  // Helper function to create a dummy XML node
  xmlWrapper::xmlNode createDummyNode(const std::string& name)
  {
    return doc.appendChild(name.c_str());
  }
};

// Test StaticParsingContext
TEST_F(InputExtensionsTest, StaticParsingContextSingleNode)
{
  auto node = createDummyNode("TestNode");
  StaticParsingContext<> context({node});

  EXPECT_EQ(context.scope(), 1);
  ASSERT_STREQ(context.getNode(0).name(), "TestNode");
}

// Test singletonParsingContext::isValid
TEST_F(InputExtensionsTest, SingletonParsingContextIsValid)
{
  auto node = createDummyNode("TestNode");
  StaticParsingContext<> context({node});

  auto isValid = singletonParsingContext::isValid<xmlWrapper::xmlNode>();
  EXPECT_TRUE(isValid(context));
}

// Test singletonParsingContext::isNamed
TEST_F(InputExtensionsTest, SingletonParsingContextIsNamed)
{
  auto node = createDummyNode("TestNode");
  StaticParsingContext<> context({node});

  auto isNamed = singletonParsingContext::isNamed<xmlWrapper::xmlNode>("TestNode");
  EXPECT_TRUE(isNamed(context));
}

// Test Attribute transformations
TEST_F(InputExtensionsTest, AttributeTransformations)
{
  auto node = createDummyNode("TestNode");
  StaticParsingContext<> context({node});

  auto identityAttr = mapAttribute::identity<>("oldName", "newName");
  EXPECT_EQ(identityAttr.transform("value"), "value");

  auto appendAttr = mapAttribute::append<>("oldName", "newName", "_suffix");
  EXPECT_EQ(appendAttr.transform("value"), "value_suffix");

  auto atomAsListAttr = mapAttribute::atomAsList<>("oldName", "newName");
  EXPECT_EQ(atomAsListAttr.transform("value"), "{ value }");
}

// Test Rule application
TEST_F(InputExtensionsTest, RuleApplication)
{
  auto node = createDummyNode("TestNode");
  StaticParsingContext<> context({node});

  Rule< > rule;
  rule.subTrees.push_back({ "SubNode" });
  rule.determineRoot = singletonParsingContext::thisNode<xmlWrapper::xmlNode>();

  InputExtender<> extender({});
  std::vector< Rule<> > extensions{ rule };
  extender.applyExtensions(context, extensions);

  ASSERT_STREQ(node.first_child().name(), "SubNode");
}

// Test dynamic attribute generation
TEST_F(InputExtensionsTest, DynamicAttributeGeneration)
{
  auto node = createDummyNode("TestNode");
  auto oldAttr = node.append_attribute("oldName");
  oldAttr.set_value("1.0");

  StaticParsingContext<> context({node});

  Rule< > rule;
  rule.determineRoot = singletonParsingContext::thisNode<xmlWrapper::xmlNode>();

  Node< > subNode;
  subNode.identifier = "SubNode";
  subNode.dynamicAttributes = [](const StaticParsingContext<>&)
  {
    // Take the oldAttribute from the context root and map it onto the subnode
    Attribute< > attr;
    attr.source = "oldName";
    attr.destination = "newName";
    attr.sourceNode = []( auto & c ){ return c.getNode(0); };
    return std::vector< Attribute<> >{ attr };
  };
  rule.subTrees.push_back( subNode );

  InputExtender<> extender({});
  std::vector< Rule< > > extensions{ rule };
  extender.applyExtensions(context, extensions );

  EXPECT_TRUE(node.first_child().attribute("newName"));
}

TEST_F(InputExtensionsTest, DynamicAttributeHandling)
{
    auto node = createDummyNode("DynamicNode");
    node.append_attribute("initial").set_value("original");

    StaticParsingContext<> context({node});

    Rule<> rule;
    rule.determineRoot = singletonParsingContext::thisNode<xmlWrapper::xmlNode>();

    Node<> subNode;
    subNode.identifier = "SubNode";
    subNode.dynamicAttributes = [](const StaticParsingContext<>& ctx) {
        Attribute<> attr;
        attr.source = "initial";
        attr.destination = "transformed";
        attr.transform = [](const std::string & value) { return value + "_modified"; };
        attr.sourceNode = [](auto & c) { return c.getNode(0); };
        return std::vector<Attribute<>>{attr};
    };
    rule.subTrees.push_back(subNode);

    InputExtender<> extender({});
    std::vector<Rule<>> extensions{rule};
    extender.applyExtensions(context, extensions);

    EXPECT_STREQ(node.first_child().attribute("transformed").value(), "original_modified");
}

TEST_F(InputExtensionsTest, MergableNodeCreation)
{
    auto node = createDummyNode("TestNode");
    node.append_child("ContainerNode");

    StaticParsingContext<> context({node});

    std::set<std::string> mergable{"ContainerNode"};
    InputExtender<> extender(mergable);

    Node<> containerNode;
    containerNode.identifier = "ContainerNode";
    Rule<> rule;
    rule.subTrees.push_back(containerNode);
    rule.determineRoot = singletonParsingContext::thisNode<xmlWrapper::xmlNode>();

    std::vector<Rule<>> extensions{rule};
    extender.applyExtensions(context, extensions);

    // Verify that no extra "ContainerNode" is created under the original "ContainerNode"
    ASSERT_EQ(std::distance(node.children().begin(), node.children().end()), 1);
}

static bool called = false;
void registerCall() { called = true; }

TEST_F(InputExtensionsTest, CallbackExecution)
{
    auto node = createDummyNode("CallbackNode");
    StaticParsingContext<> context({node});

    Rule< > rule;
    rule.determineRoot = singletonParsingContext::thisNode< xmlWrapper::xmlNode >();
    rule.flags.insert( 1 ); // Example flag to trigger callback

    InputExtender< > extender( {} );
    std::function< void( Rule< > &, StaticParsingContext< > &, xmlWrapper::xmlNode & ) > callback  = []( Rule< > & r, StaticParsingContext< > &, xmlWrapper::xmlNode & )
    {
      registerCall();
    };
    extender.registerCallback( 1, callback );

    std::vector< Rule< > > extensions{rule};
    extender.applyExtensions( context, extensions );

    EXPECT_TRUE( called );
}


// Test InputRemover
// TEST_F(InputExtensionsTest, InputRemover)
// {
//   auto node = createDummyNode("TestNode");
//   StaticParsingContext<> context({node});

//   Rule<> rule;
//   rule.subTrees.push_back({ "SubNode" });
//   rule.appliesWhen = [](const StaticParsingContext<>&) { return true; };
//   rule.determineRoot = singletonParsingContext::thisNode<>();

//   internal::InputRemover<> remover;
//   std::vector< Rule < > > extensions{ rule };
//   remover.applyRemovals(context, extensions);

//   EXPECT_FALSE(node.first_child());
// }
