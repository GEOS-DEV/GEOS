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

#include "../Patterns.hpp"

#include <gtest/gtest.h>
#include <memory>

using namespace geos;

// Mock Observer that records updates
template<typename T>
class MockObserver : public Observer<T>
{
public:
  T lastUpdate;
  bool updated = false;

  void update(const T& data) override
  {
    lastUpdate = data;
    updated = true;
  }
};

// Test fixture for Observable tests
template<typename T>
class ObservableTest : public ::testing::Test
{
protected:
  std::shared_ptr<Observable<T>> observable;
  std::shared_ptr<MockObserver<T>> observer1, observer2;

  void SetUp() override
  {
      observable = std::make_shared<Observable<T>>();
      observer1 = std::make_shared<MockObserver<T>>();
      observer2 = std::make_shared<MockObserver<T>>();
      observable->observe(observer1);
      observable->observe(observer2);
  }
};

// Define test types
using ObservableTypes = ::testing::Types<int, std::string>;
TYPED_TEST_SUITE(ObservableTest, ObservableTypes);

// Test observer notification
TYPED_TEST(ObservableTest, NotifiesObservers)
{
  TypeParam testData = TypeParam{};  // Initialize with default value
  this->observable->notify(testData);
  ASSERT_TRUE(this->observer1->updated);
  ASSERT_TRUE(this->observer2->updated);
  ASSERT_EQ(this->observer1->lastUpdate, testData);
  ASSERT_EQ(this->observer2->lastUpdate, testData);
}

// Test observer removal
TYPED_TEST(ObservableTest, RemovesObservers)
{
  this->observable->forget(this->observer1);
  TypeParam testData = TypeParam{};
  this->observable->notify(testData);
  ASSERT_FALSE(this->observer1->updated);
  ASSERT_TRUE(this->observer2->updated);
}


// Define a test fixture class for ActiveObserver tests
template<typename T>
class ActiveObserverTest : public ::testing::Test
{
protected:
  T observedValue;
  bool predicateCalled = false;
  bool callbackCalled = false;

  // Setup with default observed value
  void SetUp() override
  {
    observedValue = T{};
  }

  // Helper functions to simulate predicate and callback
  bool testPredicate(const T& val)
  {
    predicateCalled = true;
    return true;  // Always indicate a change
  }

  T testCallback(const T& val)
  {
    callbackCalled = true;
    return val;  // Return the input value
  }
};

// Define test types
using ActiveObserverTypes = ::testing::Types<int, std::string>;
TYPED_TEST_SUITE(ActiveObserverTest, ActiveObserverTypes);

// Test change detection
TYPED_TEST(ActiveObserverTest, DetectsChanges)
{
  ActiveObserver<TypeParam, TypeParam> observer(
    this->observedValue,
    [this](const TypeParam& val) { return this->testPredicate(val); },
    [this](const TypeParam& val) { return this->testCallback(val); }
  );

  observer.check();
  ASSERT_TRUE(this->predicateCalled);
  ASSERT_TRUE(this->callbackCalled);
}

// Test no change detection
TYPED_TEST(ActiveObserverTest, DoesNotDetectChangesWhenPredicateFalse)
{
  ActiveObserver<TypeParam, TypeParam> observer(
    this->observedValue,
    [this](const TypeParam& val) -> bool { this->predicateCalled = true; return false; },  // Always return false
    [this](const TypeParam& val) -> TypeParam { this->callbackCalled = true; return val; }
  );

  observer.check();
  ASSERT_TRUE(this->predicateCalled);
  ASSERT_FALSE(this->callbackCalled);  // Callback should not be called
}

// Simple Node class for testing
struct Node
{
  int m_value;
  std::list< Node > m_children;

  explicit Node(int val) : m_value(val) {}

  void addChild(const Node & child)
  {
    m_children.push_back(child);
  }

  std::list< Node > & children()
  {
    return m_children;
  }

  const std::list< Node > & children() const
  {
    return m_children;
  }

  friend bool operator<( const Node & one, const Node & two )
  {
    return one.m_value < two.m_value;
  }
};

// Test fixture for StaticTreeIteration tests
class StaticTreeIterationTest : public ::testing::Test
{
protected:
  Node root;

  StaticTreeIterationTest() : root(0) {}

  void SetUp() override
  {
    // Build a simple tree structure
    Node child1(1);
    child1.addChild(Node(11));
    child1.addChild(Node(12));

    Node child2(2);
    child2.addChild(Node(21));
    child2.addChild(Node(22));

    root.addChild(child1);
    root.addChild(child2);
  }
};

TEST_F(StaticTreeIterationTest, ProcessTreeAppliesFunctionToEveryNode)
{
  std::vector< int > values;
  StaticTreeIteration< >::processTree(root, [&](Node & node)
  {
    values.push_back(node.m_value);
  });

  std::vector< int > expected = {0, 1, 11, 12, 2, 21, 22};
  ASSERT_EQ(values, expected);
}

// Test fixture for DynamicTreeIteration tests
class DynamicTreeIterationTest : public ::testing::Test
{
protected:
  Node root;
  // we're just using the node values as identifiers for these tests
  struct ValueAsIdentity
  {
    auto operator()( Node const & node ) const
    {
      return node.m_value;
    }
  };

  DynamicTreeIterationTest() : root(0) {}

  void SetUp() override
  {
    // Initial tree setup
    root.addChild(Node(1));
    root.addChild(Node(2));
  }
};

TEST_F(DynamicTreeIterationTest, ProcessTreeHandlesDynamicAdditions)
{
  std::set<int> values;
  auto lambda = [&](Node & node)
  {
    values.insert(node.m_value);
    // Simulate dynamic addition
    if (node.m_value == 1)
    {
      node.addChild(Node(11));
    }
  };
  DynamicTreeIteration< ValueAsIdentity >::processTree( root, lambda );
  std::set<int> expected = {0, 1, 2, 11};  // Including dynamically added node 11
  ASSERT_EQ(values, expected);
}

TEST_F(DynamicTreeIterationTest, AddChildToGrandparent)
{
  DynamicTreeIteration<Node> iter;
  std::set<int> values;
  auto lambda = [&](Node & node)
  {
    values.insert(node.m_value);
    // Simulate dynamic addition to the root node when processing the first child
    if (node.m_value == 1)
    {
      root.addChild(Node(12));  // Adding a new child to the root
    }
  };
  DynamicTreeIteration< ValueAsIdentity >::processTree( root, lambda );
  std::set<int> expected = {0, 1, 2, 12};  // Including dynamically added node 12
  ASSERT_EQ(values, expected);
}

TEST_F(DynamicTreeIterationTest, MultipleLevelsOfAdditions)
{
  DynamicTreeIteration<Node> iter;
  std::set<int> values;
  auto lambda = [&](Node & node)
  {
    values.insert(node.m_value);
    if (node.m_value == 1)
    {
      node.addChild(Node(11));  // Add child
    }
    if (node.m_value == 11)
    {
      node.addChild(Node(111));  // Add child to newly added child
    }
  };
  DynamicTreeIteration< ValueAsIdentity >::processTree( root, lambda );
  std::set<int> expected = {0, 1, 2, 11, 111};
  ASSERT_EQ(values, expected);
}

TEST_F(DynamicTreeIterationTest, SimultaneousAdditions)
{
  DynamicTreeIteration<Node> iter;
  std::set<int> values;
  auto lambda = [&](Node & node)
  {
    values.insert(node.m_value);
    if (node.m_value == 1 || node.m_value == 2)
    {
      node.addChild(Node(node.m_value * 10));  // Each node adds a child
    }
  };
  DynamicTreeIteration< ValueAsIdentity >::processTree( root, lambda );
  std::set<int> expected = {0, 1, 2, 10, 20};
  ASSERT_EQ(values, expected);
}

TEST_F(DynamicTreeIterationTest, NoAdditions)
{
  DynamicTreeIteration<Node> iter;
  std::set<int> values;
  auto lambda = [&](Node & node)
  {
    values.insert(node.m_value);
  };
  DynamicTreeIteration< ValueAsIdentity >::processTree( root, lambda );
  std::set<int> expected = {0, 1, 2};  // No new nodes added
  ASSERT_EQ(values, expected);
}

TEST(FlagCallbacksTest, RegisterAndExecuteSingleCallback)
{
  FlagCallbacks< int, std::tuple< std::string, int > > callbacks;
  bool callbackCalled = false;

  std::function< void( std::string, int ) > callback = [&](std::string key, int value)
  {
    callbackCalled = true;
    EXPECT_EQ(key, "test_key");
    EXPECT_EQ(value, 42);
  };

  callbacks.registerCallback( 1, callback );

  callbacks.executeCallbacks({1}, "test_key", 42);
  EXPECT_TRUE(callbackCalled);
}

TEST(FlagCallbacksTest, ExecuteNonRegisteredCallback)
{
  FlagCallbacks< int, std::tuple< std::string, int > > callbacks;
  bool callbackCalled = false;

  std::function< void( std::string, int ) > callback = [&](std::string key, int value)
  {
    callbackCalled = true;
  };

  callbacks.registerCallback(1, callback);

  callbacks.executeCallbacks({2}, "test_key", 42);
  EXPECT_FALSE(callbackCalled);
}

TEST(FlagCallbackCoordinatorTest, RegisterAndExecuteCoordinatorCallback)
{
  FlagCallbackCoordinator< int, std::tuple< std::string, int > > coordinator;
  bool callbackCalled = false;

  std::function< void( std::string, int ) > callback = [&](std::string key, int value)
  {
    callbackCalled = true;
    EXPECT_EQ(key, "test_key");
    EXPECT_EQ(value, 42);
  };

  coordinator.registerCallback(1, callback);

  std::string keyObject = "test_key";
  coordinator.executeCallbacks({1}, keyObject, 42);
  EXPECT_TRUE(callbackCalled);
}

TEST(FlagCallbackCoordinatorTest, MultipleKeyTypes)
{
  FlagCallbackCoordinator< int, std::tuple< std::string, int >, std::tuple< int, int > > coordinator;
  bool stringCallbackCalled = false;
  bool intCallbackCalled = false;

  std::function< void( std::string, int ) > callback1 = [&](std::string key, int value)
  {
    stringCallbackCalled = true;
    EXPECT_EQ(key, "test_key");
    EXPECT_EQ(value, 42);
  };

  std::function< void( int, int ) > callback2 = [&](int key, int value)
  {
    intCallbackCalled = true;
    EXPECT_EQ(key, 7);
    EXPECT_EQ(value, 24);
  };

  coordinator.registerCallback(1, callback1);
  coordinator.registerCallback(2, callback2);

  std::string stringKeyObject = "test_key";
  int intKeyObject = 7;
  coordinator.executeCallbacks({1}, stringKeyObject, 42);
  coordinator.executeCallbacks({2}, intKeyObject, 24);
  EXPECT_TRUE(stringCallbackCalled);
  EXPECT_TRUE(intCallbackCalled);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}