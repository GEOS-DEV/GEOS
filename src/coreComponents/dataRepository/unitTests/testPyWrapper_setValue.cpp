#include <gtest/gtest.h>
#include <Python.h>
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"

#include "dataRepository/python/PyWrapper.hpp"



using namespace geos::python;
using namespace geos::dataRepository;


class PythonEnvironment : public ::testing::Environment
{
public:
  void SetUp() override { Py_Initialize(); }
  void TearDown() override { Py_Finalize(); }
};

::testing::Environment * const pythonEnv = ::testing::AddGlobalTestEnvironment( new PythonEnvironment );


template< typename T >
void testSetValue( PyObject * arg, const T & expected, const std::string & typeName )
{
  conduit::Node node;
  Group rootGroup( "root", node );
  Group dummyGroup( "dummy", &rootGroup );
  T * value = new T( expected );
  auto * wrapper = new Wrapper< T >( "test", dummyGroup, value );

  PyWrapper pyWrapper;
  pyWrapper.wrapper = wrapper;

  PyObject * result = PyWrapper_setValue( &pyWrapper, arg );
  ASSERT_NE( result, nullptr ) << "Failed to set value for type: " << typeName;

  EXPECT_EQ( wrapper->reference(), expected );

  delete wrapper;
  Py_DECREF( arg );
  Py_XDECREF( result );
}


TEST( PyWrapperSetValue, SetString )
{
  PyObject * arg = Py_BuildValue( "(s)", "hello" );
  testSetValue< std::string >( arg, "hello", "string" );
}

TEST( PyWrapperSetValue, SetDouble )
{
  PyObject * arg = Py_BuildValue( "(d)", 3.14 );
  testSetValue< double >( arg, 3.14, "double" );
}

TEST( PyWrapperSetValue, SetInt )
{
  PyObject * arg = Py_BuildValue( "(i)", 42 );
  testSetValue< int >( arg, 42, "int" );
}

TEST( PyWrapperSetValue, SetLong )
{
  PyObject * arg = Py_BuildValue( "(l)", 123456789L );
  testSetValue< long >( arg, 123456789L, "long" );
}


TEST( PyWrapperSetValue, TypeMismatch )
{
  conduit::Node node;
  Group rootGroup( "root", node );
  Group dummyGroup( "dummy", &rootGroup );
  int * value = new int(0);
  auto * wrapper = new Wrapper< int >( "test", dummyGroup, value );

  PyWrapper pyWrapper;
  pyWrapper.wrapper = wrapper;

  PyObject * arg = Py_BuildValue( "(s)", "not an int" );
  PyObject * result = PyWrapper_setValue( &pyWrapper, arg );

  EXPECT_EQ( result, nullptr );
  EXPECT_TRUE( PyErr_Occurred() );
  PyErr_Clear();

  delete wrapper;
  Py_DECREF( arg );
}


TEST( PyWrapperSetValue, UnsupportedType )
{
  conduit::Node node;
  Group rootGroup( "root", node );
  Group dummyGroup( "dummy", &rootGroup );
  auto * wrapper = new Wrapper< std::vector< int > >( "test", dummyGroup, new std::vector< int >() );

  PyWrapper pyWrapper;
  pyWrapper.wrapper = wrapper;

  PyObject * dict = PyDict_New();
  PyObject * arg = Py_BuildValue( "(O)", dict );
  PyObject * result = PyWrapper_setValue( &pyWrapper, arg );

  EXPECT_EQ( result, nullptr );
  EXPECT_TRUE( PyErr_Occurred());
  PyErr_Clear();

  delete wrapper;
  Py_DECREF( dict );
  Py_DECREF( arg );
}
