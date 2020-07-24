// Python.h must be included first.
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "managers/initialization.hpp"
#include "managers/GeosxState.hpp"

// System includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#define NPY_NO_DEPRECATED_API NPY_1_15_API_VERSION
#include <numpy/arrayobject.h>
#pragma GCC diagnostic pop
#include <chrono>

namespace geosx
{

// TODO: corbett This should go in LvArray
class PyObjectRef
{
public:

  /**
   * @brief Create an uninitialized (nullptr) reference.
   */
  PyObjectRef() = default;
  
  /**
   * @brief Take ownership of a reference to @p src.
   * @p src The object to be referenced.
   */
  explicit PyObjectRef( PyObject * src ):
    m_object( src )
  {}

  // Theese could be implemented but I don't have a use case yet.
  PyObjectRef( PyObjectRef const & ) = delete;
  PyObjectRef( PyObjectRef && ) = delete;

  /**
   * @brief Decrease the reference count to the current object.
   */
  ~PyObjectRef()
  {
    if ( m_object != nullptr )
    { Py_DECREF( m_object ); }
  }

  PyObjectRef & operator=( PyObjectRef const & ) = delete;
  PyObjectRef & operator=( PyObjectRef && ) = delete;

  /**
   * @brief Decrease the reference count to the current object and take ownership
   *   of a new reference.
   * @p src The new object to be referenced.
   * @return *this.
   */
  PyObjectRef & operator=( PyObject * src )
  {
    if ( m_object != nullptr )
    { Py_DECREF( m_object ); }

    m_object = src;
    return *this;
  }

  operator PyObject*()
  { return m_object; }

private:
  PyObject * m_object = nullptr;
};



static std::unique_ptr< GeosxState > state;

/**
 * 
 */
static PyObject * initialize( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );

  GEOSX_ERROR_IF( state != nullptr, "Already initialized." );

  long pythonMPIRank;
  PyObject * list;
  if ( !PyArg_ParseTuple( args, "lO", &pythonMPIRank, &list ) )
  { return nullptr; }

  PyObjectRef iterator{ PyObject_GetIter( list ) };
  if ( iterator == nullptr )
  { return nullptr; }

  std::vector< std::string > stringArgs;
  for( PyObjectRef item{ PyIter_Next( iterator ) }; item != nullptr; item = PyIter_Next( iterator ) )
  {
    PyObjectRef ascii { PyUnicode_AsASCIIString( item ) };
    if ( ascii == nullptr )
    { return nullptr; }

    char const * const stringValue = PyBytes_AsString( ascii );
    if ( stringValue == nullptr )
    { return nullptr; }

    stringArgs.push_back( stringValue );
  }

  if ( PyErr_Occurred() )
  { return nullptr; }

  std::vector< char * > argv( stringArgs.size() + 1 );
  for ( std::size_t i = 0; i < stringArgs.size(); ++i )
  { argv[ i ] = const_cast< char * >( stringArgs[ i ].data() ); }

  basicSetup( argv.size() - 1, argv.data(), true );

  GEOSX_ERROR_IF_NE( pythonMPIRank, MpiWrapper::Comm_rank() );

  state = std::make_unique< GeosxState >();
  state->initializeDataRepository();

  Py_RETURN_NONE;
}

/**
 * 
 */
static PyObject * get( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );

  PyObject * unicodePath;
  int modify;
  if ( !PyArg_ParseTuple( args, "Up", &unicodePath, &modify ) )
  { return nullptr; }

  PyObjectRef asciiPath { PyUnicode_AsASCIIString( unicodePath ) };
  if ( asciiPath == nullptr )
  { return nullptr; }

  char const * const charPath = PyBytes_AsString( asciiPath );
  if ( charPath == nullptr )
  { return nullptr; }

  std::string groupPath, wrapperName;
  splitPath( charPath, groupPath, wrapperName );

  dataRepository::Group * const group = state->getGroupByPath( groupPath );
  if ( group == nullptr )
  {
    GEOSX_LOG_RANK( "The group doesn't exist: " << groupPath );
    Py_RETURN_NONE;
  }

  dataRepository::WrapperBase * const wrapper = group->getWrapperBase( wrapperName );
  if ( wrapper == nullptr )
  {
    GEOSX_LOG_RANK( "Goup " << groupPath << " doesn't have a wrapper " << wrapper );
    Py_RETURN_NONE;
  }

  PyObject * const ret = wrapper->createPythonObject( modify );
  if ( ret == nullptr )
  { Py_RETURN_NONE; }
  else
  { return ret; }
}

/**
 * 
 */
static PyObject * applyInitialConditions( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );

  if ( !PyArg_ParseTuple( args, "" ) )
  { return nullptr; }

  GEOSX_ERROR_IF( state == nullptr, "State must be initialized." );
  state->applyInitialConditions();
  Py_RETURN_NONE;
}

/**
 * 
 */
static PyObject * run( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );

  if ( !PyArg_ParseTuple( args, "" ) )
  { return nullptr; }

  GEOSX_ERROR_IF( state == nullptr, "State must be initialized." );
  state->run();
  return PyLong_FromLong( static_cast< int >( state->getState() ) );
}

/**
 * 
 */
static PyObject * finalize( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );

  if ( !PyArg_ParseTuple( args, "" ) )
  { return nullptr; }

  state = nullptr;
  basicCleanup();

  Py_RETURN_NONE;
}

} // namespace geosx

/**
 * 
 */
static PyMethodDef pygeosxFuncs[] = {
  { "initialize", geosx::initialize, METH_VARARGS, 
    "" },
  { "get", geosx::get, METH_VARARGS, 
    "" },
  { "applyInitialConditions", geosx::applyInitialConditions, METH_VARARGS, 
    "" },
  { "run", geosx::run, METH_VARARGS, 
    "" },
  { "finalize", geosx::finalize, METH_VARARGS, 
    "" },
  { NULL, NULL, 0, NULL }        /* Sentinel */
};

/**
 * Initialize the module object for Python with the exported functions
 */
static struct PyModuleDef pygeosxModuleFunctions = {
  PyModuleDef_HEAD_INIT,
  "pygeosx",   /* name of module */
  "Module for testing numpy views of LvArray::Array objects", /* module documentation, may be NULL */
  -1,       /* size of per-interpreter state of the module or -1 if the module keeps state in global variables. */
  pygeosxFuncs,
  NULL,
  NULL,
  NULL,
  NULL,
};

/**
 * 
 */
PyMODINIT_FUNC
PyInit_pygeosx(void)
{
  import_array();
  return PyModule_Create(&pygeosxModuleFunctions);
}
