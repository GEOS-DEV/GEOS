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

  PyObject * get() const
  { return m_object; }

private:
  PyObject * m_object = nullptr;
};

static std::unique_ptr< GeosxState > state;

/**
 *
 */
struct Group
{
  PyObject_HEAD
  dataRepository::Group * group;
};






/**
 *
 */
static PyObject * initialize( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );

  if ( state != nullptr ){
    PyErr_SetString( PyExc_RuntimeError, "state already initialized" );
    return NULL;
  }

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

  if ( pythonMPIRank != MpiWrapper::Comm_rank() ){
    PyErr_SetString( PyExc_ValueError, "Python MPI rank does not align with GEOSX MPI rank" );
    return NULL;
  }

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
  GEOSX_UNUSED_VAR( args );

  if ( state == nullptr ){
    PyErr_SetString( PyExc_RuntimeError, "state must be initialized" );
    return NULL;
  }
  state->applyInitialConditions();
  Py_RETURN_NONE;
}

/**
 *
 */
static PyObject * run( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );
  GEOSX_UNUSED_VAR( args );

  if ( state == nullptr ){
    PyErr_SetString( PyExc_RuntimeError, "state must be initialized" );
    return NULL;
  }
  state->run();
  return PyLong_FromLong( static_cast< int >( state->getState() ) );
}

/**
 *
 */
static PyObject * finalize( PyObject * self, PyObject * args )
{
  GEOSX_UNUSED_VAR( self );
  GEOSX_UNUSED_VAR( args );

  if ( state == nullptr ){
    PyErr_SetString( PyExc_RuntimeError, "state already finalized" );
    return NULL;
  }
  state = nullptr;
  basicCleanup();

  Py_RETURN_NONE;
}

} // namespace geosx

/**
 * Add geosx::State enums to the given module. Return the module, or NULL on failure
 */
static PyObject * addConstants( PyObject * module )
{
  std::array< std::pair< long, char const * >, 4 > const constants = { {
    { static_cast< long >( geosx::State::COMPLETED ), "COMPLETED" },
    { static_cast< long >( geosx::State::INITIALIZED ), "INITIALIZED" },
    { static_cast< long >( geosx::State::UNINITIALIZED ), "UNINITIALIZED" },
    { static_cast< long >( geosx::State::READY_TO_RUN ), "READY_TO_RUN" }
  } };

  for ( std::pair< long, char const * > const & pair : constants )
  {
    if ( PyModule_AddIntConstant( module, pair.second, pair.first ) )
    {
      PyErr_SetString( PyExc_RuntimeError, "couldn't add constant" );
      return nullptr;
    }
  }

  return module;
}

/**
 * Add exit handler to a module. Register the module's 'finalize' function,
 * which should take no arguments, with the `atexit` standard library module.
 */
static bool addExitHandler( PyObject * module ){
  geosx::PyObjectRef atexit_module { PyImport_ImportModule( "atexit" ) };
  
  if ( atexit_module == nullptr )
  { return 0; }

  geosx::PyObjectRef atexit_register_pyfunc { PyObject_GetAttrString( atexit_module, "register" ) };
  if ( atexit_register_pyfunc == nullptr )
  { return 0; }

  geosx::PyObjectRef finalize_pyfunc { PyObject_GetAttrString( module, "finalize" ) };
  if ( finalize_pyfunc == nullptr )
  { return 0; }

  if ( !PyCallable_Check( atexit_register_pyfunc ) || !PyCallable_Check( finalize_pyfunc ) )
  { return 0; }

  geosx::PyObjectRef returnval { PyObject_CallFunctionObjArgs( atexit_register_pyfunc, finalize_pyfunc.get(), NULL ) };
  
  return returnval != nullptr;
}

// Allow mixing designated and non-designated initializers in the same initializer list.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-designator"

static PyTypeObject GroupType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
  .tp_name = "geosx.Group",
  .tp_basicsize = sizeof( geosx::Group ),
  .tp_itemsize = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT,
  .tp_doc = "",
  .tp_new = PyType_GenericNew,
};

/**
 *
 */
static PyMethodDef pygeosxFuncs[] = {
  { "initialize", geosx::initialize, METH_VARARGS,
    "" },
  { "get", geosx::get, METH_VARARGS,
    "" },
  { "applyInitialConditions", geosx::applyInitialConditions, METH_NOARGS,
    "" },
  { "run", geosx::run, METH_NOARGS,
    "" },
  { "finalize", geosx::finalize, METH_NOARGS,
    "" },
  { NULL, NULL, 0, NULL }        /* Sentinel */
};

/**
 * Initialize the module object for Python with the exported functions
 */
static struct PyModuleDef pygeosxModuleFunctions = {
  PyModuleDef_HEAD_INIT,
  .m_name = "pygeosx",
  .m_doc = "Module for testing numpy views of LvArray::Array objects",
  .m_size = -1,
  .m_methods = pygeosxFuncs
};

/**
 * Initialize the module with functions, constants, and exit handler
 */
PyMODINIT_FUNC
PyInit_pygeosx(void)
{
  import_array();

  geosx::PyObjectRef module{ PyModule_Create( &pygeosxModuleFunctions ) };
  if ( module == nullptr )
  { return nullptr; }

  if( PyType_Ready( &GroupType ) < 0 )
  { return nullptr; }

  Py_INCREF( &GroupType );
  if ( PyModule_AddObject( module, "Group", reinterpret_cast< PyObject * >( &GroupType ) ) < 0 )
  {
    Py_DECREF( &GroupType );
    return nullptr;
  }

  if ( !addExitHandler( module ) )
  {
    if ( PyErr_Occurred() == nullptr )
    { PyErr_SetString( PyExc_RuntimeError, "couldn't add exit handler" ); }
    return nullptr;
  }

  return addConstants( module );
}

#pragma GCC diagnostic pop
