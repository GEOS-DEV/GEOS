#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "../pygeosx.hpp"
#include "PyHistory.hpp"

#include "fileIO/timeHistory/TimeHistoryCollection.hpp"


#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->collection.size() == 0, PyExc_RuntimeError, "The PyHistory is not initialized.", nullptr )


namespace geosx
{
extern std::unique_ptr< GeosxState > g_state;

extern bool g_alreadyInitialized;

namespace python
{

struct PyHistory
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to HistoryCollection.";

  std::map< string, geosx::HistoryCollection * > collection;
  std::map< string, geosx::TimeHistoryOutput * > output;
};

static void PyHistory_dealloc( PyHistory * self )
{
  Py_TYPE( self )->tp_free((PyObject *)self );
}



static PyObject * PyHistory_new( PyTypeObject *type, PyObject *args, PyObject *kwds )
{
  GEOSX_UNUSED_VAR( args, kwds );
  PyHistory *self;
  std::map< string, geosx::HistoryCollection * > collection;
  std::map< string, geosx::TimeHistoryOutput * > output;

  self = (PyHistory *)type->tp_alloc( type, 0 );
  if( self != nullptr )
  {
    self->collection = collection;
    self->output = output;
  }

  return (PyObject *)self;
}



static int PyHistory_init( PyHistory *self, PyObject *args, PyObject *kwds )
{
  GEOSX_UNUSED_VAR( kwds, args );

  PYTHON_ERROR_IF( !g_alreadyInitialized, PyExc_RuntimeError, "GEOSX must be initialized first. Call pygeosx.initialize().", -1 );

  geosx::HistoryCollection *collection=nullptr;
  geosx::TimeHistoryOutput *output=nullptr;

  geosx::EventManager & eventManager = g_state->getProblemManager().getEventManager();

  eventManager.forSubGroups< geosx::EventBase >( [&]( geosx::EventBase & subEvent )
  {
    subEvent.getTargetReferences();
  } );

  geosx::EventBase * subEvent = nullptr;
  for( int currentSubEvent = 0; currentSubEvent<eventManager.numSubGroups(); ++currentSubEvent )
  {
    subEvent = static_cast< geosx::EventBase * >( eventManager.getSubGroups()[currentSubEvent]);
    string const subEventName = subEvent->getTargetName();

    if( subEventName.find( "Tasks" ) > 0 && subEventName.find( "Tasks" ) != std::string::npos )
    {
      collection = static_cast< geosx::HistoryCollection * >(subEvent->getEventTarget());
      int const firstChar = subEventName.find( "/", 1 ) + 1;
      int const lenName = subEventName.find( "Collection" ) - firstChar;
      self->collection.insert( {subEventName.substr( firstChar, lenName ), collection} );
    }
    else if( subEventName.find( "Outputs" ) > 0 && subEventName.find( "Outputs" ) != std::string::npos )
    {
      output = static_cast< geosx::TimeHistoryOutput * >(subEvent->getEventTarget());
      int const firstChar = subEventName.find( "/", 1 ) + 1;
      int const lenName = subEventName.find( "Output", 2 ) - firstChar;
      self->output.insert( {subEventName.substr( firstChar, lenName ), output} );
    }
    else
    {
      subEvent = nullptr;
    }
  }
  PYTHON_ERROR_IF( self->collection.size() == 0, PyExc_RuntimeError, "Target not found", -1 );
  PYTHON_ERROR_IF( self->output.size() == 0, PyExc_RuntimeError, "Target not found", -1 );

  return 0;
}



static PyObject * PyHistory_repr( PyObject * const obj ) noexcept
{
  PyHistory const * const pyHistory = LvArray::python::convert< PyHistory >( obj, getPyHistoryType() );
  if( pyHistory == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pyHistory );

  //std::map<std::string, geosx::PackCollection*>::iterator it;

  string repr = "Collections : \n";
  for( auto it = (pyHistory->collection).begin(); it != (pyHistory->collection).end(); it++ )
  {
    string const path = it->second->getPath();
    string const type = LvArray::system::demangle( typeid( *(it->second) ).name() );
    repr += path + " ( " + type + " )";
  }
  repr+="\nOutputs : \n";
  for( auto it = (pyHistory->output).begin(); it != (pyHistory->output).end(); it++ )
  {
    string const path = it->second->getPath();
    string const type = LvArray::system::demangle( typeid( *(it->second) ).name() );
    repr += path + " ( " + type + " )";
  }

  return PyUnicode_FromString( repr.c_str() );

}



static PyObject * collect( PyHistory * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  double time;
  double dt;
  PyObject * unicodePath;

  if( !PyArg_ParseTuple( args, "Udd", &unicodePath, &time, &dt ) )
  {
    return nullptr;
  }

  LvArray::python::PyObjectRef<> asciiPath { PyUnicode_AsASCIIString( unicodePath ) };
  if( asciiPath == nullptr )
  {
    return nullptr;
  }

  char const * const key = PyBytes_AsString( asciiPath );
  if( key == nullptr )
  {
    return nullptr;
  }

  geosx::DomainPartition & domain = g_state->getProblemManager().getDomainPartition();

  try
  {
    self->collection.at( key )->geosx::HistoryCollection::execute( time, dt, 0, 0, 0, domain );
  }
  catch( std::out_of_range const & e )
  {
    std::cout<<"Target \""<<key<<"\" not found. Impossible output."<<std::endl;
  }
  Py_RETURN_NONE;
}


static PyObject * output( PyHistory * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  double time;
  double dt;
  PyObject * unicodePath;

  if( !PyArg_ParseTuple( args, "Udd", &unicodePath, &time, &dt ) )
  {
    return nullptr;
  }

  LvArray::python::PyObjectRef<> asciiPath { PyUnicode_AsASCIIString( unicodePath ) };
  if( asciiPath == nullptr )
  {
    return nullptr;
  }

  char const * const key = PyBytes_AsString( asciiPath );
  if( key == nullptr )
  {
    return nullptr;
  }

  geosx::DomainPartition & domain = g_state->getProblemManager().getDomainPartition();

  try
  {
    self->output.at( key )->execute( time, dt, 0, 0, 0, domain );
  }
  catch( std::out_of_range const & e )
  {
    std::cout<<"Target \""<<key<<"\" not found. Impossible output."<<std::endl;
  }
  Py_RETURN_NONE;
}

static PyMethodDef PyHistory_methods[] = {
  { "collect", (PyCFunction) collect, METH_VARARGS, "wrapper to routine TimeHistoryCollection::execute" },
  { "output", (PyCFunction) output, METH_VARARGS, "wrapper to routine TimeHistoryOutput::execute"},
  { nullptr, nullptr, 0, nullptr }      /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PyHistoryType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pysolver.History",
  .tp_basicsize = sizeof( PyHistory ),
  .tp_itemsize = 0,
  .tp_dealloc = (destructor)PyHistory_dealloc,
  .tp_repr = PyHistory_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PyHistory::docString,
  .tp_methods = PyHistory_methods,
  .tp_init = (initproc)PyHistory_init,
  .tp_new = PyHistory_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS

static PyModuleDef pyhistorymodule = {
  PyModuleDef_HEAD_INIT,
  "pyhistory",
  "pyhistory module for HistoryCollection",
  -1,
  NULL, NULL, NULL, NULL, NULL
};


PyMODINIT_FUNC
PyInit_pyhistory( void )
{
  LvArray::python::PyObjectRef<> module{ PyModule_Create( &pyhistorymodule ) };
  if( module == nullptr )
  {
    return nullptr;
  }

  if( !LvArray::python::addTypeToModule( module, geosx::python::getPyHistoryType(), "History" ) )
  {
    return nullptr;
  }

  return module.release();
}

PyTypeObject * getPyHistoryType()
{ return &PyHistoryType; }

}
}
