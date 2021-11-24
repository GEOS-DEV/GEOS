#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "../pygeosx.hpp"
#include "PyHDF5.hpp"

#include "fileIO/timeHistory/TimeHistoryCollection.hpp"


#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->collection.size() == 0, PyExc_RuntimeError, "The PyHDF5 is not initialized.", nullptr )


namespace geosx
{
  extern std::unique_ptr< GeosxState > g_state;

  extern bool g_alreadyInitialized;

namespace python
{

struct PyHDF5
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to HDF5.";

  std::map<string, geosx::PackCollection*> collection;
  std::map<string, geosx::TimeHistoryOutput*> output;
  geosx::ProblemManager * pb_manager;
};

static void PyHDF5_dealloc(PyHDF5* self)
{
    Py_TYPE(self)->tp_free((PyObject*)self);
}



static PyObject * PyHDF5_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  GEOSX_UNUSED_VAR( args, kwds );
  PyHDF5 *self;
  geosx::ProblemManager * pb_manager = &(g_state->getProblemManager());
  std::map<string, geosx::PackCollection*> collection;
  std::map<string, geosx::TimeHistoryOutput*> output;

    self = (PyHDF5 *)type->tp_alloc(type, 0);
    if (self != nullptr)
    {
      self->collection = collection;
      self->output = output;
      self->pb_manager = pb_manager;

      if (self->pb_manager == nullptr)
      {
	Py_DECREF(self);
	return nullptr;
      }

    }

    return (PyObject *)self;
}



static int PyHDF5_init(PyHDF5 *self, PyObject *args, PyObject *kwds)
{
  GEOSX_UNUSED_VAR( kwds, args );

  PYTHON_ERROR_IF( !g_alreadyInitialized, PyExc_RuntimeError, "GEOSX must be initialized first. Call pygeosx.initialize().", -1 );

  geosx::PackCollection *collection=nullptr;
  geosx::TimeHistoryOutput *output=nullptr;

  geosx::EventManager & eventManager = self->pb_manager->getEventManager();

  eventManager.forSubGroups< geosx::EventBase >( [&]( geosx::EventBase & subEvent )
  {
    subEvent.getTargetReferences();
  } );

  geosx::EventBase * subEvent = nullptr;
  for(int currentSubEvent = 0; currentSubEvent<eventManager.numSubGroups(); ++currentSubEvent )
  {
    subEvent = static_cast< geosx::EventBase * >( eventManager.getSubGroups()[currentSubEvent]);
    string const subEventName = subEvent->getEventName();

    if (subEventName.find("Tasks") > 0 && subEventName.find("Tasks") != std::string::npos)
    {
      collection = static_cast<geosx::PackCollection*>(subEvent->getEventTarget());
      int const firstChar = subEventName.find("/", 1) + 1;
      int const lenName = subEventName.find("Collection") - firstChar;
      self->collection.insert( {subEventName.substr(firstChar, lenName), collection} );
    }
    else if (subEventName.find("Outputs") > 0 && subEventName.find("Outputs") != std::string::npos)
    {
      output = static_cast<geosx::TimeHistoryOutput*>(subEvent->getEventTarget());
      int const firstChar = subEventName.find("/", 1) + 1;
      int const lenName = subEventName.find("Output", 2) - firstChar;
      self->output.insert( {subEventName.substr(firstChar, lenName), output} );
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



static PyObject * PyHDF5_repr( PyObject * const obj ) noexcept
{
  PyHDF5 const * const pyHDF5 = LvArray::python::convert< PyHDF5 >( obj, getPyHDF5Type() );
  if( pyHDF5 == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pyHDF5 );

  //std::map<std::string, geosx::PackCollection*>::iterator it;

  string repr = "Collections : \n";
  for (auto it = (pyHDF5->collection).begin(); it != (pyHDF5->collection).end(); it++)
  {
    string const path = it->second->getPath();
    string const type = LvArray::system::demangle( typeid( *(it->second) ).name() );
    repr += path + " ( " + type + " )";
  }
  repr+="\nOutputs : \n";
  for (auto it = (pyHDF5->output).begin(); it != (pyHDF5->output).end(); it++)
  {
    string const path = it->second->getPath();
    string const type = LvArray::system::demangle( typeid( *(it->second) ).name() );
    repr += path + " ( " + type + " )";
  }

  return PyUnicode_FromString( repr.c_str() );

}



static PyObject * collect(PyHDF5 * self, PyObject * args)
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

  geosx::DomainPartition & domain = self->pb_manager->getDomainPartition();

  try
  {
    self->collection.at(key)->geosx::HistoryCollection::execute(time, dt, 0, 0, 0, domain);
  }
  catch(std::out_of_range const & e)
  {
    std::cout<<"Target \""<<key<<"\" not found. Impossible output."<<std::endl;
  }
  Py_RETURN_NONE;
}


static PyObject * output(PyHDF5 * self, PyObject * args)
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

  geosx::DomainPartition & domain = self->pb_manager->getDomainPartition();

  try
  {
    self->output.at(key)->geosx::TimeHistoryOutput::execute(time, dt, 0, 0, 0, domain);
  }
  catch(std::out_of_range const & e)
  {
    std::cout<<"Target \""<<key<<"\" not found. Impossible output."<<std::endl;
  }
  Py_RETURN_NONE;
}

static PyMethodDef PyHDF5_methods[] = {
{ "collect", (PyCFunction) collect, METH_VARARGS, "wrapper to routine TimeHistoryCollection::execute" },
{ "output", (PyCFunction) output, METH_VARARGS, "wrapper to routine TimeHistoryOutput::execute"},
{ nullptr, nullptr, 0, nullptr }        /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PyHDF5Type = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pysolver.HDF5",
  .tp_basicsize = sizeof( PyHDF5 ),
  .tp_itemsize = 0,
  .tp_dealloc = (destructor)PyHDF5_dealloc,
  .tp_repr = PyHDF5_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PyHDF5::docString,
  .tp_methods = PyHDF5_methods,
  .tp_init = (initproc)PyHDF5_init,
  .tp_new = PyHDF5_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS

static PyModuleDef pyhdf5module = {
    PyModuleDef_HEAD_INIT,
    "pyhdf5",
    "pyhdf5 module for HDF5",
    -1,
    NULL, NULL, NULL, NULL, NULL
};


PyMODINIT_FUNC
PyInit_pyhdf5(void)
{
    PyObject* module = PyModule_Create(&pyhdf5module);

    if (PyType_Ready(&PyHDF5Type) < 0)
        return nullptr;

    if (module == nullptr)
        return nullptr;

    Py_INCREF(&PyHDF5Type);
    PyModule_AddObject(module, "HDF5", (PyObject *)&PyHDF5Type);

    return module;
}

PyTypeObject * getPyHDF5Type()
{ return &PyHDF5Type; }

}
}
