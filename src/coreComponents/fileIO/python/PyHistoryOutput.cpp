#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "fileIO/Outputs/TimeHistoryOutput.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"

#include "PyHistoryOutputType.hpp"
#include "dataRepository/python/PyGroupType.hpp"


#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PyHistoryOutput is not initialized.", nullptr )


namespace geosx
{

namespace python
{

struct PyHistoryOutput
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to TimeHistoryOutput.";

  geosx::TimeHistoryOutput * group;
};


static PyObject * PyHistoryOutput_new( PyTypeObject *type, PyObject *args, PyObject *kwds )
{
  GEOSX_UNUSED_VAR( args, kwds );
  PyHistoryOutput *self;

  self = (PyHistoryOutput *)type->tp_alloc( type, 0 );
  if( self != nullptr )
  {
    self->group = nullptr;
  }

  return (PyObject *)self;
}


static PyObject * PyHistoryOutput_repr( PyObject * const obj ) noexcept
{
  PyHistoryOutput const * const pyHistoryOutput = LvArray::python::convert< PyHistoryOutput >( obj, getPyHistoryOutputType() );
  if( pyHistoryOutput == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pyHistoryOutput );

  string repr;
  string const path = pyHistoryOutput->group->getPath();
  string const type = LvArray::system::demangle( typeid( *(pyHistoryOutput->group) ).name() );
  repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );

}


static PyObject * output( PyHistoryOutput * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  double time;
  double dt;

  if( !PyArg_ParseTuple( args, "dd", &time, &dt ) )
  {
    return nullptr;
  }

  geosx::GeosxState * g_state = &getGlobalState();
  geosx::DomainPartition & domain = g_state->getProblemManager().getDomainPartition();

  int cycleNumber = int(time/dt);
  try
  {
    self->group->cleanup( time, cycleNumber, 0, 0, domain );
  }
  catch( std::out_of_range const & e )
  {
    std::cout<<"Target not found. Impossible output."<<std::endl;
  }
  Py_RETURN_NONE;
}

static PyObject * setOutputName( PyHistoryOutput * self, PyObject * args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  PyObject * unicodePath;
  if( !PyArg_ParseTuple( args, "U", &unicodePath ) )
  {
    return nullptr;
  }

  LvArray::python::PyObjectRef<> asciiPath { PyUnicode_AsASCIIString( unicodePath ) };
  if( asciiPath == nullptr )
  {
    return nullptr;
  }

  char const * const path = PyBytes_AsString( asciiPath );
  if( path == nullptr )
  {
    return nullptr;
  }

  try
  {
    self->group->setFileName( path );
  }
  catch( std::out_of_range const & e )
  {
    std::cout<<"Target not found."<<std::endl;
  }
  Py_RETURN_NONE;
}


static PyObject * reinit( PyHistoryOutput * self, PyObject *args )
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOSX_UNUSED_VAR( args );

  self->group->reinit();

  Py_RETURN_NONE;
}


static PyMethodDef PyHistoryOutput_methods[] = {
  { "output", (PyCFunction) output, METH_VARARGS, "wrapper to routine TimeHistoryOutput::execute"},
  { "setOutputName", (PyCFunction) setOutputName, METH_VARARGS, "wrapper to routine TimeHistoryOutput::setFileName()"},
  { "reinit", (PyCFunction) reinit, METH_VARARGS, "reinitialization function"},
  { nullptr, nullptr, 0, nullptr }      /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PyHistoryOutputType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pygeosx.HistoryOutput",
  .tp_basicsize = sizeof( PyHistoryOutput ),
  .tp_itemsize = 0,
  .tp_repr = PyHistoryOutput_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PyHistoryOutput::docString,
  .tp_methods = PyHistoryOutput_methods,
  .tp_base = getPyGroupType(),
  .tp_new = PyHistoryOutput_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS

PyTypeObject * getPyHistoryOutputType()
{ return &PyHistoryOutputType; }

}
}
