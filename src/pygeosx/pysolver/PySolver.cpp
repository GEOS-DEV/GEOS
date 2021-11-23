#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Source includes
#include "../pygeosx.hpp"
#include "PySolver.hpp"
#include "../PyGroup.hpp"

#define VERIFY_NON_NULL_SELF( self ) \
  PYTHON_ERROR_IF( self == nullptr, PyExc_RuntimeError, "Passed a nullptr as self.", nullptr )

#define VERIFY_INITIALIZED( self ) \
  PYTHON_ERROR_IF( self->group == nullptr, PyExc_RuntimeError, "The PySolver is not initialized.", nullptr )


namespace geosx
{
  extern std::unique_ptr< GeosxState > g_state;

  extern bool g_alreadyInitialized;

namespace python
{

struct PySolver
{
  PyObject_HEAD

  static constexpr char const * docString =
    "A Python interface to geosx::SolverBase.";

  geosx::SolverBase * group;
  geosx::ProblemManager * pb_manager;
};

static void PySolver_dealloc(PySolver* self)
{
    Py_TYPE(self)->tp_free((PyObject*)self);
}



static PyObject * PySolver_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  GEOSX_UNUSED_VAR( args, kwds );
  PySolver *self;
  geosx::ProblemManager * pb_manager = &(g_state->getProblemManager());

    self = (PySolver *)type->tp_alloc(type, 0);
    if (self != nullptr)
    {
      self->group = nullptr;
      self->pb_manager = pb_manager;

      if (self->pb_manager == nullptr)
      {
	Py_DECREF(self);
	return nullptr;
      }

    }

    return (PyObject *)self;
}



static int PySolver_init(PySolver *self, PyObject *args, PyObject *kwds)
{
  GEOSX_UNUSED_VAR( kwds );

  PYTHON_ERROR_IF( !g_alreadyInitialized, PyExc_RuntimeError, "GEOSX must be initialized first. Call pygeosx.initialize().", -1 );

  geosx::SolverBase *group=nullptr;

  PyObject * unicodePath;
  if( !PyArg_ParseTuple( args, "U", &unicodePath ) )
  {
    return -1;
  }

  LvArray::python::PyObjectRef<> asciiPath { PyUnicode_AsASCIIString( unicodePath ) };
  if( asciiPath == nullptr )
  {
    return -1;
  }

  char const * const path = PyBytes_AsString( asciiPath );
  if( path == nullptr )
  {
    return -1;
  }

  try
  {
    geosx::EventManager & eventManager = self->pb_manager->getEventManager();

    eventManager.forSubGroups< geosx::EventBase >( [&]( geosx::EventBase & subEvent )
    {
      subEvent.getTargetReferences();
    } );

    geosx::EventBase * subEvent = nullptr;
    for(int currentSubEvent = 0; currentSubEvent<eventManager.numSubGroups(); ++currentSubEvent )
    {

      subEvent = static_cast< geosx::EventBase * >( eventManager.getSubGroups()[currentSubEvent]);
      if (subEvent->getEventName() == path)
      {
	break;
      }
      else
      {
	subEvent = nullptr;
      }
    }
    PYTHON_ERROR_IF( subEvent == nullptr, PyExc_RuntimeError, "Target not found", -1 );

    group = static_cast<geosx::SolverBase*>(subEvent->getEventTarget());
    self->group = group;

    return 0;
  }
  catch( std::domain_error const & e )
  {
    // If no default return value was specified then this results in a Python exception.
    PyErr_SetString( PyExc_KeyError, e.what() );
    return -1;
  }

}


static PyObject * PySolver_repr( PyObject * const obj ) noexcept
{
  PySolver const * const pySolver = LvArray::python::convert< PySolver >( obj, getPySolverType() );
  if( pySolver == nullptr )
  {
    return nullptr;
  }

  VERIFY_INITIALIZED( pySolver );

  string const path = pySolver->group->getPath();
  string const type = LvArray::system::demangle( typeid( *(pySolver->group) ).name() );
  string const repr = path + " ( " + type + " )";

  return PyUnicode_FromString( repr.c_str() );
}



static PyObject * postProcessInput(PySolver* self, PyObject* args)
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOSX_UNUSED_VAR( args);

  self->group->postProcessInput();

  Py_RETURN_NONE;
}


static PyObject * initPostInitialConditions(PySolver* self, PyObject* args)
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );
  GEOSX_UNUSED_VAR( args);

  self->group->initializePostInitialConditionsPreSubGroups();

  Py_RETURN_NONE;
}


static PyObject * explicitStep(PySolver * self, PyObject * args)
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  double time;
  double dt;
  if( !PyArg_ParseTuple( args, "dd", &time, &dt ) )
  {
    return nullptr;
  }

  geosx::DomainPartition & domain = self->pb_manager->getDomainPartition();

  self->group->explicitStep(time, dt, 0, domain);

  Py_RETURN_NONE;
}



static PyObject * linearImplicitStep(PySolver * self, PyObject * args)
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  double time;
  double dt;
  if( !PyArg_ParseTuple( args, "dd", &time, &dt ) )
  {
    return nullptr;
  }

  geosx::DomainPartition & domain = self->pb_manager->getDomainPartition();

  self->group->linearImplicitStep(time, dt, 0, domain);

  Py_RETURN_NONE;
}


static PyObject * nonlinearImplicitStep(PySolver * self, PyObject * args)
{
  VERIFY_NON_NULL_SELF( self );
  VERIFY_INITIALIZED( self );

  double time;
  double dt;
  if( !PyArg_ParseTuple( args, "dd", &time, &dt ) )
  {
    return nullptr;
  }

  geosx::DomainPartition & domain = self->pb_manager->getDomainPartition();

  self->group->nonlinearImplicitStep(time, dt, 0, domain);

  Py_RETURN_NONE;
}


static PyMethodDef PySolver_methods[] = {
{ "explicitStep", (PyCFunction) explicitStep, METH_VARARGS, "explicit Step" },
{ "linearImplicitStep", (PyCFunction) linearImplicitStep, METH_VARARGS, "linear implicit step" },
{ "nonlinearImplicitStep", (PyCFunction) nonlinearImplicitStep, METH_VARARGS, "non linear implicit step" },
{ "postProcessInput", (PyCFunction) postProcessInput, METH_NOARGS, "post processing input"},
{ "initPostInitialConditions", (PyCFunction) initPostInitialConditions, METH_NOARGS, "call initializePostInitialConditionsPreSubGroup"},
{ "get_wrapper", (PyCFunction) PyGroup_getWrapper<PySolver>, METH_VARARGS, PyGroup_getWrapperDocString },
{ nullptr, nullptr, 0, nullptr }        /* Sentinel */
};


/**
 * Initialize the module object for Python with the exported functions
 */

BEGIN_ALLOW_DESIGNATED_INITIALIZERS

static PyTypeObject PySolverType = {
  PyVarObject_HEAD_INIT( nullptr, 0 )
    .tp_name = "pysolver.Solver",
  .tp_basicsize = sizeof( PySolver ),
  .tp_itemsize = 0,
  .tp_dealloc = (destructor)PySolver_dealloc,
  .tp_repr = PySolver_repr,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = PySolver::docString,
  .tp_methods = PySolver_methods,
  .tp_init = (initproc)PySolver_init,
  .tp_new = PySolver_new,
};

END_ALLOW_DESIGNATED_INITIALIZERS

static PyModuleDef pysolvermodule = {
    PyModuleDef_HEAD_INIT,
    "pysolver",
    "pysolver module for SolverBase",
    -1,
    NULL, NULL, NULL, NULL, NULL
};


PyMODINIT_FUNC
PyInit_pysolver(void)
{
    PyObject* module = PyModule_Create(&pysolvermodule);

    if (PyType_Ready(&PySolverType) < 0)
        return nullptr;

    if (module == nullptr)
        return nullptr;

    Py_INCREF(&PySolverType);
    PyModule_AddObject(module, "Solver", (PyObject *)&PySolverType);

    return module;
}

PyTypeObject * getPySolverType()
{ return &PySolverType; }

}
}
